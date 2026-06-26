# -*- coding: utf-8 -*-
# =============================================================================
# DEM2STL.pyt  —  ArcGIS Pro Python Toolbox
# Converts a Digital Elevation Model raster to a watertight binary STL file
# suitable for slicing and 3D printing.
#
# Tool:   DEM to STL
# Author: David A. Riggs <david.a.riggs@gmail.com>
# =============================================================================

import arcpy
import numpy as np
import os
import struct
import tempfile


# =============================================================================
class Toolbox(object):
    def __init__(self):
        self.label = "DEM to STL Toolbox"
        self.alias = "dem2stl"
        self.tools = [DEMToSTL]


# =============================================================================
class DEMToSTL(object):
    """
    Converts a DEM raster to a watertight binary STL file for 3-D printing.

    Workflow overview
    -----------------
    1.  Describe the input DEM — get real-world extent and native cell size.
    2.  Compute XY scale (mm per map unit) so the longest axis fills the bed.
    3.  Derive a target cell size from the user's Minimum Detail Size parameter
        (clamped to the native resolution so we never upsample).
    4.  Resample the DEM to the target cell size using BILINEAR interpolation.
        Simplification happens HERE, before triangulation, because triangulating
        a 1 m-resolution km-scale DEM first would waste hundreds of MB of RAM
        on triangles that are far smaller than any printer can resolve.
    5.  Load the resampled array into NumPy; delete the temp raster.
    6.  Fill / mark NoData cells according to user choice.
    7.  Build a watertight triangle mesh:
            · Terrain surface  — 2 triangles per quad (CCW, +Z normal, verified)
            · South wall       — 2 triangles per column strip (-Y normal, verified)
            · North wall       — 2 triangles per column strip (+Y normal, verified)
            · West wall        — 2 triangles per row strip    (-X normal, verified)
            · East wall        — 2 triangles per row strip    (+X normal, verified)
            · Bottom face      — 2 triangles                  (-Z normal, verified)
    8.  Compute per-triangle normals via cross product (float64 precision).
    9.  Write a binary STL in 100 K-triangle chunks using a vectorised
        byte-interleaving trick to insert the 2-byte attribute field without a
        Python loop.

    Coordinate system (model space, mm)
    ------------------------------------
        X : 0 (west edge)  →  (nc-1) x cell_mm_w  (east edge)
        Y : 0 (south edge) →  (nr-1) x cell_mm_h  (north edge)
        Z : 0 (print plate) →  base_thick + terrain_relief  (terrain peaks)

    Array orientation: row 0 = north (max Y), row nr-1 = south (Y=0).
    """

    def __init__(self):
        self.label = "DEM to STL"
        self.description = (
            "Converts a Digital Elevation Model (DEM) raster to a watertight binary "
            "STL file suitable for slicing and 3-D printing. The model is scaled to "
            "fit a specified print-bed size with configurable vertical exaggeration, "
            "minimum detail size (controls mesh density), and solid base thickness."
        )
        self.canRunInBackground = False

    # ------------------------------------------------------------------
    def getParameterInfo(self):

        # 0 — Input DEM
        p0 = arcpy.Parameter(
            displayName="Input DEM",
            name="in_dem",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input",
        )

        # 1 — Output STL file
        p1 = arcpy.Parameter(
            displayName="Output STL File",
            name="out_stl",
            datatype="DEFile",
            parameterType="Required",
            direction="Output",
        )
        p1.filter.list = ["stl"]

        # 2 — Max print-bed dimension
        p2 = arcpy.Parameter(
            displayName="Maximum Print-Bed Dimension (mm)",
            name="max_bed_dim",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p2.value = 200.0

        # 3 — Vertical exaggeration
        p3 = arcpy.Parameter(
            displayName="Vertical Exaggeration Factor",
            name="vert_exag",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p3.value = 1.0

        # 4 — Minimum detail size (controls resampling / mesh density)
        p4 = arcpy.Parameter(
            displayName="Minimum Detail Size (mm)",
            name="min_detail",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p4.value = 0.5

        # 5 — Base thickness
        p5 = arcpy.Parameter(
            displayName="Base Thickness (mm)",
            name="base_thick",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p5.value = 3.0

        # 6 — NoData handling
        p6 = arcpy.Parameter(
            displayName="NoData Handling",
            name="nodata_handling",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p6.filter.type = "ValueList"
        p6.filter.list = [
            "Set to Minimum Elevation",
            "Set to Sea Level (0)",
            "Leave as Void (open mesh)",
        ]
        p6.value = "Set to Minimum Elevation"

        return [p0, p1, p2, p3, p4, p5, p6]

    # ------------------------------------------------------------------
    def isLicensed(self):
        return True

    # ------------------------------------------------------------------
    def updateParameters(self, parameters):
        pass

    # ------------------------------------------------------------------
    def updateMessages(self, parameters):
        p_bed    = parameters[2]
        p_exag   = parameters[3]
        p_detail = parameters[4]
        p_base   = parameters[5]

        if p_bed.value is not None:
            if p_bed.value <= 0:
                p_bed.setErrorMessage("Maximum bed dimension must be greater than 0.")
            elif p_bed.value < 10:
                p_bed.setWarningMessage("Very small bed dimension — the model may lack usable surface detail.")

        if p_exag.value is not None and p_exag.value <= 0:
            p_exag.setErrorMessage("Vertical exaggeration must be greater than 0.")

        if p_detail.value is not None and p_detail.value <= 0:
            p_detail.setErrorMessage("Minimum detail size must be greater than 0.")

        if p_base.value is not None and p_base.value < 0:
            p_base.setErrorMessage("Base thickness cannot be negative.")

    # ------------------------------------------------------------------
    def execute(self, parameters, messages):

        in_dem      = parameters[0].valueAsText
        out_stl     = parameters[1].valueAsText
        max_bed     = float(parameters[2].value)
        vert_exag   = float(parameters[3].value)
        min_detail  = float(parameters[4].value)
        base_thick  = float(parameters[5].value)
        nodata_mode = parameters[6].valueAsText

        if not out_stl.lower().endswith(".stl"):
            out_stl += ".stl"

        # ── Step 1 / 8 — Describe DEM ────────────────────────────────────────
        messages.addMessage("Step 1/8 — Analyzing input DEM...")
        desc    = arcpy.Describe(in_dem)
        ext     = desc.extent
        orig_cs = desc.meanCellWidth   # map units per native cell

        rw_w = ext.width
        rw_h = ext.height
        messages.addMessage(f"  Extent    : {rw_w:.3f} x {rw_h:.3f} map units")
        messages.addMessage(f"  Cell size : {orig_cs:.4f} map units")

        # ── Step 2 / 8 — Scale and target cell size ──────────────────────────
        messages.addMessage("Step 2/8 — Computing scale factors...")
        xy_scale  = max_bed / max(rw_w, rw_h)           # mm per map unit
        target_rw = max(min_detail / xy_scale, orig_cs) # map units; never finer than native

        messages.addMessage(f"  XY scale       : 1:{1/xy_scale:,.0f}  ({xy_scale:.6f} mm/map unit)")
        messages.addMessage(f"  Target cell    : {target_rw:.4f} map units  ({target_rw * xy_scale:.3f} mm in model)")

        # ── Step 3 / 8 — Resample DEM ────────────────────────────────────────
        # Simplification happens HERE, before triangulation.
        # A 1 m / 1 km DEM has ~1 M cells; at a 1:5 000 scale each maps to
        # 0.2 mm — far below any printer resolution. Resampling first gives
        # exactly the right triangle count with no wasted intermediate memory.
        messages.addMessage("Step 3/8 — Resampling DEM to target cell size (BILINEAR)...")
        scratch    = arcpy.env.scratchFolder or tempfile.gettempdir()
        tmp_raster = os.path.join(scratch, "dem2stl_resampled.tif")

        try:
            arcpy.management.Resample(
                in_raster=in_dem,
                out_raster=tmp_raster,
                cell_size=target_rw,
                resampling_type="BILINEAR",
            )

            # ── Step 4 / 8 — Load into NumPy ─────────────────────────────────
            messages.addMessage("Step 4/8 — Loading raster into memory...")

            # RasterToNumPyArray raises ValueError if you pass nodata_to_value=nan
            # for integer pixel types (U8, S16, S32, U16, U32 ...).  Detect the
            # pixel type AFTER resampling (BILINEAR always produces a float output,
            # but the raster's stored type may still be reported as integer on some
            # ArcGIS Pro builds).  The safe strategy: for non-float types, load
            # with the raster's own NoData integer value as sentinel, upcast to
            # float64, then replace the sentinel with NaN.
            tmp_desc = arcpy.Describe(tmp_raster)
            pix_type = getattr(tmp_desc, "pixelType", "F32").upper()
            nd_raw   = tmp_desc.noDataValue   # None if no NoData is defined

            messages.addMessage(f"  Pixel type : {pix_type}")

            if pix_type.startswith("F"):
                # Float raster — NaN fills NoData directly
                arr = arcpy.RasterToNumPyArray(tmp_raster, nodata_to_value=np.nan).astype(np.float64)
            else:
                # Integer raster — load with integer sentinel then swap to NaN.
                # Use the raster's own NoData value; fall back to a large negative
                # number that is safely outside any real elevation range.
                if nd_raw is not None:
                    sentinel = int(float(nd_raw))
                else:
                    sentinel = np.iinfo(np.int32).min   # –2 147 483 648
                arr = arcpy.RasterToNumPyArray(tmp_raster, nodata_to_value=sentinel).astype(np.float64)
                if nd_raw is not None:
                    arr[arr == float(sentinel)] = np.nan

            # arr shape: (nr, nc); row 0 = north (max Y), row nr-1 = south (Y=0)

        finally:
            # Always clean up temp raster, even on error
            if arcpy.Exists(tmp_raster):
                arcpy.management.Delete(tmp_raster)

        nr, nc = arr.shape
        messages.addMessage(f"  Grid : {nc} cols x {nr} rows  ({nc * nr:,} vertices)")

        if nr < 2 or nc < 2:
            raise RuntimeError("Resampled grid is too small (< 2x2). Increase the bed size or reduce the Minimum Detail Size.")

        # Vertex spacing in mm
        # Vertices span the full real-world extent, so:
        #   cell_mm_w = total_model_width  / (num_cols - 1)
        #   cell_mm_h = total_model_height / (num_rows - 1)
        model_w_mm = rw_w * xy_scale
        model_h_mm = rw_h * xy_scale
        cell_mm_w  = model_w_mm / (nc - 1)
        cell_mm_h  = model_h_mm / (nr - 1)
        messages.addMessage(f"  Cell spacing : {cell_mm_w:.3f} x {cell_mm_h:.3f} mm")

        # ── Step 5 / 8 — NoData handling and Z array ─────────────────────────
        messages.addMessage("Step 5/8 — Processing elevation values...")
        nan_mask = np.isnan(arr)
        elev_min_raw = float(np.nanmin(arr))   # original map-unit values
        elev_max_raw = float(np.nanmax(arr))   # before any NoData fill

        if nodata_mode == "Set to Minimum Elevation":
            if nan_mask.any():
                arr[nan_mask] = np.nanmin(arr)
            z_ref = arr.min()

        elif nodata_mode == "Set to Sea Level (0)":
            if nan_mask.any():
                arr[nan_mask] = 0.0
            z_ref = min(np.nanmin(arr), 0.0)   # floor at 0 so sea level maps to base

        else:  # "Leave as Void (open mesh)"
            z_ref = np.nanmin(arr)

        # Convert to mm above the print plate.
        # z=0  → bottom face (print plate contact)
        # z=base_thick → lowest terrain point
        # z=base_thick + relief_mm → highest terrain point
        z_mm = (arr - z_ref) * xy_scale * vert_exag + base_thick

        if nodata_mode == "Leave as Void (open mesh)":
            z_mm[nan_mask] = np.nan   # holes remain open; walls seal the perimeter

        relief_mm  = float(np.nanmax(z_mm)) - base_thick
        total_z_mm = float(np.nanmax(z_mm))
        elev_range_raw = elev_max_raw - elev_min_raw
        messages.addMessage(f"  Elev range     : {elev_min_raw:.2f} – {elev_max_raw:.2f}  (Δ{elev_range_raw:.2f} map units)")
        messages.addMessage(f"  Terrain relief : {relief_mm:.2f} mm  (x{vert_exag:.2f} exag)")
        messages.addMessage(f"  Total height   : {total_z_mm:.2f} mm  (terrain + {base_thick:.1f} mm base)")

        # ── Step 6 / 8 — Build triangle mesh ─────────────────────────────────
        messages.addMessage("Step 6/8 — Building watertight triangle mesh...")

        # Model-space vertex coordinates for the whole grid
        #   X[j] = j x cell_mm_w           (west = 0)
        #   Y[i] = (nr-1-i) x cell_mm_h    (row 0 = north = max Y; row nr-1 = south = 0)
        cols_x = np.arange(nc, dtype=np.float64) * cell_mm_w
        rows_y = (nr - 1 - np.arange(nr, dtype=np.float64)) * cell_mm_h

        X_grid, Y_grid = np.meshgrid(cols_x, rows_y)   # (nr, nc)
        Z_grid = z_mm                                  # (nr, nc), NaN for voids

        Xmax = (nc - 1) * cell_mm_w
        Ymax = (nr - 1) * cell_mm_h

        all_v0, all_v1, all_v2 = [], [], []

        # ── 6a. Terrain surface ──────────────────────────────────────────────
        # For each quad (i,j) the four corners are:
        #   A = (i,   j  )  NW — larger Y, smaller X
        #   B = (i,   j+1)  NE — larger Y, larger X
        #   C = (i+1, j+1)  SE — smaller Y, larger X
        #   D = (i+1, j  )  SW — smaller Y, smaller X
        #
        # Winding (CCW from above → outward normal +Z, verified):
        #   Tri1: A, D, C   (NW → SW → SE)
        #   Tri2: A, C, B   (NW → SE → NE)
        #
        # Cross-product verification (Tri1):
        #   e1 = D-A = (0, -ch, zD-zA);  e2 = C-A = (cw, -ch, zC-zA)
        #   n_z = cw·ch > 0  ✓
        # Cross-product verification (Tri2):
        #   e1 = C-A = (cw, -ch, zC-zA);  e2 = B-A = (cw, 0, zB-zA)
        #   n_z = cw·ch > 0  ✓

        ri0, ri1 = slice(None, nr - 1), slice(1, None)
        ci0, ci1 = slice(None, nc - 1), slice(1, None)

        def _v(rs, cs):
            """Return (N, 3) float64 array of vertex positions for a grid slice."""
            return np.stack(
                [X_grid[rs, cs].ravel(), Y_grid[rs, cs].ravel(), Z_grid[rs, cs].ravel()],
                axis=1,
            )

        A = _v(ri0, ci0)
        B = _v(ri0, ci1)
        C = _v(ri1, ci1)
        D = _v(ri1, ci0)

        if nodata_mode == "Leave as Void (open mesh)":
            ok1 = ~(np.isnan(A[:, 2]) | np.isnan(D[:, 2]) | np.isnan(C[:, 2]))
            ok2 = ~(np.isnan(A[:, 2]) | np.isnan(C[:, 2]) | np.isnan(B[:, 2]))
            all_v0 += [A[ok1], A[ok2]]
            all_v1 += [D[ok1], C[ok2]]
            all_v2 += [C[ok1], B[ok2]]
        else:
            all_v0 += [A, A]
            all_v1 += [D, C]
            all_v2 += [C, B]

        del A, B, C, D

        # Helper — replace NaN edge elevations with base_thick so walls close
        def _safe(z_arr):
            v = z_arr.copy()
            v[np.isnan(v)] = base_thick
            return v

        # Shared column strip index arrays (length nc-1)
        j_idx = np.arange(nc - 1, dtype=np.float64)
        sx    = j_idx * cell_mm_w           # west X of each column strip
        sx1   = (j_idx + 1) * cell_mm_w     # east X

        # Shared row strip Y arrays (length nr-1)
        i_idx = np.arange(nr - 1, dtype=np.float64)
        Yi    = (nr - 1 - i_idx) * cell_mm_h   # north Y of each row strip
        Yi1   = (nr - 2 - i_idx) * cell_mm_h   # south Y

        # ── 6b. South wall — row nr-1, Y=0, outward normal -Y ────────────────
        # Quad per column strip (viewed from south, looking +Y):
        #   TL=(j·cw,  0, Zj)   TR=((j+1)·cw, 0, Zj+1)
        #   BL=(j·cw,  0, 0 )   BR=((j+1)·cw, 0, 0   )
        #
        # Tri1: TL, BL, BR  →  n=(0, -Zj·cw, 0) → -Y ✓
        # Tri2: TL, BR, TR  →  n=(0, -Zj+1·cw, 0) → -Y ✓

        Zs   = _safe(Z_grid[nr - 1, :])
        szj  = Zs[:nc - 1]
        szj1 = Zs[1:]
        sy   = np.zeros(nc - 1)
        sbz  = np.zeros(nc - 1)

        all_v0 += [np.c_[sx,  sy, szj],  np.c_[sx,  sy, szj]]
        all_v1 += [np.c_[sx,  sy, sbz],  np.c_[sx1, sy, sbz]]
        all_v2 += [np.c_[sx1, sy, sbz],  np.c_[sx1, sy, szj1]]

        # ── 6c. North wall — row 0, Y=Ymax, outward normal +Y ────────────────
        # Winding reversed vs. south wall:
        # Tri1: TL, BR, BL  →  n=(0, +Zj·cw, 0) → +Y ✓
        # Tri2: TL, TR, BR  →  n=(0, +Zj+1·cw, 0) → +Y ✓

        Zn   = _safe(Z_grid[0, :])
        nzj  = Zn[:nc - 1]
        nzj1 = Zn[1:]
        ny   = np.full(nc - 1, Ymax)
        nbz  = np.zeros(nc - 1)

        all_v0 += [np.c_[sx,  ny, nzj],  np.c_[sx,  ny, nzj]]
        all_v1 += [np.c_[sx1, ny, nbz],  np.c_[sx1, ny, nzj1]]
        all_v2 += [np.c_[sx,  ny, nbz],  np.c_[sx1, ny, nbz]]

        # ── 6d. West wall — col 0, X=0, outward normal -X ────────────────────
        # For each row strip i:
        #   TN=(0, Yi,  Zi)   TS=(0, Yi1, Zi+1)   BN=(0, Yi, 0)   BS=(0, Yi1, 0)
        #
        # Tri1: TN, BN, BS  →  e1=(0,0,-Zi), e2=(0,-ch,0)
        #                       n_x = -Zi·ch < 0 → -X ✓
        # Tri2: TN, BS, TS  →  n_x = -Zi+1·ch < 0 → -X ✓

        Zw   = _safe(Z_grid[:, 0])
        wzi  = Zw[:nr - 1]
        wzi1 = Zw[1:]
        wx   = np.zeros(nr - 1)
        wbz  = np.zeros(nr - 1)

        all_v0 += [np.c_[wx, Yi,  wzi],  np.c_[wx, Yi,  wzi]]
        all_v1 += [np.c_[wx, Yi,  wbz],  np.c_[wx, Yi1, wbz]]
        all_v2 += [np.c_[wx, Yi1, wbz],  np.c_[wx, Yi1, wzi1]]

        # ── 6e. East wall — col nc-1, X=Xmax, outward normal +X ──────────────
        # Winding reversed vs. west wall:
        # Tri1: TN, BS, BN  →  n_x = +Zi·ch > 0 → +X ✓
        # Tri2: TN, TS, BS  →  n_x = +Zi+1·ch > 0 → +X ✓

        Ze   = _safe(Z_grid[:, nc - 1])
        ezi  = Ze[:nr - 1]
        ezi1 = Ze[1:]
        ex   = np.full(nr - 1, Xmax)
        ebz  = np.zeros(nr - 1)

        all_v0 += [np.c_[ex, Yi,  ezi],  np.c_[ex, Yi,  ezi]]
        all_v1 += [np.c_[ex, Yi1, ebz],  np.c_[ex, Yi1, ezi1]]
        all_v2 += [np.c_[ex, Yi,  ebz],  np.c_[ex, Yi1, ebz]]

        # ── 6f. Bottom face — Z=0, outward normal -Z ─────────────────────────
        # Corners: SW=(0,0,0) NW=(0,Ymax,0) NE=(Xmax,Ymax,0) SE=(Xmax,0,0)
        #
        # Tri1: SW, NW, NE  →  n_z = -Ymax·Xmax < 0 → -Z ✓
        # Tri2: SW, NE, SE  →  n_z = -Ymax·Xmax < 0 → -Z ✓

        SW = np.array([[0.0,  0.0,  0.0]])
        NW = np.array([[0.0,  Ymax, 0.0]])
        NE = np.array([[Xmax, Ymax, 0.0]])
        SE = np.array([[Xmax, 0.0,  0.0]])

        all_v0 += [SW, SW]
        all_v1 += [NW, NE]
        all_v2 += [NE, SE]

        # ── Step 7 / 8 — Compute per-triangle normals ────────────────────────
        messages.addMessage("Step 7/8 — Computing normals...")

        V0 = np.concatenate(all_v0).astype(np.float32)
        V1 = np.concatenate(all_v1).astype(np.float32)
        V2 = np.concatenate(all_v2).astype(np.float32)
        del all_v0, all_v1, all_v2

        total_tris = len(V0)
        messages.addMessage(f"  Total triangles : {total_tris:,}")

        # Cross products in float64 for numerical precision, then normalise
        e1  = V1.astype(np.float64) - V0.astype(np.float64)
        e2  = V2.astype(np.float64) - V0.astype(np.float64)
        nrm = np.cross(e1, e2)
        mag = np.linalg.norm(nrm, axis=1, keepdims=True)
        mag[mag == 0.0] = 1.0   # guard against degenerate (zero-area) triangles
        nrm = (nrm / mag).astype(np.float32)
        del e1, e2, mag

        # ── Step 8 / 8 — Write binary STL ────────────────────────────────────
        # Binary STL record (50 bytes per triangle):
        #   normal    : 3 x float32  (12 bytes)
        #   vertex 0  : 3 x float32  (12 bytes)
        #   vertex 1  : 3 x float32  (12 bytes)
        #   vertex 2  : 3 x float32  (12 bytes)
        #   attribute : uint16        ( 2 bytes) — always 0
        #
        # Written in chunks of CHUNK triangles to limit peak RAM.
        # Byte interleaving is fully vectorised (no Python loop per triangle).
        messages.addMessage(f"Step 8/8 — Writing binary STL: {out_stl}")
        CHUNK = 100_000

        with open(out_stl, "wb") as f:
            # 80-byte ASCII header (informational only)
            hdr = b"DEM to STL | ArcGIS Pro Python Toolbox"
            f.write(hdr.ljust(80, b" "))

            # Triangle count as little-endian uint32
            f.write(struct.pack("<I", total_tris))

            for start in range(0, total_tris, CHUNK):
                end = min(start + CHUNK, total_tris)
                n   = end - start

                # Stack into (n, 12) float32:
                #   columns 0–2   = normal, 3–5 = V0, 6–8 = V1, 9–11 = V2
                block = np.hstack([
                    nrm[start:end],
                    V0[start:end],
                    V1[start:end],
                    V2[start:end],
                ]).astype(np.float32)   # n x 12 floats = n x 48 bytes

                # Build (n, 50) uint8 output buffer; place float bytes in cols
                # 0–47; cols 48–49 stay 0 (attribute = 0).
                buf = np.zeros((n, 50), dtype=np.uint8)
                buf[:, :48] = np.frombuffer(block.tobytes(), dtype=np.uint8).reshape(n, 48)
                f.write(buf.tobytes())

        fsize_mb = os.path.getsize(out_stl) / (1024 ** 2)
        messages.addMessage("")
        messages.addMessage("✓  STL written successfully!")
        messages.addMessage(f"  Model dimensions : {model_w_mm:.1f} x {model_h_mm:.1f} x {total_z_mm:.2f} mm  (W x D x H)")
        messages.addMessage(f"  Triangles        : {total_tris:,}")
        messages.addMessage(f"  File size        : {fsize_mb:.2f} MB")
        messages.addMessage(f"  Output           : {out_stl}")

    # ------------------------------------------------------------------
    def postExecute(self, parameters):
        return
