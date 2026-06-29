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
    6.  Establish the Z floor (minimum elevation or sea level). In Rectangular
        mode fill NoData cells with the floor value. In Tight mode keep NoData
        as NaN so the boundary can be traced.
    7.  Build a watertight triangle mesh:
            · Terrain surface  — 2 triangles per valid quad (+Z normal)
            · Boundary walls   — 2 triangles per boundary edge, all 4 directions
                                 (each wall faces outward from the valid data)
            · Bottom face      — 2 triangles per valid quad (-Z normal, reversed
                                 winding vs. terrain; exactly mirrors the terrain
                                 footprint so every wall-bottom edge is shared)
    8.  Compute per-triangle normals via cross product (float64 precision).
    9.  Write a binary STL in 100 K-triangle chunks using a vectorised
        byte-interleaving trick to insert the 2-byte attribute field without a
        Python loop.

    Footprint modes
    ---------------
    Rectangular : NoData cells are filled with the floor elevation. The model
                  has a flat rectangular base and four straight perimeter walls.
    Tight       : NoData cells (always at the raster edge, never interior) are
                  treated as void. Walls follow the irregular data boundary,
                  hugging the actual DEM outline. Produces a smaller, more
                  printable model for non-rectangular survey areas.

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

        # 6 — Z floor (vertical reference for the base of the model)
        p6 = arcpy.Parameter(
            displayName="Z Floor Reference",
            name="z_floor",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p6.filter.type = "ValueList"
        p6.filter.list = [
            "Minimum Elevation",
            "Sea Level (0)",
        ]
        p6.value = "Minimum Elevation"

        # 7 — Model footprint (rectangular vs. tight boundary)
        p7 = arcpy.Parameter(
            displayName="Model Footprint",
            name="model_footprint",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p7.filter.type = "ValueList"
        p7.filter.list = [
            "Rectangular",
            "Tight (Follows DEM Boundary)",
        ]
        p7.value = "Rectangular"

        return [p0, p1, p2, p3, p4, p5, p6, p7]

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

        in_dem          = parameters[0].valueAsText
        out_stl         = parameters[1].valueAsText
        max_bed         = float(parameters[2].value)
        vert_exag       = float(parameters[3].value)
        min_detail      = float(parameters[4].value)
        base_thick      = float(parameters[5].value)
        z_floor_mode    = parameters[6].valueAsText
        footprint_mode  = parameters[7].valueAsText

        tight = (footprint_mode == "Tight (Follows DEM Boundary)")

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

        # ── Step 5 / 8 — Z floor and elevation array ─────────────────────────
        messages.addMessage("Step 5/8 — Processing elevation values...")
        nan_mask     = np.isnan(arr)
        elev_min_raw = float(np.nanmin(arr))
        elev_max_raw = float(np.nanmax(arr))

        if z_floor_mode == "Sea Level (0)":
            z_ref = min(elev_min_raw, 0.0)   # sea level maps to base_thick
        else:  # "Minimum Elevation"
            z_ref = elev_min_raw

        if not tight:
            # Rectangular: fill NoData with the floor value so every cell is valid
            if nan_mask.any():
                arr[nan_mask] = z_ref

        # Convert elevations to mm above the print plate.
        # z=0          → bottom face (print plate contact)
        # z=base_thick → floor elevation (minimum terrain or sea level)
        # z=base_thick + relief_mm → terrain peaks
        # NaN cells survive the linear transform as NaN (tight mode only).
        z_mm = (arr - z_ref) * xy_scale * vert_exag + base_thick

        relief_mm      = float(np.nanmax(z_mm)) - base_thick
        total_z_mm     = float(np.nanmax(z_mm))
        elev_range_raw = elev_max_raw - elev_min_raw
        messages.addMessage(f"  Elev range     : {elev_min_raw:.2f} – {elev_max_raw:.2f}  (Δ{elev_range_raw:.2f} map units)")
        messages.addMessage(f"  Terrain relief : {relief_mm:.2f} mm  (x{vert_exag:.2f} exag)")
        messages.addMessage(f"  Total height   : {total_z_mm:.2f} mm  (terrain + {base_thick:.1f} mm base)")

        # ── Step 6 / 8 — Build triangle mesh ─────────────────────────────────
        messages.addMessage("Step 6/8 — Building watertight triangle mesh...")

        # Model-space vertex coordinates
        #   X[j] = j × cell_mm_w           (west = 0)
        #   Y[i] = (nr-1-i) × cell_mm_h    (row 0 = north = max Y; row nr-1 = south = 0)
        cols_x = np.arange(nc, dtype=np.float64) * cell_mm_w
        rows_y = (nr - 1 - np.arange(nr, dtype=np.float64)) * cell_mm_h

        X_grid, Y_grid = np.meshgrid(cols_x, rows_y)   # (nr, nc)
        Z_grid = z_mm                                   # (nr, nc); NaN in tight mode

        # valid_cell[i,j] = True iff all 4 corners of quad (i,j) are non-NaN.
        # Rectangular mode: all-True (NoData filled above).
        # Tight mode: True only where the DEM has data on all four corners.
        valid_cell = (
            ~np.isnan(Z_grid[:-1, :-1]) &
            ~np.isnan(Z_grid[:-1, 1:])  &
            ~np.isnan(Z_grid[1:,  :-1]) &
            ~np.isnan(Z_grid[1:,  1:])
        )

        if not valid_cell.any():
            raise RuntimeError(
                "No valid terrain quads found after resampling. "
                "Check the input DEM for NoData coverage."
            )

        all_v0, all_v1, all_v2 = [], [], []

        # Flat integer index arrays for all valid quads — used by terrain and bottom.
        ii, jj = np.where(valid_cell)   # quad row, quad col

        def _gv(ri, ci):
            """Vertex array (N,3) from grid row indices ri and col indices ci."""
            return np.stack([X_grid[ri, ci], Y_grid[ri, ci], Z_grid[ri, ci]], axis=1)

        # ── 6a. Terrain surface ──────────────────────────────────────────────
        # Quad corners for quad (i,j):
        #   A = (i,   j  )  NW    B = (i,   j+1)  NE
        #   D = (i+1, j  )  SW    C = (i+1, j+1)  SE
        #
        # Winding (CCW from above → outward +Z normal):
        #   Tri1: A, D, C   (NW→SW→SE)   n_z = cw·ch > 0 ✓
        #   Tri2: A, C, B   (NW→SE→NE)   n_z = cw·ch > 0 ✓

        A = _gv(ii,     jj)
        B = _gv(ii,     jj + 1)
        C = _gv(ii + 1, jj + 1)
        D = _gv(ii + 1, jj)

        all_v0 += [A, A]
        all_v1 += [D, C]
        all_v2 += [C, B]
        del A, B, C, D

        # ── 6b–6e. Boundary walls ────────────────────────────────────────────
        # A wall is emitted on each face of a valid quad that borders an invalid
        # or absent quad. In rectangular mode (all quads valid) this produces
        # exactly the four perimeter walls. In tight mode it traces the
        # irregular DEM boundary.
        #
        # Neighbour masks (False where no neighbour exists → border ⟹ wall needed)
        above = np.zeros_like(valid_cell); above[1:,  :] = valid_cell[:-1, :]
        below = np.zeros_like(valid_cell); below[:-1, :] = valid_cell[1:,  :]
        left  = np.zeros_like(valid_cell); left[:,  1:]  = valid_cell[:, :-1]
        right = np.zeros_like(valid_cell); right[:, :-1] = valid_cell[:, 1:]

        need_n = valid_cell & ~above   # north face (+Y outward)
        need_s = valid_cell & ~below   # south face (-Y outward)
        need_w = valid_cell & ~left    # west  face (-X outward)
        need_e = valid_cell & ~right   # east  face (+X outward)

        def _ns_wall(mask, row_of_edge, plus_y):
            # Wall on a horizontal grid edge (connects col j to col j+1).
            # row_of_edge(i) gives the grid row index of the edge for quad row i.
            # Winding verified by cross product:
            #   +Y (north): Tri1 TL,BR,BL  Tri2 TL,TR,BR
            #   -Y (south): Tri1 TL,BL,BR  Tri2 TL,BR,TR
            nonlocal all_v0, all_v1, all_v2
            wi, wj = np.where(mask)
            if len(wi) == 0:
                return
            ri  = row_of_edge(wi)
            TL  = np.stack([X_grid[ri, wj],   Y_grid[ri, wj],   Z_grid[ri, wj]],   axis=1)
            TR  = np.stack([X_grid[ri, wj+1], Y_grid[ri, wj+1], Z_grid[ri, wj+1]], axis=1)
            BL  = np.stack([X_grid[ri, wj],   Y_grid[ri, wj],   np.zeros(len(wi))], axis=1)
            BR  = np.stack([X_grid[ri, wj+1], Y_grid[ri, wj+1], np.zeros(len(wi))], axis=1)
            if plus_y:
                all_v0 += [TL, TL]; all_v1 += [BR, TR]; all_v2 += [BL, BR]
            else:
                all_v0 += [TL, TL]; all_v1 += [BL, BR]; all_v2 += [BR, TR]

        def _ew_wall(mask, col_of_edge, plus_x):
            # Wall on a vertical grid edge (connects row i to row i+1).
            # col_of_edge(j) gives the grid col index of the edge for quad col j.
            # Winding verified by cross product:
            #   -X (west): Tri1 TN,BN,BS  Tri2 TN,BS,TS
            #   +X (east): Tri1 TN,BS,BN  Tri2 TN,TS,BS
            nonlocal all_v0, all_v1, all_v2
            wi, wj = np.where(mask)
            if len(wi) == 0:
                return
            ci  = col_of_edge(wj)
            TN  = np.stack([X_grid[wi,   ci], Y_grid[wi,   ci], Z_grid[wi,   ci]], axis=1)
            TS  = np.stack([X_grid[wi+1, ci], Y_grid[wi+1, ci], Z_grid[wi+1, ci]], axis=1)
            BN  = np.stack([X_grid[wi,   ci], Y_grid[wi,   ci], np.zeros(len(wi))], axis=1)
            BS  = np.stack([X_grid[wi+1, ci], Y_grid[wi+1, ci], np.zeros(len(wi))], axis=1)
            if plus_x:
                all_v0 += [TN, TN]; all_v1 += [BS, TS]; all_v2 += [BN, BS]
            else:
                all_v0 += [TN, TN]; all_v1 += [BN, BS]; all_v2 += [BS, TS]

        _ns_wall(need_n, lambda i: i,     plus_y=True)   # north edge at grid row i
        _ns_wall(need_s, lambda i: i + 1, plus_y=False)  # south edge at grid row i+1
        _ew_wall(need_w, lambda j: j,     plus_x=False)  # west  edge at grid col j
        _ew_wall(need_e, lambda j: j + 1, plus_x=True)   # east  edge at grid col j+1

        # ── 6f. Bottom face — Z=0, outward normal -Z ─────────────────────────
        # Mirrors the terrain: same valid quads, same XY positions, all Z=0,
        # reversed winding. Every wall-bottom edge is thereby shared by exactly
        # one bottom-face triangle — no T-junctions, fully manifold.
        #
        # Winding (reversed from terrain → outward -Z normal):
        #   Tri1: A, C, D   (NW→SE→SW)   n_z = −cw·ch < 0 ✓
        #   Tri2: A, B, C   (NW→NE→SE)   n_z = −cw·ch < 0 ✓

        def _bv(ri, ci):
            return np.stack([X_grid[ri, ci], Y_grid[ri, ci], np.zeros(len(ri))], axis=1)

        bA = _bv(ii,     jj)
        bB = _bv(ii,     jj + 1)
        bC = _bv(ii + 1, jj + 1)
        bD = _bv(ii + 1, jj)

        all_v0 += [bA, bA]
        all_v1 += [bC, bB]
        all_v2 += [bD, bC]
        del bA, bB, bC, bD

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
