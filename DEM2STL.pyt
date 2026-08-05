# -*- coding: utf-8 -*-
# =============================================================================
# DEM2STL.pyt  -  ArcGIS Pro Python Toolbox
# Converts a Digital Elevation Model raster to a mesh file suitable for 3D printing.
#
# Tools:
#  - DEM to STL
#  - Split DEM to STL
#  - DEM to 3MF
#  - Split DEM to 3MF
#
# Author: David A. Riggs <david.a.riggs@gmail.com>
# =============================================================================

import arcpy
import numpy as np
import os
import struct
import tempfile
import zipfile
from typing import Optional


# Common filament colors, name → uppercase hex (Bambu's parser requires uppercase).
PAINT_COLORS = {
    "White":  "#FFFFFF",
    "Black":  "#000000",
    "Grey":   "#808080",
    "Red":    "#FF0000",
    "Orange": "#FF8000",
    "Yellow": "#FFFF00",
    "Green":  "#00AA00",
    "Blue":   "#0000FF",
    "Purple": "#800080",
    "Brown":  "#8B4513",
    "Pink":   "#FF69B4",
    "Cyan":   "#00FFFF",
}


# =============================================================================
# Module-level helpers shared by all tools
# =============================================================================

def _detect_elev_units(dem_path: str) -> Optional[str]:
    """Returns 'Meters', 'Feet', or None if no VCS is defined on the DEM.

    The arcpy VerticalCoordinateSystem object exposes .name (not .linearUnitName,
    which belongs to SpatialReference for XY units).  VCS names encode the unit
    when feet are used, e.g. 'NAVD_1988_height_ftUS', 'NGVD_1929_Height_US_Ft'.
    A VCS that exists but contains no foot token is assumed to be in meters.
    """
    try:
        sr = arcpy.Describe(dem_path).spatialReference
        if sr and sr.VCS:
            name = sr.VCS.name.lower()
            if any(tok in name for tok in ("foot", "feet", "_ft", "ftus")):
                return "Feet"
            return "Meters"
    except Exception:
        pass
    return None


def _load_raster_array(raster_path: str) -> np.ndarray:
    """Load a raster into a float64 NumPy array with NoData replaced by NaN.

    RasterToNumPyArray raises ValueError if nodata_to_value=nan is used with
    integer pixel types.  This function detects the pixel type and uses an
    integer sentinel for non-float rasters, then swaps to NaN afterward.
    """
    desc     = arcpy.Describe(raster_path)
    pix_type = getattr(desc, "pixelType", "F32").upper()
    nd_raw   = desc.noDataValue

    if pix_type.startswith("F"):
        return arcpy.RasterToNumPyArray(raster_path, nodata_to_value=np.nan).astype(np.float64)

    sentinel = int(float(nd_raw)) if nd_raw is not None else np.iinfo(np.int32).min
    arr = arcpy.RasterToNumPyArray(raster_path, nodata_to_value=sentinel).astype(np.float64)
    if nd_raw is not None:
        arr[arr == float(sentinel)] = np.nan
    return arr


def _build_mesh(z_mm: np.ndarray, cell_mm_w: float, cell_mm_h: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a watertight triangle mesh from a z_mm elevation array.

    Returns three float64 (N, 3) arrays (V0, V1, V2) - the first, second, and
    third vertex of each of the N triangles, in that order.  Winding is
    counter-clockwise as seen from outside the solid, so the right-hand-rule
    normal points outward.  This convention is exactly what both binary STL and
    3MF expect, so the same mesh feeds either writer unchanged.

    z_mm    : float64 (nr, nc) array; NaN cells are treated as void (no geometry).
              In Rectangular mode the caller fills NaN with z_ref before calling here,
              so all cells are valid and only perimeter walls are generated.
              In Tight mode NaN cells remain and walls follow the irregular boundary.
    cell_mm_w / cell_mm_h : vertex spacing in mm (X / Y axes respectively).

    Because every shared grid corner is computed from the same X_grid/Y_grid/Z_grid
    entry, coincident vertices are bit-identical across triangles - so a downstream
    writer can weld them into an indexed manifold mesh with an exact (round-free)
    np.unique.
    """
    nr, nc = z_mm.shape

    # Model-space vertex coordinates
    #   X[j] = j x cell_mm_w           (west = 0)
    #   Y[i] = (nr-1-i) x cell_mm_h    (row 0 = north = max Y; row nr-1 = south = 0)
    cols_x = np.arange(nc, dtype=np.float64) * cell_mm_w
    rows_y = (nr - 1 - np.arange(nr, dtype=np.float64)) * cell_mm_h

    X_grid, Y_grid = np.meshgrid(cols_x, rows_y)
    Z_grid = z_mm

    # valid_cell[i,j] = True iff all 4 corners of quad (i,j) are non-NaN.
    valid_cell = (
        ~np.isnan(Z_grid[:-1, :-1]) &
        ~np.isnan(Z_grid[:-1, 1:])  &
        ~np.isnan(Z_grid[1:,  :-1]) &
        ~np.isnan(Z_grid[1:,  1:])
    )

    if not valid_cell.any():
        raise RuntimeError("No valid terrain quads found. Check the input DEM for NoData coverage.")

    all_v0, all_v1, all_v2 = [], [], []
    ii, jj = np.where(valid_cell)

    def _gv(ri, ci):
        return np.stack([X_grid[ri, ci], Y_grid[ri, ci], Z_grid[ri, ci]], axis=1)

    # ── Terrain surface ──────────────────────────────────────────────────────
    # Quad corners for quad (i,j):
    #   A = (i,   j  )  NW    B = (i,   j+1)  NE
    #   D = (i+1, j  )  SW    C = (i+1, j+1)  SE
    #
    # Winding (CCW from above → outward +Z normal):
    #   Tri1: A, D, C   (NW→SW→SE)   n_z = cw·ch > 0 ✓
    #   Tri2: A, C, B   (NW→SE→NE)   n_z = cw·ch > 0 ✓
    A = _gv(ii,     jj);     B = _gv(ii,     jj + 1)
    C = _gv(ii + 1, jj + 1); D = _gv(ii + 1, jj)
    all_v0 += [A, A]; all_v1 += [D, C]; all_v2 += [C, B]
    del A, B, C, D

    # ── Boundary walls ───────────────────────────────────────────────────────
    # A wall is emitted on each face of a valid quad that borders an invalid
    # or absent quad. In rectangular mode (all quads valid) this produces
    # exactly the four perimeter walls. In tight mode it traces the
    # irregular DEM boundary.
    above = np.zeros_like(valid_cell); above[1:,  :] = valid_cell[:-1, :]
    below = np.zeros_like(valid_cell); below[:-1, :] = valid_cell[1:,  :]
    left  = np.zeros_like(valid_cell); left[:,  1:]  = valid_cell[:, :-1]
    right = np.zeros_like(valid_cell); right[:, :-1] = valid_cell[:, 1:]

    need_n = valid_cell & ~above
    need_s = valid_cell & ~below
    need_w = valid_cell & ~left
    need_e = valid_cell & ~right

    def _ns_wall(mask, row_of_edge, plus_y):
        # Wall on a horizontal grid edge (connects col j to col j+1).
        # Winding: +Y (north): Tri1 TL,BR,BL  Tri2 TL,TR,BR
        #          -Y (south): Tri1 TL,BL,BR  Tri2 TL,BR,TR
        nonlocal all_v0, all_v1, all_v2
        wi, wj = np.where(mask)
        if len(wi) == 0:
            return
        ri = row_of_edge(wi)
        TL = np.stack([X_grid[ri, wj],   Y_grid[ri, wj],   Z_grid[ri, wj]],    axis=1)
        TR = np.stack([X_grid[ri, wj+1], Y_grid[ri, wj+1], Z_grid[ri, wj+1]],  axis=1)
        BL = np.stack([X_grid[ri, wj],   Y_grid[ri, wj],   np.zeros(len(wi))], axis=1)
        BR = np.stack([X_grid[ri, wj+1], Y_grid[ri, wj+1], np.zeros(len(wi))], axis=1)
        if plus_y:
            all_v0 += [TL, TL]; all_v1 += [BR, TR]; all_v2 += [BL, BR]
        else:
            all_v0 += [TL, TL]; all_v1 += [BL, BR]; all_v2 += [BR, TR]

    def _ew_wall(mask, col_of_edge, plus_x):
        # Wall on a vertical grid edge (connects row i to row i+1).
        # Winding: -X (west): Tri1 TN,BN,BS  Tri2 TN,BS,TS
        #          +X (east): Tri1 TN,BS,BN  Tri2 TN,TS,BS
        nonlocal all_v0, all_v1, all_v2
        wi, wj = np.where(mask)
        if len(wi) == 0:
            return
        ci = col_of_edge(wj)
        TN = np.stack([X_grid[wi,   ci], Y_grid[wi,   ci], Z_grid[wi,   ci]],  axis=1)
        TS = np.stack([X_grid[wi+1, ci], Y_grid[wi+1, ci], Z_grid[wi+1, ci]],  axis=1)
        BN = np.stack([X_grid[wi,   ci], Y_grid[wi,   ci], np.zeros(len(wi))], axis=1)
        BS = np.stack([X_grid[wi+1, ci], Y_grid[wi+1, ci], np.zeros(len(wi))], axis=1)
        if plus_x:
            all_v0 += [TN, TN]; all_v1 += [BS, TS]; all_v2 += [BN, BS]
        else:
            all_v0 += [TN, TN]; all_v1 += [BN, BS]; all_v2 += [BS, TS]

    _ns_wall(need_n, lambda i: i,     plus_y=True)
    _ns_wall(need_s, lambda i: i + 1, plus_y=False)
    _ew_wall(need_w, lambda j: j,     plus_x=False)
    _ew_wall(need_e, lambda j: j + 1, plus_x=True)

    # ── Bottom face - Z=0, outward normal -Z ────────────────────────────────
    # Mirrors the terrain: same valid quads, same XY positions, all Z=0,
    # reversed winding. Every wall-bottom edge is shared by exactly one
    # bottom triangle - no T-junctions, fully manifold.
    #
    # Winding (reversed from terrain → outward -Z normal):
    #   Tri1: A, C, D   (NW→SE→SW)   n_z = −cw·ch < 0 ✓
    #   Tri2: A, B, C   (NW→NE→SE)   n_z = −cw·ch < 0 ✓
    def _bv(ri, ci):
        return np.stack([X_grid[ri, ci], Y_grid[ri, ci], np.zeros(len(ri))], axis=1)

    bA = _bv(ii,     jj);     bB = _bv(ii,     jj + 1)
    bC = _bv(ii + 1, jj + 1); bD = _bv(ii + 1, jj)
    all_v0 += [bA, bA]; all_v1 += [bC, bB]; all_v2 += [bD, bC]
    del bA, bB, bC, bD

    V0 = np.concatenate(all_v0)
    V1 = np.concatenate(all_v1)
    V2 = np.concatenate(all_v2)
    return V0, V1, V2


def _write_stl(V0: np.ndarray, V1: np.ndarray, V2: np.ndarray, out_stl: str) -> tuple[int, float]:
    """Write a triangle mesh (V0/V1/V2 vertex arrays) as a binary STL file.

    Returns (total_triangles, file_size_mb).

    Binary STL record (50 bytes per triangle):
        normal   : 3 x float32  (12 bytes)
        vertex 0 : 3 x float32  (12 bytes)
        vertex 1 : 3 x float32  (12 bytes)
        vertex 2 : 3 x float32  (12 bytes)
        attribute: uint16       ( 2 bytes) - always 0
    """
    V0 = V0.astype(np.float32)
    V1 = V1.astype(np.float32)
    V2 = V2.astype(np.float32)
    total_tris = len(V0)

    # Cross products in float64 for numerical precision, then normalise
    e1  = V1.astype(np.float64) - V0.astype(np.float64)
    e2  = V2.astype(np.float64) - V0.astype(np.float64)
    nrm = np.cross(e1, e2)
    mag = np.linalg.norm(nrm, axis=1, keepdims=True)
    mag[mag == 0.0] = 1.0
    nrm = (nrm / mag).astype(np.float32)
    del e1, e2, mag

    CHUNK = 100_000
    with open(out_stl, "wb") as f:
        f.write(b"DEM to STL | ArcGIS Pro Python Toolbox".ljust(80, b" "))
        f.write(struct.pack("<I", total_tris))
        for start in range(0, total_tris, CHUNK):
            end = min(start + CHUNK, total_tris)
            n   = end - start
            block = np.hstack([
                nrm[start:end], V0[start:end], V1[start:end], V2[start:end],
            ]).astype(np.float32)
            buf = np.zeros((n, 50), dtype=np.uint8)
            buf[:, :48] = np.frombuffer(block.tobytes(), dtype=np.uint8).reshape(n, 48)
            f.write(buf.tobytes())

    fsize_mb = os.path.getsize(out_stl) / (1024 ** 2)
    return total_tris, fsize_mb


def _array_to_stl(z_mm: np.ndarray, cell_mm_w: float, cell_mm_h: float, out_stl: str) -> tuple[int, float]:
    """Build a watertight mesh from an elevation array and write it as binary STL.

    Convenience wrapper: _build_mesh followed by _write_stl.
    Returns (total_triangles, file_size_mb).
    """
    V0, V1, V2 = _build_mesh(z_mm, cell_mm_w, cell_mm_h)
    return _write_stl(V0, V1, V2, out_stl)


def _prepare_paint_src(paint_fc: str, scratch: str, buffer_map_units: Optional[float],
                       suffix: str = "") -> tuple:
    """Buffer/prepare a paint feature class geometry for rasterization.

    Returns (src_fc_path, tmp_fc_list).  src_fc_path is paint_fc itself when no
    buffering is needed, or a buffered copy in the memory workspace.  tmp_fc_list
    contains any temp FCs created - the caller must delete them when done.

    suffix makes names unique when multiple paint layers are pre-buffered concurrently
    (e.g. Split DEM to 3MF pre-buffers all layers before its piece loop).
    """
    buf_name      = f"memory/dem23mf_paintbuf{suffix}"
    boundary_name = f"memory/dem23mf_paintboundary{suffix}"
    shape_type    = arcpy.Describe(paint_fc).shapeType
    tmp_fcs       = []
    try:
        if buffer_map_units and shape_type in ("Polyline", "Point", "Multipoint"):
            if arcpy.Exists(buf_name):
                arcpy.management.Delete(buf_name)
            arcpy.analysis.Buffer(paint_fc, buf_name, buffer_map_units, dissolve_option="ALL")
            tmp_fcs.append(buf_name)
            return buf_name, tmp_fcs
        elif buffer_map_units and shape_type == "Polygon":
            # Width specified → paint the buffered boundary band, not the full interior.
            # Extract polygon rings as lines, then buffer those lines.
            if arcpy.Exists(boundary_name):
                arcpy.management.Delete(boundary_name)
            arcpy.management.PolygonToLine(paint_fc, boundary_name)
            tmp_fcs.append(boundary_name)
            if arcpy.Exists(buf_name):
                arcpy.management.Delete(buf_name)
            arcpy.analysis.Buffer(boundary_name, buf_name, buffer_map_units, dissolve_option="ALL")
            tmp_fcs.append(buf_name)
            return buf_name, tmp_fcs
        else:
            return paint_fc, tmp_fcs
    except Exception:
        for fc in tmp_fcs:
            try:
                arcpy.management.Delete(fc)
            except Exception:
                pass
        raise


def _rasterize_src(src_fc: str, snap_raster: str, arr_shape: tuple,
                   scratch: str) -> np.ndarray:
    """Rasterize a pre-prepared feature class to a boolean mask matching the snap raster.

    arcpy.env snap/extent/cellSize are set to align output pixels exactly to the snap
    raster so indices correspond 1-to-1 with the resampled DEM array.  Environment
    variables are restored on exit.  Returns a bool (nr, nc) array matching arr_shape.
    """
    tmp_ras   = os.path.join(scratch, "dem23mf_paintras.tif")
    snap_desc = arcpy.Describe(snap_raster)

    old_snap   = arcpy.env.snapRaster
    old_extent = arcpy.env.extent
    old_cs     = arcpy.env.cellSize
    old_crs    = arcpy.env.outputCoordinateSystem
    try:
        arcpy.env.snapRaster             = snap_raster
        arcpy.env.extent                 = snap_desc.extent
        arcpy.env.cellSize               = snap_desc.meanCellWidth
        arcpy.env.outputCoordinateSystem = snap_desc.spatialReference

        oid_field = arcpy.Describe(src_fc).OIDFieldName
        if arcpy.Exists(tmp_ras):
            arcpy.management.Delete(tmp_ras)
        arcpy.conversion.FeatureToRaster(src_fc, oid_field, tmp_ras, snap_desc.meanCellWidth)
    finally:
        arcpy.env.snapRaster             = old_snap
        arcpy.env.extent                 = old_extent
        arcpy.env.cellSize               = old_cs
        arcpy.env.outputCoordinateSystem = old_crs

    try:
        # Use -1 as nodata sentinel (not 0) because shapefile FID starts at 0,
        # which would otherwise be indistinguishable from the nodata fill value.
        paint_arr = arcpy.RasterToNumPyArray(tmp_ras, nodata_to_value=-1).astype(int)
        mask = (paint_arr != -1)

        # Guard against off-by-one from raster snapping
        nr, nc = arr_shape
        mr, mc = mask.shape
        if (mr, mc) != (nr, nc):
            result = np.zeros((nr, nc), dtype=bool)
            result[:min(nr, mr), :min(nc, mc)] = mask[:min(nr, mr), :min(nc, mc)]
            mask = result

        return mask
    finally:
        if arcpy.Exists(tmp_ras):
            try:
                arcpy.management.Delete(tmp_ras)
            except Exception:
                pass


def _rasterize_paint_layer(paint_fc: str, snap_raster: str, arr_shape: tuple,
                           scratch: str, buffer_map_units: Optional[float]) -> np.ndarray:
    """Rasterize a vector paint layer to a boolean mask matching the snap raster grid.

    Lines and points are pre-buffered by buffer_map_units (in the raster's map-unit CRS)
    before rasterization.  Polygons are rasterized directly when buffer_map_units is None;
    when buffer_map_units is provided, only a band of that width along the polygon boundary
    is painted (PolygonToLine then Buffer, matching the line-layer path).  Returns a bool
    (nr, nc) array - True where paint features are present - with the same shape as arr_shape.

    Chains _prepare_paint_src (buffering) with _rasterize_src (FeatureToRaster + mask
    conversion) and cleans up temp FCs on exit.
    """
    src, tmp_fcs = _prepare_paint_src(paint_fc, scratch, buffer_map_units)
    try:
        return _rasterize_src(src, snap_raster, arr_shape, scratch)
    finally:
        for fc in tmp_fcs:
            try:
                arcpy.management.Delete(fc)
            except Exception:
                pass


def _write_3mf(V0: np.ndarray, V1: np.ndarray, V2: np.ndarray, out_3mf: str,
               object_name: str, base_color: str = "#FFFFFF",
               tri_color_indices: Optional[np.ndarray] = None,
               paint_colors: Optional[list] = None) -> tuple[int, int, float]:
    """Write a triangle mesh (V0/V1/V2 vertex arrays) as a standards-compliant 3MF file.

    A 3MF file is an OPC (ZIP) package with three required members:
        [Content_Types].xml   - declares .rels and .model content types
        _rels/.rels           - points the package root at the 3D model part
        3D/3dmodel.model      - model geometry plus optional color metadata

    Unlike STL's triangle soup, 3MF uses an *indexed* mesh: a shared vertex list
    plus triangles referencing vertices by index.  Coincident vertices are
    bit-identical (see _build_mesh), so np.unique produces a fully-welded manifold.

    Color support (3MF Materials and Properties extension, namespace "m"):
        base_color        - hex color (#RRGGBB) assigned to the whole object.  Always
                            emitted; BambuStudio 1.9+ maps this to an AMS filament slot
                            via its "Standard 3MF Color Parsing" dialog.
        tri_color_indices - optional int array length total_tris; 0 = base_color,
                            1..N index into paint_colors.  Non-zero terrain triangles
                            receive a per-triangle color override; walls and bottom
                            always use index 0 (base_color).
        paint_colors      - list of N hex strings matching indices 1..N.

    Uppercase hex digits are required - Bambu's parser silently drops lowercase.

    Returns (total_triangles, vertex_count, file_size_mb).
    """
    total_tris = len(V0)

    # Weld the per-triangle vertices into a shared, indexed vertex list.
    # stacked = [ all v0 ; all v1 ; all v2 ]  →  triangle i uses rows i, N+i, 2N+i.
    stacked = np.concatenate([V0, V1, V2], axis=0)
    verts, inv = np.unique(stacked, axis=0, return_inverse=True)
    inv = inv.ravel()
    tris = np.stack([inv[:total_tris], inv[total_tris:2 * total_tris], inv[2 * total_tris:]], axis=1)
    vert_count = len(verts)

    safe_name = (object_name or "DEM").replace("&", "&amp;").replace("<", "&lt;").replace('"', "&quot;")

    # ── Vertex XML ───────────────────────────────────────────────────────────
    vert_xml = "".join(
        '<vertex x="%.4f" y="%.4f" z="%.4f"/>' % (x, y, z) for x, y, z in verts
    )

    # ── Triangle XML ─────────────────────────────────────────────────────────
    # Painted triangles carry an explicit pid/p1 index; unpainted triangles
    # inherit the object-level default (index 0 = base_color).
    use_multicolor = tri_color_indices is not None and paint_colors
    mat_ns = ' xmlns:m="http://schemas.microsoft.com/3dmanufacturing/material/2015/02"'

    if use_multicolor:
        indices = tri_color_indices
        tri_xml = "".join(
            ('<triangle v1="%d" v2="%d" v3="%d" pid="2" p1="%d"/>' %
             (int(tris[k, 0]), int(tris[k, 1]), int(tris[k, 2]), int(indices[k]))
             if indices[k] > 0
             else '<triangle v1="%d" v2="%d" v3="%d"/>' %
             (int(tris[k, 0]), int(tris[k, 1]), int(tris[k, 2])))
            for k in range(total_tris)
        )
    else:
        tri_xml = "".join(
            '<triangle v1="%d" v2="%d" v3="%d"/>' % (a, b, c) for a, b, c in tris
        )

    # ── Assemble the 3D/3dmodel.model XML ────────────────────────────────────
    # The 3MF Materials extension (xmlns:m) is always included so the model
    # carries a base color.  When paint_colors is provided, one entry per layer
    # is appended after the base color entry.
    if use_multicolor:
        color_entries = '<m:color color="' + base_color + '"/>'
        for pc in paint_colors:
            color_entries += '<m:color color="' + pc + '"/>'
        color_group_xml = '<m:colorgroup id="2">' + color_entries + '</m:colorgroup>'
    else:
        color_group_xml = (
            '<m:colorgroup id="2">'
            '<m:color color="' + base_color + '"/>'
            '</m:colorgroup>'
        )

    # Object inherits base color from the group; individual tris may override.
    obj_attrs = ' pid="2" p1="0"'

    model_xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\r\n'
        '<model unit="millimeter" xml:lang="en-US"'
        ' xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02"'
        + mat_ns + '>'
        '<metadata name="Application">DEM to 3MF | ArcGIS Pro Python Toolbox</metadata>'
        '<metadata name="Title">' + safe_name + '</metadata>'
        '<resources>'
        + color_group_xml +
        '<object id="1" type="model" name="' + safe_name + '"' + obj_attrs + '>'
        '<mesh>'
        '<vertices>' + vert_xml + '</vertices>'
        '<triangles>' + tri_xml + '</triangles>'
        '</mesh>'
        '</object>'
        '</resources>'
        '<build><item objectid="1"/></build>'
        '</model>'
    )

    # ── OPC package ──────────────────────────────────────────────────────────
    content_types = (
        '<?xml version="1.0" encoding="UTF-8"?>\r\n'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>'
        '</Types>'
    )

    rels = (
        '<?xml version="1.0" encoding="UTF-8"?>\r\n'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rel0" Target="/3D/3dmodel.model"'
        ' Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>'
        '</Relationships>'
    )

    with zipfile.ZipFile(out_3mf, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", rels)
        z.writestr("3D/3dmodel.model", model_xml)

    fsize_mb = os.path.getsize(out_3mf) / (1024 ** 2)
    return total_tris, vert_count, fsize_mb


def _array_to_3mf(z_mm: np.ndarray, cell_mm_w: float, cell_mm_h: float, out_3mf: str,
                  object_name: str, quad_paint_masks: Optional[list] = None,
                  base_color: str = "#FFFFFF",
                  paint_colors: Optional[list] = None) -> tuple[int, int, float]:
    """Build a watertight mesh from an elevation array and write it as a 3MF file.

    Convenience wrapper: _build_mesh followed by _write_3mf.

    quad_paint_masks : optional list of bool (nr-1, nc-1) arrays, one per paint layer.
                       Layers are applied in order; later entries win on overlap.
                       Walls and bottom always use base_color.  When None or empty,
                       the model is written in base_color only.
    paint_colors     : list of hex strings (#RRGGBB), one per entry in quad_paint_masks.

    Returns (total_triangles, vertex_count, file_size_mb).
    """
    V0, V1, V2 = _build_mesh(z_mm, cell_mm_w, cell_mm_h)

    tri_color_indices = None
    if quad_paint_masks and paint_colors:
        # valid_cell mirrors the quad validity check in _build_mesh.
        valid_cell = (
            ~np.isnan(z_mm[:-1, :-1]) & ~np.isnan(z_mm[:-1, 1:]) &
            ~np.isnan(z_mm[1:,  :-1]) & ~np.isnan(z_mm[1:,  1:])
        )
        total_tris        = len(V0)
        n_terrain         = 2 * int(valid_cell.sum())   # two tris per valid quad
        tri_color_indices = np.zeros(total_tris, dtype=np.int32)
        # Apply each paint layer in order; later layers win on overlap.
        # _build_mesh lays out terrain tris as [all Tri1s, all Tri2s] - np.tile
        # replicates the per-quad mask into both halves correctly.
        for k, qmask in enumerate(quad_paint_masks, start=1):
            quad_colors = qmask[valid_cell]
            layer_tris  = np.tile(quad_colors, 2)
            tri_color_indices[:n_terrain][layer_tris] = k

    return _write_3mf(V0, V1, V2, out_3mf, object_name, base_color,
                      tri_color_indices, paint_colors)


# =============================================================================
class Toolbox(object):
    def __init__(self):
        self.label = "DEM to STL Toolbox"
        self.alias = "dem2stl"
        self.tools = [DEMToSTL, SplitDEMToSTL, DEMTo3MF, SplitDEMTo3MF]


# =============================================================================
class DEMToSTL(object):
    """
    Converts a DEM raster to a watertight binary STL file for 3-D printing.

    Workflow overview
    -----------------
    1.  Describe the input DEM - get real-world extent and native cell size.
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
    7.  Build a watertight triangle mesh and write binary STL via _array_to_stl.

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

        # 0 - Input DEM
        p0 = arcpy.Parameter(
            displayName="Input DEM",
            name="in_dem",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input",
        )

        # 1 - Elevation Units
        p1 = arcpy.Parameter(
            displayName="Elevation Units",
            name="elev_units",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p1.filter.type = "ValueList"
        p1.filter.list = ["Meters", "Feet"]
        # No default - forces an explicit choice when auto-detection fails.

        # 2 - Output STL file
        p2 = arcpy.Parameter(
            displayName="Output STL File",
            name="out_stl",
            datatype="DEFile",
            parameterType="Required",
            direction="Output",
        )
        p2.filter.list = ["stl"]

        # 3 - Max print-bed dimension
        p3 = arcpy.Parameter(
            displayName="Maximum Print-Bed Dimension (mm)",
            name="max_bed_dim",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p3.value = 180.0

        # 4 - Vertical exaggeration
        p4 = arcpy.Parameter(
            displayName="Vertical Exaggeration Factor",
            name="vert_exag",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p4.value = 1.0

        # 5 - Minimum detail size (controls resampling / mesh density)
        p5 = arcpy.Parameter(
            displayName="Minimum Detail Size (mm)",
            name="min_detail",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p5.value = 0.2

        # 6 - Base thickness
        p6 = arcpy.Parameter(
            displayName="Base Thickness (mm)",
            name="base_thick",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p6.value = 3.0

        # 7 - Z floor (vertical reference for the base of the model)
        p7 = arcpy.Parameter(
            displayName="Z Floor Reference",
            name="z_floor",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p7.filter.type = "ValueList"
        p7.filter.list = ["Sea Level (0)", "Minimum Elevation"]
        p7.value = "Sea Level (0)"

        # 8 - Model footprint (rectangular vs. tight boundary)
        p8 = arcpy.Parameter(
            displayName="Model Footprint",
            name="model_footprint",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p8.filter.type = "ValueList"
        p8.filter.list = ["Tight (Follows DEM Boundary)", "Rectangular"]
        p8.value = "Tight (Follows DEM Boundary)"

        return [p0, p1, p2, p3, p4, p5, p6, p7, p8]

    # ------------------------------------------------------------------
    def isLicensed(self):
        return True

    # ------------------------------------------------------------------
    def updateParameters(self, parameters):
        dem_param   = parameters[0]
        units_param = parameters[1]

        # Re-detect elevation units every time the DEM changes.
        # altered=True + hasBeenValidated=False means the DEM was just modified
        # by the user, so we overwrite whatever was previously in the units field
        # (intentionally non-sticky: re-selecting a different DEM always re-detects).
        if dem_param.value and dem_param.altered and not dem_param.hasBeenValidated:
            units_param.value = _detect_elev_units(dem_param.valueAsText)

    # ------------------------------------------------------------------
    def updateMessages(self, parameters):
        dem_param   = parameters[0]
        units_param = parameters[1]
        p_bed    = parameters[3]
        p_exag   = parameters[4]
        p_detail = parameters[5]
        p_base   = parameters[6]

        if dem_param.value and not units_param.value:
            units_param.setWarningMessage(
                "Elevation units could not be auto-detected from this DEM's spatial "
                "reference (no Vertical Coordinate System defined). Please select."
            )

        if p_bed.value is not None:
            if p_bed.value <= 0:
                p_bed.setErrorMessage("Maximum bed dimension must be greater than 0.")
            elif p_bed.value < 10:
                p_bed.setWarningMessage("Very small bed dimension - the model may lack usable surface detail.")

        if p_exag.value is not None and p_exag.value <= 0:
            p_exag.setErrorMessage("Vertical exaggeration must be greater than 0.")

        if p_detail.value is not None and p_detail.value <= 0:
            p_detail.setErrorMessage("Minimum detail size must be greater than 0.")

        if p_base.value is not None and p_base.value < 0:
            p_base.setErrorMessage("Base thickness cannot be negative.")

    # ------------------------------------------------------------------
    def execute(self, parameters, messages):

        in_dem         = parameters[0].valueAsText
        elev_units     = parameters[1].valueAsText
        out_stl        = parameters[2].valueAsText
        max_bed        = float(parameters[3].value)
        vert_exag      = float(parameters[4].value)
        min_detail     = float(parameters[5].value)
        base_thick     = float(parameters[6].value)
        z_floor_mode   = parameters[7].valueAsText
        footprint_mode = parameters[8].valueAsText

        tight = (footprint_mode == "Tight (Follows DEM Boundary)")

        if not out_stl.lower().endswith(".stl"):
            out_stl += ".stl"

        # ── Step 1/6 - Describe DEM ──────────────────────────────────────────
        messages.addMessage("Step 1/6 - Analyzing input DEM...")
        desc    = arcpy.Describe(in_dem)
        sr      = desc.spatialReference
        ext     = desc.extent
        orig_cs = desc.meanCellWidth

        rw_w = ext.width
        rw_h = ext.height
        messages.addMessage(f"  Extent    : {rw_w:.3f} x {rw_h:.3f} map units")
        messages.addMessage(f"  Cell size : {orig_cs:.4f} map units")

        # z_to_xy converts elevation values to XY map units before xy_scale is
        # applied - handles the common case where Z is in feet but XY is meters.
        xy_unit_name = sr.linearUnitName.lower() if sr else "meter"
        xy_is_feet   = "foot" in xy_unit_name or "feet" in xy_unit_name

        if elev_units is None:
            elev_units = "Feet" if xy_is_feet else "Meters"
            messages.addMessage(f"  Elevation units : not specified - assuming {elev_units} (same as XY units)")

        z_is_feet = (elev_units == "Feet")
        if z_is_feet and not xy_is_feet:
            z_to_xy = 0.3048
        elif not z_is_feet and xy_is_feet:
            z_to_xy = 1.0 / 0.3048
        else:
            z_to_xy = 1.0

        if z_to_xy != 1.0:
            messages.addMessage(
                f"  Elevation units : {elev_units} → XY units are "
                f"{'feet' if xy_is_feet else 'meters'}  "
                f"(z_to_xy factor: {z_to_xy:.7f})"
            )
        else:
            messages.addMessage(f"  Elevation units : {elev_units} (matches XY units, no conversion)")

        # ── Step 2/6 - Scale and target cell size ────────────────────────────
        messages.addMessage("Step 2/6 - Computing scale factors...")
        xy_scale  = max_bed / max(rw_w, rw_h)
        target_rw = max(min_detail / xy_scale, orig_cs)

        messages.addMessage(f"  XY scale    : 1:{1/xy_scale:,.0f}  ({xy_scale:.6f} mm/map unit)")
        messages.addMessage(f"  Target cell : {target_rw:.4f} map units  ({target_rw * xy_scale:.3f} mm in model)")

        # ── Step 3/6 - Resample DEM ──────────────────────────────────────────
        # Simplification happens HERE, before triangulation.
        # A 1 m / 1 km DEM has ~1 M cells; at a 1:5 000 scale each maps to
        # 0.2 mm - far below any printer resolution. Resampling first gives
        # exactly the right triangle count with no wasted intermediate memory.
        messages.addMessage("Step 3/6 - Resampling DEM to target cell size (BILINEAR)...")
        scratch    = arcpy.env.scratchFolder or tempfile.gettempdir()
        tmp_raster = os.path.join(scratch, "dem2stl_resampled.tif")

        try:
            arcpy.management.Resample(
                in_raster=in_dem,
                out_raster=tmp_raster,
                cell_size=target_rw,
                resampling_type="BILINEAR",
            )

            # ── Step 4/6 - Load into NumPy ───────────────────────────────────
            messages.addMessage("Step 4/6 - Loading raster into memory...")
            arr = _load_raster_array(tmp_raster)

        finally:
            if arcpy.Exists(tmp_raster):
                arcpy.management.Delete(tmp_raster)

        nr, nc = arr.shape
        messages.addMessage(f"  Grid : {nc} cols x {nr} rows  ({nc * nr:,} vertices)")

        if nr < 2 or nc < 2:
            raise RuntimeError("Resampled grid is too small (< 2x2). Increase the bed size or reduce the Minimum Detail Size.")

        model_w_mm = rw_w * xy_scale
        model_h_mm = rw_h * xy_scale
        cell_mm_w  = model_w_mm / (nc - 1)
        cell_mm_h  = model_h_mm / (nr - 1)
        messages.addMessage(f"  Cell spacing : {cell_mm_w:.3f} x {cell_mm_h:.3f} mm")

        # ── Step 5/6 - Z floor and elevation array ───────────────────────────
        messages.addMessage("Step 5/6 - Processing elevation values...")
        nan_mask     = np.isnan(arr)
        elev_min_raw = float(np.nanmin(arr))
        elev_max_raw = float(np.nanmax(arr))

        if z_floor_mode == "Sea Level (0)":
            z_ref = min(elev_min_raw, 0.0)
        else:
            z_ref = elev_min_raw

        if not tight:
            if nan_mask.any():
                arr[nan_mask] = z_ref

        # z_to_xy normalises elevation values into XY map units before xy_scale
        # converts map units → mm.  When Z and XY share the same unit, z_to_xy=1.
        z_mm = (arr - z_ref) * z_to_xy * xy_scale * vert_exag + base_thick

        relief_mm      = float(np.nanmax(z_mm)) - base_thick
        total_z_mm     = float(np.nanmax(z_mm))
        elev_range_raw = elev_max_raw - elev_min_raw
        messages.addMessage(f"  Elev range     : {elev_min_raw:.2f} – {elev_max_raw:.2f}  (Δ{elev_range_raw:.2f} {elev_units.lower()})")
        messages.addMessage(f"  Terrain relief : {relief_mm:.2f} mm  (x{vert_exag:.2f} exag)")
        messages.addMessage(f"  Total height   : {total_z_mm:.2f} mm  (terrain + {base_thick:.1f} mm base)")

        # ── Step 6/6 - Build mesh and write STL ──────────────────────────────
        messages.addMessage(f"Step 6/6 - Building mesh and writing STL: {out_stl}")
        total_tris, fsize_mb = _array_to_stl(z_mm, cell_mm_w, cell_mm_h, out_stl)

        messages.addMessage("")
        messages.addMessage("✓  STL written successfully!")
        messages.addMessage(f"  Model dimensions : {model_w_mm:.1f} x {model_h_mm:.1f} x {total_z_mm:.2f} mm  (W x D x H)")
        messages.addMessage(f"  Triangles        : {total_tris:,}")
        messages.addMessage(f"  File size        : {fsize_mb:.2f} MB")
        messages.addMessage(f"  Output           : {out_stl}")

    # ------------------------------------------------------------------
    def postExecute(self, parameters):
        return


# =============================================================================
class SplitDEMToSTL(object):
    """
    Splits a DEM by a polygon feature class and converts each piece to its own
    watertight binary STL file for 3-D printing.

    Scale is computed globally: the largest polygon piece is measured first,
    then a single XY scale is derived so that piece fits within the Maximum
    Print-Bed Dimension.  All other pieces use the same scale, so they print
    at the correct relative size and can be physically arranged together.

    Each output file is named:  <Base Name>_<Name Attribute value>.stl
    or, when no Name Attribute is chosen:  <Base Name>_1.stl, _2.stl, …

    Workflow
    --------
    1.  Describe the DEM - establish extent, cell size, and z_to_xy factor.
    2.  Scan all polygon extents (intersected with the DEM extent) to find the
        largest dimension; compute a single global XY scale and target cell size.
    3.  For each polygon:
          a. Clip the DEM to the polygon boundary (arcpy.management.Clip).
          b. Resample to the target cell size.
          c. Load into NumPy, compute z_mm, build mesh, write STL.

    Z Floor Reference
    -----------------
    "Minimum Elevation (per piece)" : each piece independently stands up from
        its own local minimum.  Pieces may differ in absolute model height but
        every model uses its full base thickness over its own lowest terrain.
    "Sea Level (0)" : z_ref = min(local_min, 0) for every piece - pieces that
        share the same datum will have consistent relative heights and can be
        arranged side-by-side to form a continuous landscape.
    """

    def __init__(self):
        self.label = "Split DEM to STL"
        self.description = (
            "Splits a DEM raster by a polygon feature class and converts each piece "
            "to a separate watertight binary STL file. A single consistent scale is "
            "chosen so the largest piece fits within the specified print-bed size; "
            "all pieces are produced at that same scale."
        )
        self.canRunInBackground = False

    # ------------------------------------------------------------------
    def getParameterInfo(self):

        # 0 - Input DEM
        p0 = arcpy.Parameter(
            displayName="Input DEM",
            name="in_dem",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input",
        )

        # 1 - Elevation Units
        p1 = arcpy.Parameter(
            displayName="Elevation Units",
            name="elev_units",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p1.filter.type = "ValueList"
        p1.filter.list = ["Meters", "Feet"]

        # 2 - Split Polygon Feature Class
        p2 = arcpy.Parameter(
            displayName="Split Polygon Feature Class",
            name="split_polys",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input",
        )
        p2.filter.list = ["Polygon"]

        # 3 - Name Attribute (optional field from p2 used in output filenames)
        p3 = arcpy.Parameter(
            displayName="Name Attribute",
            name="name_attr",
            datatype="Field",
            parameterType="Optional",
            direction="Input",
        )
        p3.parameterDependencies = [p2.name]
        p3.filter.list = ["Short", "Long", "BigInteger", "Double", "Single", "Text"]

        # 4 - Output Folder
        p4 = arcpy.Parameter(
            displayName="Output Folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input",
        )

        # 5 - Output Base Name
        p5 = arcpy.Parameter(
            displayName="Output Base Name",
            name="base_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )

        # 6 - Max print-bed dimension
        p6 = arcpy.Parameter(
            displayName="Maximum Print-Bed Dimension (mm)",
            name="max_bed_dim",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p6.value = 180.0

        # 7 - Vertical exaggeration
        p7 = arcpy.Parameter(
            displayName="Vertical Exaggeration Factor",
            name="vert_exag",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p7.value = 1.0

        # 8 - Minimum detail size
        p8 = arcpy.Parameter(
            displayName="Minimum Detail Size (mm)",
            name="min_detail",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p8.value = 0.2

        # 9 - Base thickness
        p9 = arcpy.Parameter(
            displayName="Base Thickness (mm)",
            name="base_thick",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p9.value = 3.0

        # 10 - Z floor reference
        p10 = arcpy.Parameter(
            displayName="Z Floor Reference",
            name="z_floor",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p10.filter.type = "ValueList"
        p10.filter.list = ["Sea Level (0)", "Minimum Elevation (per piece)"]
        p10.value = "Sea Level (0)"

        # 11 - Model footprint
        p11 = arcpy.Parameter(
            displayName="Model Footprint",
            name="model_footprint",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p11.filter.type = "ValueList"
        p11.filter.list = ["Tight (Follows Polygon Boundary)", "Rectangular"]
        p11.value = "Tight (Follows Polygon Boundary)"

        return [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]

    # ------------------------------------------------------------------
    def isLicensed(self):
        return True

    # ------------------------------------------------------------------
    def updateParameters(self, parameters):
        dem_param   = parameters[0]
        units_param = parameters[1]
        poly_param  = parameters[2]
        attr_param  = parameters[3]

        if dem_param.value and dem_param.altered and not dem_param.hasBeenValidated:
            units_param.value = _detect_elev_units(dem_param.valueAsText)

        # Clear stale name attribute selection when the polygon FC changes.
        if poly_param.altered and not poly_param.hasBeenValidated:
            attr_param.value = None

    # ------------------------------------------------------------------
    def updateMessages(self, parameters):
        dem_param   = parameters[0]
        units_param = parameters[1]
        p_bed    = parameters[6]
        p_exag   = parameters[7]
        p_detail = parameters[8]
        p_base   = parameters[9]

        if dem_param.value and not units_param.value:
            units_param.setWarningMessage(
                "Elevation units could not be auto-detected from this DEM's spatial "
                "reference (no Vertical Coordinate System defined). Please select."
            )

        if p_bed.value is not None:
            if p_bed.value <= 0:
                p_bed.setErrorMessage("Maximum bed dimension must be greater than 0.")
            elif p_bed.value < 10:
                p_bed.setWarningMessage("Very small bed dimension - the model may lack usable surface detail.")

        if p_exag.value is not None and p_exag.value <= 0:
            p_exag.setErrorMessage("Vertical exaggeration must be greater than 0.")

        if p_detail.value is not None and p_detail.value <= 0:
            p_detail.setErrorMessage("Minimum detail size must be greater than 0.")

        if p_base.value is not None and p_base.value < 0:
            p_base.setErrorMessage("Base thickness cannot be negative.")

    # ------------------------------------------------------------------
    def execute(self, parameters, messages):

        in_dem         = parameters[0].valueAsText
        elev_units     = parameters[1].valueAsText
        poly_fc        = parameters[2].valueAsText
        name_attr      = parameters[3].valueAsText   # None if not chosen
        out_folder     = parameters[4].valueAsText
        base_name      = parameters[5].valueAsText
        max_bed        = float(parameters[6].value)
        vert_exag      = float(parameters[7].value)
        min_detail     = float(parameters[8].value)
        base_thick     = float(parameters[9].value)
        z_floor_mode   = parameters[10].valueAsText
        footprint_mode = parameters[11].valueAsText

        tight = footprint_mode.startswith("Tight")

        os.makedirs(out_folder, exist_ok=True)

        # ── Step 1/3 - Describe DEM ──────────────────────────────────────────
        messages.addMessage("Step 1/3 - Analyzing input DEM...")
        dem_desc = arcpy.Describe(in_dem)
        sr       = dem_desc.spatialReference
        dem_ext  = dem_desc.extent
        orig_cs  = dem_desc.meanCellWidth

        messages.addMessage(f"  Extent    : {dem_ext.width:.3f} x {dem_ext.height:.3f} map units")
        messages.addMessage(f"  Cell size : {orig_cs:.4f} map units")

        xy_unit_name = sr.linearUnitName.lower() if sr else "meter"
        xy_is_feet   = "foot" in xy_unit_name or "feet" in xy_unit_name

        if elev_units is None:
            elev_units = "Feet" if xy_is_feet else "Meters"
            messages.addMessage(f"  Elevation units : not specified - assuming {elev_units} (same as XY units)")

        z_is_feet = (elev_units == "Feet")
        if z_is_feet and not xy_is_feet:
            z_to_xy = 0.3048
        elif not z_is_feet and xy_is_feet:
            z_to_xy = 1.0 / 0.3048
        else:
            z_to_xy = 1.0

        if z_to_xy != 1.0:
            messages.addMessage(
                f"  Elevation units : {elev_units} → XY units are "
                f"{'feet' if xy_is_feet else 'meters'}  "
                f"(z_to_xy factor: {z_to_xy:.7f})"
            )
        else:
            messages.addMessage(f"  Elevation units : {elev_units} (matches XY units, no conversion)")

        # ── Step 2/3 - Scan polygons, compute global scale ───────────────────
        messages.addMessage("Step 2/3 - Scanning split polygons for consistent scale...")

        cursor_fields = ["OID@", "SHAPE@"]
        if name_attr:
            cursor_fields.append(name_attr)

        oid_field  = arcpy.Describe(poly_fc).OIDFieldName
        scratch    = arcpy.env.scratchFolder or tempfile.gettempdir()
        poly_layer = arcpy.management.MakeFeatureLayer(poly_fc, "split_stl_polys")[0]

        # Pre-filter to polygons that intersect the DEM bounding box.
        # For large continent-wide layers this avoids iterating thousands of
        # irrelevant features. SelectLayerByLocation handles CRS differences.
        total_poly_count = int(arcpy.management.GetCount(poly_layer)[0])
        dem_extent_fc    = "memory/dem2stl_extent_filter"
        try:
            if arcpy.Exists(dem_extent_fc):
                arcpy.management.Delete(dem_extent_fc)
            ring = arcpy.Array([
                arcpy.Point(dem_ext.XMin, dem_ext.YMin),
                arcpy.Point(dem_ext.XMin, dem_ext.YMax),
                arcpy.Point(dem_ext.XMax, dem_ext.YMax),
                arcpy.Point(dem_ext.XMax, dem_ext.YMin),
            ])
            arcpy.management.CopyFeatures([arcpy.Polygon(ring, sr)], dem_extent_fc)
            arcpy.management.SelectLayerByLocation(poly_layer, "INTERSECT", dem_extent_fc)
        except Exception:
            arcpy.management.SelectLayerByAttribute(poly_layer, "CLEAR_SELECTION")
        finally:
            if arcpy.Exists(dem_extent_fc):
                arcpy.management.Delete(dem_extent_fc)

        selected_poly_count = int(arcpy.management.GetCount(poly_layer)[0])
        if selected_poly_count < total_poly_count:
            messages.addMessage(
                f"  Pre-filtered from {total_poly_count:,} to {selected_poly_count:,} "
                f"polygon(s) intersecting DEM bounding box"
            )

        pieces = []
        with arcpy.da.SearchCursor(poly_layer, cursor_fields) as cur:
            for row in cur:
                oid  = row[0]
                geom = row[1]
                nval = str(row[2]) if name_attr else None
                # Project to DEM CRS for extent comparison and overlap detection.
                # Clip uses OID-based selection and reprojects automatically,
                # so only this stored geometry (used for extent math) needs projection.
                try:
                    if sr and sr.name != "Unknown":
                        geom = geom.projectAs(sr)
                except Exception:
                    pass
                pieces.append((oid, geom, nval))

        if not pieces:
            raise RuntimeError("No polygon features overlap the DEM extent.")

        # Find the largest piece dimension (polygon extent clipped to DEM extent).
        # This determines the global scale; all pieces are produced at that scale.
        global_max_dim = 0.0
        for _, geom, _ in pieces:
            pe     = geom.extent
            clip_w = min(pe.XMax, dem_ext.XMax) - max(pe.XMin, dem_ext.XMin)
            clip_h = min(pe.YMax, dem_ext.YMax) - max(pe.YMin, dem_ext.YMin)
            global_max_dim = max(global_max_dim, max(0.0, clip_w), max(0.0, clip_h))

        if global_max_dim == 0.0:
            raise RuntimeError("No polygon features overlap the DEM extent.")

        xy_scale  = max_bed / global_max_dim
        target_rw = max(min_detail / xy_scale, orig_cs)

        messages.addMessage(f"  {len(pieces)} polygon(s) to process")
        messages.addMessage(f"  Largest piece dimension : {global_max_dim:.1f} map units")
        messages.addMessage(f"  XY scale    : 1:{1/xy_scale:,.0f}  ({xy_scale:.6f} mm/map unit)")
        messages.addMessage(f"  Target cell : {target_rw:.4f} map units  ({target_rw * xy_scale:.3f} mm)")

        # ── Step 3/3 - Process each polygon piece ────────────────────────────
        messages.addMessage(f"Step 3/3 - Processing {len(pieces)} piece(s)...")

        skipped = 0

        try:
            for idx, (oid, geom, nval) in enumerate(pieces):
                label    = nval if nval is not None else str(idx + 1)
                safe     = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(label)).strip("_")
                out_stl  = os.path.join(out_folder, f"{base_name}_{safe}.stl")

                messages.addMessage(f"\n  [{idx+1}/{len(pieces)}] {label}")

                # Skip polygons that don't overlap the DEM
                pe      = geom.extent
                clip_x0 = max(pe.XMin, dem_ext.XMin)
                clip_y0 = max(pe.YMin, dem_ext.YMin)
                clip_x1 = min(pe.XMax, dem_ext.XMax)
                clip_y1 = min(pe.YMax, dem_ext.YMax)
                if clip_x1 <= clip_x0 or clip_y1 <= clip_y0:
                    messages.addWarningMessage("    Polygon does not overlap the DEM - skipping.")
                    skipped += 1
                    continue

                rect       = f"{clip_x0} {clip_y0} {clip_x1} {clip_y1}"
                tmp_clip   = os.path.join(scratch, f"dem2stl_clip_{idx}.tif")
                tmp_resamp = os.path.join(scratch, f"dem2stl_resamp_{idx}.tif")

                try:
                    arcpy.management.SelectLayerByAttribute(
                        poly_layer, "NEW_SELECTION", f"{oid_field} = {oid}"
                    )
                    arcpy.management.Clip(
                        in_raster=in_dem,
                        rectangle=rect,
                        out_raster=tmp_clip,
                        in_template_dataset=poly_layer,
                        clipping_geometry="ClippingGeometry",
                    )
                    arcpy.management.Resample(
                        in_raster=tmp_clip,
                        out_raster=tmp_resamp,
                        cell_size=target_rw,
                        resampling_type="BILINEAR",
                    )

                    arr    = _load_raster_array(tmp_resamp)
                    nr, nc = arr.shape

                    if nr < 2 or nc < 2:
                        messages.addWarningMessage("    Piece is too small after resampling (<2x2 grid) - skipping.")
                        skipped += 1
                        continue

                    # Use actual raster extent for accurate cell spacing
                    clip_ext  = arcpy.Describe(tmp_resamp).extent
                    rw_w      = clip_ext.width
                    rw_h      = clip_ext.height
                    cell_mm_w = rw_w * xy_scale / (nc - 1)
                    cell_mm_h = rw_h * xy_scale / (nr - 1)

                    # Z floor
                    nan_mask     = np.isnan(arr)
                    elev_min_raw = float(np.nanmin(arr))
                    z_ref        = min(elev_min_raw, 0.0) if z_floor_mode == "Sea Level (0)" else elev_min_raw

                    if not tight:
                        arr[nan_mask] = z_ref

                    z_mm       = (arr - z_ref) * z_to_xy * xy_scale * vert_exag + base_thick
                    total_z_mm = float(np.nanmax(z_mm))
                    model_w_mm = rw_w * xy_scale
                    model_h_mm = rw_h * xy_scale

                    num_tris, fsize_mb = _array_to_stl(z_mm, cell_mm_w, cell_mm_h, out_stl)

                    messages.addMessage(
                        f"    ✓  {os.path.basename(out_stl)}"
                        f"  ({model_w_mm:.1f}x{model_h_mm:.1f}x{total_z_mm:.2f} mm,"
                        f"  {num_tris:,} tris,  {fsize_mb:.2f} MB)"
                    )

                except Exception as piece_err:
                    messages.addWarningMessage(f"    Failed: {piece_err}")
                    skipped += 1

                finally:
                    for tmp in (tmp_clip, tmp_resamp):
                        if arcpy.Exists(tmp):
                            arcpy.management.Delete(tmp)

        finally:
            arcpy.management.Delete(poly_layer)

        processed = len(pieces) - skipped
        messages.addMessage(f"\n✓  {processed} of {len(pieces)} piece(s) written to: {out_folder}")
        if skipped:
            messages.addMessage(f"  ({skipped} skipped - see warnings above)")

    # ------------------------------------------------------------------
    def postExecute(self, parameters):
        return


# =============================================================================
class DEMTo3MF(object):
    """
    Converts a DEM raster to a watertight 3MF file for 3-D printing.

    This tool is the 3MF counterpart of "DEM to STL". The DEM analysis, scaling,
    resampling, elevation handling, and mesh construction are identical - the two
    tools share _build_mesh - so a given DEM produces the same geometry either way.
    Only the output container differs: STL writes a binary triangle soup, while
    this tool writes a standards-compliant 3MF (an OPC/ZIP package holding an
    indexed, welded, manifold mesh in 3D/3dmodel.model), which imports directly
    into BambuStudio for printing on a Bambu Lab printer.

    Multicolor support (optional):
        A Paint Layer (polygon, polyline, or point feature class) can be draped onto
        the terrain surface.  The tool rasterizes the layer onto the same grid as the
        resampled DEM, then assigns the paint color to every terrain triangle whose
        quad overlaps the paint features.  Lines and points are pre-buffered to give
        them area on the raster grid.  Wall and bottom triangles always use the base
        color.

        Colors are written using the 3MF Materials and Properties extension
        (xmlns:m="http://schemas.microsoft.com/3dmanufacturing/material/2015/02").
        BambuStudio 1.9+ imports them via its "Standard 3MF Color Parsing" dialog,
        which maps each color group to an AMS filament slot.

    See DEMToSTL for the full workflow, footprint-mode, and coordinate-system
    documentation.
    """

    def __init__(self):
        self.label = "DEM to 3MF"
        self.description = (
            "Converts a Digital Elevation Model (DEM) raster to a watertight, standards-"
            "compliant 3MF file suitable for slicing and 3-D printing (e.g. in BambuStudio). "
            "Supports optional per-triangle multicolor painting by draping a vector feature "
            "layer onto the terrain surface.  Colors are encoded using the 3MF Materials "
            "extension and imported by BambuStudio via its Standard 3MF Color Parsing dialog."
        )
        self.canRunInBackground = False

    # ------------------------------------------------------------------
    def getParameterInfo(self):

        # 0 - Input DEM
        p0 = arcpy.Parameter(
            displayName="Input DEM",
            name="in_dem",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input",
        )

        # 1 - Elevation Units
        p1 = arcpy.Parameter(
            displayName="Elevation Units",
            name="elev_units",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p1.filter.type = "ValueList"
        p1.filter.list = ["Meters", "Feet"]
        # No default - forces an explicit choice when auto-detection fails.

        # 2 - Output 3MF file
        p2 = arcpy.Parameter(
            displayName="Output 3MF File",
            name="out_3mf",
            datatype="DEFile",
            parameterType="Required",
            direction="Output",
        )
        p2.filter.list = ["3mf"]

        # 3 - Max print-bed dimension
        p3 = arcpy.Parameter(
            displayName="Maximum Print-Bed Dimension (mm)",
            name="max_bed_dim",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p3.value = 180.0

        # 4 - Vertical exaggeration
        p4 = arcpy.Parameter(
            displayName="Vertical Exaggeration Factor",
            name="vert_exag",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p4.value = 1.0

        # 5 - Minimum detail size (controls resampling / mesh density)
        p5 = arcpy.Parameter(
            displayName="Minimum Detail Size (mm)",
            name="min_detail",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p5.value = 0.2

        # 6 - Base thickness
        p6 = arcpy.Parameter(
            displayName="Base Thickness (mm)",
            name="base_thick",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p6.value = 3.0

        # 7 - Z floor (vertical reference for the base of the model)
        p7 = arcpy.Parameter(
            displayName="Z Floor Reference",
            name="z_floor",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p7.filter.type = "ValueList"
        p7.filter.list = ["Sea Level (0)", "Minimum Elevation"]
        p7.value = "Sea Level (0)"

        # 8 - Model footprint (rectangular vs. tight boundary)
        p8 = arcpy.Parameter(
            displayName="Model Footprint",
            name="model_footprint",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p8.filter.type = "ValueList"
        p8.filter.list = ["Tight (Follows DEM Boundary)", "Rectangular"]
        p8.value = "Tight (Follows DEM Boundary)"

        # ── Colors group ─────────────────────────────────────────────────────

        # 9 - Base Color (whole-model default; always written to the 3MF color group)
        p9 = arcpy.Parameter(
            displayName="Base Color",
            name="base_color",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Colors",
        )
        p9.filter.type = "ValueList"
        p9.filter.list = list(PAINT_COLORS.keys())
        p9.value = "White"

        # 10 - Paint Layers (value table: one row per paint layer)
        # Columns: Layer (feature class draped onto terrain), Color (named color),
        # Width mm (buffer radius for Line/Point layers; ignored for Polygon layers).
        p10 = arcpy.Parameter(
            displayName="Paint Layers",
            name="paint_layers",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input",
            category="Colors",
        )
        p10.columns = [
            ["GPFeatureLayer", "Layer"],
            ["GPString",       "Color"],
            ["GPDouble",       "Width (mm)"],
        ]
        p10.filters[1].type = "ValueList"
        p10.filters[1].list = list(PAINT_COLORS.keys())

        return [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]

    # ------------------------------------------------------------------
    def isLicensed(self):
        return True

    # ------------------------------------------------------------------
    def updateParameters(self, parameters):
        dem_param   = parameters[0]
        units_param = parameters[1]

        # Re-detect elevation units every time the DEM changes (see DEMToSTL).
        if dem_param.value and dem_param.altered and not dem_param.hasBeenValidated:
            units_param.value = _detect_elev_units(dem_param.valueAsText)

    # ------------------------------------------------------------------
    def updateMessages(self, parameters):
        dem_param   = parameters[0]
        units_param = parameters[1]
        p_bed    = parameters[3]
        p_exag   = parameters[4]
        p_detail = parameters[5]
        p_base   = parameters[6]

        if dem_param.value and not units_param.value:
            units_param.setWarningMessage(
                "Elevation units could not be auto-detected from this DEM's spatial "
                "reference (no Vertical Coordinate System defined). Please select."
            )

        if p_bed.value is not None:
            if p_bed.value <= 0:
                p_bed.setErrorMessage("Maximum bed dimension must be greater than 0.")
            elif p_bed.value < 10:
                p_bed.setWarningMessage("Very small bed dimension - the model may lack usable surface detail.")

        if p_exag.value is not None and p_exag.value <= 0:
            p_exag.setErrorMessage("Vertical exaggeration must be greater than 0.")

        if p_detail.value is not None and p_detail.value <= 0:
            p_detail.setErrorMessage("Minimum detail size must be greater than 0.")

        if p_base.value is not None and p_base.value < 0:
            p_base.setErrorMessage("Base thickness cannot be negative.")

    # ------------------------------------------------------------------
    def execute(self, parameters, messages):

        in_dem          = parameters[0].valueAsText
        elev_units      = parameters[1].valueAsText
        out_3mf         = parameters[2].valueAsText
        max_bed         = float(parameters[3].value)
        vert_exag       = float(parameters[4].value)
        min_detail      = float(parameters[5].value)
        base_thick      = float(parameters[6].value)
        z_floor_mode    = parameters[7].valueAsText
        footprint_mode  = parameters[8].valueAsText
        base_color_name = parameters[9].valueAsText   # e.g. "White"

        # Parse the Paint Layers value table.
        # Each row: [GPFeatureLayer, GPString color name, GPDouble width mm]
        paint_layer_entries = []   # list of (layer_path, color_name, color_hex, width_mm)
        vtab = parameters[10]
        if vtab.value:
            for row in vtab.value:
                layer_path = str(row[0]) if row[0] is not None else None
                if not layer_path:
                    continue
                color_name = str(row[1]) if row[1] is not None else list(PAINT_COLORS.keys())[0]
                color_hex  = PAINT_COLORS.get(color_name, "#FF0000")
                width_mm   = float(row[2]) if row[2] else None
                paint_layer_entries.append((layer_path, color_name, color_hex, width_mm))

        tight = (footprint_mode == "Tight (Follows DEM Boundary)")

        if not out_3mf.lower().endswith(".3mf"):
            out_3mf += ".3mf"

        base_color_hex = PAINT_COLORS.get(base_color_name or "White", "#FFFFFF")

        # ── Step 1/6 - Describe DEM ──────────────────────────────────────────
        messages.addMessage("Step 1/6 - Analyzing input DEM...")
        desc    = arcpy.Describe(in_dem)
        sr      = desc.spatialReference
        ext     = desc.extent
        orig_cs = desc.meanCellWidth

        rw_w = ext.width
        rw_h = ext.height
        messages.addMessage(f"  Extent    : {rw_w:.3f} x {rw_h:.3f} map units")
        messages.addMessage(f"  Cell size : {orig_cs:.4f} map units")

        # z_to_xy converts elevation values to XY map units before xy_scale is
        # applied - handles the common case where Z is in feet but XY is meters.
        xy_unit_name = sr.linearUnitName.lower() if sr else "meter"
        xy_is_feet   = "foot" in xy_unit_name or "feet" in xy_unit_name

        if elev_units is None:
            elev_units = "Feet" if xy_is_feet else "Meters"
            messages.addMessage(f"  Elevation units : not specified - assuming {elev_units} (same as XY units)")

        z_is_feet = (elev_units == "Feet")
        if z_is_feet and not xy_is_feet:
            z_to_xy = 0.3048
        elif not z_is_feet and xy_is_feet:
            z_to_xy = 1.0 / 0.3048
        else:
            z_to_xy = 1.0

        if z_to_xy != 1.0:
            messages.addMessage(
                f"  Elevation units : {elev_units} → XY units are "
                f"{'feet' if xy_is_feet else 'meters'}  "
                f"(z_to_xy factor: {z_to_xy:.7f})"
            )
        else:
            messages.addMessage(f"  Elevation units : {elev_units} (matches XY units, no conversion)")

        # ── Step 2/6 - Scale and target cell size ────────────────────────────
        messages.addMessage("Step 2/6 - Computing scale factors...")
        xy_scale  = max_bed / max(rw_w, rw_h)
        target_rw = max(min_detail / xy_scale, orig_cs)

        messages.addMessage(f"  XY scale    : 1:{1/xy_scale:,.0f}  ({xy_scale:.6f} mm/map unit)")
        messages.addMessage(f"  Target cell : {target_rw:.4f} map units  ({target_rw * xy_scale:.3f} mm in model)")

        # ── Step 3/6 - Resample DEM ──────────────────────────────────────────
        messages.addMessage("Step 3/6 - Resampling DEM to target cell size (BILINEAR)...")
        scratch    = arcpy.env.scratchFolder or tempfile.gettempdir()
        tmp_raster = os.path.join(scratch, "dem23mf_resampled.tif")

        quad_paint_masks  = []   # list of bool (nr-1, nc-1) arrays, one per layer
        paint_colors_hex  = []   # parallel list of hex strings
        paint_layer_names = []   # parallel list of display names (for final summary)
        try:
            arcpy.management.Resample(
                in_raster=in_dem,
                out_raster=tmp_raster,
                cell_size=target_rw,
                resampling_type="BILINEAR",
            )

            # ── Step 4/6 - Load into NumPy ───────────────────────────────────
            messages.addMessage("Step 4/6 - Loading raster into memory...")
            arr = _load_raster_array(tmp_raster)

            # Rasterize each paint layer while the snap raster still exists on disk.
            # _rasterize_paint_layer aligns its output to tmp_raster's grid exactly,
            # so the resulting boolean mask indexes 1-to-1 with arr.
            if paint_layer_entries:
                n = len(paint_layer_entries)
                messages.addMessage(f"  Rasterizing {n} paint layer(s)...")
                total_quads = (arr.shape[0] - 1) * (arr.shape[1] - 1)

                dem_paint_extent_fc = "memory/dem2stl_paint_extent"
                try:
                    if arcpy.Exists(dem_paint_extent_fc):
                        arcpy.management.Delete(dem_paint_extent_fc)
                    ring = arcpy.Array([
                        arcpy.Point(ext.XMin, ext.YMin),
                        arcpy.Point(ext.XMin, ext.YMax),
                        arcpy.Point(ext.XMax, ext.YMax),
                        arcpy.Point(ext.XMax, ext.YMin),
                    ])
                    arcpy.management.CopyFeatures([arcpy.Polygon(ring, sr)], dem_paint_extent_fc)
                except Exception:
                    dem_paint_extent_fc = None

                try:
                    for i, (layer_path, color_name, color_hex, width_mm) in enumerate(paint_layer_entries):
                        try:
                            shape_type = arcpy.Describe(layer_path).shapeType
                        except Exception:
                            shape_type = ""
                        buf_units    = None
                        display_width = None
                        if shape_type in ("Polyline", "Point", "Multipoint"):
                            effective_width = width_mm if width_mm is not None else 1.0
                            buf_units     = (effective_width / 2.0) / xy_scale
                            display_width = effective_width   # always shown for line/point
                        elif shape_type == "Polygon" and width_mm is not None and width_mm > 0:
                            buf_units     = (width_mm / 2.0) / xy_scale
                            display_width = width_mm          # only shown when explicitly set
                        layer_basename = os.path.basename(str(layer_path))
                        buf_info = f", {display_width:.4g} mm width" if display_width is not None else ""
                        messages.addMessage(
                            f"  Layer {i+1}/{n}: {layer_basename}"
                            f"  ({shape_type or 'Unknown'} → {color_name}{buf_info})"
                        )
                        paint_src       = layer_path
                        tmp_filtered_fc = None
                        if dem_paint_extent_fc:
                            try:
                                paint_fl = arcpy.management.MakeFeatureLayer(
                                    layer_path, f"dem2stl_paintfl_{i}"
                                )[0]
                                try:
                                    total_feat    = int(arcpy.management.GetCount(paint_fl)[0])
                                    arcpy.management.SelectLayerByLocation(
                                        paint_fl, "INTERSECT", dem_paint_extent_fc
                                    )
                                    selected_feat = int(arcpy.management.GetCount(paint_fl)[0])
                                    if selected_feat < total_feat:
                                        tmp_filtered_fc = f"memory/dem2stl_paintfilt_{i}"
                                        if arcpy.Exists(tmp_filtered_fc):
                                            arcpy.management.Delete(tmp_filtered_fc)
                                        arcpy.management.CopyFeatures(paint_fl, tmp_filtered_fc)
                                        paint_src = tmp_filtered_fc
                                        messages.addMessage(
                                            f"    Pre-filtered from {total_feat:,} to "
                                            f"{selected_feat:,} features intersecting DEM extent"
                                        )
                                finally:
                                    try:
                                        arcpy.management.Delete(paint_fl)
                                    except Exception:
                                        pass
                            except Exception:
                                if tmp_filtered_fc and arcpy.Exists(tmp_filtered_fc):
                                    try:
                                        arcpy.management.Delete(tmp_filtered_fc)
                                    except Exception:
                                        pass
                                tmp_filtered_fc = None
                                paint_src = layer_path
                        try:
                            paint_grid = _rasterize_paint_layer(
                                paint_src, tmp_raster, arr.shape, scratch, buf_units
                            )
                        finally:
                            if tmp_filtered_fc and arcpy.Exists(tmp_filtered_fc):
                                try:
                                    arcpy.management.Delete(tmp_filtered_fc)
                                except Exception:
                                    pass
                        messages.addMessage(f"    Paint cells : {int(paint_grid.sum()):,} of {paint_grid.size:,}")
                        quad_mask = (
                            paint_grid[:-1, :-1] | paint_grid[:-1, 1:] |
                            paint_grid[1:,  :-1] | paint_grid[1:,  1:]
                        )
                        painted_quads = int(quad_mask.sum())
                        messages.addMessage(
                            f"    Paint quads : {painted_quads:,} of {total_quads:,}"
                            f"  ({painted_quads * 2:,} terrain triangles)"
                        )
                        if painted_quads > 0:
                            quad_paint_masks.append(quad_mask)
                            paint_colors_hex.append(color_hex)
                            paint_layer_names.append(layer_basename)
                finally:
                    if dem_paint_extent_fc and arcpy.Exists(dem_paint_extent_fc):
                        try:
                            arcpy.management.Delete(dem_paint_extent_fc)
                        except Exception:
                            pass

        finally:
            if arcpy.Exists(tmp_raster):
                arcpy.management.Delete(tmp_raster)

        nr, nc = arr.shape
        messages.addMessage(f"  Grid : {nc} cols x {nr} rows  ({nc * nr:,} vertices)")

        if nr < 2 or nc < 2:
            raise RuntimeError("Resampled grid is too small (< 2x2). Increase the bed size or reduce the Minimum Detail Size.")

        model_w_mm = rw_w * xy_scale
        model_h_mm = rw_h * xy_scale
        cell_mm_w  = model_w_mm / (nc - 1)
        cell_mm_h  = model_h_mm / (nr - 1)
        messages.addMessage(f"  Cell spacing : {cell_mm_w:.3f} x {cell_mm_h:.3f} mm")

        # ── Step 5/6 - Z floor and elevation array ───────────────────────────
        messages.addMessage("Step 5/6 - Processing elevation values...")
        nan_mask     = np.isnan(arr)
        elev_min_raw = float(np.nanmin(arr))
        elev_max_raw = float(np.nanmax(arr))

        if z_floor_mode == "Sea Level (0)":
            z_ref = min(elev_min_raw, 0.0)
        else:
            z_ref = elev_min_raw

        if not tight:
            if nan_mask.any():
                arr[nan_mask] = z_ref

        # z_to_xy normalises elevation values into XY map units before xy_scale
        # converts map units → mm.  When Z and XY share the same unit, z_to_xy=1.
        z_mm = (arr - z_ref) * z_to_xy * xy_scale * vert_exag + base_thick

        relief_mm      = float(np.nanmax(z_mm)) - base_thick
        total_z_mm     = float(np.nanmax(z_mm))
        elev_range_raw = elev_max_raw - elev_min_raw
        messages.addMessage(f"  Elev range     : {elev_min_raw:.2f} – {elev_max_raw:.2f}  (Δ{elev_range_raw:.2f} {elev_units.lower()})")
        messages.addMessage(f"  Terrain relief : {relief_mm:.2f} mm  (x{vert_exag:.2f} exag)")
        messages.addMessage(f"  Total height   : {total_z_mm:.2f} mm  (terrain + {base_thick:.1f} mm base)")

        # ── Step 6/6 - Build mesh and write 3MF ──────────────────────────────
        messages.addMessage(f"Step 6/6 - Building mesh and writing 3MF: {out_3mf}")
        object_name = os.path.splitext(os.path.basename(out_3mf))[0]
        total_tris, vert_count, fsize_mb = _array_to_3mf(
            z_mm, cell_mm_w, cell_mm_h, out_3mf, object_name,
            quad_paint_masks or None,
            base_color_hex,
            paint_colors_hex or None,
        )

        messages.addMessage("")
        messages.addMessage("✓  3MF written successfully!")
        messages.addMessage(f"  Model dimensions : {model_w_mm:.1f} x {model_h_mm:.1f} x {total_z_mm:.2f} mm  (W x D x H)")
        messages.addMessage(f"  Vertices         : {vert_count:,}")
        messages.addMessage(f"  Triangles        : {total_tris:,}")
        if paint_layer_names:
            messages.addMessage(f"  Base color       : {base_color_name}")
            layer_parts = []
            for name, hex_val in zip(paint_layer_names, paint_colors_hex):
                cname = next((k for k, v in PAINT_COLORS.items() if v == hex_val), hex_val)
                layer_parts.append(f"{name} ({cname})")
            messages.addMessage(f"  Paint layers     : {', '.join(layer_parts)}")
        else:
            messages.addMessage(f"  Color            : {base_color_name}")
        messages.addMessage(f"  File size        : {fsize_mb:.2f} MB")
        messages.addMessage(f"  Output           : {out_3mf}")

    # ------------------------------------------------------------------
    def postExecute(self, parameters):
        return


# =============================================================================
class SplitDEMTo3MF(object):
    """
    Splits a DEM raster by a polygon feature class and converts each piece to a
    separate watertight 3MF file for 3-D printing.

    This tool is the 3MF counterpart of "Split DEM to STL".  A single consistent
    scale is chosen so the largest polygon piece fits within the specified print-bed
    size; all pieces are produced at that scale, preserving their relative proportions.
    Multicolor paint layers from "DEM to 3MF" are fully supported - paint geometries
    are buffered once (outside the piece loop) and then rasterized per piece against
    each piece's own resampled DEM, guaranteeing pixel-perfect alignment.

    Paint layer strategy
    --------------------
    - Buffering (Buffer, PolygonToLine) depends only on the paint FC geometry, not on
      any piece extent.  It runs once per paint layer before the piece loop.
    - Rasterization (FeatureToRaster) must snap to each piece's own resampled raster
      so that the resulting boolean mask aligns 1-to-1 with the piece's elevation array.
      It runs once per (piece, paint layer) pair, but each rasterization is cheap
      because it covers only the piece's smaller extent.
    """

    def __init__(self):
        self.label = "Split DEM to 3MF"
        self.description = (
            "Splits a DEM raster by a polygon feature class and converts each piece "
            "to a separate watertight 3MF file. A single consistent scale is chosen "
            "so the largest piece fits within the specified print-bed size; all pieces "
            "are produced at that same scale. Supports optional multicolor paint layers "
            "identical to those in the DEM to 3MF tool."
        )
        self.canRunInBackground = False

    # ------------------------------------------------------------------
    def getParameterInfo(self):

        # 0 - Input DEM
        p0 = arcpy.Parameter(
            displayName="Input DEM",
            name="in_dem",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input",
        )

        # 1 - Elevation Units
        p1 = arcpy.Parameter(
            displayName="Elevation Units",
            name="elev_units",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p1.filter.type = "ValueList"
        p1.filter.list = ["Meters", "Feet"]

        # 2 - Split Polygon Feature Class
        p2 = arcpy.Parameter(
            displayName="Split Polygon Feature Class",
            name="split_polys",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input",
        )
        p2.filter.list = ["Polygon"]

        # 3 - Name Attribute (optional field from p2 used in output filenames)
        p3 = arcpy.Parameter(
            displayName="Name Attribute",
            name="name_attr",
            datatype="Field",
            parameterType="Optional",
            direction="Input",
        )
        p3.parameterDependencies = [p2.name]
        p3.filter.list = ["Short", "Long", "BigInteger", "Double", "Single", "Text"]

        # 4 - Output Folder
        p4 = arcpy.Parameter(
            displayName="Output Folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input",
        )

        # 5 - Output Base Name
        p5 = arcpy.Parameter(
            displayName="Output Base Name",
            name="base_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )

        # 6 - Max print-bed dimension
        p6 = arcpy.Parameter(
            displayName="Maximum Print-Bed Dimension (mm)",
            name="max_bed_dim",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p6.value = 180.0

        # 7 - Vertical exaggeration
        p7 = arcpy.Parameter(
            displayName="Vertical Exaggeration Factor",
            name="vert_exag",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p7.value = 1.0

        # 8 - Minimum detail size
        p8 = arcpy.Parameter(
            displayName="Minimum Detail Size (mm)",
            name="min_detail",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p8.value = 0.2

        # 9 - Base thickness
        p9 = arcpy.Parameter(
            displayName="Base Thickness (mm)",
            name="base_thick",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        p9.value = 3.0

        # 10 - Z floor reference
        p10 = arcpy.Parameter(
            displayName="Z Floor Reference",
            name="z_floor",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p10.filter.type = "ValueList"
        p10.filter.list = ["Sea Level (0)", "Minimum Elevation (per piece)"]
        p10.value = "Sea Level (0)"

        # 11 - Model footprint
        p11 = arcpy.Parameter(
            displayName="Model Footprint",
            name="model_footprint",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        p11.filter.type = "ValueList"
        p11.filter.list = ["Tight (Follows Polygon Boundary)", "Rectangular"]
        p11.value = "Tight (Follows Polygon Boundary)"

        # ── Colors group ─────────────────────────────────────────────────────

        # 12 - Base Color
        p12 = arcpy.Parameter(
            displayName="Base Color",
            name="base_color",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Colors",
        )
        p12.filter.type = "ValueList"
        p12.filter.list = list(PAINT_COLORS.keys())
        p12.value = "White"

        # 13 - Paint Layers
        p13 = arcpy.Parameter(
            displayName="Paint Layers",
            name="paint_layers",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input",
            category="Colors",
        )
        p13.columns = [
            ["GPFeatureLayer", "Layer"],
            ["GPString",       "Color"],
            ["GPDouble",       "Width (mm)"],
        ]
        p13.filters[1].type = "ValueList"
        p13.filters[1].list = list(PAINT_COLORS.keys())

        return [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13]

    # ------------------------------------------------------------------
    def isLicensed(self):
        return True

    # ------------------------------------------------------------------
    def updateParameters(self, parameters):
        dem_param   = parameters[0]
        units_param = parameters[1]
        poly_param  = parameters[2]
        attr_param  = parameters[3]

        if dem_param.value and dem_param.altered and not dem_param.hasBeenValidated:
            units_param.value = _detect_elev_units(dem_param.valueAsText)

        if poly_param.altered and not poly_param.hasBeenValidated:
            attr_param.value = None

    # ------------------------------------------------------------------
    def updateMessages(self, parameters):
        dem_param   = parameters[0]
        units_param = parameters[1]
        p_bed    = parameters[6]
        p_exag   = parameters[7]
        p_detail = parameters[8]
        p_base   = parameters[9]

        if dem_param.value and not units_param.value:
            units_param.setWarningMessage(
                "Elevation units could not be auto-detected from this DEM's spatial "
                "reference (no Vertical Coordinate System defined). Please select."
            )

        if p_bed.value is not None:
            if p_bed.value <= 0:
                p_bed.setErrorMessage("Maximum bed dimension must be greater than 0.")
            elif p_bed.value < 10:
                p_bed.setWarningMessage("Very small bed dimension - the model may lack usable surface detail.")

        if p_exag.value is not None and p_exag.value <= 0:
            p_exag.setErrorMessage("Vertical exaggeration must be greater than 0.")

        if p_detail.value is not None and p_detail.value <= 0:
            p_detail.setErrorMessage("Minimum detail size must be greater than 0.")

        if p_base.value is not None and p_base.value < 0:
            p_base.setErrorMessage("Base thickness cannot be negative.")

    # ------------------------------------------------------------------
    def execute(self, parameters, messages):

        in_dem          = parameters[0].valueAsText
        elev_units      = parameters[1].valueAsText
        poly_fc         = parameters[2].valueAsText
        name_attr       = parameters[3].valueAsText   # None if not chosen
        out_folder      = parameters[4].valueAsText
        base_name       = parameters[5].valueAsText
        max_bed         = float(parameters[6].value)
        vert_exag       = float(parameters[7].value)
        min_detail      = float(parameters[8].value)
        base_thick      = float(parameters[9].value)
        z_floor_mode    = parameters[10].valueAsText
        footprint_mode  = parameters[11].valueAsText
        base_color_name = parameters[12].valueAsText

        # Parse the Paint Layers value table.
        paint_layer_entries = []   # list of (layer_path, color_name, color_hex, width_mm)
        vtab = parameters[13]
        if vtab.value:
            for row in vtab.value:
                layer_path = str(row[0]) if row[0] is not None else None
                if not layer_path:
                    continue
                color_name = str(row[1]) if row[1] is not None else list(PAINT_COLORS.keys())[0]
                color_hex  = PAINT_COLORS.get(color_name, "#FF0000")
                width_mm   = float(row[2]) if row[2] else None
                paint_layer_entries.append((layer_path, color_name, color_hex, width_mm))

        tight = footprint_mode.startswith("Tight")

        base_color_hex = PAINT_COLORS.get(base_color_name or "White", "#FFFFFF")

        os.makedirs(out_folder, exist_ok=True)

        # ── Step 1/3 - Describe DEM ──────────────────────────────────────────
        messages.addMessage("Step 1/3 - Analyzing input DEM...")
        dem_desc = arcpy.Describe(in_dem)
        sr       = dem_desc.spatialReference
        dem_ext  = dem_desc.extent
        orig_cs  = dem_desc.meanCellWidth

        messages.addMessage(f"  Extent    : {dem_ext.width:.3f} x {dem_ext.height:.3f} map units")
        messages.addMessage(f"  Cell size : {orig_cs:.4f} map units")

        xy_unit_name = sr.linearUnitName.lower() if sr else "meter"
        xy_is_feet   = "foot" in xy_unit_name or "feet" in xy_unit_name

        if elev_units is None:
            elev_units = "Feet" if xy_is_feet else "Meters"
            messages.addMessage(f"  Elevation units : not specified - assuming {elev_units} (same as XY units)")

        z_is_feet = (elev_units == "Feet")
        if z_is_feet and not xy_is_feet:
            z_to_xy = 0.3048
        elif not z_is_feet and xy_is_feet:
            z_to_xy = 1.0 / 0.3048
        else:
            z_to_xy = 1.0

        if z_to_xy != 1.0:
            messages.addMessage(
                f"  Elevation units : {elev_units} → XY units are "
                f"{'feet' if xy_is_feet else 'meters'}  "
                f"(z_to_xy factor: {z_to_xy:.7f})"
            )
        else:
            messages.addMessage(f"  Elevation units : {elev_units} (matches XY units, no conversion)")

        # ── Step 2/3 - Scan polygons, compute global scale ───────────────────
        messages.addMessage("Step 2/3 - Scanning split polygons for consistent scale...")

        cursor_fields = ["OID@", "SHAPE@"]
        if name_attr:
            cursor_fields.append(name_attr)

        oid_field  = arcpy.Describe(poly_fc).OIDFieldName
        scratch    = arcpy.env.scratchFolder or tempfile.gettempdir()
        poly_layer = arcpy.management.MakeFeatureLayer(poly_fc, "split_3mf_polys")[0]

        # Pre-filter to polygons that intersect the DEM bounding box.
        # For large continent-wide layers this avoids iterating thousands of
        # irrelevant features. SelectLayerByLocation handles CRS differences.
        total_poly_count = int(arcpy.management.GetCount(poly_layer)[0])
        dem_extent_fc    = "memory/dem2stl_extent_filter"
        try:
            if arcpy.Exists(dem_extent_fc):
                arcpy.management.Delete(dem_extent_fc)
            ring = arcpy.Array([
                arcpy.Point(dem_ext.XMin, dem_ext.YMin),
                arcpy.Point(dem_ext.XMin, dem_ext.YMax),
                arcpy.Point(dem_ext.XMax, dem_ext.YMax),
                arcpy.Point(dem_ext.XMax, dem_ext.YMin),
            ])
            arcpy.management.CopyFeatures([arcpy.Polygon(ring, sr)], dem_extent_fc)
            arcpy.management.SelectLayerByLocation(poly_layer, "INTERSECT", dem_extent_fc)
        except Exception:
            arcpy.management.SelectLayerByAttribute(poly_layer, "CLEAR_SELECTION")
        finally:
            if arcpy.Exists(dem_extent_fc):
                arcpy.management.Delete(dem_extent_fc)

        selected_poly_count = int(arcpy.management.GetCount(poly_layer)[0])
        if selected_poly_count < total_poly_count:
            messages.addMessage(
                f"  Pre-filtered from {total_poly_count:,} to {selected_poly_count:,} "
                f"polygon(s) intersecting DEM bounding box"
            )

        pieces = []
        with arcpy.da.SearchCursor(poly_layer, cursor_fields) as cur:
            for row in cur:
                oid  = row[0]
                geom = row[1]
                nval = str(row[2]) if name_attr else None
                # Project to DEM CRS for extent comparison and overlap detection.
                # Clip uses OID-based selection and reprojects automatically,
                # so only this stored geometry (used for extent math) needs projection.
                try:
                    if sr and sr.name != "Unknown":
                        geom = geom.projectAs(sr)
                except Exception:
                    pass
                pieces.append((oid, geom, nval))

        if not pieces:
            raise RuntimeError("No polygon features overlap the DEM extent.")

        global_max_dim = 0.0
        for _, geom, _ in pieces:
            pe     = geom.extent
            clip_w = min(pe.XMax, dem_ext.XMax) - max(pe.XMin, dem_ext.XMin)
            clip_h = min(pe.YMax, dem_ext.YMax) - max(pe.YMin, dem_ext.YMin)
            global_max_dim = max(global_max_dim, max(0.0, clip_w), max(0.0, clip_h))

        if global_max_dim == 0.0:
            raise RuntimeError("No polygon features overlap the DEM extent.")

        xy_scale  = max_bed / global_max_dim
        target_rw = max(min_detail / xy_scale, orig_cs)

        messages.addMessage(f"  {len(pieces)} polygon(s) to process")
        messages.addMessage(f"  Largest piece dimension : {global_max_dim:.1f} map units")
        messages.addMessage(f"  XY scale    : 1:{1/xy_scale:,.0f}  ({xy_scale:.6f} mm/map unit)")
        messages.addMessage(f"  Target cell : {target_rw:.4f} map units  ({target_rw * xy_scale:.3f} mm)")

        # ── Step 3/3 - Process each polygon piece ────────────────────────────
        messages.addMessage(f"Step 3/3 - Processing {len(pieces)} piece(s)...")

        skipped = 0

        # Pre-buffer all paint layer geometries once before the piece loop.
        # Buffering is independent of piece extent; rasterization snaps to each
        # piece's own DEM raster, so it must happen per piece.
        prepared_srcs = []   # list of (src_fc, tmp_fcs, color_hex, display_width, basename)
        if paint_layer_entries:
            n = len(paint_layer_entries)
            messages.addMessage(f"  Preparing {n} paint layer(s)...")

            dem_paint_extent_fc = "memory/dem2stl_paint_extent"
            try:
                if arcpy.Exists(dem_paint_extent_fc):
                    arcpy.management.Delete(dem_paint_extent_fc)
                ring = arcpy.Array([
                    arcpy.Point(dem_ext.XMin, dem_ext.YMin),
                    arcpy.Point(dem_ext.XMin, dem_ext.YMax),
                    arcpy.Point(dem_ext.XMax, dem_ext.YMax),
                    arcpy.Point(dem_ext.XMax, dem_ext.YMin),
                ])
                arcpy.management.CopyFeatures([arcpy.Polygon(ring, sr)], dem_paint_extent_fc)
            except Exception:
                dem_paint_extent_fc = None

            for i, (layer_path, color_name, color_hex, width_mm) in enumerate(paint_layer_entries):
                try:
                    shape_type = arcpy.Describe(layer_path).shapeType
                except Exception:
                    shape_type = ""
                buf_units     = None
                display_width = None
                if shape_type in ("Polyline", "Point", "Multipoint"):
                    effective_width = width_mm if width_mm is not None else 1.0
                    buf_units       = (effective_width / 2.0) / xy_scale
                    display_width   = effective_width
                elif shape_type == "Polygon" and width_mm is not None and width_mm > 0:
                    buf_units     = (width_mm / 2.0) / xy_scale
                    display_width = width_mm
                layer_basename = os.path.basename(str(layer_path))
                buf_info = f", {display_width:.4g} mm width" if display_width is not None else ""
                messages.addMessage(
                    f"  Layer {i+1}/{n}: {layer_basename}"
                    f"  ({shape_type or 'Unknown'} → {color_name}{buf_info})"
                )
                paint_src       = layer_path
                tmp_filtered_fc = None
                if dem_paint_extent_fc:
                    try:
                        paint_fl = arcpy.management.MakeFeatureLayer(
                            layer_path, f"dem2stl_paintfl_{i}"
                        )[0]
                        try:
                            total_feat    = int(arcpy.management.GetCount(paint_fl)[0])
                            arcpy.management.SelectLayerByLocation(
                                paint_fl, "INTERSECT", dem_paint_extent_fc
                            )
                            selected_feat = int(arcpy.management.GetCount(paint_fl)[0])
                            if selected_feat < total_feat:
                                tmp_filtered_fc = f"memory/dem2stl_paintfilt_{i}"
                                if arcpy.Exists(tmp_filtered_fc):
                                    arcpy.management.Delete(tmp_filtered_fc)
                                arcpy.management.CopyFeatures(paint_fl, tmp_filtered_fc)
                                paint_src = tmp_filtered_fc
                                messages.addMessage(
                                    f"    Pre-filtered from {total_feat:,} to "
                                    f"{selected_feat:,} features intersecting DEM extent"
                                )
                        finally:
                            try:
                                arcpy.management.Delete(paint_fl)
                            except Exception:
                                pass
                    except Exception:
                        if tmp_filtered_fc and arcpy.Exists(tmp_filtered_fc):
                            try:
                                arcpy.management.Delete(tmp_filtered_fc)
                            except Exception:
                                pass
                        tmp_filtered_fc = None
                        paint_src = layer_path
                src_fc, tmp_fcs = _prepare_paint_src(paint_src, scratch, buf_units, suffix=f"_{i}")
                if tmp_filtered_fc:
                    tmp_fcs.append(tmp_filtered_fc)
                prepared_srcs.append((src_fc, tmp_fcs, color_hex, display_width, layer_basename))
            if dem_paint_extent_fc and arcpy.Exists(dem_paint_extent_fc):
                try:
                    arcpy.management.Delete(dem_paint_extent_fc)
                except Exception:
                    pass

        try:
            for idx, (oid, geom, nval) in enumerate(pieces):
                label   = nval if nval is not None else str(idx + 1)
                safe    = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(label)).strip("_")
                out_3mf = os.path.join(out_folder, f"{base_name}_{safe}.3mf")

                messages.addMessage(f"\n  [{idx+1}/{len(pieces)}] {label}")

                pe      = geom.extent
                clip_x0 = max(pe.XMin, dem_ext.XMin)
                clip_y0 = max(pe.YMin, dem_ext.YMin)
                clip_x1 = min(pe.XMax, dem_ext.XMax)
                clip_y1 = min(pe.YMax, dem_ext.YMax)
                if clip_x1 <= clip_x0 or clip_y1 <= clip_y0:
                    messages.addWarningMessage("    Polygon does not overlap the DEM - skipping.")
                    skipped += 1
                    continue

                rect       = f"{clip_x0} {clip_y0} {clip_x1} {clip_y1}"
                tmp_clip   = os.path.join(scratch, f"dem23mf_clip_{idx}.tif")
                tmp_resamp = os.path.join(scratch, f"dem23mf_resamp_{idx}.tif")

                try:
                    arcpy.management.SelectLayerByAttribute(
                        poly_layer, "NEW_SELECTION", f"{oid_field} = {oid}"
                    )
                    arcpy.management.Clip(
                        in_raster=in_dem,
                        rectangle=rect,
                        out_raster=tmp_clip,
                        in_template_dataset=poly_layer,
                        clipping_geometry="ClippingGeometry",
                    )
                    arcpy.management.Resample(
                        in_raster=tmp_clip,
                        out_raster=tmp_resamp,
                        cell_size=target_rw,
                        resampling_type="BILINEAR",
                    )

                    arr    = _load_raster_array(tmp_resamp)
                    nr, nc = arr.shape

                    if nr < 2 or nc < 2:
                        messages.addWarningMessage("    Piece is too small after resampling (<2x2 grid) - skipping.")
                        skipped += 1
                        continue

                    clip_ext  = arcpy.Describe(tmp_resamp).extent
                    rw_w      = clip_ext.width
                    rw_h      = clip_ext.height
                    cell_mm_w = rw_w * xy_scale / (nc - 1)
                    cell_mm_h = rw_h * xy_scale / (nr - 1)

                    nan_mask     = np.isnan(arr)
                    elev_min_raw = float(np.nanmin(arr))
                    z_ref        = min(elev_min_raw, 0.0) if z_floor_mode == "Sea Level (0)" else elev_min_raw

                    if not tight:
                        arr[nan_mask] = z_ref

                    z_mm       = (arr - z_ref) * z_to_xy * xy_scale * vert_exag + base_thick
                    total_z_mm = float(np.nanmax(z_mm))
                    model_w_mm = rw_w * xy_scale
                    model_h_mm = rw_h * xy_scale

                    # Rasterize each pre-buffered paint layer against this piece's snap raster.
                    quad_paint_masks  = []
                    paint_colors_hex  = []
                    paint_layer_names = []
                    for src_fc, _, color_hex, _, layer_basename in prepared_srcs:
                        paint_grid = _rasterize_src(src_fc, tmp_resamp, arr.shape, scratch)
                        quad_mask  = (
                            paint_grid[:-1, :-1] | paint_grid[:-1, 1:] |
                            paint_grid[1:,  :-1] | paint_grid[1:,  1:]
                        )
                        if quad_mask.any():
                            quad_paint_masks.append(quad_mask)
                            paint_colors_hex.append(color_hex)
                            paint_layer_names.append(layer_basename)

                    object_name = f"{base_name}_{safe}"
                    total_tris, vert_count, fsize_mb = _array_to_3mf(
                        z_mm, cell_mm_w, cell_mm_h, out_3mf, object_name,
                        quad_paint_masks or None,
                        base_color_hex,
                        paint_colors_hex or None,
                    )

                    color_info = ""
                    if paint_layer_names:
                        color_info = f",  {len(paint_layer_names)} paint layer(s)"
                    messages.addMessage(
                        f"    ✓  {os.path.basename(out_3mf)}"
                        f"  ({model_w_mm:.1f}x{model_h_mm:.1f}x{total_z_mm:.2f} mm,"
                        f"  {total_tris:,} tris,  {fsize_mb:.2f} MB{color_info})"
                    )

                except Exception as piece_err:
                    messages.addWarningMessage(f"    Failed: {piece_err}")
                    skipped += 1

                finally:
                    for tmp in (tmp_clip, tmp_resamp):
                        if arcpy.Exists(tmp):
                            arcpy.management.Delete(tmp)

        finally:
            arcpy.management.Delete(poly_layer)
            for _, tmp_fcs, _, _, _ in prepared_srcs:
                for fc in tmp_fcs:
                    try:
                        arcpy.management.Delete(fc)
                    except Exception:
                        pass

        processed = len(pieces) - skipped
        messages.addMessage(f"\n✓  {processed} of {len(pieces)} piece(s) written to: {out_folder}")
        if skipped:
            messages.addMessage(f"  ({skipped} skipped - see warnings above)")

    # ------------------------------------------------------------------
    def postExecute(self, parameters):
        return
