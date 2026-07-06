# DEM to STL Toolbox

An ArcGIS Pro Python Toolbox (`.pyt`) that converts Digital Elevation Model rasters to binary STL files for 3D printing.

**Tools**
- [DEM to STL](#dem-to-stl) — converts a single DEM to one STL file
- [Split DEM to STL](#split-dem-to-stl) — divides a DEM by a polygon feature class and produces a separate STL file for each piece at a consistent scale

**Requirements:** ArcGIS Pro 3.x with a Standard or Advanced license. No Spatial Analyst extension required.

---

## General Notes

### Elevation Units and Vertical Exaggeration

Many DEMs store elevation in feet while their horizontal coordinate system uses meters (e.g., UTM with NAVD 88 feet). If the tool treated the Z values as meters it would introduce an unintended ~3.28x vertical exaggeration. The **Elevation Units** parameter corrects this by applying the appropriate conversion factor before scaling.

When a DEM's spatial reference includes a Vertical Coordinate System (VCS), the elevation unit is auto-detected on input. Most real-world DEMs do not carry a formal VCS; in those cases the field is left blank and must be set manually. The default when no value is to assume Z units are the same as X and Y units.

### Output Format

Both tools write standard binary STL files compatible with any 3D printing slicer. The mesh is fully watertight: terrain surface, perimeter walls, and a flat bottom are all closed with consistent outward-facing normals.

### Mesh Density

The DEM is resampled to a target resolution *before* triangulation. This means triangle count is proportional to the chosen **Minimum Detail Size**, not the raw DEM resolution, so even high-resolution DEMs produce a manageable file at print scale.

---

## DEM to STL

Converts a single DEM raster to one STL file. The model is scaled so its longest horizontal axis fills the specified print-bed dimension.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Input DEM** | — | Single-band elevation raster. Both floating-point and integer pixel types are supported. |
| **Elevation Units** | *(auto-detect)* | Vertical unit of the elevation values: Meters or Feet. Auto-detected from the DEM's Vertical Coordinate System when available; otherwise left blank and must be set manually. |
| **Output STL File** | — | Output path for the binary STL file. The `.stl` extension is added if omitted. |
| **Maximum Print-Bed Dimension (mm)** | 180 | Longest dimension of the printer's build plate in millimeters. The model is uniformly scaled so its longest horizontal axis equals this value. |
| **Vertical Exaggeration Factor** | 1.0 | Multiplier applied to terrain relief after unit conversion. Use 1.0 for a geometrically accurate model. Higher values amplify relief, which is useful for flat terrain. Does not affect base thickness. |
| **Minimum Detail Size (mm)** | 0.2 | Smallest surface feature to resolve, in millimeters. Controls the resampling resolution before triangulation. Smaller values produce more triangles and larger files. The cell size is never made finer than the DEM's native resolution. |
| **Base Thickness (mm)** | 3.0 | Thickness of the solid flat base beneath the lowest terrain point. A minimum of 1–2 mm is recommended to prevent warping during printing. |
| **Z Floor Reference** | Sea Level (0) | Which elevation maps to the bottom of the terrain surface (the top of the base layer). **Sea Level (0)** maps elevation 0 to the terrain floor — useful when the datum relationship matters or the area includes below-sea-level terrain. **Minimum Elevation** maps the lowest cell in the DEM to the floor, maximizing visible relief. |
| **Model Footprint** | Tight (Follows DEM Boundary) | How non-rectangular DEM boundaries are handled. **Tight** leaves NoData cells as void and traces walls along the actual data boundary. **Rectangular** fills NoData cells with the floor elevation, producing a flat-based rectangular block. |

---

## Split DEM to STL

Splits a DEM by a polygon feature class and produces one STL file per polygon. Scale is computed globally: the largest polygon piece is measured first, a single XY scale is derived so that piece fits within the print-bed dimension, and all other pieces are produced at the same scale. This allows pieces to be physically arranged together at the correct relative size.

Polygons that do not overlap the DEM extent are skipped with a warning. Processing continues if a single piece fails.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Input DEM** | — | Single-band elevation raster. Both floating-point and integer pixel types are supported. |
| **Elevation Units** | *(auto-detect)* | Vertical unit of the elevation values: Meters or Feet. Applied consistently to all output pieces. See the [General Notes](#elevation-units-and-vertical-exaggeration) above. |
| **Split Polygon Feature Class** | — | Polygon feature class whose features define the DEM pieces to produce. One STL file is generated per feature that overlaps the DEM. The feature class may be in a different coordinate system from the DEM. |
| **Name Attribute** | *(optional)* | A field from the polygon feature class whose values are appended to the Output Base Name to form each filename (e.g., `Region_ValleyFloor.stl`). Characters invalid in filenames are replaced with underscores. If omitted, pieces are numbered sequentially (`Region_1.stl`, `Region_2.stl`, …). |
| **Output Folder** | — | Folder where all STL files are written. Created automatically if it does not exist. |
| **Output Base Name** | — | Filename prefix for all output files. Each file is named `<Base Name>_<piece name or number>.stl`. |
| **Maximum Print-Bed Dimension (mm)** | 180 | Longest build-plate dimension in millimeters. The largest polygon piece is scaled to fit this dimension; all other pieces use the same scale and print proportionally smaller. |
| **Vertical Exaggeration Factor** | 1.0 | Relief multiplier applied consistently to every piece. |
| **Minimum Detail Size (mm)** | 0.2 | Smallest surface feature to resolve, in millimeters. A single target cell size is derived from this value and applied to all pieces, ensuring consistent mesh density across the set. |
| **Base Thickness (mm)** | 3.0 | Solid base thickness beneath each piece's terrain floor, in millimeters. Applied to every piece. |
| **Z Floor Reference** | Sea Level (0) | Controls the terrain floor elevation for each piece. **Sea Level (0)** uses elevation 0 as the floor for every piece — because the datum reference is shared, adjacent pieces will have consistent relative heights and can be arranged side-by-side to form a continuous landscape. **Minimum Elevation (per piece)** uses each piece's own local minimum as its floor; pieces stand independently but may differ in absolute model height. |
| **Model Footprint** | Tight (Follows Polygon Boundary) | How the clipped boundary of each piece is handled. **Tight** traces walls along the actual polygon boundary and is recommended for split mode. **Rectangular** fills areas outside the polygon with the floor elevation, producing a rectangular block aligned to the polygon's bounding box. |
