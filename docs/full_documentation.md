# ncartifactgenerator Full Documentation

## Table of Contents

1. [Overview](#1-overview)
2. [Serverless Access Architecture](#2-serverless-access-architecture)
3. [Repository Structure](#3-repository-structure)
4. [Artefact Types and File Layout](#4-artefact-types-and-file-layout)
5. [Processing Pipeline](#5-processing-pipeline)
   - 5.1 [Step 0: Coordinate Ordering](#51-step-0-coordinate-ordering)
   - 5.2 [Step 1: Time-Series Rechunking (`-t.nc`)](#52-step-1-time-series-rechunking--tnc)
   - 5.3 [Step 2: Time-Series Binary Index (`-t.bin`)](#53-step-2-time-series-binary-index--tbin)
   - 5.4 [Step 3: Map Rechunking (`-xy.nc`)](#54-step-3-map-rechunking--xync)
   - 5.5 [Step 4: Map Binary Index (`-xy.bin`)](#55-step-4-map-binary-index--xybin)
   - 5.6 [Step 5: Metadata Accumulation](#56-step-5-metadata-accumulation)
6. [times.json Schema](#6-timesjson-schema)
7. [Function Reference](#7-function-reference)
   - 7.1 [generate_artifacts()](#71-generate_artifacts)
   - 7.2 [writeJson()](#72-writejson)
   - 7.3 [config_web()](#73-config_web)
   - 7.4 [write_nc_chunk_t()](#74-write_nc_chunk_t)
   - 7.5 [write_nc_chunk_xy()](#75-write_nc_chunk_xy)
   - 7.6 [write_nc_t_chunk_dir_iter()](#76-write_nc_t_chunk_dir_iter)
   - 7.7 [write_nc_xy_chunk_dir_iter()](#77-write_nc_xy_chunk_dir_iter)
   - 7.8 [fusion_pen_can()](#78-fusion_pen_can)
   - 7.9 [Utility Functions](#79-utility-functions)
8. [Input NetCDF Requirements](#8-input-netcdf-requirements)
9. [Usage Examples](#9-usage-examples)


---

## 1. Overview

`ncartifactgenerator` is an R package that converts 3D gridded NetCDF files (lon × lat × time) into a set of optimised files that can be served to web visualisation clients without a server-side GeoServer, WCS, or similar middleware. The web viewer accesses data using plain HTTP range requests against static files hosted on any web server or object store.

The package is used as the final step in the PTI+ Clima pipeline, applied to the NetCDF output of the [`indexCalc`](https://github.com/PTI-Clima/indexCalc) and [`data_flow`](https://github.com/PTI-Clima/data_flow) repositories. It is invoked by the `viewer_*_web_data` Airflow DAGs in the [`lcsc-dags`](https://github.com/PTI-Clima/lcsc-dags) orchestration repository.


---

## 2. Serverless Access Architecture

A standard NetCDF file (in HDF5 format) stores its data in chunks. The chunk layout determines which access patterns are fast. Typical climate NetCDF files use chunk dimensions that match the full spatial extent of a single time step (X × Y × 1) or some other layout suited to the producing application — but not necessarily to the web viewer's needs.

The PTI+ Clima web viewers support two access patterns:

1. **Pixel time-series**: the user clicks a grid cell and requests the full time series for that location. This requires reading a slice `[x, y, :]` across all time steps.
2. **Date map**: the user selects a date and requests the spatial field for that instant. This requires reading a slice `[:, :, t]` for a single time step.

These two patterns are in direct conflict in terms of optimal chunk layout. The package solves this by maintaining **two rechunked copies** of each NetCDF file:

| File suffix | HDF5 chunk size | Fast for |
|---|---|---|
| `-t.nc` | `[1, 1, T]` | Time-series queries (one pixel, all times) |
| `-xy.nc` | `[X, Y, 1]` | Map queries (all pixels, one time) |

Since HDF5 chunk byte addresses are not contiguous and are not directly addressable from a plain HTTP client, the package also generates a **binary index file** for each rechunked NetCDF. This index maps logical positions (pixel index or time index) to the byte offset and byte size of the corresponding HDF5 chunk within the file. The web viewer reads the index first (a tiny HTTP request), then fetches exactly the right bytes with a single `Range` header.

```
Browser                       Static file server
  │                                │
  │── GET times.json ─────────────►│ (metadata: variables, dates, extents)
  │◄── times.json ────────────────│
  │                                │
  │── GET ffd_pen-t.bin ──────────►│ (12 × lon_num × lat_num bytes)
  │◄── full index file ───────────│
  │                                │
  │── GET ffd_pen-t.nc             │
  │   Range: bytes=<offset>-<offset+size-1> ►│
  │◄── single HDF5 chunk ─────────│ (one pixel's time series)
```

This approach eliminates server-side processing entirely: no tiling service, no dynamic subsetting, no session state. The only server requirement is support for the HTTP `Range` header, which all modern web servers and object stores (including MinIO) provide.


---

## 3. Repository Structure

```
ncartifactgenerator/
├── R/
│   ├── functions.R      ← Utility functions + generate_artifacts() top-level pipeline
│   ├── chunk.R          ← Rechunking (write_nc_chunk_t/xy) and index generation
│   └── config-web.R     ← Metadata accumulation (config_web) and JSON output (writeJson)
├── man/                 ← Roxygen2-generated .Rd documentation files
├── DESCRIPTION          ← Package metadata, dependencies
├── NAMESPACE            ← Exported symbols
└── LICENSE.md           ← GPL-3 licence
```

Source files:

| File | Key exports |
|---|---|
| `R/functions.R` | `generate_artifacts`, `read_epsg`, `returnXYNames`, `read_times`, `read_coords`, `readMinMax`, `mergeArrays`, `timePosition`, `getVarName`, `arrayRtojs`, `get_struct_typecode` |
| `R/chunk.R` | `write_nc_chunk_t`, `write_nc_chunk_xy`, `write_nc_t_chunk_dir`, `write_nc_t_chunk_dir_iter`, `write_nc_xy_chunk_dir`, `write_nc_xy_chunk_dir_iter`, `fusion_pen_can`, `create_nc_name` |
| `R/config-web.R` | `config_web`, `writeJson` |


---

## 4. Artefact Types and File Layout

For a single call to `generate_artifacts()` with `var_id = "ffd"` and `portion = "_pen"`, the following files are created under `out_root`:

```
<out_root>/
├── nc/
│   ├── ffd_pen-t.nc    ← rechunked NC: chunk = [1, 1, T]; WGS84 CRS variable added
│   ├── ffd_pen-t.bin   ← pixel index: lon_num × lat_num entries, each 12 bytes
│   ├── ffd_pen-xy.nc   ← rechunked NC: chunk = [X, Y, 1]; WGS84 CRS variable added
│   └── ffd_pen-xy.bin  ← date index: T entries, each 12 bytes
└── times.json          ← written once after all generate_artifacts() calls complete
```

**Binary index entry format** (both `-t.bin` and `-xy.bin`):

| Bytes | Type | Description |
|---|---|---|
| 0–7 | `int64`, little-endian | Byte offset of the HDF5 chunk in the `.nc` file |
| 8–11 | `int32`, little-endian | Byte size of the HDF5 chunk |

For a `-t.bin` file, entries are ordered row-major: entry at position `lat_index × lon_num + lon_index` corresponds to the pixel at `(lon_index, lat_index)`. For a `-xy.bin` file, entry at position `t_index` corresponds to time step `t_index`.

All rechunked NetCDF files use:
- NetCDF-4 format (`force_v4 = TRUE`)
- Compression level 9 (zlib)
- `NaN` as fill value
- A `crs` scalar variable with full WGS84 (EPSG:4326) `crs_wkt` attribute


---

## 5. Processing Pipeline

`generate_artifacts()` executes the following five steps in sequence for each `(nc_filename, portion)` pair. Steps 1–4 are skipped if the output file already exists and has a modification time newer than the source NetCDF.

### 5.1 Step 0: Coordinate Ordering

Opens the source NetCDF file in write mode and checks whether the longitude and latitude dimension arrays are in ascending order.

- If `lon[1] > lon[end]`, the data array is reversed along the longitude axis (in batches of 100 columns to avoid memory overflow) and the `lon` coordinate vector is reversed.
- If `lat[1] > lat[end]`, the same is done along the latitude axis.

This normalisation is required because the rechunking and index-building steps assume ascending coordinates.

### 5.2 Step 1: Time-Series Rechunking (`-t.nc`)

Calls `write_nc_chunk_t()` to create a new NetCDF-4 file with:
- Chunk dimensions `[1, 1, T]` — one pixel's full time series fits in a single HDF5 chunk.
- All global attributes, dimension attributes, and variable attributes copied from the source.
- A `crs` scalar variable with WGS84 attributes added (so the file is CF-conformant).

Data is read and written in spatial batches of `lon_by × lat_by` (default 100 × 100) to limit peak memory use.

### 5.3 Step 2: Time-Series Binary Index (`-t.bin`)

Calls `write_nc_t_chunk_dir_iter()`, which:
1. Opens the `-t.nc` file with `hdf5r` to access HDF5 chunk metadata directly.
2. Iterates over all chunks using `$chunk_iter()`, retrieving each chunk's byte offset and size.
3. For each chunk at HDF5 offset `(0, lat_index, lon_index)`, writes the 8-byte address and 4-byte size at position `(OFFSET_TYPE_SIZE + SIZE_TYPE_SIZE) × (lon_index + lat_index × lon_num)` in the output binary file.

The result is a flat binary array of 12-byte records, indexed by `lat × lon` position, that the browser can use to locate any pixel's time series with a single seek.

A legacy non-iterator variant (`write_nc_t_chunk_dir`) uses `get_chunk_info_by_coord()` instead and may be slower on large files.

### 5.4 Step 3: Map Rechunking (`-xy.nc`)

Calls `write_nc_chunk_xy()` to create a new NetCDF-4 file with:
- Chunk dimensions `[X, Y, 1]` — one complete spatial map fits in a single HDF5 chunk.
- Same attribute copying and CRS variable as `-t.nc`.

Data is read and written in time batches of `time_by` (default 100) to limit peak memory.

### 5.5 Step 4: Map Binary Index (`-xy.bin`)

Calls `write_nc_xy_chunk_dir_iter()`, which:
1. Opens the `-xy.nc` file with `hdf5r`.
2. Iterates over all chunks; each chunk corresponds to one time step.
3. For each chunk at offset `(time_index, 0, 0)`, writes address + size at byte position `(OFFSET_TYPE_SIZE + SIZE_TYPE_SIZE) × time_index`.

### 5.6 Step 5: Metadata Accumulation

Calls `config_web()` with the source NetCDF file and the accumulating `info_js` object. `config_web()`:
- Reads spatial extents (lon/lat min, max, count).
- Reads the time coordinate and converts it to `Date` objects.
- Optionally calculates per-time-step min/max values (`readMinMax()`), or accepts pre-computed values (`varmin=-1, varmax=-1` to skip).
- Merges per-portion extents and per-time-step value ranges into the shared `info_js` accumulator using `mergeArrays()`.

The `info_js` object is returned and passed to the next `generate_artifacts()` call, building up a complete picture of all variables and both spatial portions.


---

## 6. `times.json` Schema

After all variables are processed, `writeJson()` serialises `info_js` to `<out_root>/times.json`. The JSON structure is:

```json
{
  "center": { "lat": 40.0, "lng": -3.5 },
  "times": {
    "ffd": ["1961-01-01", "1962-01-01", ...],
    "lfd": [...]
  },
  "varMin": {
    "ffd": [15.0, 12.0, ...],
    ...
  },
  "varMax": {
    "ffd": [320.0, 318.0, ...],
    ...
  },
  "varTitle": {
    "ffd": ["First frost day"],
    "lfd": ["Last frost day"]
  },
  "legendTitle": {
    "ffd": ["days"],
    "lfd": ["days"]
  },
  "portions": {
    "ffd": ["_pen", "_can"],
    "lfd": ["_pen", "_can"]
  },
  "lonMin": { "ffd_pen": -9.5, "ffd_can": -18.2, ... },
  "lonMax": { "ffd_pen": 4.3, "ffd_can": -13.4, ... },
  "lonNum": { "ffd_pen": 545, "ffd_can": 189, ... },
  "latMin": { "ffd_pen": 35.8, ... },
  "latMax": { "ffd_pen": 43.8, ... },
  "latNum": { "ffd_pen": 341, ... },
  "minVal": { "ffd": 0.0, ... },
  "maxVal": { "ffd": 365.0, ... },
  "varType": "f",
  "offsetType": "Q",
  "sizeType": "I",
  "projection": "EPSG:4326"
}
```

| Field | Description |
|---|---|
| `center` | Geographic centroid of all processed domains, in EPSG:4326 |
| `times.<var>` | Array of date strings (one per time step) for each variable |
| `varMin.<var>` / `varMax.<var>` | Per-time-step min/max value arrays (after merging across portions) |
| `varTitle.<var>` | Human-readable display title for each variable |
| `legendTitle.<var>` | Legend unit string for each variable |
| `portions.<var>` | Array of portion suffixes available for each variable (e.g., `["_pen", "_can"]`) |
| `lonMin/lonMax/lonNum.<var><portion>` | Spatial extent and resolution per variable × portion |
| `latMin/latMax/latNum.<var><portion>` | Spatial extent and resolution per variable × portion |
| `minVal/maxVal.<var>` | Global min/max over all time steps and portions |
| `varType` | Struct format code for the data type (`"f"` = float32, `"d"` = float64, `"i"` = int32) |
| `offsetType` | Struct format code for chunk offsets in `.bin` files (`"Q"` = uint64) |
| `sizeType` | Struct format code for chunk sizes in `.bin` files (`"I"` = uint32) |
| `projection` | CRS identifier string |


---

## 7. Function Reference

### 7.1 `generate_artifacts()`

```r
generate_artifacts(
  nc_root, out_root, nc_filename, portion, var_id,
  epsg, info_js,
  lon_name = NA, lat_name = NA, time_name = NA,
  write = FALSE, calcMaxMin = FALSE, nc_dims = 3
)
```

Main entry point. Runs the full 5-step pipeline for one `(nc_filename, portion)` pair.

**Parameters**:

| Parameter | Description |
|---|---|
| `nc_root` | Directory containing the source NetCDF file |
| `out_root` | Output directory; an `nc/` subdirectory is created inside it |
| `nc_filename` | Path relative to `nc_root`, without the `portion` suffix and `.nc` extension |
| `portion` | Spatial domain suffix appended to filenames, e.g. `"_pen"` or `"_can"` |
| `var_id` | Output variable identifier; used as the base name for output files |
| `epsg` | EPSG code of the input CRS (e.g., `"4326"`) |
| `info_js` | Accumulated metadata object from previous calls; pass `NA` on first call |
| `lon_name`, `lat_name` | Names of the longitude and latitude dimensions; auto-detected if `NA` |
| `time_name` | Name of the time dimension; auto-detected if `NA` |
| `write` | Reserved parameter; currently unused |
| `calcMaxMin` | If `TRUE`, computes per-time-step min/max from data; if `FALSE`, sets them to `-1` (skipped) |
| `nc_dims` | Number of dimensions: `3` for standard lon×lat×time (default); `2` for 2D variables without time |

**Returns**: Updated `info_js` list with this variable's metadata merged in.

**Idempotency**: Each of the four output files (steps 1–4) is regenerated only if the source NetCDF has a more recent modification time. Re-running is safe and efficient.

### 7.2 `writeJson()`

```r
writeJson(folder, infoJs, varTitle, legendTitle = "Legend",
          offsetType = "Q", sizeType = "I", minify = TRUE)
```

Serialises the accumulated `infoJs` metadata to `<folder>/times.json`.

**Parameters**:

| Parameter | Description |
|---|---|
| `folder` | Output directory (same as `out_root` passed to `generate_artifacts`) |
| `infoJs` | Accumulated metadata object returned by the last `generate_artifacts()` call |
| `varTitle` | Named list mapping `var_id` → display title string |
| `legendTitle` | Named list mapping `var_id` → legend unit string; or a single string if shared |
| `offsetType` | Struct format code for chunk offsets; `"Q"` = unsigned 64-bit integer |
| `sizeType` | Struct format code for chunk sizes; `"I"` = unsigned 32-bit integer |
| `minify` | If `TRUE` (default), outputs compact single-line JSON |

**Returns**: The JSON string (also written to file).

### 7.3 `config_web()`

```r
config_web(file, folder, epsg, formatdates, varmin, varmax,
           varName, infoJs = NA, lon_name, lat_name,
           time_name = "time", portion)
```

Accumulates metadata for one variable × portion into `infoJs`. Called internally by `generate_artifacts()`.

For each call, it:
- Reads lon/lat extents and dimension counts.
- Reads the time coordinate and converts to `Date` (or formatted strings if `formatdates` is supplied).
- Computes or receives per-time-step min/max values.
- Merges value ranges with any previously accumulated ranges for the same `varName` using `mergeArrays()` (taking the minimum of per-step minima and maximum of per-step maxima across portions).

**Returns**: Updated `infoJs` list.

### 7.4 `write_nc_chunk_t()`

```r
write_nc_chunk_t(in_file, out_file, lon_by = -1, lat_by = -1,
                 lon_name, lat_name, time_name, signif_digits)
```

Creates a rechunked copy of `in_file` optimised for pixel time-series access.

- Output chunk size: `[1, 1, time_num]` — one pixel's full time series.
- Processes data in `lon_by × lat_by` spatial tiles (default: all at once if `-1`).
- Adds a `crs` scalar variable with WGS84 grid mapping attributes.
- Copies all global and variable attributes from the source.
- Optional `signif_digits` parameter rounds data values to reduce file size.

### 7.5 `write_nc_chunk_xy()`

```r
write_nc_chunk_xy(in_file, out_file, time_by = -1,
                  lon_name, lat_name, time_name, signif_digits)
```

Creates a rechunked copy of `in_file` optimised for date-map access.

- Output chunk size: `[lon_num, lat_num, 1]` — one full spatial map per chunk.
- Processes data in `time_by` temporal batches.
- Same CRS variable and attribute copying as `write_nc_chunk_t()`.

### 7.6 `write_nc_t_chunk_dir_iter()`

```r
write_nc_t_chunk_dir_iter(in_file, out_file)
```

Creates the binary pixel index for a `-t.nc` file using `hdf5r`'s `$chunk_iter()` API.

Each record in the output binary file is 12 bytes: 8-byte int64 chunk address + 4-byte int32 chunk size. Records are ordered `lat × lon` (row-major), so the record for pixel `(lon_i, lat_j)` is at byte offset `12 × (lon_i + lat_j × lon_num)`.

A legacy variant `write_nc_t_chunk_dir()` uses `$get_chunk_info_by_coord()` and iterates over pixels explicitly; it produces identical output but may be slower on large files.

### 7.7 `write_nc_xy_chunk_dir_iter()`

```r
write_nc_xy_chunk_dir_iter(in_file, out_file)
```

Creates the binary date index for a `-xy.nc` file.

Each record is 12 bytes: 8-byte int64 address + 4-byte int32 size. Record at position `12 × t_index` corresponds to the spatial map for time step `t_index`.

A legacy variant `write_nc_xy_chunk_dir()` iterates over time steps explicitly.

### 7.8 `fusion_pen_can()`

```r
fusion_pen_can(can_filename, pen_filename, fusion_filename,
               nc_chunk_num = 16, signif_digits,
               lon_name = "lon", lat_name = "lat",
               time_name = "time", grid_size)
```

Merges two spatially disjoint NetCDF files — one covering the Iberian Peninsula + Balearics (`pen`) and one covering the Canary Islands (`can`) — into a single NetCDF covering the full geographic bounding box of both domains. Pixels in the bounding box not covered by either domain are left as `NA`.

This is useful when a web viewer needs to display a single combined map of all Spanish territory. The merged file shares the time dimension from the peninsular file; both input files must have the same time coordinate.

**Parameters**:

| Parameter | Description |
|---|---|
| `can_filename` | Path to the Canary Islands NetCDF |
| `pen_filename` | Path to the Peninsular NetCDF |
| `fusion_filename` | Path for the output merged NetCDF |
| `nc_chunk_num` | Number of chunks per dimension in the output (default 16; adjusts automatically if smaller than `time_num`) |
| `signif_digits` | Optional: round values to N significant digits |
| `grid_size` | Grid resolution in degrees; auto-detected from the peninsular file if not supplied |

The output uses CF-1.8 conventions, WGS84 CRS, and compression level 9.

### 7.9 Utility Functions

| Function | Signature | Description |
|---|---|---|
| `returnXYNames()` | `(nc)` | Returns the names and dimension positions of the X and Y (lon/lat) dimensions; handles `X/Y`, `x/y`, `longitude/latitude`, and `lon/lat` naming conventions |
| `read_times()` | `(nc, formatdates)` | Reads the time dimension, parses numeric values using the `units` attribute (e.g., `"days since 1961-01-01"`), and returns a `Date` vector or formatted strings |
| `read_epsg()` | `(nc)` | Extracts the EPSG code from global NetCDF attributes by searching for `epsg:NNNN` patterns |
| `read_coords()` | `(nc, epsg)` | Returns a matrix of lon/lat coordinates in EPSG:4326 for all grid cells |
| `readMinMax()` | `(nc)` | Computes per-time-step minimum and maximum values across all spatial cells, reading in blocks of `T/100` time steps to limit memory use |
| `mergeArrays()` | `(arr1, arr2, positions, aggregation_func)` | Merges two arrays at specified positions using an aggregation function (e.g., `min` or `max`); used to combine value ranges across spatial portions |
| `timePosition()` | `(nc)` | Returns the index of the dimension named `"time"` (case-insensitive) in the NetCDF dimension list |
| `getVarName()` | `(nc)` | Returns the name of the primary data variable, excluding `crs` and `time_bounds` |
| `get_struct_typecode()` | `(nc_type)` | Converts NetCDF precision strings (`"float"`, `"double"`, `"int"`) to struct format codes (`"f"`, `"d"`, `"i"`) |
| `create_nc_name()` | `(file_name, sufix, ext)` | Generates a derived filename by inserting a suffix before the `.nc` extension |
| `arrayRtojs()` | `(name, value, type, digits, value_array)` | Converts an R named list to a JavaScript variable assignment string (legacy; used in an older JS output mode) |


---

## 8. Input NetCDF Requirements

Input NetCDF files must satisfy:

- **Format**: NetCDF-4 (HDF5 backend); the rechunking functions use `hdf5r` for chunk metadata and require the file to be in chunked (not contiguous) HDF5 storage.
- **Dimensions**: three dimensions ordered `(lon, lat, time)` (or `(X, Y, time)`, `(longitude, latitude, time)`, `(x, y, time)`). Two-dimensional files (without a time axis) are supported via `nc_dims = 2`.
- **Time units**: CF-convention `"days since YYYY-MM-DD"` format is required for `read_times()` to parse dates correctly.
- **CRS**: any projected or geographic CRS is accepted; the EPSG code is passed explicitly and used for coordinate transformation to EPSG:4326 in `config_web()`.
- **Variable**: one primary data variable per file (plus optional `crs` and `time_bounds`). `getVarName()` selects the first variable that is not `crs` or `time_bounds`.

Files produced by `data_flow` and `indexCalc` conform to CF-1.11 conventions and meet all these requirements.


---

## 9. Usage Examples

### Minimal single-variable example

```r
library('ncartifactgenerator')

info_js <- generate_artifacts(
  nc_root     = "/genoma/aemet-sc/nc/amm",
  out_root    = "/genoma/aemet-sc/web-dev/amm",
  nc_filename = "primera_helada/ff",
  portion     = "_pen",
  var_id      = "ffd",
  epsg        = "4326",
  info_js     = NA
)

writeJson(
  folder      = "/genoma/aemet-sc/web-dev/amm",
  infoJs      = info_js,
  varTitle    = list("ffd" = "First frost day"),
  legendTitle = list("ffd" = "days")
)
```

### Multi-variable, two-domain example (standard viewer pattern)

```r
library('ncartifactgenerator')

nc_root  <- "/genoma/aemet-sc/nc/amm"
out_root <- "/genoma/aemet-sc/web-dev/amm"
epsg     <- "4326"
portions <- c("_pen", "_can")

ncVars <- data.frame(
  var  = c("ffd", "lfd", "gdd_corn", "gdd_wine", "gdd_winter_cereal"),
  file = c("primera_helada/ff", "ultima_helada/lf",
           "acum_calor/gdd_corn", "acum_calor/gdd_wine",
           "acum_calor/gdd_winter_cereal"),
  stringsAsFactors = FALSE
)

info_js <- NA
for (i in seq_len(nrow(ncVars))) {
  for (portion in portions) {
    info_js <- generate_artifacts(
      nc_root     = nc_root,
      out_root    = out_root,
      nc_filename = ncVars$file[i],
      portion     = portion,
      var_id      = ncVars$var[i],
      epsg        = epsg,
      info_js     = info_js,
      write       = FALSE
    )
  }
}

var_title <- list(
  "ffd"               = "First frost day",
  "lfd"               = "Last frost day",
  "gdd_corn"          = "Growing degree-days for corn",
  "gdd_wine"          = "Growing degree-days for wine",
  "gdd_winter_cereal" = "Growing degree-days for winter cereal"
)
legend_title <- list(
  "ffd"               = "days",
  "lfd"               = "days",
  "gdd_corn"          = "°C",
  "gdd_wine"          = "°C",
  "gdd_winter_cereal" = "°C"
)

writeJson(
  folder      = out_root,
  infoJs      = info_js,
  varTitle    = var_title,
  legendTitle = legend_title,
  minify      = TRUE
)
```

### Merging peninsular and Canary Islands domains

```r
fusion_pen_can(
  can_filename    = "/genoma/aemet-sc/nc/amm/primera_helada/ff_can.nc",
  pen_filename    = "/genoma/aemet-sc/nc/amm/primera_helada/ff_pen.nc",
  fusion_filename = "/genoma/aemet-sc/nc/amm/primera_helada/ff_all.nc",
  nc_chunk_num    = 16,
  signif_digits   = 4
)
```

### Rechunking only (without the full pipeline)

```r
# Time-series rechunking
write_nc_chunk_t(
  in_file  = "/data/ffd_pen.nc",
  out_file = "/data/nc/ffd_pen-t.nc",
  lon_by   = 100,
  lat_by   = 100,
  lon_name = "lon",
  lat_name = "lat",
  time_name = "time"
)

# Time-series binary index
write_nc_t_chunk_dir_iter(
  in_file  = "/data/nc/ffd_pen-t.nc",
  out_file = "/data/nc/ffd_pen-t.bin"
)
```
