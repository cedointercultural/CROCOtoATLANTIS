# croco2atlantis

[Español](README.md) | **English**

A Python tool to convert output from the [CROCO](https://www.croco-ocean.org/) ocean model into the hydrographic forcing files of the [ATLANTIS](https://github.com/Atlantis-Ecosystem-Model) ecosystem model: temperature, salinity and water exchange between boxes and layers.

**Author:** Ricardo Cavieses.

## Status

**Under development.** The full conversion does not work yet, and the CLI only prints the version.

| Module | Status | What it does |
|---|---|---|
| `readers/croco.py` | Done and tested | Opens a CROCO file; reads the grid (with the ghost ring masked out and the open boundary kept apart), the vertical parameters and the dates; computes the depth of the s-levels (`Vtransform` 1 and 2) |
| `geometry/bgm.py` | Done and tested | Reads the BGM: boxes, faces, normals, projection |
| `vertical.py` | Done and tested | Layers of each box, overlap of s-levels with layers, and thickness-weighted remapping |
| `horizontal.py` | Done and tested | Assigns cells to boxes, computes volume-weighted means per box and layer, and fills boxes without data |
| `fluxes.py` | Done and tested | Water transport between pairs of boxes, per layer, on the faces of the CROCO cells |
| `exchanges.py` | Pending | Source and destination of each exchange, vertical flux, hyperdiffusion correction, conversion to m³ per time step and ATLANTIS layer order |
| `writers/atlantis_nc.py` | Pending | Writes `exchange`, `dest_b`, `dest_k`, temperature and salinity |
| `validate.py` | Pending | Volume, continuity, physical ranges and layer order |
| `cli.py` | Pending | CLI with YAML configuration |
| `r/` (R package) | Done and tested | R interface to the finished modules, through `reticulate` (see [r/README.md](r/README.md), in Spanish) |

## Installation and tests

```
pip install -e ".[dev]"
pytest
```

Without installing: `PYTHONPATH=src python3 -m pytest`. A single test: `pytest tests/test_fluxes.py::test_nondivergent_flow_closes_every_box`.

## Problem

CROCO uses terrain-following vertical levels (s-coordinates) on a curvilinear grid. ATLANTIS uses polygons (boxes) with fixed-thickness layers and faces shared between boxes. Values cannot be copied from one grid to the other: they have to be remapped while conserving heat, salt and water volume.

## Method

### Input data: averaged output (`avg`)

The tool uses CROCO's average files (`avg`), not the instantaneous ones (`his`). In a run with tides, one snapshot every 12 h undersamples the tide: M2 (12.42 h) shows up as a spurious oscillation of about 15 days, and S2 (exactly 12 h) is always sampled at the same phase and becomes a constant bias in sea level. Both kinds of file have the same layout, so the reader works on either. In `avg` files, the time of each record is the centre of its averaging window.

A 12 h average removes almost all of the semidiurnal tides, but leaves a good part of the diurnal ones (K1, O1).

### Vertical

- ATLANTIS layers are defined by the depth of their bottoms, for example `(18, 54, 270, 684, 918, 1782, 3100)` m.
- Each box keeps the layers whose top lies above its mean depth (`botz` in the BGM), and the last one is truncated at that depth.
- For averaging, the deepest layer of each box takes the whole water column below its top (option `deepest_to_bottom`), because some model columns are deeper than the box mean.
- The first layer starts at the free surface, so it includes the elevation `zeta`.
- CROCO levels are split among the layers by their actual overlap in depth.
- Internal convention: in CROCO index 0 is the bottom; for ATLANTIS layers index 0 is the shallowest until the files are written, where layer 0 becomes the bottom of each box.

### Horizontal

- Each CROCO water cell goes to the box that contains its centre. If the centre falls exactly on the edge between two boxes, the first one wins.
- Water cells whose centre lies outside every box, or inside a BGM island that the CROCO grid does not resolve, are moved to the nearest box with water. This way no water or flow is lost. An optional `max_distance` leaves farther cells without a box.
- The value of a box in a layer is the volume-weighted mean of its cells (cell area × thickness of its water inside the layer).
- A missing value (NaN) in a water cell is an error and stops the computation. In CROCO land is 0, not NaN, which is why `mask_rho` is always applied.
- Boxes or layers without data (boxes with no cells, or layers deeper than all of their columns) are filled with the profile of the nearest box that has data in **all** of the missing layers. The filled column is thus a real profile and not layers from different boxes. An optional `max_distance` avoids copying values from too far away, for example into boxes outside the CROCO domain.

### Transport between boxes

- Transport is computed on the faces of the cells of CROCO's C-grid, where `u` and `v` live. Each face between two cells of different boxes contributes velocity × face length × water thickness, integrated within each layer. The levels of a face are the mean of those of its two cells, as in CROCO.
- The water budget of each box is therefore the model's own divergence, consistent with the volumes used in the means.
- Velocity is not sampled along the BGM faces: tested on real data, that method gives an error in the net budget of each box of the same order as the budget itself.
- The result is per pair of boxes (`pairs`), positive from `pairs[:, 0]` to `pairs[:, 1]`.
- Box −1 stands for outside the boxes: the open boundaries of CROCO (water cells of the ghost ring) and any water cell left without a box.
- Transport layers are the standard ones, from the surface to the bottom, without truncation. Mapping them to the layers of each box is done when the exchanges are built.
- With `avg` files, the transport is the mean velocity times the thickness from the mean `zeta`. It leaves out the correlation between tidal current and tidal elevation (tidal Stokes transport), which CROCO does not store.

### Current usage (API)

```python
from croco2atlantis.readers import croco
from croco2atlantis.geometry.bgm import read_bgm
from croco2atlantis import horizontal as hz, fluxes

bottoms = [18, 54, 270, 684, 918, 1782, 3100]
ds = croco.open_output("croco_avg.nc")
grid = croco.read_grid(ds)
params = croco.read_vertical_params(ds)
bgm = read_bgm("model.bgm")

cellmap = hz.assign_cells_from_grid(grid, bgm)   # once
faces = fluxes.cell_faces(cellmap, grid)         # once

t = 0
z_w = croco.s_depths(grid.h, ds.zeta.isel(time=t).values, params, "w")
weights = hz.box_layer_weights(z_w, cellmap, grid.cell_area, bgm.depth, bottoms)
temp = weights.apply(ds.temp.isel(time=t).values)            # (box, layer)
temp, filled = hz.fill_missing(temp, bgm.depth, bottoms, centroids)
q = fluxes.face_transport(faces, ds.u.isel(time=t).values, ds.v.isel(time=t).values, z_w, bottoms)
inflow = fluxes.net_inflow(faces, q)                          # budget per box and layer
```

`centroids` is a (box, 2) array in metres in the BGM projection, for example the interior points of the boxes: `numpy.array([b.inside for b in bgm.boxes])`.

### Usage from R

The R package in `r/` calls these same modules through `reticulate` and returns data frames and `sf` objects; the computation stays in Python. Installation and an example are in [r/README.md](r/README.md):

```r
remotes::install_github("rcavieses/CROCOtoATLANTIS", subdir = "r")
croco2atlantis::install_croco2atlantis()   # Python virtualenv with the package
```

## Pending decisions

- **Open boundary:** which ATLANTIS box receives the water that enters or leaves through CROCO's open boundaries (box −1): the neighbouring box outside the domain or a boundary box.
- **Pairs without a BGM face:** since cells are assigned by their centre, two boxes that only touch at a vertex, or are separated by a thin box, can exchange water without sharing a face (`fluxes.pairs_without_face`). It remains to decide whether they stay as a direct exchange or are routed through the box in between. The hyperdiffusion correction needs the geometry of a face.
- **Diurnal tides:** whether to average further to 24 h or more (for example, the running mean of two 12 h records) before building the fluxes.
- **Vertical flux:** from continuity of the horizontal fluxes or from `w`.
- **Hyperdiffusion correction:** divide by the distance across the box along the face angle, by the length, or by the area.
- **Boxes outside the CROCO domain:** whether to fill them (and with what `max_distance`) or leave them without data.

## Change log

### 2026-10-22: R package

- New R package in `r/`: reads the BGM (boxes and faces as `sf`) and CROCO, assigns cells to boxes, computes box-layer means and the transport between boxes by calling the Python code through `reticulate`. `install_croco2atlantis()` creates the Python virtualenv.
- In R, `layer` and `time` count from 1; `box` keeps the BGM box number, from 0.

### 2026-10-22: code and logic review

- Input data: switched from `his` to `avg` because of tidal undersampling. `open_history` was renamed `open_output`.
- Transport: sampling the velocity of the nearest cell along the BGM faces was replaced by the transport on the C-grid cell faces, which closes the budget of each box. Added `cell_faces`, `net_inflow` and `pairs_without_face`; removed `rho_velocity` and `face_samples`.
- Open boundary: `CrocoGrid.open_boundary` keeps the water cells of the ghost ring; their transport goes to box −1.
- Cells without a box or on BGM islands: now moved to the nearest box with water (`move_to_water_boxes`); before, they were dropped.
- Filling: each box takes all of its missing layers from a single box; before, each layer could come from a different box.
- Cell assignment: a centre on the edge between two boxes goes to the first one; before, it went to none.
- `BoxLayerWeights.apply` accepts any array and fails on NaN in water; `s_depths` only accepts `"w"` or `"rho"`.

## Scope

- CROCO only; it does not aim to cover other ROMS models.
- Meant to be reused by other researchers working with CROCO, with their own domain and their own ATLANTIS geometry (BGM file). Everything that depends on a region (BGM, layer thicknesses, file names) will go in a YAML configuration file, not in the code.
- Internal use for now; the code is published on GitHub.
- It does not use Hydro or HydroConstruct. The method follows existing projects from ROMS and other models (NEMOW Climate Data Workshop, ROMS_Atlantis_translation, SSMtoAtlantis).

## Layout

```
src/croco2atlantis/
  readers/croco.py      # the only module that knows about CROCO
  geometry/bgm.py       # BGM reader: boxes, faces, normals
  vertical.py           # ATLANTIS layers and vertical remapping
  horizontal.py         # cell→box, box means and filling
  fluxes.py             # transport between boxes on the cell faces
  exchanges.py          # source/destination, vertical flux, hyperdiffusion and layer order
  writers/atlantis_nc.py
  validate.py
  cli.py
tests/
r/                      # R package (interface through reticulate)
examples/
```

## License

GPL-3.0 (see [LICENSE](LICENSE)).
