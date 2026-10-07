# urban-fire

Code and figure-ready data for:

**Extreme wildfires disproportionately reach the wildland-urban interface in California**

This repository includes **MATLAB**, **Python**, and **Google Earth Engine (GEE)** scripts for data processing, analysis, and figure generation.

## Repository layout

```text
urban-fire/
├── DtPro.m, Fig1.m … Fig5.m   # Main-text MATLAB (run with cwd = repo root)
├── matlab/                    # Full MATLAB set (main + SI helpers)
├── python/                    # Revision / SI analyses and figure scripts
│   ├── *.py                   # Analysis pipelines
│   └── Figs/                  # Python figure scripts
├── GEEcode/                   # Earth Engine extraction scripts
└── dataFig/                   # Intermediate tables used for plotting
```

## MATLAB

Primary processing and main-text figures.

| Script | Role |
|--------|------|
| `DtPro.m` | Fire-size distributions, β estimates, stratified summaries |
| `Fig1.m`–`Fig5.m` | Main figures |
| `matlab/Fig6.m`, `FigSI_*.m`, `elevation.m`, `z-score.m`, … | SI / supporting MATLAB |

**How to run:** set the MATLAB current folder to the **repository root** (so paths like `dataFig/...` resolve). Full reprocessing also needs fire perimeter inputs under `dataPrc/firePrmt/` (not shipped here; see Data below).

## Python

Additional analyses and revised / SI figures (response-letter updates).

| Script | Role |
|--------|------|
| `python/calfire_fwi_gridmet.py` | Canadian FWI moisture codes (FFMC, DMC, DC) from gridMET |
| `python/fit_size_distributions_alt.py` | Size-distribution fits (power law / alternatives) |
| `python/size_dist_AIC_LR_comparison.py` | AIC and likelihood-ratio / Vuong tests |
| `python/socal_sawri_fire_flags.py` | Santa Ana wind-event flags (SoCal) |
| `python/fire_perimeter_forest_SG_composition.py` | Within-perimeter forest vs shrub/grass cover |
| `python/wui_vs_wildland_fuel_abundance.py` | Fuel / vegetation abundance contrasts |
| `python/wui_vs_wildland_lc_fractions_2020.py` | Land-cover fractions |
| `python/westcoast_wui_wildland_beta.py` | Western U.S. MTBS WUI vs wildland β |
| `python/vpdmax_landscape_fireseason.py` | Landscape VPDmax / fire-season helpers |
| `python/response_revision_figures.py` | Entry point to regenerate revision figures |
| `python/Figs/*.py` | Individual figure scripts (e.g. Fig. S5, S6, S9, S10) |

Paths are relative to the repository root (`Path(__file__).resolve().parents[...]`).

**How to run (examples):**

```bash
# from repository root
python python/Figs/FigS6.py
python python/Figs/Fig_elev_dist_Forest_SG.py
python python/response_revision_figures.py list
python python/response_revision_figures.py S1 S2
```

Typical Python dependencies: `numpy`, `pandas`, `matplotlib`, `scipy`, `geopandas` (and `ee` / `pyFWI` for FWI recompute).

## Google Earth Engine

| Script | Role |
|--------|------|
| `GEEcode/NDVI.js` | NDVI extraction |
| `GEEcode/elevation.js` | Elevation |
| `GEEcode/landCover-IGBP.js` | Land cover |
| `GEEcode/LANDFIRE_fuel_veg_CA.js` | LANDFIRE fuel / vegetation |

Run in the [Earth Engine Code Editor](https://code.earthengine.google.com/) with your project assets / fire perimeters.

## Data

- **`dataFig/`** — CSV tables used to redraw most published figures (included).
- **`dataPrc/`** — larger processed products (fire perimeters, daily weather caches, etc.) are **not** fully included. Reproduce those steps with GEE + the MATLAB/Python pipelines, or obtain companion deposits if available (e.g. Zenodo).

Outputs from Python figure scripts are written under `Fig/` by default (created on first run).

## Citation

If you use this code, please cite the paper and this repository.
