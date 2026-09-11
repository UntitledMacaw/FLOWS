# FLOWS

**FLOWS: Flood Likelihood from pre-compiled Watershed Analysis**

> A lightweight flood-risk estimation approach designed to run on modest hardware — using a topology-based basins pre-compilation instead of full hydrodynamic simulation.

[![Status](https://img.shields.io/badge/status-in%20development-yellow)]()
[![Python](https://img.shields.io/badge/python-3.14-blue)]()
[![License](https://img.shields.io/badge/license-MIT-green)]()

---

## About the project

Flood-risk estimation can require the integration of several different types of data, including terrain,
soil properties, land cover, precipitation, etc. This can become particularly challenging for intense, short-duration rainfall events - where conditions may change quickly and useful estimates may need to be produced within
a limited time window.

FLOWS starts from a simple question: can the computationally expensive analysis of terrain topology and storage capacity be precomputed beforehand, so that flood-risk estimation can be generated
quickly as soon as new rainfall data becomes available, while accepting some loss in physical fidelity?

Rather than attempting to deal with the full complexity of a hydrodynamic model, FLOWS precomputes
the terrain structure and its potential storage relationships, then combines this information with
rainfall-derived runoff estimates. Under an intentionally conservative scenario of high antecedent moisture, the model uses AMC III soil conditions to estimate water accumulation and depth distribution. This preprocessing works by finding bowl-like depressions and their contributing areas, where runoff tends to accumulate, and computes:

- Area
- Elevation – volume relationship
- Geographical positioning
- Composite Curve Number

Once this data is compiled, rainfall forecasts can be converted into water volumes and handled by our pre-calculated values — shifting most of the computational cost from runtime to a one-time preprocessing stage.

## Goal

To investigate whether a precompiled representation of terrain storage can provide useful flood-risk estimates while remaining lightweight and quick enough to run on everyday laptops — making flood-awareness tools more accessible to small municipalities, schools, researchers, and citizen initiatives.

## Technical approach

1. Obtain a Digital Elevation Model (DEM) with acceptable resolution
2. Identify bowl-like topographic depressions ("basins")
3. Calculate the storage–volume relationship for each basin
4. Incorporate rainfall forecasts and infiltration estimates
5. Estimate water accumulation and flood risk via the already calculated basin properties.
6. Test on real-life scenarios.

### Main libraries and tools

| Tool | Purpose |
|---|---|
| **Python** | Main project language |
| **GDAL / `gdalbuildvrt`** | Virtual raster mosaic (VRT) of DEM tiles, without duplicating data on disk |
| **WhiteboxTools** | Hydrological processing (pit breaching, D8 pointer, sinks, watersheds) |
| **rasterio** | Raster (.tif) reading/writing |
| **geopandas / shapely** | Geometry handling and coordinate reprojection |
| **NumPy / SciPy** | Masking, centroid calculation, and array operations |
| **Matplotlib** | Visual diagnostics of the results |

### Elevation data

The project uses ANADEM, a DEM produced by Brazil's National Water and Basic Sanitation Agency (ANA) in collaboration with UFRGS — a refined version of Copernicus GLO-30 with reduced vegetation bias, developed specifically for hydrological analysis in South America.

> ⚠️ The data sources and parameters used here were chosen for **Brazilian territory**. To adapt this project to another region, use an equivalent DEM and parameters for your area of study. But relax: documentation will be provided.

## How it works (pipeline overview)

Since the full South American DEM exceeds 60 GB in storage, processing is done in dynamic windows ("core boxes") to avoid loading everything into memory:

1. Define a central analysis box (core box)
2. Clip the terrain with an initial margin (buffer)
3. Identify basins within that margin
4. Check whether any basin "belonging" to the core box touches the buffer's edge (which would indicate an incorrect cut)
5. If so, increase the buffer and repeat, up to a safety limit
6. If not, save the result — free from tile boundary artifacts

Whether a basin "belongs" to the core box or "leaked" to the edge is determined by comparing each basin's centroid against the original core box bounds — ensuring every basin is counted exactly once.

Next step in development: Handling the entire dataset.

## Current status

Actively in development. Already implemented:

- [x] VRT generation from DEM tiles
- [x] Bowl-like depression identification
- [x] Final mask removing basins outside the area of interest
- [x] Visual diagnostics of the process
- [x] Infiltration calculations via the Curve Number method
- [x] Pre-calculate results

In progress / next steps:

- [ ] Final flood-risk estimation for the entire dataset
- [ ] Testing with real-life events

## Reference hardware

The project is developed and tested on low-cost hardware, to validate its accessibility goal:

- Lenovo ThinkPad E470
- Intel Core i5-7200U
- 8 GB DDR4 RAM
- 120 GB SATA SSD
- Fedora Workstation 44


## 📄 License

This project is licensed under the [MIT License](LICENSE).

## Author

**Thiago Borges** aka UntitledMacaw

---

*This README documents a work-in-progress project, built by an enthusiast for study and experimentation in low-cost computational hydrology.*