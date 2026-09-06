# FLOWS

**Flood Likelihood from pre-compiled Topological Analysis**

> A lightweight, fast flood-risk estimation approach designed to run on modest hardware — using a topology-based basins instead of full hydrodynamic simulation.

[![Status](https://img.shields.io/badge/status-in%20development-yellow)]()
[![Python](https://img.shields.io/badge/python-3.14-blue)]()
[![License](https://img.shields.io/badge/license-MIT-green)]()

---

## 📖 About the project

Traditional hydraulic models (HEC-RAS, WRF-Hydro, national forecasting systems, etc) achieve high accuracy and are well-established, valuable tools. For most use cases — 1D or 2D modeling over small/medium domains — they run fine on ordinary hardware. What demands heavier resources (substantial RAM and CPU, sometimes cluster infrastructure) are specific scenarios: high-resolution 2D simulations over large domains, or large-scale probabilistic studies. Even though that heavier tier is scenario-specific, it's still a real barrier for many schools, local governments, independent researchers, and developers.

**FLOWS** starts from a simple question: *can a useful flood-risk estimate be produced, trading off some physical complexity, in order to run on consumer-grade hardware?*

The core idea: instead of solving fluid dynamics equations over the entire terrain in real time, the terrain is **pre-processed once** into a network of natural storage basins (topographic depressions). Each basin is described by:

- Area
- Storage capacity
- Elevation–volume relationship
- Overflow point
- Connections to neighboring basins

Once this network is compiled, rainfall forecasts can be converted into water volumes and propagated through the basin graph — shifting most of the computational cost from runtime to a one-time preprocessing stage.

## 🎯 Goal

To investigate whether a graph-based representation of terrain storage can provide useful flood-risk estimates while remaining lightweight enough to run on everyday laptops — making flood-awareness tools more accessible to small municipalities, schools, researchers, and citizen initiatives.

## 🛠️ Technical approach

1. Obtain a Digital Elevation Model (DEM) with acceptable resolution
2. Identify topographic depressions ("basins")
3. Calculate the storage–volume relationship for each basin
4. Build a connectivity graph between basins
5. Incorporate rainfall forecasts and infiltration estimates
6. Estimate water accumulation and flood risk via graph traversal, rather than full hydrodynamic simulation

### Tech stack

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

The project uses **ANADEM**, a DEM produced by Brazil's National Water and Basic Sanitation Agency (ANA) in collaboration with UFRGS — a refined version of Copernicus GLO-30 with reduced vegetation bias, developed specifically for hydrological analysis in South America.

> ⚠️ The data sources and parameters used here were chosen for **Brazilian territory**. To adapt this project to another region, use an equivalent DEM and parameters for your area of study. But relax: documentation will be provided.

## 🧩 How it works (pipeline overview)

Since the full South American DEM exceeds 60 GB, processing is done in **dynamic windows** ("core boxes") to avoid loading everything into memory:

1. Define a central analysis box (*core box*)
2. Clip the terrain with an initial margin (*buffer*)
3. Identify basins within that margin
4. Check whether any basin "belonging" to the core box touches the buffer's edge (which would indicate an incorrect cut)
5. If so, increase the buffer and repeat, up to a safety limit
6. If not, save the result — free from tile boundary artifacts

Whether a basin "belongs" to the core box or "leaked" to the edge is determined by comparing each basin's **centroid** against the original core box bounds — ensuring every basin is counted exactly once.

Next step in development: incorporating soil water infiltration via the **Curve Number (CN) method**, developed by the NRCS/USDA, to estimate precipitation excess as a function of cumulative rainfall, land cover, and antecedent soil moisture

## 📌 Current status

Actively in development. Already implemented:

- [x] VRT generation from DEM tiles
- [x] Bowl-like depression identification
- [x] Final mask removing basins outside the area of interest
- [x] Visual diagnostics of the process

In progress / next steps:

- [ ] Infiltration calculation via the Curve Number method
- [ ] Pre-calculate results
- [ ] Final flood-risk estimation
- [ ] Testing with real-life events

## 💻 Reference hardware

The project is developed and tested on **low-cost hardware**, to validate its accessibility goal:

- Lenovo ThinkPad E470
- Intel Core i5-7200U
- 8 GB DDR4 RAM
- 120 GB SATA SSD
- Fedora Workstation 44


## 📄 License

This project is licensed under the [MIT License](LICENSE).

## 👤 Author

**Thiago Ricardo Borges** aka UntitledMacaw

---

*This README documents a work-in-progress project, built for study and experimentation in low-cost computational hydrology.*
