# FLOWS
(Flood Risk Likelihood from pre-compiled Watershed Analysis)

**A python tool for generating and evaluating flood-risk data without relying on constant hydrodynamic simulation and high-end hardware**

[![Status](https://img.shields.io/badge/status-in%20development-yellow)]()
[![Python](https://img.shields.io/badge/python-3.14-blue)]()
[![License](https://img.shields.io/badge/license-MIT-green)]()

[![alt screenshot](images/running.png)]()

Github web page: https://untitledmacaw.github.io/FLOWS/

## About
Rather than dealing with complex hydrodynamic models, FLOWS analyses the terrain and finds bowl-like depressions - where water
tends to accumulate in the terrain - and their contributing areas (where rainwater flows down until it gets trapped in a respective
bowl-like depression) so it can analyse them. Why is that good: For a fixed amount n of rainfall (in mm), the maximum height that water will reach in a respective  bowl-like depression + 
contributing area (let's call this union a "basin") will stay the same. This means that we can take basins and simulate beforehand various rainfall and antecedent soil moisture scenarios
and store those values so we can simply access them later - thus saving up a lot of computational power.

## How it works

[![alt finding](images/diagnostic.png)]()

FLOWS finds basins (bowl-like depressions + contributing areas) and for each one stores:

- Area
- Water height for each rainfall and soil moisture scenarios
- Geographical positioning (so we can actually know where basins are in the map)
- Composite Curve Number (for rainfall infiltration)

> IMPORTANT NOTE: The SCS Curve Number method can be faulty on some certain types of terrain, and therefore data needs adjustments before it can be safely used with it. Be aware when using this tool for your specific region. For the author's needs (Brazilian territory), adjustments and notes were used from: Sartori et al (2005).

For actually finding basins without exceeding RAM usage limits, FLOWS cuts the entire dataset into multiple analysis boxes (core boxes), and for each one:

1. Inserts an initial small buffer of surrounding terrain
2. Runs WhiteboxTools to identify basins within the total boundaries (initial core box boundaries + buffer)
3. Runs a collision test where it checks whether any basins that belong to the core box touch the edge (defined by the total boundaries)
4. If a basin touch the edge, that indicates an incorrect cut, so we increase the buffer and go back to step 2
5. If not, we save our results (making sure we will be free from boundary artifacts)
6. We will keep increasing the edge until a safety limit (so the computer doesn't crash)

> For checking whether a basin belongs in a respective core box, we compare its centroid (from the 2D top-down view) and check if it is contained within the initial core box boundaries (excluding buffer). Since the centroid is unique, we make sure that basins are counted once.

## Usage and requirements
The FLOWS command tool can be found on the `terrain_cut_preprocessing.py` file (other alternatives exist for backup purposes).
FLOWS will ask you to input:
- Path for VRT (must be in EPSG:4326 coordinate system, so far no other systems are supported. Make sure to convert your dataset before proceeding)
- Path for Curve Number raster (must also be in EPSG:4326 coordinate system) containing CN values for each pixel.
- A folder for storing output

It will return:
- A .csv table containing water height data from each scenario from each basin from each core box
- A geopackage file (.gpkg) that serves as a map of where all of our basins are in the real world
- A bunch of images from individual core boxes showing how basin finding turned out visually

Here is an image of the .csv table after a test:

[![alt table](images/data.png)]()

As of requirements, Python 3.14 is used and all libraries can be found in requirements.txt, and can be installed with pip by running:
`pip install -r requirements.txt`
> NOTE: you must have GDAL installed before running pip install -r requirements.txt, or else you'll get an error. Follow your distro specific instructions for getting it.

## Acknowledgements
FLOWS development heavily depended on:
- Brazil's National Water and Basic Sanitation Agency, that not only provided a DEM with reduced vegetation bias and a CN raster but also answered my emails and explained doubts!
- Python: Main programming language
- GDAL (`gdalbuildvrt`): Creating our virtual raster mosaic so we don't have to merge our dataset (that was provided in tiles. If your data is provided as a single big DEM, you probably can skip that step)
- WhiteboxTools: Terrain analysis for finding basins (pit breaching, D8 pointer, sink, and other fancy geographical tools).
- rasterio: Raster (.tif) reading/writing
- geopandas / shapely: Geometry handling and coordinates reprojection
- NumPy / SciPy: Masking (removing unwanted basins), centroid calculation and array stuff.
- Matplotlib: Visual diagnostics
- Terminal CSS: Making a pretty github page!

and many other libraries that were very important to the project.

## Hardware for this project (call this minimal requirements if you want to)
To prove FLOWS is able to be run on low-end hardware, here are the specs in which the project has been built from the ground up:

- Lenovo Thinkpad E470
- Intel Core i5-7200U (4) @ 3.10 GHz
- 8 GB of DDR4 RAM
- 120GB SATA SSD
> NOTE: 120 GB of storage is acceptable here since the dataset used by me (the author) weighted around 62GB. If your DEM and stuff weights more than that, consider revising storage for your own needs.

## Current compiling process (Brazil)
The FLOWS tool present in `terrain_cut_preprocessing.py` has currently been running on my hardware and using ANADEM as dataset. ANADEM is a DEM provided by Brazil's National Water and Basic Sanitation Agency (ANA) in collaboration with the Federal University of Rio Grande do Sul (UFGRS). It is a refined version of Copernicus GLO-30, with reduced vegetation bias, made just for hydrological analysis like we're doing here.

The FLOWS tool cut the dataset into 2820 core boxes, and so far it still hasn't finished all of them.

Core boxes compiled so far: 0 / 2820 (last updated September 26, 2026)

When runtime is over, that data will be made open and available to everyone that wants to use it (usage rules determined by the project's license -> MIT License)

## AI usage
AI was used for learning how to use WhiteboxTools and debugging / verifying code and its output
AI did not develop the project from scratch based on a "please do something that does that" prompt, neither it developed the main idea of finding bowl-like depressions and finding them.

## License
This project is licensed under the [MIT License](LICENCE).

## Author

Thiago Borges aka UntitledMacaw

*Note: this is a work-in-progress project, built by an enthusiast for study and experimentation in low-cost computational hydrology. Take that into consideration when applying it to your specific use case*
