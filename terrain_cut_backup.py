import os
import geopandas as gpd
from shapely.geometry import box
from osgeo import gdal
from whitebox import WhiteboxTools
import rasterio
import numpy as np
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap

gdal.UseExceptions()

def plot_tile_diagnostics(
    basins_data,
    transform,
    core_bounds_5800,
    buf_deg_bounds,
    centers,
    unique_ids,
    owned_by_core_ids,
    tile_id,
):
    """
    Takes our results from the dynamic cutting and outputs them as a visual .png file so we can have visual
    proof that all centroids are contained inside the core box
    It is advised to read the dynamic_cut function before reading this one (down below).

    :param basins_data: Our basins raster array
    :param transform: Our transformation matrix (provided by rasterio)
    :param core_bounds_5800: Tuple (minx, miny, maxx, maxy) that defines our core box limits
    :param buf_deg_bounds: Tuple (buf_minx, buf_miny, buf_maxx, buf_maxy) that deines our core box + buffer limits
    :param centers: List with all coordinates (line and column) from each of the basins centroids
    :param unique_ids: Array containing all of the unique basin IDs present in the image (removing 0)
    :param tile_id: Our tile ID
    
    """

    # WARNING: DO NOT REMOVE FIG. 
    fig, ax = plt.subplots(figsize=(12, 10), facecolor='white')

    # 1. Generating unique colors for each ID
    # Generates random colors for each ID to try to avoid visual problems with neighbor basins
    max_id = int(basins_data.max())
    np.random.seed(
        42
    )
    color_array = np.random.rand(max_id + 1, 4)
    color_array[
        0
    ] = [  # Defines background color (ID 0 = noData)
        0,
        0,
        0,
        0,
    ]
    custom_cmap = ListedColormap(color_array)

    # 2. Plots basins raster with custom colors
    im = ax.imshow(
        basins_data,
        cmap=custom_cmap,
        interpolation='nearest',  # Avoids blurring
        extent=[
            transform[2],
            transform[2] + transform[0] * basins_data.shape[1],
            transform[5] + transform[4] * basins_data.shape[0],
            transform[5],
        ],
    )

    buf_minx, buf_miny, buf_maxx, buf_maxy = buf_deg_bounds
    df_buf = gpd.GeoDataFrame(
        {'geometry': [box(buf_minx, buf_miny, buf_maxx, buf_maxy)]},
        crs='EPSG:4326',
    )
    b_minx, b_miny, b_maxx, b_maxy = df_buf.to_crs('EPSG:5880').total_bounds
    c_minx, c_miny, c_maxx, c_maxy = core_bounds_5800

    # Draws Core Box limits
    core_rect = patches.Rectangle(
        (c_minx, c_miny),
        c_maxx - c_minx,
        c_maxy - c_miny,
        linewidth=2.5,
        edgecolor='red',
        facecolor='none',
        linestyle='--',
        label='Core Tile',
    )
    ax.add_patch(core_rect)

    # Draws Core Box + Buffer tile
    buf_rect = patches.Rectangle(
        (b_minx, b_miny),
        b_maxx - b_minx,
        b_maxy - b_miny,
        linewidth=2,
        edgecolor='blue',
        facecolor='none',
        label='Buffer Tile (Final)',
    )
    ax.add_patch(buf_rect)

    # Plots centroids
    first_core = True
    for basin_id, (row_center, col_center) in zip(unique_ids, centers):
        x_coord, y_coord = transform * (col_center, row_center)

        if basin_id in owned_by_core_ids:
            label_text = 'Centroid Core' if first_core else ''
            ax.plot(
                x_coord,
                y_coord,
                marker='*',  
                color='black',
                markersize=5,  
                markeredgecolor='white',
                markeredgewidth=0.5,
                label=label_text,
            )
            first_core = False

        else:
            ax.plot(
                x_coord,
                y_coord,
                marker='o',
                color='dimgray',
                markersize=2,
                alpha=0.4,
            )

    # Ajust zoom
    ax.set_xlim(b_minx, b_maxx)
    ax.set_ylim(b_miny, b_maxy)

    # Graphic aesthetics
    ax.set_facecolor('#f4f4f4')
    ax.grid(True, linestyle=':', alpha=0.5, color='gray')
    plt.title(
        f'Basins check - Tile: {tile_id}',
        fontsize=14,
        fontweight='bold',
    )
    plt.xlabel('Easting (Meters - EPSG:5880)', fontsize=11)
    plt.ylabel('Northing (Meters - EPSG:5880)', fontsize=11)
    plt.legend(loc='upper right', frameon=True, facecolor='white', shadow=True)

    plt.tight_layout()
    plt.savefig(f'diagnostic_{tile_id}.png', dpi=300)
    plt.close()
    print(
        f'    -> Image stored as: diagnostic_{tile_id}.png'
    )

def dynamic_cut(path_vrt, out_dir, core_bounds_deg, tile_id):
    """
    This will cut the DEM according to topography, preserving basin boundaries (AKA no basins get cut)
    Criteria for cutting: All basins with a centroid that belongs to core tile A will be grouped together on the same TIF

    :param path_vrt: Path to our .vrt file
    :param out_dir: Path for our output
    :param code_bounds_deg: Tuple (min_x, min_y, max_x, max_y)
    :param tile_id: Unique tile ID

    """

    os.makedirs(out_dir, exist_ok=True)
    minx, miny, maxx, maxy = core_bounds_deg

    # Projects Core Tile limist from 4326 (degress) to 5800 (meters)
    df_core = gpd.GeoDataFrame({'geometry': [box(minx, miny, maxx, maxy)]}, crs="EPSG:4326")
    df_core_5800 = df_core.to_crs("EPSG:5880")
    core_bounds_5800 = df_core_5800.total_bounds # (minx, miny, maxx, maxy) em metros

    # Dynamic buffer initial setup
    buffer_current = 0.1
    buffer_steps = 0.1
    buffer_limit = 4
    sucess = False

    wbt = WhiteboxTools()
    wbt.set_working_dir(os.path.abspath(out_dir))
    wbt.verbose = False

    while not sucess and buffer_current <= buffer_limit:
        print(f"[{tile_id}] Trying to process with a margin of {buffer_current}")

        # Define limits
        buf_minx = minx - buffer_current
        buf_miny = miny - buffer_current
        buf_maxx = maxx + buffer_current
        buf_maxy = maxy + buffer_current

        # Names of temporary files
        filename_dem = f"temp_{tile_id}_dem.tif"
        filename_dem_clean = f"temp_{tile_id}_dem_clean.tif"
        filename_d8 = f"temp_{tile_id}_d8.tif"
        filename_sinks = f"temp_{tile_id}_sinks.tif"
        filename_basins = f"basins_final_{tile_id}.tif"
        filename_depth = f"depth_final_{tile_id}.tif"

        path_dem = os.path.join(out_dir, filename_dem)
        path_basins = os.path.join(out_dir, filename_basins)

        try:
            # Cuts OG VRT
            gdal.Warp(
                destNameOrDestDS=path_dem,
                srcDSOrSrcDSTab=path_vrt,
                format='GTiff',
                outputBounds=(buf_minx, buf_miny, buf_maxx, buf_maxy),
                outputBoundsSRS='EPSG:4326',
                dstSRS='EPSG:5880',
                xRes=30, yRes=30,
                creationOptions=['COMPRESS=DEFLATE', 'TILED=YES']
            )

            # WhiteboxTools pipeline here
            print(f"    -> Running WhiteboxTools ...")
            wbt.breach_single_cell_pits(dem=filename_dem, output=filename_dem_clean)
            wbt.d8_pointer(dem=filename_dem_clean, output=filename_d8)
            wbt.sink(i=filename_dem_clean, output=filename_sinks, zero_background=True)
            wbt.watershed(d8_pntr=filename_d8, pour_pts=filename_sinks, output=filename_basins)

            # Colision Test pipeline here

            print(f"    -> Running colision test ...")

            with rasterio.open(path_basins) as src:
                basins_data = src.read(1)
                transform = src.transform

            # a) Grabs IDS from all basins that touch the bordes of the current matrix (core box + current buffer)
            top_edge = basins_data[0, :]
            bottom_edge = basins_data[-1, :]
            left_edge = basins_data[:, 0]
            right_edge = basins_data[:, -1]

            edge_ids = set(np.unique(np.concatenate([top_edge, bottom_edge, left_edge, right_edge])))
            edge_ids.discard(0) # Removes noData

            # b) Calculates the centroid of each basin present inside the raster file
            unique_ids = np.unique(basins_data)
            unique_ids = unique_ids[unique_ids > 0] # Removes noData (I hope)

            centers = ndimage.center_of_mass(basins_data > 0, basins_data, unique_ids)

            owned_by_core_ids = set()
            c_minx, c_miny, c_maxx, c_maxy = core_bounds_5800

            # c) Filters which basins have their centroid INSIDE core tile
            for basin_id, (row_center, col_center) in zip(unique_ids, centers):
                # Converts pixel (line, column) to geo cords (x, y in meters)
                x_coord, y_coord = transform * (col_center, row_center)

                if (c_minx <= x_coord <= c_maxx) and (c_miny <= y_coord <= c_maxy):
                    owned_by_core_ids.add(basin_id)

            # d) Is there any basin with centroid inside core tile that leaks outside of the buffer?
            leaking_basins = owned_by_core_ids.intersection(edge_ids)

            if len(leaking_basins) > 0:
                print(f"[FAIL] {len(leaking_basins)} basin(s) leaked outside current buffer")
                buffer_current += buffer_steps

            else:
                print(f"[SUCESS] All interesting basins fit within a buffer of {buffer_current}")

                wbt.depth_in_sink(dem=filename_dem_clean, output=filename_depth, zero_background=True)

                # testing quick visuals
                print(f"    -> Generating image ...")
                buf_bounds = (buf_minx, buf_miny, buf_maxx, buf_maxy)
                plot_tile_diagnostics(
                    basins_data,
                    transform,
                    core_bounds_5800,
                    buf_bounds,
                    centers,
                    unique_ids,
                    owned_by_core_ids,
                    tile_id,
                )
                sucess = True



        except Exception as e:
            print(f"[ERROR]: {e}")
            print(f"    -> Ending runtime ...")
            break

        # 
        temp_files = [filename_dem, filename_dem_clean, filename_d8, filename_sinks]
        if not sucess:
            temp_files.extend([filename_basins, filename_depth]) # Deletes if failed

        for f in temp_files:
            path_f = os.path.join(out_dir, f)
            if os.path.exists(path_f):
                os.remove(path_f)

    if not sucess and buffer_current >= buffer_limit:
        print(f"[WARNING]: Buffer limit reached. Consider revising terrain.")

if True:
    PATH_VRT = '/home/untitledmacaw/Downloads/output.vrt'
    OUT_DIR = 'test_basins_dynamic'
    BOUNDS_TESTING = (-37.7, -9.7, -35.2, -9.2)

    dynamic_cut(
        path_vrt=PATH_VRT,
        out_dir=OUT_DIR,
        core_bounds_deg=BOUNDS_TESTING,
        tile_id="test_01"
    )
