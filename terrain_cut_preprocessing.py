import os
import geopandas as gpd
from shapely.geometry import box, shape
from osgeo import gdal
from whitebox import WhiteboxTools
import rasterio
import rasterio.features
import numpy as np
import pandas as pd
import scipy.ndimage as ndimage
import matplotlib
matplotlib.use('Agg') # avoiding some Rich + Tkinter issues
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
from rich.console import Console
import time
from datetime import datetime
from contextlib import contextmanager

console = Console()

gdal.UseExceptions()

LOG_FILE = None

def log_noprint(*args, **kwargs):
    """
    Little brother of log_print function just below this one. Saves status to log without printing it to the user.
    """

    if LOG_FILE:
        try:
            message = " ".joim(str(a) for a in args)
            with open(LOG_FILE, "a", encoding="utf-8") as f:
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                f.write(f"[{timestamp}] {message}\n")
                
        except Exception:
            pass

def log_print(*args, **kwargs):
    """
    Normal print, but this time it also stores stuff to log!
    """

    print(*args, **kwargs)
    if LOG_FILE:
        try:
            message = " ".join(str(a) for a in args)
            with open(LOG_FILE, "a", encoding="utf-8") as f:
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                f.write(f"[{timestamp}] {message}\n")

        except Exception:
            pass

@contextmanager
def logged_status(msg, spinner="dots"):
    """
    Wrapper for console.status that also saves stuff in log
    """
    log_noprint(f"[STATUS] {msg.strip()}")
    with console.status(msg, spinner="dots") as status:
        yield status

def generate_core_boxes(vrt_path, width_deg, height_deg):
    """
    Will take the dimentions of the vrt_path and cut it into core boxes for our other functions

    :param vrt_path: Path to VRT file
    :param width_deg: Desired core box width
    :param height_deg: Desired core box height
    """
    try:
        with rasterio.open(vrt_path) as src:
            bounds = src.bounds
            vrt_minx, vrt_miny, vrt_maxx, vrt_maxy = bounds.left, bounds.bottom, bounds.right, bounds.top
            log_print(f"    -> Total VRT extension:\n   Lon: {vrt_minx} - {vrt_maxx}\n   Lat: {vrt_miny} - {vrt_maxy}")

        boxes = []

        x_coords = np.arange(vrt_minx, vrt_maxx, width_deg)
        y_coords = np.arange(vrt_miny, vrt_maxy, height_deg)

        tile_count = 1

        for y in y_coords:
            for x in x_coords:
                box_maxx = min(x + width_deg, vrt_maxx)
                box_maxy = min(y + height_deg, vrt_maxy)

                bbox = (x, y, box_maxx, box_maxy)

                boxes.append({
                    "tile_id": f"tile_{tile_count:04d}",
                    "bounds": bbox
                })

                tile_count += 1

        log_print(f"    -> Total of generated core boxes: {len(boxes)}")

        return boxes
    except Exception as e:
        raise RuntimeError(f"Error while generating core boxes: {e}")

def simulate_basin_volumes_worst_case(path_basins, path_depth, path_cn, owned_ids, out_dir, tile_id):
    """
    Analyses individual basins, and for each one calculates:
    -> Composite CN (after we convert each one from AMC II to AMC III)
    -> Volume / Water height relations
    -> Other important stuff

    and then stores the values in a cute table!

    :param path_basins: Path to our basins_final .tif file
    :param path_depth: Path to our depth_final .tif file
    :param path_cn: Path to our CN table .tif file (observation: their CN values are AMC II standard. We must convert them first)
    :param owned_ids: Our owned_by_core_ids variable (go to the main function - dynamic_cut - if in doubt)
    :param out_dir: Our output directory
    :param tile_id: Our tile id (go to the main function if in doubt)

    """
    try:
        # saves different storm strenghs to test (in mm)
        # according to opticalscientific.com:
        # light rain: < 60mm / 24hrs ( < 2.5mm / hour, and we converted it to 24hrs)
        # moderate rain: 60 to 180mm / 24hrs (2.5 to 7.5 mm / hour)
        # heavy rain: 180mm to 1200mm / 24hrs (7.5 to 50mm / hour)
        # Violent rain: > 1200 / 24hrs ( > 50mm / hour)
        # here, we will add some in-between values for calculations' sake
        scenarios = np.array([30, 60, 100, 150, 250, 400, 700], dtype=float)

        all_dfs = []

        with rasterio.open(path_basins) as src_b, \
            rasterio.open(path_depth) as src_d, \
            rasterio.open(path_cn) as src_cn:   

            # Validates dimentions (between our CN and DEM)
            if not (src_b.shape == src_d.shape == src_cn.shape):
                raise ValueError(
                    f"[ALIGNMENT ERROR] Dimentions are not compatible\n"
                    f"      -> Basins: {src_b.shape}\n"
                    f"      -> Depth:  {src_d.shape}\n"
                    f"      -> CN:     {src_cn.shape}"
                )

            # 2. Validates resolution and geographical positioning (Geotransform)
            # Allows a tiny tiny tolerance for floating point Warp imperfections
            if not np.allclose(src_b.transform, src_cn.transform, atol=1e-5):
                raise ValueError(
                    "[ALIGNMENT ERROR] The matrixes cover entirely different areas or have different spacial resolution"
                )

            basins_data = src_b.read(1)
            depth_data = src_d.read(1)
            cn_data = src_cn.read(1)

            pixel_width = abs(src_b.transform[0])
            pixel_height = abs(src_b.transform[4])
            pixel_area = pixel_width * pixel_height

            for basin_id in owned_ids:
                mask = (basins_data == basin_id)
                if not np.any(mask):
                    continue

                basin_area = np.sum(mask) * pixel_area

                valid_cn = cn_data[mask]
                valid_cn = valid_cn[valid_cn > 0]
                # Converts array to float to prevent uint8 overflow

                if len(valid_cn) == 0:
                    log_print(f"    -> Skipping basin {basin_id}: No valid CN data.")
                    continue

                valid_cn_float = valid_cn.astype(float)

                # converts entire matrix from AMC II to AMC III pixel by pixel
                valid_cn_iii = (23 * valid_cn_float) / (10 + 0.13 * valid_cn_float)

                # calculates composite already on the worst case scenario
                cn_iii = np.mean(valid_cn_iii)
                cn_ii_mean_original = np.mean(valid_cn)

                # SCS CN with saturated soil
                S_iii = (25400 / cn_iii) - 254
                Q_iii = np.zeros_like(scenarios)

                valid_p_iii = scenarios > (0.2 * S_iii)
                Q_iii[valid_p_iii] = ((scenarios[valid_p_iii] - 0.2 * S_iii)**2) / (scenarios[valid_p_iii] + 0.8 * S_iii)

                volume_rain_m3_iii = (Q_iii / 1000) * basin_area

                # SCS CN on regular moisture conditions
                S_ii = (25400 / cn_ii_mean_original) - 254
                Q_ii = np.zeros_like(scenarios)

                valid_p_ii = scenarios > (0.2 * S_ii)
                Q_ii[valid_p_ii] = ((scenarios[valid_p_ii] - 0.2 * S_ii)**2 ) / (scenarios[valid_p_ii] + 0.8 * S_ii)

                volume_rain_m3_ii = (Q_ii / 1000) * basin_area

                #height interpolation

                depths = depth_data[mask]
                max_depth = depths.max()

                h_range = np.linspace(0, max_depth, 50)
                v_range = []

                for h in h_range:
                    water_depth = np.maximum(0, h - (max_depth - depths))
                    v_range.append(np.sum(water_depth) * pixel_area)

                final_heights_iii = np.interp(volume_rain_m3_iii, v_range, h_range)
                final_heights_ii = np.interp(volume_rain_m3_ii, v_range, h_range)

                #exporting to final file

                df = pd.DataFrame({
                    'Tile ID': tile_id,
                    'Basin ID': int(basin_id),
                    'Precipitation_mm': scenarios,
                    'AMC_II_Original': np.round(cn_ii_mean_original, 2),
                    'AMC_III_Saturated': np.round(cn_iii, 2),
                    'Runoff_Q_mm_AMC_II': np.round(Q_ii, 2),
                    'Volume_m3_AMC_II': np.round(volume_rain_m3_ii, 2),
                    'Max_height_m_AMC_II': np.round(final_heights_ii, 2),
                    'Runoff_Q_mm_AMC_III': np.round(Q_iii, 2),
                    'Volume_m3_AMC_III': np.round(volume_rain_m3_iii, 2),
                    'Max_height_m_AMC_III': np.round(final_heights_iii, 2)
                })

                all_dfs.append(df)

        if all_dfs:
            combined_df = pd.concat(all_dfs, ignore_index=True)
            return combined_df

        return None # returns nothing if we don't find anything
    except Exception as e:
        raise RuntimeError(f"Error on basin rainfall simulaton: {e}")


def plot_tile_diagnostics(
    basins_data,
    transform,
    core_bounds_5800,
    buf_deg_bounds,
    centers,
    unique_ids,
    owned_by_core_ids,
    tile_id,
    out_dir
):
    """
    Takes our results from the dynamic cutting and outputs them as a visual .png file so we can have visual
    proof that all centroids are contained inside the core box
    It is advised to read the dynamic_cut function before reading this one (down below).
    """
    try:
        # WARNING: DO NOT REMOVE FIG. 
        fig, ax = plt.subplots(figsize=(12, 10), facecolor='white')

        # 1. Generating unique colors for each ID
        # Generates random colors for each ID to try to avoid visual problems with neighbor basins
        max_id = int(basins_data.max())
        np.random.seed(42)
        color_array = np.random.rand(max_id + 1, 4)
        color_array[0] = [0, 0, 0, 0]  # Defines background color (ID 0 = noData)
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
                    marker='o',  # <-- Alterado para o ponto preto aparecer!
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
                    marker='o',  # <-- Alterado para o ponto cinza aparecer!
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
        img_path = os.path.join(out_dir, f'diagnostic_{tile_id}.png')
        plt.savefig(img_path, dpi=300)
        plt.close()
        log_print(
            f'    -> Image stored as: diagnostic_{tile_id}.png'
        )

    except Exception as e:
        raise RuntimeError(f"Erro while ploting image diagnostics: {e}")

def dynamic_cut(path_vrt, out_dir, core_bounds_deg, tile_id, path_cn_input):
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

    df_tile = None

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
        log_print(f"[{tile_id}] Trying to process with a margin of {buffer_current}")

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
        filename_cn = f"temp_{tile_id}_cn.tif"

        path_dem = os.path.join(out_dir, filename_dem)
        path_basins = os.path.join(out_dir, filename_basins)
        path_cn_temp = os.path.join(out_dir, filename_cn)

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

            # Cuts our CN 
            path_cn_original = path_cn_input
            gdal.Warp(
                destNameOrDestDS=path_cn_temp,
                srcDSOrSrcDSTab=path_cn_original,
                format='GTiff',
                outputBounds=(buf_minx, buf_miny, buf_maxx, buf_maxy),
                outputBoundsSRS='EPSG:4326',
                dstSRS='EPSG:5880',
                xRes=30, yRes=30,
                dstNodata=0, # <-- makes sure to have 0's on noDATA zones
                creationOptions=['COMPRESS=DEFLATE', 'TILED=YES']
            )

            # WhiteboxTools pipeline here
            with logged_status("Running WhiteboxTools ", spinner="dots"):
                wbt.breach_single_cell_pits(dem=filename_dem, output=filename_dem_clean)
                wbt.d8_pointer(dem=filename_dem_clean, output=filename_d8)
                wbt.sink(i=filename_dem_clean, output=filename_sinks, zero_background=True)
                wbt.watershed(d8_pntr=filename_d8, pour_pts=filename_sinks, output=filename_basins)


                # Colision Test pipeline here
            with logged_status("Running colision test ", spinner="dots"):
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

                # c.1) Some basins may be basically just ocean, thus won't containg
                # any basins. Let's skip the core boxes where that happens
                if len(owned_by_core_ids) == 0:
                    log_print(f"[SKIP] No valid basins found (likely pure ocean or Nodata).")

                    temp_files = [filename_dem, filename_dem_clean, filename_d8, filename_sinks, filename_cn, filename_basins]
                    for f in temp_files:
                        path_f = os.path.join(out_dir, f)
                        if os.path.exists(path_f):
                            os.remove(path_f)
                            
                    return None

                # d) Is there any basin with centroid inside core tile that leaks outside of the buffer?
                leaking_basins = owned_by_core_ids.intersection(edge_ids)

            if len(leaking_basins) > 0:
                log_print(f"[FAIL] {len(leaking_basins)} basin(s) leaked outside current buffer")
                buffer_current += buffer_steps

            else:
                log_print(f"[SUCESS] All interesting basins fit within a buffer of {buffer_current}")

                wbt.depth_in_sink(dem=filename_dem_clean, output=filename_depth, zero_background=True)

                # Removing excess basins
                with logged_status("Applying mask to keep only basins of interest ", spinner="dots"):

                    # a) Transforms the valid IDs set into a numpy array
                    core_ids_array = np.array(list(owned_by_core_ids))

                    # b) Creates a boolean mask (True if the pixel belongs to a basin of interest)
                    mask = np.isin(basins_data, core_ids_array)

                    # c) Checks which noData value is used by DEM and applies mask over basin_data
                    with rasterio.open(path_basins, 'r+') as dst:
                        lines, columns = basins_data.shape
                        area_m_square = (lines * 30) * (columns * 30)
                        log_print(f"    -> Total area processed: {area_m_square / 1000000:.2f} km2")
                        nodata_val = dst.nodata if dst.nodata is not None else 0

                        filtered_basins = np.where(mask, basins_data, nodata_val)
                        dst.write(filtered_basins, 1)

                    # d) Open, filers and overwrites depth TIF
                    path_depth = os.path.join(out_dir, filename_depth)
                    with rasterio.open(path_depth, 'r+') as dst:
                        depth_data = dst.read(1)
                        nodata_val_depth = dst.nodata if dst.nodata is not None else 0

                        filtered_depth = np.where(mask, depth_data, nodata_val_depth)
                        dst.write(filtered_depth, 1)

                log_print(f"    -> Mask applied!")
                # testing quick visuals
                with logged_status("Generating image ", spinner="dots"):
                    buf_bounds = (buf_minx, buf_miny, buf_maxx, buf_maxy)
                    plot_tile_diagnostics(
                        filtered_basins,
                        transform,
                        core_bounds_5800,
                        buf_bounds,
                        centers,
                        unique_ids,
                        owned_by_core_ids,
                        tile_id,
                        out_dir
                    )

                # calculating basin geometry and volume / water height table, alongside some other
                # important stuff
                with logged_status("Calculating worst-case hydrological scenarios\n  Alonside some other important stuff ", spinner="dots"):
                    shapes_gen = rasterio.features.shapes(
                        filtered_basins.astype(np.int32),
                        transform=transform
                    )

                    records = []

                    for geom, value in shapes_gen:
                        if value > 0 and int (value) in owned_by_core_ids:
                            poly = shape(geom)
                            records.append({
                                'basin_id': int(value),
                                'geometry': poly
                            })

                    if records:
                        gdf_basins = gpd.GeoDataFrame(records, crs="EPSG:5880")
                        gdf_basins = gdf_basins.dissolve(by='basin_id', as_index=False)
                        gdf_basins['area_m2'] = gdf_basins.geometry.area
                        gdf_basins['Tile_ID'] = tile_id

                        gpkg_path = os.path.join(out_dir, f'core_basins_{tile_id}.gpkg')
                        gdf_basins.to_file(gpkg_path, driver="GPKG")
                        log_print(f"    -> Core basin vector stored as: core_basins_{tile_id}.gpkg")

                    df_tile = simulate_basin_volumes_worst_case(
                        path_basins=path_basins,
                        path_depth=path_depth,
                        path_cn=path_cn_temp,
                        owned_ids=owned_by_core_ids,
                        out_dir=out_dir,
                        tile_id=tile_id
                    )
                sucess = True

        except Exception as e:
            raise RuntimeError(f"Critical Failure on {tile_id}: {e}")

        finally:
            # Limpeza incluindo o filename_cn!
            temp_files = [filename_dem, filename_dem_clean, filename_d8, filename_sinks, filename_cn]
            if not sucess:
                temp_files.extend([filename_basins, filename_depth]) # Deletes if failed

            for f in temp_files:
                path_f = os.path.join(out_dir, f)
                if os.path.exists(path_f):
                    os.remove(path_f)

    if not sucess and buffer_current >= buffer_limit:
        raise RuntimeError(f"Buffer limit reached on {tile_id}. Consider revising terrain.")

    return df_tile

if __name__ == "__main__":
    log_print("Hello! Welcome to FLOWS command tool")
    time.sleep(0.5)
    
    PATH_VRT = str(input("Path for VRT: ")).strip('\"')
    PATH_CN_INPUT = str(input("Path for CN raster: ")).strip('\"')
    OUT_DIR = str(input("Output folder name (e.g., output_dir): ")).strip('\"')
    
    TILE_WIDTH = float(input("Tile width in degrees (e.g., 2.5): "))
    TILE_HEIGHT = float(input("Tile height in degrees (e.g., 0.5): "))

    os.makedirs(OUT_DIR, exist_ok=True)
    LOG_FILE = os.path.join(OUT_DIR, "log.txt")
    
    time.sleep(0.5)
    log_print(f"Log will be saved in {OUT_DIR}/log.txt")
    time.sleep(0.5)

    start_time = time.time()

    with logged_status("Generating grid", spinner="dots"):
        all_core_boxes = generate_core_boxes(PATH_VRT, TILE_WIDTH, TILE_HEIGHT)

    all_tile_dfs = []

    for box_info in all_core_boxes:
        tile_name = box_info["tile_id"]
        tile_bounds = box_info["bounds"]

        try:
            with logged_status(f"Processing {tile_name} | Bounds: {tile_bounds}"):
                df_result = dynamic_cut(
                    path_vrt=PATH_VRT,
                    out_dir=OUT_DIR,
                    core_bounds_deg=tile_bounds,
                    tile_id=tile_name,
                    path_cn_input=PATH_CN_INPUT
                )

                if df_result is not None:
                    all_tile_dfs.append(df_result)

        except Exception as e:
            log_print(f"[ERROR]: {e}")
            time.sleep(0.5)
            log_print(f"    -> Aborting all subsequent tiles")
            break

    with logged_status("Saving CSV ", spinner="dots"):
        if all_tile_dfs:
            final_df = pd.concat(all_tile_dfs, ignore_index=True)
            final_csv_path = os.path.join(OUT_DIR, 'FLOWS_all_basins.csv')
            final_df.to_csv(final_csv_path, index=False)
            log_print(f"    -> CSV file saved as FLOWS_all_basins.csv")
        else:
            log_print(f"[WARNING] No basins were processed. No CSV generated")

    with logged_status("Generating main GeoPackage and cleaning auxiliaries"):
            gpkg_files = [
                os.path.join(OUT_DIR, f)
                for f in os.listdir(OUT_DIR)
                if f.startswith("core_basins_tile_") and f.endswith(".gpkg")
            ]

            if gpkg_files:
                gdf_list = [gpd.read_file(f) for f in gpkg_files]

                merge_gdf = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True), crs=gdf_list[0].crs)

                final_gpkg_path = os.path.join(OUT_DIR, 'FLOWS_basins_map.gpkg')
                merge_gdf.to_file(final_gpkg_path, driver="GPKG")
                log_print(f"    -> Merged map saves as FLOWS_basins_map.gpkg")

                # cleaning up temp files
                for f in gpkg_files:
                    os.remove(f)

                log_print(f"    -> Temporary files removed sucessfully")
            
    end_time = time.time()
    elapsed_time = end_time - start_time
    mins = int(elapsed_time // 60)
    secs = int(elapsed_time % 60)
    log_print("[END OF RUNTIME] Program finished")
    log_print(f"Total elapsed time: {mins}m {secs}s")