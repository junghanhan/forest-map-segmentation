import warnings
import os
import sys
from shapely.ops import polygonize_full
import matplotlib.pyplot as plt
from line_proc import get_dot_points, get_geojson_list, get_shapely_geom, plot_line, get_image_bbox
from dot_dash import extract_dot_dashed_lines
from txt_proc import recognize_texts
from settings import MODEL_PATH, TARGET_ALPHABETS, OUTPUT_DIR
import rasterio
from shapely.geometry import Polygon
from line_proc import affine_transform
from gis_io import write_line_shapefile, write_poly_shapefile
import logging
from traceback import print_exc
from multiprocessing import Pool, cpu_count, current_process


def main(p_idx: int, image_file: str, input_dir: str, output_dir: str) -> None:
    image_path = os.path.join(input_dir, image_file)
    shapefile_file = image_file.split(".")[0]
    line_shapefile_path = os.path.join(os.path.join(output_dir, 'line'), shapefile_file)
    poly_shapefile_path = os.path.join(os.path.join(output_dir, 'poly'), shapefile_file)
    logfile_path = os.path.join(output_dir, f'{image_file.split(".")[0]}.txt')

    logging.basicConfig(filename=logfile_path,
                        level=logging.INFO,
                        format='%(asctime)s | %(name)s | %(levelname)s | %(message)s',
                        filemode='w')

    try:
        # ----- Line Extraction
        print(f'{p_idx}:Extracting dot-dashed lines')
        # dot detection
        logging.info('Detecting dots')
        blob_points = get_dot_points(image_path)

        # get polygons' geojsons from raster image
        logging.info('Getting GeoJSON from image')
        geojson_list = get_geojson_list(image_path, 0)  # masking black
        logging.info(f'Number of found polygons : {len(geojson_list)}')

        # get all polygons on the image
        logging.info('Getting Shapely polygons from GeoJSONs')
        geoms = get_shapely_geom(geojson_list)

        logging.info('Extracting dot dashed lines')
        image_bbox = get_image_bbox(image_path)
        lines = extract_dot_dashed_lines(blob_points, geoms, image_bbox)

        logging.info('Writing extracted dot dashed lines to Shapefile')
        write_line_shapefile(list(lines), line_shapefile_path)

        print(f'{p_idx}:Polygonizing the extracted dot-dashed lines')
        logging.info('Polygonizing the extracted dot dashed lines')
        result_polys, dangles, cuts, invalids = polygonize_full(lines)
        result_polys = list(result_polys)
        logging.info(f'Number of created final polygons: {len(result_polys)}')

        # ----- Label Extraction
        print(f'{p_idx}:Extracting labels')
        try:
            logging.info('Extracting labels on the map')
            ocr_result = recognize_texts(image_path, MODEL_PATH, TARGET_ALPHABETS)
            # plot_prediction_result(image_path, ocr_result)

            # affine transform for recognized labels
            labels = []
            transform = rasterio.open(image_path).transform
            for word, bbox in ocr_result:
                geo_box_coords = []
                for px_coord in bbox:
                    geo_coord = affine_transform(transform, px_coord)
                    geo_box_coords.append(geo_coord)
                labels.append((word, Polygon(geo_box_coords).centroid))
        except Exception as err:
            logging.error(err, exc_info=True)
            print_exc()
            print(err)
            labels = []

        # ----- Writing Shapefile
        print(f'{p_idx}:Writing Shapefile')
        logging.info('Writing extracted polygons and labels to Shapefile')
        write_poly_shapefile(result_polys, labels, poly_shapefile_path)

        print(f'{p_idx}:Done Conversion')
        logging.info('End of main')
    except Exception as err:
        logging.error(err, exc_info=True)
        print_exc()
        print(err)
        print(f'Please check the log ({logfile_path})')


# suppress warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


if __name__ == '__main__':
    # command line argument check
    if len(sys.argv) < 2:
        print('usage: main.py input_gtiff_file1 input_gtiff_file2 input_gtiff_file3 ...')
        sys.exit(2)
    # if the number of input images is more than the hardware capability, terminate
    elif len(sys.argv) > cpu_count() + 1:
        print('error: Too many input files')
        print(f'cpu count: {cpu_count()}')

    args = []
    # images are processed in parallel
    for i in range(1, len(sys.argv)):
        input_dir, image_file = os.path.split(sys.argv[i])
        args.append((i, image_file, input_dir, OUTPUT_DIR))

    if len(args) == 1:
        p_idx, image_file, input_dir, output_dir = args[0]
        main(p_idx, image_file, input_dir, output_dir)
    else:
        with Pool() as pool:
            pool.starmap(main, args)

    print(f'All output files are stored in: {OUTPUT_DIR}')





