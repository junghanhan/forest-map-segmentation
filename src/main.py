from shapely.ops import polygonize_full
import matplotlib.pyplot as plt
from line_proc import get_dot_points, get_geojson_list, get_shapely_geom, plot_line
from dot_dash import extract_dot_dashed_lines
from txt_proc import recognize_texts
from settings import MODEL_PATH, TARGET_ALPHABETS, MULTITHREAD, PLOT_RESULT, \
    OUTPUT_DIR, INPUT_DIR, IMAGE_FILES, IMAGE_FILE
import rasterio
from shapely.geometry import Polygon
from line_proc import affine_transform
from gis_io import write_line_shapefile, write_poly_shapefile
import os
import logging


def main(image_file, input_dir, output_dir, model_path=MODEL_PATH, target_alphabet=TARGET_ALPHABETS):
    image_path = os.path.join(input_dir, image_file)
    shapefile_file = image_file.split(".")[0]
    line_shapefile_path = os.path.join(os.path.join(output_dir, 'line'), shapefile_file)
    poly_shapefile_path = os.path.join(os.path.join(output_dir, 'poly'), shapefile_file)
    logfile_path = os.path.join(output_dir, f'{image_file.split(".")[0]}.txt')

    logging.basicConfig(filename=logfile_path,
                        level=logging.DEBUG,
                        format='%(asctime)s | %(name)s | %(levelname)s | %(message)s',
                        filemode='w')

    # TODO: temporary exception handling for debugging purpose
    try:
        # ----- Line Extraction
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
        lines = extract_dot_dashed_lines(blob_points, geoms)

        logging.info('Writing extracted dot dashed lines to Shapefile')
        write_line_shapefile(list(lines), line_shapefile_path)

        logging.info('Polygonizing the extracted dot dashed lines')
        result_polys, dangles, cuts, invalids = polygonize_full(lines)
        result_polys = list(result_polys)
        logging.info(f'Number of created final polygons: {len(result_polys)}')

        # ----- Label Extraction
        try:
            logging.info('Extracting labels on the map')
            ocr_result = recognize_texts(image_path, model_path, target_alphabet)
            # plot_prediction_result(IMAGE_PATH, prediction_result)

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
            logging.debug(err)
            labels = []

        # ----- Writing Shapefile
        logging.info('Writing extracted polygons and labels to Shapefile')
        write_poly_shapefile(result_polys, labels, poly_shapefile_path)

        logging.info('End of main')
    except Exception as err:
        logging.debug(err)


if __name__ == '__main__':
    if MULTITHREAD:
        # multi thread
        from multiprocessing import Pool

        with Pool() as pool:
            pool.starmap(main,
                         [(image_file, INPUT_DIR, OUTPUT_DIR)
                          for image_file in IMAGE_FILES])
    else:
        # single thread
        main(IMAGE_FILE, INPUT_DIR, OUTPUT_DIR)




