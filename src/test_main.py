from shapely.ops import polygonize_full
import matplotlib.pyplot as plt
from line_proc import get_dot_points, get_geojson_list, get_shapely_geom, plot_line, get_image_bbox
from dot_dash import extract_dot_dashed_lines
from txt_proc import recognize_texts
from settings import MODEL_PATH, TARGET_ALPHABETS, MULTITHREAD, \
    OUTPUT_DIR, INPUT_DIR, IMAGE_FILES, IMAGE_FILE, OUTER_IMAGE_BBOX_OFFSET, \
    INNER_IMAGE_BBOX_OFFSET
import rasterio
from shapely.geometry import Polygon
from line_proc import affine_transform
from gis_io import write_line_shapefile, write_poly_shapefile
import cv2
import numpy as np
from osgeo import gdal
import logging
import os
from traceback import print_exc


def main(image_file, input_dir, output_dir):
    plt.figure(figsize=(21, 19))

    image_path = os.path.join(input_dir, image_file)
    shapefile_file = image_file.split(".")[0]
    line_shapefile_path = os.path.join(os.path.join(output_dir, 'line'), shapefile_file)
    poly_shapefile_path = os.path.join(os.path.join(output_dir, 'poly'), shapefile_file)
    logfile_path = os.path.join(output_dir, f'{image_file.split(".")[0]}.txt')

    logging.basicConfig(filename=logfile_path,
                        level=logging.DEBUG,
                        format='%(asctime)s | %(name)s | %(levelname)s | %(message)s',
                        filemode='w')

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
    outer_image_bbox = get_image_bbox(image_path, offset=OUTER_IMAGE_BBOX_OFFSET)
    inner_image_bbox = get_image_bbox(image_path, offset=INNER_IMAGE_BBOX_OFFSET)
    lines = extract_dot_dashed_lines(blob_points, geoms, outer_image_bbox, inner_image_bbox)

    logging.info('Writing extracted dot dashed lines to Shapefile')
    write_line_shapefile(list(lines), line_shapefile_path)

    logging.info('Polygonizing the extracted dot dashed lines')
    result_polys, dangles, cuts, invalids = polygonize_full(lines)
    result_polys = list(result_polys)
    logging.info(f'Number of created final polygons: {len(result_polys)}')

    # ----- Label Extraction
    try:
        logging.info('Extracting labels on the map')
        ocr_result = recognize_texts(image_path, MODEL_PATH, TARGET_ALPHABETS)
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
        logging.error(err, exc_info=True)
        print_exc()
        print(err)
        labels = []

    # ----- Writing Shapefile
    logging.info('Writing extracted polygons and labels to Shapefile')
    write_poly_shapefile(result_polys, labels, poly_shapefile_path)

    logging.info('End of main')


    # plt.figure(figsize=(21, 19))

    # # plot GeoTiff
    # ds = gdal.Open(image_path)
    # gt = ds.GetGeoTransform()
    # ulx, xres, _, uly, _, yres = gt
    # lrx = ulx + (ds.RasterXSize * xres)
    # lry = uly + (ds.RasterYSize * yres)
    # extent = [ulx, lrx, lry, uly]
    # array = ds.ReadAsArray()
    #
    # plt.imshow(array, extent=extent)
    #
    # plot dots
    # for p in blob_points:
    #     plt.plot(p.x, p.y, marker="o", markerfacecolor="red")
    #
    # # plot lines
    # # plot_line(lines)
    #
    # # plot polygons
    # for poly in result_polys:
    #     plt.plot(*poly.exterior.xy)

    # plt.show()

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




