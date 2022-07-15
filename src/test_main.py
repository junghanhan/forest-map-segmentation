from shapely.ops import polygonize_full
import matplotlib.pyplot as plt
from line_proc import get_dot_points, get_geojson_list, get_shapely_geom, plot_line
from dot_dash import draw_dot_dashed_lines
from txt_proc import recognize_texts
from settings import IMAGE_PATH, MODEL_PATH, LINE_SHAPEFILE_PATH, POLY_SHAPEFILE_PATH, TARGET_ALPHABETS, \
    IMAGE_PATHS, LINE_SHAPEFILE_PATHS, POLY_SHAPEFILE_PATHS, MULTITHREAD, PLOT_RESULT
import rasterio
from shapely.geometry import Polygon
from line_proc import affine_transform
from gis_io import write_line_shapefile, write_poly_shapefile
import cv2
import numpy as np
from osgeo import gdal

def main(image_path, line_shapefile_path, poly_shapefile_path, model_path=MODEL_PATH, target_alphabet=TARGET_ALPHABETS):
    # ----- Line Extraction
    # dot detection
    blob_points = get_dot_points(image_path)

    # get polygons' geojsons from raster image
    geojson_list = get_geojson_list(image_path, 0)  # masking black
    print(f'Number of found polygons : {len(geojson_list)}')

    # get all polygons on the image
    geoms = get_shapely_geom(geojson_list)

    lines = draw_dot_dashed_lines(blob_points, geoms)
    write_line_shapefile(list(lines), line_shapefile_path)

    # result_polys, dangles, cuts, invalids = polygonize_full(lines)
    # result_polys = list(result_polys)
    #
    # print(f'Created final polygons: {len(result_polys)}')
    # for poly in result_polys:
    #     plt.plot(*poly.exterior.xy)
    #
    # # ----- Label Extraction
    # ocr_result = recognize_texts(image_path, model_path, target_alphabet)
    # # plot_prediction_result(IMAGE_PATH, prediction_result)
    #
    # # affine transform for recognized labels
    # labels = []
    # transform = rasterio.open(image_path).transform
    # for word, bbox in ocr_result:
    #     geo_box_coords = []
    #     for px_coord in bbox:
    #         geo_coord = affine_transform(transform, px_coord)
    #         geo_box_coords.append(geo_coord)
    #     labels.append((word, Polygon(geo_box_coords).centroid))
    #
    # # ----- Writing Shapefile
    # write_poly_shapefile(result_polys, labels, poly_shapefile_path)

    if PLOT_RESULT:
        plt.figure(figsize=(21, 19))

        # plot GeoTiff
        ds = gdal.Open(image_path)
        gt = ds.GetGeoTransform()
        ulx, xres, _, uly, _, yres = gt
        lrx = ulx + (ds.RasterXSize * xres)
        lry = uly + (ds.RasterYSize * yres)
        extent = [ulx, lrx, lry, uly]
        array = ds.ReadAsArray()

        plt.imshow(array, extent=extent)

        # plot dots
        for p in blob_points:
            plt.plot(p.x, p.y, marker="o", markerfacecolor="red")


        # plot lines
        plot_line(lines)

        plt.show()


if __name__ == '__main__':
    if MULTITHREAD:
        # multi thread
        from multiprocessing import Pool

        with Pool() as pool:
            pool.starmap(main,
                         [(image_path, line_shapefile_path, poly_shapefile_path)
                          for (image_path, line_shapefile_path, poly_shapefile_path)
                          in zip(IMAGE_PATHS, LINE_SHAPEFILE_PATHS, POLY_SHAPEFILE_PATHS)])
    else:
        # single thread
        main(IMAGE_PATH, LINE_SHAPEFILE_PATH, POLY_SHAPEFILE_PATH)




