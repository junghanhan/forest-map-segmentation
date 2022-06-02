import cv2 as cv
import gis_io
from settings import RESOURCE_DIR
import os

TEST_FILE_NAME = 'test_map'
RASTER_FILE_PATH = f'{RESOURCE_DIR}/{TEST_FILE_NAME}.gtiff'
VECTOR_FILE_PATH = f'{RESOURCE_DIR}/{TEST_FILE_NAME}.shp'
TIFF_FILE_PATH = f'{RESOURCE_DIR}/{TEST_FILE_NAME}.tiff'
DRIVER = 'ESRI Shapefile'

def test_read_geo_tiff():
    try:
        img, profile = gis_io.read_geotiff(f'{RESOURCE_DIR}/test_map.gtiff')

        # RGB to BGR
        if len(img.shape) == 3:
            img = img[:,:,::-1]

        assert profile['driver'] == 'GTiff'

        cv.imshow('img', img)
        # cv.waitKey(0)
    except Exception as e:
        print(e)
        assert False


def test_polygonize():
    try:
        gis_io.polygonize(RASTER_FILE_PATH, VECTOR_FILE_PATH, DRIVER, 0)
        if not os.path.exists(VECTOR_FILE_PATH):
            assert False

    except Exception as e:
        print(e)
        assert False


def test_plot_shape_file():
    try:
        gis_io.polygonize(RASTER_FILE_PATH, VECTOR_FILE_PATH, DRIVER, 0)
        gis_io.plot_shape_file(VECTOR_FILE_PATH)

    except Exception as e:
        print(e)
        assert False


def test_save_image_tiff():
    try:
        img, profile = gis_io.read_geotiff(RASTER_FILE_PATH)
        gis_io.save_image_tiff(img, profile, RESOURCE_DIR, TEST_FILE_NAME)

        if not os.path.exists(TIFF_FILE_PATH):
            assert False

    except Exception as e:
        print(e)
        assert False


