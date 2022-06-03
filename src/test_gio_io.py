import cv2 as cv
import gis_io
from settings import ROOT_DIR, DRIVER
import os

TEST_FILE_NAME = 'test_map'
TEST_RESOURCE_DIR = f'{ROOT_DIR}/test_res'
TEST_RASTER_FILE_PATH = f'{TEST_RESOURCE_DIR}/{TEST_FILE_NAME}.gtiff'
TEST_VECTOR_FILE_PATH = f'{TEST_RESOURCE_DIR}/{TEST_FILE_NAME}.shp'
TEST_TIFF_FILE_PATH = f'{TEST_RESOURCE_DIR}/{TEST_FILE_NAME}.tiff'

def test_read_geo_tiff():
    try:
        img, profile = gis_io.read_geotiff(TEST_RASTER_FILE_PATH)

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
        gis_io.polygonize(TEST_RASTER_FILE_PATH, TEST_VECTOR_FILE_PATH, DRIVER, 0)
        assert os.path.exists(TEST_VECTOR_FILE_PATH)

    except Exception as e:
        print(e)
        assert False


def test_plot_shape_file():
    try:
        gis_io.polygonize(TEST_RASTER_FILE_PATH, TEST_VECTOR_FILE_PATH, DRIVER, 0)
        gis_io.plot_shape_file(TEST_VECTOR_FILE_PATH)

    except Exception as e:
        print(e)
        assert False


def test_save_image_tiff():
    try:
        img, profile = gis_io.read_geotiff(TEST_RASTER_FILE_PATH)
        gis_io.save_image_tiff(img, profile, TEST_RESOURCE_DIR, TEST_FILE_NAME)

        assert os.path.exists(TEST_TIFF_FILE_PATH)

    except Exception as e:
        print(e)
        assert False


