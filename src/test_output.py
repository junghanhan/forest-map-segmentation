import output
import input
import os
from settings import RESOURCE_DIR

TEST_FILE_NAME = 'test_map'
RASTER_FILE_PATH = f'{RESOURCE_DIR}/{TEST_FILE_NAME}.gtiff'
VECTOR_FILE_PATH = f'{RESOURCE_DIR}/{TEST_FILE_NAME}.shp'
TIFF_FILE_PATH = f'{RESOURCE_DIR}/{TEST_FILE_NAME}.tiff'
DRIVER = 'ESRI Shapefile'

def test_polygonize():
    try:
        output.polygonize(RASTER_FILE_PATH, VECTOR_FILE_PATH, DRIVER, 0)
        if not os.path.exists(VECTOR_FILE_PATH):
            assert False

    except Exception as e:
        print(e)
        assert False


def test_plot_shape_file():
    try:
        output.polygonize(RASTER_FILE_PATH, VECTOR_FILE_PATH, DRIVER, 0)
        output.plot_shape_file(VECTOR_FILE_PATH)

    except Exception as e:
        print(e)
        assert False


def test_save_image_tiff():
    try:
        img, profile = input.read_geotiff(RASTER_FILE_PATH)
        output.save_image_tiff(img, profile, RESOURCE_DIR, TEST_FILE_NAME)

        if not os.path.exists(TIFF_FILE_PATH):
            assert False

    except Exception as e:
        print(e)
        assert False


