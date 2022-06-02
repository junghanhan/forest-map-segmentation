import cv2 as cv
import input
from settings import RESOURCE_DIR

def test_read_geo_tiff():
    try:
        img, profile = input.read_geotiff(f'{RESOURCE_DIR}/test_map.gtiff')

        # RGB to BGR
        if len(img.shape) == 3:
            img = img[:,:,::-1]

        assert profile['driver'] == 'GTiff'

        cv.imshow('img', img)
        # cv.waitKey(0)
    except Exception as e:
        print(e)
        assert False