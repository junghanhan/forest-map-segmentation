import cv2 as cv
import numpy as np
import rasterio
from rasterio.plot import show

def readGeoTiff(imageFilePath):
    """Read a georeferenced tiff file and return its raster and profile data.

    Args:
        imageFilePath (str): Absolute or relative file path of the geotiff file
            to be read

    Returns:
        imageRGB (np uint8 array): 3-channel RGB raster image data of geotiff
        profile (rasterio profile object): Geo-referenced profile data
    """

    # read in RGB or intensity data from file
    dataset = rasterio.open(imageFilePath)
    # show(dataset)
    # print(dataset.bounds)
    # print(dataset.crs)
    # print(dataset.transform)

    # save the geodata and coordinates of the imported file
    profile = dataset.profile
    numChannels = dataset.read().shape[0]

    # read in RGB
    try:
        if numChannels == 3:
            imageR = dataset.read(1)
            imageG = dataset.read(2)
            imageB = dataset.read(3)

            # extract the dimensions of the image to determine the number of rows and
            # in the image RGB channels.
            shape = imageR.shape

            # create variables to store row, col, and channel sizes
            numRows = shape[0]
            numColumns = shape[1]

            # create empty black image (np.zeros makes it black)
            imageRGB = np.zeros((numRows, numColumns, numChannels), dtype='uint8')

            # replace the empty image channels with the imported RGB channels to create
            # a single RGB image
            imageRGB[:,:,0] = imageR
            imageRGB[:,:,1] = imageG
            imageRGB[:,:,2] = imageB

            return imageRGB, profile
        # read in grayscale
        elif numChannels == 1:
            imageGray = dataset.read(1)
            return imageGray, profile
        else:
            raise Exception('invalid number of channel of input file')
    except Exception as exc:
        print(exc)


img, profile = readGeoTiff('/home/student/PycharmProjects/forest-map-segmentation/res/map.gtiff')
#img, profile = readGeoTiff('/home/student/PycharmProjects/forest-map-segmentation/res/girl.jpg')

# RGB to BGR
if len(img.shape) == 3:
    img = img[:,:,::-1]

print(profile)

cv.imshow('img', img)
cv.waitKey(0)

