import numpy as np
import rasterio

def read_geotiff(image_file_path):
    """Read a georeferenced tiff file and return its raster and profile data.

    Args:
        image_file_path (str): Absolute or relative file path of the geotiff file
            to be read

    Returns:
        image (np uint8 array): 3-channel RGB raster image data or grayscale intensity
            of geotiff
        profile (rasterio profile object): Geo-referenced profile data
    """

    # read in RGB or intensity data from file
    dataset = rasterio.open(image_file_path)
    # show(dataset)
    # print(dataset.bounds)
    # print(dataset.crs)
    # print(dataset.transform)

    # save the geodata and coordinates of the imported file
    profile = dataset.profile
    num_channels = dataset.read().shape[0]

    # read in RGB
    try:
        if num_channels == 3:
            image_r = dataset.read(1)
            image_g = dataset.read(2)
            image_b = dataset.read(3)

            # extract the dimensions of the image to determine the number of rows and
            # in the image RGB channels.
            shape = image_r.shape

            # create variables to store row, col, and channel sizes
            num_rows = shape[0]
            num_columns = shape[1]

            # create empty black image (np.zeros makes it black)
            image = np.zeros((num_rows, num_columns, num_channels), dtype='uint8')

            # replace the empty image channels with the imported RGB channels to create
            # a single RGB image
            image[:,:,0] = image_r
            image[:,:,1] = image_g
            image[:,:,2] = image_b

            return image, profile
        # read in grayscale
        elif num_channels == 1:
            image = dataset.read(1)
            return image, profile
        else:
            raise Exception('invalid number of channel of input file')
    except Exception as exc:
        print(exc)