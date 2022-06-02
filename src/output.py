import os
import rasterio
from rasterio.features import shapes
import matplotlib.pyplot as plt
from pathlib import Path
import fiona
import shapefile as shp  # pyshp
import numpy as np

def plot_shape_file(shp_file_path):
    """Read a shapefile and plot using matplotlib

    :param
        shp_file_path: File path of the shapefile to read
    """

    sf = shp.Reader(shp_file_path)

    plt.figure()
    for shape in sf.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        y = [i[1] for i in shape.shape.points[:]]
        plt.plot(x, y)
    plt.show()


def polygonize(raster_file_path, vector_file_path, driver, mask_value):
    """Convert a geo-referenced raster file of a single-colour 3-channel image
        segment or a grayscale image to a vectorized polygon shapefile.
        Each polygon region of the shape file is also given a 'CATEGORY'
        attribute with a value equivalent to the raster file's name.
        This attribute is visible from within GIS software.

        This function emulates GDAL's gdal_polygonize.py and is adapted from:
        https://github.com/sgillies/rasterio/blob/master/examples/rasterio_polygonize.py

        Args:
            raster_file_path (str): File path of the TIFF raster image to be read
            vector_file_path (str): File path of the shapefile to be written
            driver (str): Name of the driver to use to create the shapefile
            mask_value (int): Integer raster value to omit from the shapefile
        """

    # get image file name and remove file extension
    raster_file_path_name_no_extension = Path(raster_file_path).stem
    vector_file_path_with_extension = os.path.basename(vector_file_path)

    with rasterio.Env():
        with rasterio.open(raster_file_path) as src:
            intensity = src.read(1)  # red in case of RGB

        if mask_value is not None:
            mask = intensity != mask_value
        else:
            mask = None

        # create geojson file structure with geometry data
        # a CATEGORY attribute is also written to each geometry regions
        # so that the map cover category can be stored with the shapefile
        results = (
            {'properties': {'CATEGORY': raster_file_path_name_no_extension}, 'geometry': s}
            for i, (s, v) in enumerate(
                shapes(intensity, mask=mask, transform=src.transform)))

        # use fiona to convert the geojson to a shapefile
        try:
            with fiona.open(
                    vector_file_path, 'w',
                    driver=driver,
                    crs=src.crs,
                    schema={'properties': [('CATEGORY', 'str')],
                            'geometry': 'Polygon'}) as dst:

                dst.writerecords(results)
        except:
            print("Could not create the shape file '"
                  + vector_file_path_with_extension + "'.\n")
            print('The file may be open in another program. Please close it '
                  + 'and re-start the tool.')

    return dst.name

def save_image_tiff(image_array, profile, destination_folder, file_name):
    """Save an image and profile data as a TIFF with a user-specified file name
    and to a user-specified directory. If the directory does not yet exist it
    will be created.

    Args:
        image_array (np uint8 array): 3-channel array with image data
        profile (rasterio profile object): Geo-referenced profile data
        destination_folder (str): Folder in which to save the image
        file_name (str): Name of the image file (without the .png extension)
    """

    # create new output directory if it does not yet exist
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)

    # determine shape of output image
    shape = image_array.shape
    height = int(shape[0])
    width = int(shape[1])

    # rearrange the array from height-width-channel to channel-height-width
    # which is the format expected by dst.write()
    tiff_array = np.moveaxis(image_array,-1,0)

    # open a gdal environment to use with rasterio
    with rasterio.Env():
        # modify the width and height to the newly cropped image so no errors
        # will occur
        profile.update(
            width=width,
            height=height,
            dtype=rasterio.ubyte,
            count=3,
            compress='lzw')

        # create a new file with the specific name
        save_path = destination_folder + os.path.sep + file_name + '.tiff'
        with rasterio.open(save_path,'w',**profile) as dst:
            dst.write(tiff_array)



