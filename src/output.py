def polygonize(rasterFile, vectorFile, driver, maskValue):
    """Convert a geo-referenced raster file of a single-colour 3-channel image
    segment to a vectorized polygon shapefile. Each polygon region of the shape
    file is also given a 'CATEGORY' attribute with a value equivalent to the
    raster file's name. This attribute is visible from within GIS software.

    This function emulates GDAL's gdal_polygonize.py and is adapted from:
    https://github.com/sgillies/rasterio/blob/master/examples/rasterio_polygonize.py

    Args:
        rasterFile (str): File path of the TIFF raster image to be read
        vectorFile (str): File path of the shapefile to be written
        driver (str): Name of the driver to use to create the shapefile
        maskValue (int): Integer raster value to omit from the shapefile
    """
    # get image file name and remove file extension
    rasterFileNameNoExtension = Path(rasterFile).stem
    vectorFileWithExtension = os.path.basename(vectorFile)

    with rasterio.Env():
        with rasterio.open(rasterFile) as src:
            imageR = src.read(1)

        if maskValue is not None:
            mask = imageR != maskValue
        else:
            mask = None

        # create geojson file structure with geometry data
        # a CATEGORY attribute is also written to each geometry regions
        # so that the map cover category can be stored with the shapefile
        results = (
            {'properties': {'CATEGORY': rasterFileNameNoExtension}, 'geometry': s}
            for i, (s, v) in enumerate(
                shapes(imageR, mask=mask, transform=src.transform)))

        # use fiona to convert the geojson to a shapefile
        try:
            with fiona.open(
                    vectorFile, 'w',
                    driver=driver,
                    crs=src.crs,
                    schema={'properties': [('CATEGORY', 'str')],
                            'geometry': 'Polygon'}) as dst:

                    dst.writerecords(results)
        except:
            print("Could not create the shape file '"
                + vectorFileWithExtension + "'.\n")
            print('The file may be open in another program. Please close it '
                  + 'and re-start the tool.')

    return dst.name
