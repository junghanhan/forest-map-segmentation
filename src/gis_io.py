import geopandas as gpd
from fiona.crs import from_epsg
from shapely.geometry import Point
from shapely.geometry.base import BaseGeometry
from shapely.strtree import STRtree
from typing import List, Tuple
from settings import EPSG_CODE


def write_line_shapefile(lines: List[BaseGeometry], shapefile_path: str, epsg_code: int = EPSG_CODE) -> None:
    """
    Make shapefile features with lines
    and write all the features in a shapefile

    :param lines: a list of Shapely LineString objects
    :param shapefile_path: a string that represents the path that the shapefile will
        be written
    :param epsg_code: 4-5 digit numbers that represent CRS definitions
    """

    df = gpd.GeoDataFrame()

    # make feature records for all polygons
    # create record
    df['geometry'] = lines
    df.crs = from_epsg(epsg_code)
    df.to_file(shapefile_path)


def write_poly_shapefile(polygons: List[BaseGeometry], labels: List[Tuple[str, Point]], shapefile_path: str,
                         epsg_code: int = EPSG_CODE) -> None:
    """
    Make shapefile features by associating polygons and labels
    and write all the features in a shapefile

    :param polygons: a list of Shapely Polygon objects
    :param labels: a list of tuples; the result of label extraction (word, Point)
        word is a string for a label and Point is a shapely Point object that contains
        the coordinate of the label
    :param shapefile_path: a string that represents the path that the shapefile will
        be written
    :param epsg_code: 4-5 digit numbers that represent CRS definitions
    """

    result_polys_tree = STRtree(polygons)

    poly_record_dict = {}  # id(poly): (feature_id, [(label_text, label_center_p), ]
    records = gpd.GeoDataFrame()

    # make feature records for all polygons
    for poly in polygons:
        # create record
        feature_id = str(len(poly_record_dict) + 1)
        poly_record_dict[id(poly)] = (feature_id, [])
        records.loc[feature_id, 'geometry'] = poly

    # associate labels with corresponding polygons
    for label in labels:
        label_text = label[0]
        label_center_p = label[1]
        assoc_poly = result_polys_tree.query(label_center_p)

        # if the searched polygon contains the centroid of label
        # create a record or update the existing record
        for poly in assoc_poly:
            if label_center_p.within(poly):
                # update record
                feature_id, feature_labels = poly_record_dict[id(poly)]
                feature_labels.append((label_text, label_center_p))
                records.loc[feature_id, f'LBL{len(feature_labels)}_TEXT'] = label_text
                records.loc[feature_id, f'LBL{len(feature_labels)}_COORD'] = str(label_center_p.coords[0])

    records.crs = from_epsg(epsg_code)
    records.to_file(shapefile_path)

