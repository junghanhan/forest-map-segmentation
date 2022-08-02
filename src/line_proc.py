from traceback import print_exc
import logging
import affine
import cv2
import numpy as np
import rasterio
from rasterio.features import shapes
from centerline.geometry import Centerline
from shapely.geometry import shape, MultiLineString, box, LineString, Polygon, Point, LinearRing, MultiPoint
from shapely.ops import linemerge, nearest_points
from itertools import combinations
import math
import networkx as nx
from settings import MIN_BRANCH_LEN, DASH_DOT_DIST, INTERPOLATION_DIST, CENTERLINE_BUFFER, DASH_SEARCH_BOX_W, \
    MIN_BLOB_AREA, MAX_BLOB_AREA, SIMPLIFY_TOLERANCE
import matplotlib.pyplot as plt
from typing import Tuple, List, Dict
from shapely.geometry.base import BaseGeometry


def get_blobs(image: np.ndarray, min_blob_area: float = MIN_BLOB_AREA, max_blob_area: float = MAX_BLOB_AREA) -> Tuple:
    """
    Get a tuple of the detected blobs

    :param image: image read as openCV2 Numpy ndarray
    :param min_blob_area: minimum area for detecting blobs
    :param max_blob_area: maximum area for detecting blobs
    :return: a tuple that contains the detected blobs as OpenCV2 KeyPoint objects
        KeyPoint objects contains the pixel coordinates of blobs
    """
    params = cv2.SimpleBlobDetector_Params()
    params.filterByArea = True
    params.minArea = min_blob_area
    params.maxArea = max_blob_area

    detector = cv2.SimpleBlobDetector_create(params)
    keypoints = detector.detect(image)  # blobs

    return keypoints


def get_dot_points(image_path: str) -> List[Point]:
    """
    Get dots as Shapely Point objects on input forest cover map.

    :param image_path: Input GeoTiff image file path
    :return: a list of Shapely Point objects representing dots on a forest cover map. The Point objects represents
        geospatial coordinates.
    """

    img = cv2.imread(image_path)

    # detecting blobs
    blob_keypoints = get_blobs(img)

    # converting the pixel coordinates of blob into geospatial coordinates
    blob_px_coords = []
    for blob in blob_keypoints:
        blob_px_coords.append(blob.pt)

    transform = rasterio.open(image_path).transform
    blob_geo_coords = [affine_transform(transform, blob_px_coord) for blob_px_coord in blob_px_coords]

    # make blobs shapely points by using their geo coordinates
    blob_points = [Point(blob) for blob in blob_geo_coords]

    return blob_points


def get_geojson_list(image_path: str, mask_value: int = None) -> List[Dict]:
    """
    Get geojson list extracted from raster image
    This function assumes the raster file is in grayscale format.

    :param image_path: Input GeoTiff image file path
    :param mask_value: grayscale intensity value of pixels that will be extracted as geojson
        e.g., 0: black pixels will be extracted as geojson
    :return: a list of geojson values represented as dictionaries
    """

    with rasterio.Env():
        with rasterio.open(image_path) as src:
            image = src.read(1)  # Band 1
            ret, inv_bin_image = cv2.threshold(image, 120, 255, cv2.THRESH_BINARY_INV)

            if mask_value is not None:
                mask = inv_bin_image != mask_value
            else:
                mask = None

            results = (
                {'properties': {'raster_val': v}, 'geometry': s}
                for i, (s, v)
                in enumerate(
                shapes(inv_bin_image, mask=mask, transform=src.transform)))
    return list(results)


def get_shapely_geom(geojson_list: List[Dict]) -> List[BaseGeometry]:
    """
    Get shapely geometries extracted from geojson

    :param geojson_list: a list of geojson values represented as dictionaries
    :return: a list of Shapely Polygon objects created by input geojson list
    """
    results = []

    for geojson in geojson_list:
        poly = shape(geojson['geometry'])
        results.append(poly)

    return results


def create_centerline(poly: Polygon, poly_line_dict: Dict[int, BaseGeometry] = None, buffer: float = CENTERLINE_BUFFER,
                      interpolation_distance: float = INTERPOLATION_DIST, tolerance: float = SIMPLIFY_TOLERANCE) -> BaseGeometry:
    """
    Create a centerline for a polygon

    :param poly: target polygon that will be converted into line
    :param poly_line_dict: Dictionary that contains polygon and its created line. {id(polygon):line}.
            This is to prevent redundant calls for the function that creates centerline of polygon.
    :param buffer: value to expand the target polygon. This value is to simplify the polygon before the centerline is
        created from the polygon. The higher the value, the bigger the polygon.
    :param interpolation_distance: the higher value, the fewer branches from the centerline
    :param tolerance: value to simplify the created centerline. The higher value, the simpler centerline.
    :return: Returns a merged LineString if lines are continuous from tip to tail.
        Otherwise, MultiLineString will be returned.
    """

    if not isinstance(poly, Polygon):
        raise TypeError(
            f'Inappropriate type: {type(poly)} for poly whereas a Polygon is expected')

    if id(poly) in poly_line_dict:
        result_line = poly_line_dict[id(poly)]
    else:
        try:
            # simplify the polygon to prevent unnecessary branches
            geom = poly.buffer(buffer, join_style=1)

            centerline_obj = Centerline(geom, interpolation_distance)

            # merge the centerlines into a single centerline
            line = MultiLineString(centerline_obj)
            result_line = linemerge(line)
        except Exception as err:  # possibly TooFewRidgesError
            logging.error(err, exc_info=True)
            print_exc()
            print(err)
            result_line = None

        result_line = result_line.simplify(tolerance)
        poly_line_dict[id(poly)] = result_line

    return result_line


def get_search_box(point: Point, width: float, height: float) -> Polygon:
    """
    Get box polygon around the point

    :param point: Shapely Point object
    :param width: width of the returned box
    :param height: height of the returned box
    :return: a square-shaped Shapely Polygon object
    """
    if not isinstance(point, Point):
        raise TypeError(
            f'Inappropriate type: {type(point)} for point whereas a Point is expected')

    return box(point.x - width / 2, point.y - height / 2, point.x + width / 2, point.y + height / 2)


def affine_transform(transform: affine.Affine, pixel_coordinate: Tuple[float, float]) -> Tuple[float, float]:
    """
    Transform a pixel coordinate to a geospatial coordinate based on input Affine

    :param transform: Transform suitable for input to AffineTransformer, GCPTransformer, or RPCTransformer.
    :param pixel_coordinate: pixel coordinate on the image; can be a tuple or an array representing a coordinate of a point
    :return: a tuple of a geospatial coordinate
    """

    g_x, g_y = rasterio.transform.xy(transform = transform,
                                       rows = pixel_coordinate[1],
                                       cols = pixel_coordinate[0])

    return (g_x, g_y)


def get_line_endpoints(line: BaseGeometry, line_ep_dict: Dict[int, List] = None) -> List[Point]:
    """
    Get endpoints of MultiLineString except for intersecting points between its LineStrings
    Get endpoints of LineString

    :param line: a Shapely MultiLineString or LineString object
    :param line_ep_dict: Dictionary that contains line and its endpoints. {id(line):[endpoints]}.
            This is to prevent redundant calls of get_line_endpoints.
    :return: a list of Shapely Point objects
    """

    if not isinstance(line, MultiLineString) and not isinstance(line, LineString):
        raise TypeError(
            f'Inappropriate type: {type(line)} for line whereas a MultiLineString or LineString is expected')

    if id(line) in line_ep_dict:
        endpoints = line_ep_dict[id(line)]
    else:
        if isinstance(line, MultiLineString):
            intersect_points = []
            for line1, line2 in combinations([l for l in line], 2):
                if line1.intersects(line2):
                    intersections = line1.intersection(line2)
                    if isinstance(intersections, list):
                        intersect_points.extend()
                    else:
                        intersect_points.append(intersections)
            endpoints = [p for p in line.boundary if p not in intersect_points]
        else:  # LineString
            endpoints = list(line.boundary)

        line_ep_dict[id(line)] = endpoints

    return endpoints


def get_common_endpoints(endpoints1: List[Point], endpoints2: List[Point]) -> List[Point]:
    """
    Get common Shapely Points between two lists of Shapely Points

    :param endpoints1: a list of Shapely Point objects
    :param endpoints2: a list of Shapely Point objects
    :return: a list of Point objects that are within both input lists
    """
    endpoints_dict1 = {}
    endpoints_dict2 = {}
    for p in endpoints1:
        endpoints_dict1[id(p)] = p
    for p in endpoints2:
        endpoints_dict2[id(p)] = p

    common_endpoints = {k: v for k, v in endpoints_dict1.items() if k in endpoints_dict2}

    return list(common_endpoints.values())


def get_extrapolated_point(line: BaseGeometry, endpoint: Point, dist: float = DASH_DOT_DIST) -> Point:
    """
    Creates an extrapolated point from a line.

    :param line: a Shapely LineString or MultiLineString object that will be used to extrapolate
    :param endpoint: a Shapely Point object that is an endpoint of the input line.
        This endpoint will be the previous point of the extrapolated point.
    :param dist: the distance from the endpoint to the extrapolated point
    :return: a Shapely Point object. The extrapolated point.
    """

    def get_prev_point(line: BaseGeometry, endpoint: Point) -> Point:
        """
        Get the previous point of endpoint in the line

        :param line: a Shapely LineString object that contains endpoint
        :param endpoint: a Shapely Point object that is an endpoint of the input line
        :return: a Shapely Point object. It is the point on the line next to the input endpoint.
        """

        l_coords = list(line.coords)
        ep_index = l_coords.index(endpoint.coords[0])
        prev_point = l_coords[ep_index + 1] if ep_index == 0 else l_coords[ep_index - 1]
        return prev_point

    if isinstance(line, LineString):
        if not endpoint.coords[0] in line.coords:
            raise Exception(f'Endpoint {endpoint} is not included in the line {line}')
        prev_point = get_prev_point(line, endpoint)
    elif isinstance(line, MultiLineString):
        is_p_in_line = False
        for l in line:
            if endpoint.coords[0] in l.coords:
                prev_point = get_prev_point(l, endpoint)
                is_p_in_line = True
        if not is_p_in_line:
            raise Exception(f'Endpoint {endpoint} is not included in the line {line}')
    else:
        raise TypeError(
            f'Inappropriate type: {type(line)} for line whereas a MultiLineString or LineString is expected')

    if not isinstance(endpoint, Point):
        raise TypeError(f'Inappropriate type: {type(endpoint)} for endpoint whereas a Point is expected')

    diff = (endpoint.x - prev_point[0], endpoint.y - prev_point[1])
    norm = math.sqrt(diff[0] ** 2 + diff[1] ** 2)
    direction = (diff[0] / norm, diff[1] / norm)
    result = Point(endpoint.x + direction[0] * dist, endpoint.y + direction[1] * dist)
    return result


def get_shortest_path_line(path: BaseGeometry, start_p: Point, end_p: Point) -> LineString:
    """
    Get a line that connects between two points by using the input line as available paths.

    :param path: a Shapely LineString or MultiLineString object that will be used as available paths
    :param start_p: a Shapely Point object. This Point should be included on the input line.
    :param end_p: a Shapely Point object. This Point should be included on the input line.
    :return: a Shapely LineString object representing the shortest path between two input points
    """

    if not isinstance(path, MultiLineString) and not isinstance(path, LineString):
        raise TypeError(
            f'Inappropriate type: {type(path)} for path whereas a MultiLineString or LineString is expected')

    if isinstance(path, MultiLineString):
        # unwrap Point to coord
        if isinstance(start_p, Point):
            start_p = start_p.coords[0]
        if isinstance(end_p, Point):
            end_p = end_p.coords[0]

        try:
            G = nx.Graph()

            for l in path:
                prev_p = None
                for curr_p in l.coords:
                    if prev_p is not None:
                        G.add_edge(prev_p, curr_p, weight=1)
                    prev_p = curr_p

            # a list of coordinates representing the shortest path
            shortest_path = nx.shortest_path(G, source=start_p, target=end_p)
            result = LineString(shortest_path) if len(shortest_path) > 1 else None
        except Exception as err:
            logging.error(err, exc_info=True)
            print_exc()
            print(err)
            result = None
    else:
        result = path

    return result


def find_nearest_geom(target: BaseGeometry, geom_list: List[BaseGeometry]) -> BaseGeometry:
    """
    Find the nearest geometric object to the target geometric object from list

    :param target: Shapely geometry object
    :param geom_list: a list of Shapely geometry objects
    :return: the nearest Shapely geometric object
    """

    if len(geom_list) == 0:
        raise Exception("geom_list is empty")

    min_dist = float('inf')
    for geom in geom_list:
        dist = target.distance(geom)
        if dist < min_dist:
            min_dist = dist
            min_dist_geom = geom

    return min_dist_geom


def get_image_bbox(image_path: str, offset: float = 0) -> LinearRing:
    """
    Get image bounding box as a Shapely LinearRing object

    :param image_path: Input image path
    :param offset: pixel offset value that will be used to reduce the size of bounding box close to center
    :return: a Shapely LinearRing object representing image bounding box
    """

    left_x, lower_y, right_x, upper_y = rasterio.open(image_path).bounds
    if offset * 2 > right_x - left_x or offset * 2 > upper_y - lower_y:
        raise ValueError("offset value is too big for the image.")

    left_x += offset
    lower_y += offset
    right_x -= offset
    upper_y -= offset

    image_bbox = LinearRing([(left_x, upper_y), (right_x, upper_y), (right_x, lower_y), (left_x, lower_y)])

    return image_bbox


def get_close_points(line1: BaseGeometry, line2: BaseGeometry) -> List[Point]:
    """
    Get points on line1 which are close to the actual intersection points between line1 and line2

    :param line1: Shapely LineString or MultiLineString object
    :param line2: Shapely LineString or MultiLineString object such as LineString
    :return: a list of Shapely Point objects
    """

    if not isinstance(line1, MultiLineString) and not isinstance(line1, LineString):
        raise TypeError(
            f'Inappropriate type: {type(line1)} for line1 whereas a MultiLineString or LineString is expected')
    if not isinstance(line2, MultiLineString) and not isinstance(line2, LineString):
        raise TypeError(
            f'Inappropriate type: {type(line2)} for line2 whereas a MultiLineString or LineString is expected')

    close_points = []
    intersect_points = line2.intersection(line1)
    if not intersect_points.is_empty:
        if isinstance(intersect_points, Point):
            intersect_points = [intersect_points]
        for intersect_p in intersect_points:
            all_points_on_line = get_points_on_line(line1)
            target_p = nearest_points(all_points_on_line, intersect_p)[0]
            # to prevent duplicate points
            if target_p not in close_points:
                close_points.append(target_p)

    return close_points


def filter_geoms(center_point: Point, geoms: List[BaseGeometry], radius: float) -> List[BaseGeometry]:
    """
    Get Shapely geometry objects that are within a circle with a certain radius centering a point.

    :param center_point: a Shapely Point object that is used as a center of a filtering circle
    :param geoms: a list of Shapely geometry objects
    :param radius: the radius of the search circle that determines whether a point is near this Dot object or not
    :return: a list of filtered Shapely geometry objects
    """

    if not isinstance(center_point, Point):
        TypeError(
            f'Inappropriate type: {type(center_point)} for center_point whereas a Point is expected')

    filtered_geoms = geoms.copy()
    ep_scircle = center_point.buffer(radius)  # endpoint search circle
    # plt.plot(*ep_scircle.exterior.xy)
    filtered_geoms = [fg for fg in filtered_geoms
                          if not ep_scircle.contains(fg)]

    return filtered_geoms


def get_points_on_line(line: BaseGeometry) -> MultiPoint:
    """
    Get all points on a line or lines as a list.

    :param line: a Shapely MultiLineString or LineString object
    :return: a list of Shapely Point objects (MultiPoint)
    """

    if not isinstance(line, MultiLineString) and not isinstance(line, LineString):
        raise TypeError(
            f'Inappropriate type: {type(line)} for line whereas a MultiLineString or LineString is expected')

    lines = [line] if isinstance(line, LineString) else line
    all_points_on_line = MultiPoint([p for line in lines for p in line.coords])

    return all_points_on_line


def plot_line(line: BaseGeometry) -> None:
    """
    This is a function to facilitate plotting lines using matplotlib.

    :param line: a Shapely LineString or MultiLineString object that will be plotted
    """

    # plot the final line for dash
    if line.type == "MultiLineString":
        multilinestr = line
        for linestr in multilinestr:
            x, y = linestr.xy
            plt.plot(x, y)
    else:
        x, y = line.xy
        plt.plot(x, y)
