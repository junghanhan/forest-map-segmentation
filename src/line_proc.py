from traceback import print_exc
import logging
import cv2
import rasterio
from rasterio.features import shapes
from centerline.geometry import Centerline
from shapely.geometry import shape, MultiLineString, box, LineString, Polygon, Point, LinearRing, MultiPoint
from shapely.ops import linemerge, nearest_points
from itertools import combinations
import math
import networkx as nx
from settings import MIN_BRANCH_LEN, DASH_DOT_DIST, INTERPOLATION_DIST, CENTERLINE_BUFFER
import matplotlib.pyplot as plt


line_endpoints = {} # id(line): [endpoints], line can be LineString, MultiLineString


# get blobs as KeyPoint obj from image
# cv::KeyPoint attributes : pt, size, etc.
# image : cv2 image object; don't need to be a binarized image
# because cv2.SimpleBlobDetector itself binarzes the input image
def get_blobs(image, min_blob_area=20, max_blob_area=100):
    params = cv2.SimpleBlobDetector_Params()
    params.filterByArea = True
    params.minArea = min_blob_area
    params.maxArea = max_blob_area

    detector = cv2.SimpleBlobDetector_create(params)
    keypoints = detector.detect(image)  # blobs

    return keypoints

def get_dot_points(image_path):
    """
    Get dots as Shapely Point objects on input forest cover map
    :param image_path: Input image path
    :return: Shapely Point objects representing dots on a forest cover map
    """

    img = cv2.imread(image_path)

    # detecting blobs
    blob_keypoints = get_blobs(img, max_blob_area=60)

    # converting the pixel coordinates of blob into geospatial coordinates
    blob_px_coords = []
    for blob in blob_keypoints:
        blob_px_coords.append(blob.pt)

    transform = rasterio.open(image_path).transform
    blob_geo_coords = [affine_transform(transform, blob_px_coord) for blob_px_coord in blob_px_coords]

    # make blobs shapely points by using their geo coordinates
    blob_points = [Point(blob) for blob in blob_geo_coords]

    return blob_points


# get geojson list extracted from raster image
# This function assumes the raster file is in grayscale format.
# raster_file_path:  GeoTiff file path
def get_geojson_list(raster_file_path, mask_value=None):
    with rasterio.Env():
        with rasterio.open(raster_file_path) as src:
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


# get shapely geometries extracted from geojson
def get_shapely_geom(geojson_list):
    results = []

    for geojson in geojson_list:
        poly = shape(geojson['geometry'])
        results.append(poly)

    return results


def create_centerline(geom, buffer=CENTERLINE_BUFFER, interpolation_distance=INTERPOLATION_DIST, min_branch_len=MIN_BRANCH_LEN):
    """
    Create a centerline for a polygon

    :param geom: target polygon that will be converted into line
    :param buffer: value to simplify polygon
    :param interpolation_distance: the higher value, the less branches from the centerline
    :param min_branch_len: value to prune unnecessary short branches when centerline is created
    :return: Returns a merged LineString if lines are continuous from tip to tail.
        Otherwise, MultiLineString will be returned.
    """

    try:
        # simplify the polygon to prevent unnecessary branches
        geom = geom.buffer(buffer, join_style=1)

        centerline_obj = Centerline(geom, interpolation_distance)

        # merge the centerlines into a single centerline
        mline = MultiLineString(centerline_obj)
        result_line = linemerge(mline)
    except Exception as err:  # possibly TooFewRidgesError
        logging.error(err, exc_info=True)
        print_exc()
        print(err)
        result_line = None

    return result_line


def get_search_box(point, distance):
    return box(point.x - distance, point.y - distance, point.x + distance, point.y + distance)


def affine_transform (transform, pixel_coordinate):
    """
    transform: (Affine or sequence of GroundControlPoint or RPC)
    – Transform suitable for input to AffineTransformer, GCPTransformer, or RPCTransformer.
    pixel_coordinate: pixel coordinate on the image; can be a tuple or an array representing a coordinate of a point
    """
    g_x, g_y = rasterio.transform.xy(transform = transform,
                                       rows = pixel_coordinate[1],
                                       cols = pixel_coordinate[0])
    return (g_x, g_y)


def get_line_endpoints(mline):
    """
    Get endpoints of MultiLineString except for intersecting points between its LineStrings
    Get endpoints of LineString

    mline: merged line, it can be MultiLineString or LineString
    """

    if not isinstance(mline, MultiLineString) and not isinstance(mline, LineString):
        raise TypeError(
            f'Inappropriate type: {type(mline)} for mline whereas a MultiLineString or LineString is expected')

    global line_endpoints  # TODO: if using global is more efficient, then other dictionaries should also be global ?

    if id(mline) in line_endpoints:
        endpoints = line_endpoints[id(mline)]
    else:
        #plot_line(mline)
        if isinstance(mline, MultiLineString):
            intersect_points = []
            for line1, line2 in combinations([line for line in mline], 2):
                if line1.intersects(line2):
                    intersections = line1.intersection(line2)
                    if isinstance(intersections, list):
                        intersect_points.extend()
                    else:
                        intersect_points.append(intersections)
            endpoints = [p for p in mline.boundary if p not in intersect_points]
        else:  # LineString
            endpoints = list(mline.boundary)

        line_endpoints[id(mline)] = endpoints

    return endpoints

# deprecated
def is_endpoint_inside(mline, sbox):
    """
    Check if at least one endpoint of the line derived from the polygon is within
    the search box.

    mline : target LineString or MultiLineString
    sbox : search box
    """

    if not isinstance(mline, MultiLineString) and not isinstance(mline, LineString):
        raise TypeError(
            f'Inappropriate type: {type(mline)} for mline whereas a MultiLineString or LineString is expected')
    if not isinstance(sbox, Polygon):
        raise TypeError(f'Inappropriate type: {type(sbox)} for sbox whereas a Polygon is expected')

    endpoints = get_line_endpoints(mline)
    for p in endpoints:
        if sbox.contains(p):
            return True
    return False


def get_extrapolated_point(line, endpoint, dist=DASH_DOT_DIST):
    """
    Creates a point extrapoled in line's direction to endpoint
    line: LineString or MultiLineString that will be used to extrapolate
    endpoint: Shapely Point object that is the endpoint of the line
    dist: distance from end point to the extrapolated point
    """

    def get_prev_point(line, endpoint):
        """
        Get the previous point of endpoint in the line
        line: LineString that contains endpoint
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


def get_path_line(mline, start_p, end_p):
    """
    Get the LineString that connects start point to end point in a MultiLineString or LineString
    start_p and end_p should be contained in mline
    """

    def add_edge_to_graph(G, e1, e2, w):
        G.add_edge(e1, e2, weight=w)

    if not isinstance(mline, MultiLineString) and not isinstance(mline, LineString):
        raise TypeError(
            f'Inappropriate type: {type(mline)} for mline whereas a MultiLineString or LineString is expected')

    if isinstance(mline, MultiLineString):
        # unwrap Point to coord
        if isinstance(start_p, Point):
            start_p = start_p.coords[0]
        if isinstance(end_p, Point):
            end_p = end_p.coords[0]

        try:
            G = nx.Graph()

            for line in mline:
                prev_p = None
                for curr_p in line.coords:
                    if prev_p is not None:
                        add_edge_to_graph(G, prev_p, curr_p, 1)
                    prev_p = curr_p

            # a list of coordinates representing the shortest path
            path = nx.shortest_path(G, source=start_p, target=end_p)
            result = LineString(path) if len(path) > 1 else None
        except Exception as err:
            logging.error(err, exc_info=True)
            print_exc()
            print(err)
            result = None

    else:
        result = mline

    return result


def find_nearest_geom(target, geom_list):
    """
    Find the nearest geometric object to the target geometric object from list

    :param target: Shapely geometric object
    :param geom_list: a list of Shapely geometric objects
    :return: the nearest Shapely geometric object
    """

    min_dist = float('inf')
    for geom in geom_list:
        dist = target.distance(geom)
        if dist < min_dist:
            min_dist = dist
            min_dist_geom = geom

    return min_dist_geom


def get_image_bbox(image_path, offset=0):
    """
    Get image bounding box as a Shapely LinearRing object

    :param image_path: Input image path
    :param offset: pixel offset value that will be used to reduce the size of bounding box close to center
    :return: a Shapely LinearRing object
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


def get_close_points(line1, line2):
    """
    Get points on geom1 which are close to the actual intersection points between geom1 and geom2
    TODO: this is a workaround of shapely.ops.split since shapely.ops.split does not work due to precision issues.

    :param line1: Shapely LineString or MultiLineString object
    :param line2: Shapely LineString or MultiLineString object such as LineString
    :return: a list of Shapely Point objects
    """

    if not isinstance(line1, MultiLineString) and not isinstance(line1, LineString):
        raise TypeError(
            f'Inappropriate type: {type(line1)} for mline whereas a MultiLineString or LineString is expected')
    if not isinstance(line2, MultiLineString) and not isinstance(line2, LineString):
        raise TypeError(
            f'Inappropriate type: {type(line2)} for mline whereas a MultiLineString or LineString is expected')

    close_points = []
    intersect_points = line2.intersection(line1)
    if not intersect_points.is_empty:
        if isinstance(intersect_points, Point):
            intersect_points = [intersect_points]
        for intersect_p in intersect_points:
            lines = [line1] if isinstance(line1, LineString) else line1
            all_points_in_line = MultiPoint([p for line in lines for p in line.coords])
            target_p = nearest_points(all_points_in_line, intersect_p)[0]
            # to prevent duplicate points
            if target_p not in close_points:
                close_points.append(target_p)

    return close_points


# matplotlib line plot helper function
def plot_line(line):
    # plot the final line for dash
    if line.type == "MultiLineString":
        multilinestr = line
        for linestr in multilinestr:
            x, y = linestr.xy
            plt.plot(x, y)
    else:
        x, y = line.xy
        plt.plot(x, y)
