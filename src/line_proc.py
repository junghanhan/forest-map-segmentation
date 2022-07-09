import cv2
import rasterio
from rasterio.features import shapes
from centerline.geometry import Centerline
from shapely.geometry import shape, MultiLineString, box, LineString, Polygon, Point
from shapely.ops import nearest_points, linemerge
from shapely.strtree import STRtree
from itertools import combinations
import math
import networkx as nx
from settings import MIN_BRANCH_LEN

# get blobs as KeyPoint obj from image
# cv::KeyPoint attributes : pt, size, etc.
# image : cv2 image object; don't need to be a binarized image
# because cv2.SimpleBlobDetector itself binarzes the input image
def get_blobs (image, min_blob_area=20, max_blob_area=100):
  params = cv2.SimpleBlobDetector_Params()
  params.filterByArea = True
  params.minArea = min_blob_area
  params.maxArea = max_blob_area

  detector = cv2.SimpleBlobDetector_create(params)
  keypoints = detector.detect(image) # blobs

  return keypoints

# get geojson list extracted from raster image
# This function assumes the raster file is in grayscale format.
# raster_file_path:  GeoTiff file path
def get_geojson_list (raster_file_path, mask_value=None):
  with rasterio.Env():
    with rasterio.open(raster_file_path) as src:
      image = src.read(1) # Band 1
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
def get_shapely_geom (geojson_list):
  results = []

  for geojson in geojson_list:
    poly = shape(geojson['geometry'])
    results.append(poly)

  return results

# create a centerline object for a polygon
# returns a merged LineString if lines are continuous from tip to tail
# otherwise a MultiLineString consists of multiple continuous LineStrings is returned
# interpolation_distance: the higher value, the less branches from the centerline
def create_centerline(geom, interpolation_distance=0.00006):
  try:
    # simplify the polygon to prevent unnecessary branches
    geom = geom.buffer(0.00001, join_style=1)

    centerline_obj = Centerline(geom, 0.00006)

    # merge the centerlines into a single centerline
    multi_line = MultiLineString(centerline_obj)
    merged_line = linemerge(multi_line)
  except Exception as inst: # possibly TooFewRidgesError
    print(f'Exception {type(inst)} occurred.')
    print(inst)
    merged_line = None


  return merged_line

# get shapely geometries extracted from geojson
def get_shapely_geom (geojson_list):
  results = []

  for geojson in geojson_list:
    poly = shape(geojson['geometry'])
    results.append(poly)

  return results

def get_search_box(point, distance):
  return box(point.x - distance, point.y - distance, point.x + distance, point.y + distance)

def affine_transform (tif_file_path, pixel_coordinates):
  results = []
  with rasterio.open(tif_file_path) as src:
    for p_x, p_y in pixel_coordinates:
      g_x, g_y = rasterio.transform.xy(transform = src.transform,
                                       rows = p_y,
                                       cols = p_x)
      results.append((g_x,g_y))
  return results

def get_line_endpoints(mline, min_branch_len=MIN_BRANCH_LEN):
  """
  Get endpoints of MultiLineString except for intersecting points between its LineStrings
  Get endpoints of LineString

  mline: merged line, it can be MultiLineString or LineString
  """
  if not isinstance(mline, MultiLineString) and not isinstance(mline, LineString):
    raise TypeError(
      f'Inappropriate type: {type(mline)} for mline whereas a MultiLineString or LineString is expected')

  # plot_line(mline)
  if isinstance(mline, MultiLineString):
    mline = MultiLineString([line for line in mline if line.length > min_branch_len])
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

  return endpoints

def is_endpoint_inside(poly, sbox):
  """
  Check if at least one endpoint of the line derived from the polygon is within
  the search box.

  poly : target polygon
  sbox : search box
  """
  if not isinstance(poly, Polygon):
    raise TypeError(f'Inappropriate type: {type(poly)} for poly whereas a Polygon is expected')
  if not isinstance(sbox, Polygon):
    raise TypeError(f'Inappropriate type: {type(sbox)} for sbox whereas a Polygon is expected')

  line = create_centerline(poly)
  if line is not None:
    endpoints = get_line_endpoints(line)
    for p in endpoints:
      if sbox.contains(p):
        return True
  return False


def get_extrapolated_point(line, endpoint, dist=0.0001):
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
    raise TypeError(f'Inappropriate type: {type(mline)} for mline whereas a MultiLineString or LineString is expected')

  if isinstance(mline, MultiLineString):
    # unwrap Point to coord
    if isinstance(start_p, Point):
      start_p = start_p.coords[0]
    if isinstance(end_p, Point):
      end_p = end_p.coords[0]

    G = nx.Graph()

    for line in mline:
      prev_p = None
      for curr_p in line.coords:
        if prev_p is not None:
          add_edge_to_graph(G, prev_p, curr_p, 1)
        prev_p = curr_p

    result = LineString(nx.shortest_path(G, source=start_p, target=end_p))
  else:
    result = mline

  return result