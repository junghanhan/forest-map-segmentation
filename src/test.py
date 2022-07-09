######## endpoint test
# import matplotlib.pyplot as plt
# from shapely.geometry import MultiLineString, Point, LineString
# from line_proc import get_line_endpoints
#
# mline = MultiLineString([((0, 0), (1, 1), (1.5, 1)), ((0, 2), (1, 1.5), (1.5, 1), (2, 0))])
# line = LineString([(0, 0), (1, 1), (1.5, 1)])
# endpoints = get_line_endpoints(line)
#
# for p in endpoints:
#   print(p)
#   plt.plot(*p.xy, marker="o")
#
# # plotting the line
# plt.plot(*line.xy)
# plt.show()

######### Dot, Dash class test
##### check whether the nearby polygons are dashes #####
# import matplotlib.pyplot as plt
# from shapely.geometry import Point, LineString
# from shapely.strtree import STRtree
# from line_proc import create_centerline
# from dot_dash import Dot, Dash
#
# # find the polygons around each dot
# dot = Dot(Point(0.005, 0.005))
# plt.plot(dot.point.x, dot.point.y, color='green', marker='o')
#
# poly1 = LineString([(0.00425, 0.005), (0.00475, 0.005)]).buffer(0.00002)
# # poly2 = LineString([(0.00525,0.005), (0.00575, 0.005)]).buffer(0.00002)
# poly2 = LineString([(0.00525, 0.00425), (0.00525, 0.00575)]).buffer(0.00002)
# poly3 = LineString([(0.005, 0.00575), (0.005, 0.00525)]).buffer(0.00002)
# poly4 = LineString([(0.005, 0.00425), (0.005, 0.00475)]).buffer(0.00002)
# poly5 = LineString([(0.00512, 0.00575), (0.00512, 0.00525)]).buffer(0.00002)
# poly6 = LineString([(0.00513, 0.00575), (0.00513, 0.00525)]).buffer(0.00002)
# poly7 = LineString([(0.00425, 0.00475), (0.00492, 0.00475)]).buffer(0.00002)
# poly8 = LineString([(0.00425, 0.00475), (0.00492, 0.00475)]).buffer(0.00002)
# plt.plot(*poly1.exterior.xy)
# plt.plot(*poly2.exterior.xy)
# plt.plot(*poly3.exterior.xy)
# plt.plot(*poly4.exterior.xy)
# plt.plot(*poly5.exterior.xy)
# plt.plot(*poly6.exterior.xy)
# plt.plot(*poly7.exterior.xy)
# plt.plot(*poly8.exterior.xy)
#
# tree = STRtree([poly1, poly2, poly3, poly4, poly5, poly6, poly7, poly8])
# dash_polys = dot.search_dash_polygons(tree)
# print(dash_polys)
# dash_pairs = [(Dash(dash1), Dash(dash2)) for dash1, dash2 in dash_polys]
# print(dash_pairs)
# dot.save_dash_pairs(dash_pairs)
# print(dot.dash_pairs)
#
# for dash1, dash2 in dot.dash_pairs:
#     plt.plot(*create_centerline(dash1.poly).xy)
#     plt.plot(*create_centerline(dash2.poly).xy)
#
# ax = plt.gca()
# ax.set_xlim([0.004, 0.006])
# ax.set_ylim([0.004, 0.006])
# plt.show()



############## Extrapolation test
# import matplotlib.pyplot as plt
# from shapely.ops import linemerge
# from shapely.geometry import LineString, Point
# from line_proc import get_extrapolated_point
#
# line1 = LineString([(0, 0), (0, 1), (1, 1), (2,1),(3,2)])
# line2 = LineString([(2,1),(3,1),(4,1),(5,1)])
# line = linemerge([line1,line2])
#
# l_coords = list(line2.coords)
# endpoint = Point(l_coords[-1])
# #endpoint = Point(l_coords[0])
#
# ext_p = get_extrapolated_point(line, endpoint)
# print(ext_p)
#
# plt.plot(ext_p.x, ext_p.y, marker="*")
# for l in line:
#   plt.plot(*l.xy)
#
# plt.show()

##### Shortest path test
# import matplotlib.pyplot as plt
# from shapely.geometry import MultiLineString
# from shapely.geometry import Point
# from line_proc import get_path_line
#
# mline = MultiLineString([((0, 0), (1, 1), (1.5, 1)), ((0, 2), (1, 1.5), (1.5, 1), (2, 0))])
#
# path = get_path_line(mline, Point(0,0), Point(2,0))
# print(path)
#
# plt.plot(*path.xy)
#
# # plotting the multiline
# # for line in mline:
# #     if line.type == "MultiLineString":
# #         multilinestr = line
# #         for linestr in multilinestr:
# #             x, y = linestr.xy
# #             plt.plot(x, y)
# #     else:
# #         x, y = line.xy
# #         plt.plot(x, y)
#
# plt.show()

################# dot-dash line drawing test #################################
from itertools import combinations
import matplotlib.pyplot as plt
from shapely.ops import linemerge
from settings import RESOURCE_DIR
import cv2
from shapely.geometry import Point, LineString, MultiLineString
import numpy as np
from line_proc import get_blobs, affine_transform, get_geojson_list, get_shapely_geom, \
  get_line_endpoints, create_centerline, get_extrapolated_point, get_path_line
from dot_dash import Dot, Dash
from shapely.strtree import STRtree

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



#INPUT_FILENAME = '093C10g_1965_D_1.gtiff'
INPUT_FILENAME = '092L3g_1969_D_1_clipped_4263x5137.gtiff'
#INPUT_FILENAME = '092L3g_1969_D_1_clipped_1000x500.gtiff'
INPUT_PATH = f'{RESOURCE_DIR}/{INPUT_FILENAME}'

plt.figure(figsize=(35,33))
img = cv2.imread(INPUT_PATH)

# detecting blobs
blob_keypoints = get_blobs (img, max_blob_area=60)
# plot blobs as a separate image
blank = np.zeros((1, 1))
blobs_img = cv2.drawKeypoints(img, blob_keypoints, blank, (0, 0, 255),cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
blobs_resized = cv2.resize(blobs_img,(960, 720))
#cv2.imshow('original image',img)

# converting the pixel coordinates of blob into geospatial coordinates
blob_px_coords = []
for blob in blob_keypoints:
  #print(f'pt: {blob.pt}, size: {blob.size}')
  blob_px_coords.append(blob.pt)
blob_geo_coords = affine_transform(INPUT_PATH, blob_px_coords)

# plotting the blobs using geospatial coordinates
# for x, y in blob_geo_coords:
#     plt.plot(x, y, marker="o", markeredgecolor="red", markerfacecolor="green")


# make blobs shapely points by using their geo coordinates
blob_points = [Point(blob) for blob in blob_geo_coords]

# get polygons' geojsons from raster image
geojson_list = get_geojson_list(INPUT_PATH, 0) # masking black
print(f'geojson num(found polygons) : {len(geojson_list)}')

geoms = get_shapely_geom(geojson_list)

# plotting map polygons
# for geom in geoms:
#   x,y = geom.exterior.xy # x,y are arrays
#   plt.plot(x,y)

def create_dash_pairs(dot, dash_poly_pairs):
  """
  Create dash pairs from dash polygon pairs
  When Dash object is already created from the same polygon, the dot is stored in that Dash object.
  When Dash object is newly created, the input dot is stored in that Dash object

  dot: Dot object
  dash_poly_pairs: dash polygon pairs around the Dot object
  """
  dash_pairs = []
  for dp1, dp2 in dash_poly_pairs:
    if Dash.is_already_created(dp1):
      dash1 = Dash.get_dash_obj(dp1)
      dash1.save_dots([dot])
    else:
      dash1 = Dash(dp1, [dot])

    if Dash.is_already_created(dp2):
      dash2 = Dash.get_dash_obj(dp2)
      dash2.save_dots([dot])
    else:
      dash2 = Dash(dp2, [dot])

    dash_pairs.append((dash1, dash2))

  return dash_pairs


def find_nearest_geom(target, geom_list):
  """
  Find the nearest geometric object to the target geometric object from list
  """
  min_dist = float('inf')
  for geom in geom_list:
    dist = target.distance(geom)
    if dist < min_dist:
      min_dist = dist
      min_dist_geom = geom

  return min_dist_geom


# dots: dots as shapely Points
# polygons: all the polygons on the map
# dash_range: the range where dashes will be detected around dots
# max_dot_length: maximum length of dot polygon; used to filter dashes from nearby polygons
# max_dash_length: maximum length of dash polygon; used to filter dashes from nearby polygons
# returns drawn lines (lines enclosing dot-dashed lines) merged as MultiLineString or LineString
# tentative default value :  dash_range=0.00035, max_dot_length=0.0005, max_dash_length=0.008
def draw_dot_dashed_lines(dots, polygons, dash_range=0.00030, max_dot_length=0.0005, max_dash_length=0.024):
  all_drawn_lines = []
  dots_tree = STRtree(dots)
  polygons_wo_dots = [poly for poly in polygons if len([dot for dot in dots_tree.query(poly) if
                                                        poly.intersects(dot) and poly.length <= max_dot_length]) == 0]

  # plotting non dot polygons
  # for geom in polygons_wo_dots:
  #   x,y = geom.exterior.xy # x,y are arrays
  #   plt.plot(x,y)

  polygons_wo_dots_tree = STRtree(polygons_wo_dots)

  for p in dots:
    dot = Dot(p)
    #plt.plot(p.x, p.y, marker="o")
    dash_poly_pairs = dot.search_dash_polygons(polygons_wo_dots_tree)
    if len(dash_poly_pairs) == 0:
      del Dot.all_dots[id(p)]
    else:
      dash_pairs = create_dash_pairs(dot, dash_poly_pairs)  # Dash object is created
      dot.save_dash_pairs(dash_pairs)

  # make virtual dots (dots not detected) for dashes having dots less than two
  dashes_copy = list(Dash.all_dashes.values())
  for dash in dashes_copy:
    if len(dash.conn_dots) < 2:
      dash_body_line = create_centerline(dash.poly)
      if dash_body_line is not None:
        endpoints = get_line_endpoints(dash_body_line)

        # filter out the endpoints with the dots
        for dot in dash.conn_dots:
          min_dist_ep = find_nearest_geom(dot.point, endpoints)
          endpoints.remove(min_dist_ep)

        # make virtual dots at the extrapolated location
        # find dashes around the virtual dots
        for ep in endpoints:
          target_p = get_extrapolated_point(dash_body_line, ep)
          #plt.plot(target_p.x, target_p.y, marker="*")
          dot = Dot(target_p)
          dash_poly_pairs = dot.search_dash_polygons(polygons_wo_dots_tree)
          if len(dash_poly_pairs) == 0:
            del Dot.all_dots[id(target_p)]
          else:
            dash_pairs = create_dash_pairs(dot, dash_poly_pairs)  # Dash object is created
            dot.save_dash_pairs(dash_pairs)

  # obtain dash lines
  for dash in Dash.all_dashes.values():
    dash_body_line = create_centerline(dash.poly)
    if dash_body_line is not None:
      endpoints = get_line_endpoints(dash_body_line)

      # prepare dash body line to be merged
      if isinstance(dash_body_line, MultiLineString):
        lines = list(dash_body_line.geoms)
      else:
        lines = [dash_body_line]

      # make the connecting line from dash's endpoint to dot
      for dot in dash.conn_dots:
        min_dist_ep = find_nearest_geom(dot.point, endpoints)
        lines.append(LineString([min_dist_ep, dot.point]))

        # merged line with the dash body line and connecting line to dots
      mline = linemerge(lines)

      # find shortest paths between dots
      path_lines = [get_path_line(mline, dot1.point, dot2.point) for dot1, dot2 in combinations(dash.conn_dots, 2)]
      if len(path_lines) > 0:
        dash_line = linemerge(path_lines)
        dash.dash_line = dash_line

        # plot the final line for dash
        #plot_line(dash_line)

        if isinstance(dash_line, MultiLineString):
          all_drawn_lines.extend(list(dash_line.geoms))
        else: # LineString
          all_drawn_lines.append(dash_line)

  return linemerge(all_drawn_lines)

# drawing extracted lines
Dash.all_dashes.clear()
lines = draw_dot_dashed_lines(blob_points, geoms)

#print(lines)
#plot_line(lines)
#plt.show()

################# Test for writing the result to a shapefile
import fiona
import os
from shapely.geometry import mapping
from settings import OUTPUT_DIR

SHAPEFILE_NAME = f'{INPUT_FILENAME.split(".")[0]}'
SHAPEFILE_PATH = os.path.join(OUTPUT_DIR, SHAPEFILE_NAME)

# assumes CRS as EPSG:4326
def write_shapefile(geo_obj, shapefile_path):
  # define schema
  schema = {
    'geometry': geo_obj.geom_type
  }

  # open a fiona object
  with fiona.open(shapefile_path, mode='w', driver='ESRI Shapefile', schema=schema, crs="EPSG:4326") as shapefile_w:
    # save record and close shapefile
    record = {
      'geometry': mapping(geo_obj)
    }

    shapefile_w.write(record)

print(SHAPEFILE_PATH)
write_shapefile(lines, SHAPEFILE_PATH) # lines : the result of draw_dot_dashed_lines