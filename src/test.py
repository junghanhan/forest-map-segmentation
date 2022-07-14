from shapely.ops import unary_union, polygonize_full
from itertools import combinations
import matplotlib.pyplot as plt
from shapely.ops import linemerge
from settings import RESOURCE_DIR
import cv2
from shapely.geometry import Point, LineString, MultiLineString
import numpy as np
from line_proc import get_blobs, affine_transform, get_geojson_list, get_shapely_geom, \
    get_line_endpoints, create_centerline, get_extrapolated_point, get_path_line, \
    plot_line
from dot_dash import Dot, Dash
from shapely.strtree import STRtree
import os
import rasterio



#IMAGE_FILE = '092L3g_1969_D_1_clipped.gtiff'
#IMAGE_FILE = '092L3g_1969_D_1_clipped_1000x500.gtiff'
IMAGE_FILE = '093C10f_1975_D_1_clipped_small.tif'
#IMAGE_FILE = '093C10f_1975_D_1_clipped.tif'
IMAGE_PATH = os.path.join(RESOURCE_DIR, IMAGE_FILE)

MODEL_DIR = '/content/drive/MyDrive/Capstone Project/trained_model'
MODEL_FILE = 'map_labels_v2.h5'
MODEL_PATH = os.path.join(MODEL_DIR, MODEL_FILE)


plt.figure(figsize=(35, 33))

# TODO: need to standardize the variable name to use folder or dir
input_path = IMAGE_PATH

img = cv2.imread(input_path)

# detecting blobs
blob_keypoints = get_blobs(img, max_blob_area=60)
# plot blobs as a separate image
blank = np.zeros((1, 1))
blobs_img = cv2.drawKeypoints(img, blob_keypoints, blank, (0, 0, 255), cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
blobs_resized = cv2.resize(blobs_img, (960, 720))
# cv2.imshow('original image',img)

# converting the pixel coordinates of blob into geospatial coordinates
blob_px_coords = []
for blob in blob_keypoints:
    # print(f'pt: {blob.pt}, size: {blob.size}')
    blob_px_coords.append(blob.pt)

transform = rasterio.open(IMAGE_PATH).transform
blob_geo_coords = [affine_transform(transform, blob_px_coord) for blob_px_coord in blob_px_coords]

# plotting the blobs using geospatial coordinates
# for x, y in blob_geo_coords:
#     plt.plot(x, y, marker="o", markeredgecolor="red", markerfacecolor="green")


# make blobs shapely points by using their geo coordinates
blob_points = [Point(blob) for blob in blob_geo_coords]

# get polygons' geojsons from raster image
geojson_list = get_geojson_list(input_path, 0)  # masking black
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
    poly_line_dict = {}  # id(polygon):centerline
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
        # plt.plot(p.x, p.y, marker="o")
        dash_poly_pairs = dot.search_dash_polygons(polygons_wo_dots_tree, poly_line_dict)
        if len(dash_poly_pairs) == 0:
            del Dot.all_dots[id(p)]
        else:
            dash_pairs = create_dash_pairs(dot, dash_poly_pairs)  # Dash object is created
            dot.save_dash_pairs(dash_pairs)

    # make virtual dots (dots not detected) for dashes having dots less than two
    dashes_copy = list(Dash.all_dashes.values())
    for dash in dashes_copy:
        if len(dash.conn_dots) < 2:
            if id(dash.poly) in poly_line_dict:
                dash_body_line = poly_line_dict[id(dash.poly)]
            else:
                dash_body_line = create_centerline(dash.poly)
                poly_line_dict[id(dash.poly)] = dash_body_line

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
                    # plt.plot(target_p.x, target_p.y, marker="*")
                    dot = Dot(target_p)
                    dash_poly_pairs = dot.search_dash_polygons(polygons_wo_dots_tree, poly_line_dict)
                    if len(dash_poly_pairs) == 0:
                        del Dot.all_dots[id(target_p)]
                    else:
                        dash_pairs = create_dash_pairs(dot, dash_poly_pairs)  # Dash object is created
                        dot.save_dash_pairs(dash_pairs)

    # obtain dash lines
    for dash in Dash.all_dashes.values():
        if id(dash.poly) in poly_line_dict:
            dash_body_line = poly_line_dict[id(dash.poly)]
        else:
            dash_body_line = create_centerline(dash.poly)
            poly_line_dict[id(dash.poly)] = dash_body_line

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
            path_lines = [get_path_line(mline, dot1.point, dot2.point) for dot1, dot2 in
                          combinations(dash.conn_dots, 2)]

            if len(path_lines) > 0:
                dash_line = unary_union(path_lines)
                dash.dash_line = dash_line

                # plot the final line for dash
                # plot_line(dash_line)

                if isinstance(dash_line, MultiLineString):
                    all_drawn_lines.extend(list(dash_line.geoms))
                else:  # LineString
                    all_drawn_lines.append(dash_line)

    return linemerge(all_drawn_lines)



# clear prev results
Dash.all_dashes.clear()
Dot.all_dots.clear()
lines = draw_dot_dashed_lines(blob_points, geoms)

result_polys, dangles, cuts, invalids = polygonize_full(lines)

for poly in result_polys:
    plt.plot(*poly.exterior.xy)

plt.show()



# # INPUT_FILE = '093C10g_1965_D_1.gtiff'
# # INPUT_FILE = '092L3g_1969_D_1_clipped_4263x5137.gtiff'
# # INPUT_FILE = '092L3g_1969_D_1_clipped_1000x500.gtiff'
# # INPUT_PATH = f'{RESOURCE_DIR}/{INPUT_FILE}'
#
# INPUT_FILE1 = '093C10g_1965_D_1.gtiff'
# INPUT_FILE2 = '093C10e_1975_D_1.gtiff'
# INPUT_FILE3 = '092B5h_1970_D_1.gtiff'
#
# if __name__ == '__main__':
#     # import cProfile
#     # import pstats
#     # from pstats import SortKey
#     # profile = cProfile.Profile()
#     # profile.runcall(temp_func, INPUT_PATH)
#     # ps = pstats.Stats(profile)
#     # ps.sort_stats(SortKey.CUMULATIVE).print_stats()
#
#     from multiprocessing import Pool
#     with Pool() as pool:
#         pool.starmap(temp_func,
#                      [(RESOURCE_DIR, INPUT_FILE1),
#                       (RESOURCE_DIR, INPUT_FILE2),
#                       (RESOURCE_DIR, INPUT_FILE3)])
#
#     #temp_func(INPUT_PATH)