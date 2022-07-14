from shapely.ops import unary_union, polygonize_full
from itertools import combinations
import matplotlib.pyplot as plt
from shapely.ops import linemerge
from settings import RESOURCE_DIR
from shapely.geometry import LineString, MultiLineString
from line_proc import get_dot_points, get_geojson_list, get_shapely_geom, \
    get_line_endpoints, create_centerline, get_extrapolated_point, \
    get_path_line, find_nearest_geom, plot_line
from dot_dash import Dot, Dash, create_dash_pairs
from shapely.strtree import STRtree
import os


#IMAGE_FILE = '092L3g_1969_D_1_clipped.gtiff'
#IMAGE_FILE = '092L3g_1969_D_1_clipped_1000x500.gtiff'
IMAGE_FILE = '093C10f_1975_D_1_clipped_small.tif'
#IMAGE_FILE = '093C10f_1975_D_1_clipped.tif'
IMAGE_PATH = os.path.join(RESOURCE_DIR, IMAGE_FILE)

plt.figure(figsize=(21, 19))


blob_points = get_dot_points(IMAGE_PATH)
for p in blob_points:
    plt.plot(p.x, p.y, marker="o")

# get polygons' geojsons from raster image
geojson_list = get_geojson_list(IMAGE_PATH, 0)  # masking black
print(f'geojson num(found polygons) : {len(geojson_list)}')

# get all polygons on the image
geoms = get_shapely_geom(geojson_list)


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
    for geom in polygons_wo_dots:
      x,y = geom.exterior.xy # x,y are arrays
      plt.plot(x,y)

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
                    plt.plot(target_p.x, target_p.y, marker="*")
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
                plot_line(dash_line)

                if isinstance(dash_line, MultiLineString):
                    all_drawn_lines.extend(list(dash_line.geoms))
                else:  # LineString
                    all_drawn_lines.append(dash_line)

    return linemerge(all_drawn_lines)


lines = draw_dot_dashed_lines(blob_points, geoms)

result_polys, dangles, cuts, invalids = polygonize_full(lines)

print(f'Created final polygons: {len(result_polys)}')
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