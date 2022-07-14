from shapely.ops import polygonize_full
import matplotlib.pyplot as plt
from line_proc import get_dot_points, get_geojson_list, get_shapely_geom
from settings import IMAGE_PATH
from dot_dash import draw_dot_dashed_lines

plt.figure(figsize=(21, 19))

blob_points = get_dot_points(IMAGE_PATH)

# get polygons' geojsons from raster image
geojson_list = get_geojson_list(IMAGE_PATH, 0)  # masking black
print(f'geojson num(found polygons) : {len(geojson_list)}')

# get all polygons on the image
geoms = get_shapely_geom(geojson_list)

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