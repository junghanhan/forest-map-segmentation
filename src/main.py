from shapely.ops import polygonize_full
import matplotlib.pyplot as plt
from line_proc import get_dot_points, get_geojson_list, get_shapely_geom
from dot_dash import draw_dot_dashed_lines
from txt_proc import recognize_texts
from settings import IMAGE_PATH, MODEL_PATH, SHAPEFILE_PATH, TARGET_ALPHABETS
import rasterio
from shapely.geometry import Polygon
from line_proc import affine_transform
from gis_io import write_shapefile


def main():
    # ----- Line Extraction
    # dot detection
    blob_points = get_dot_points(IMAGE_PATH)

    # get polygons' geojsons from raster image
    geojson_list = get_geojson_list(IMAGE_PATH, 0)  # masking black
    print(f'Number of found polygons : {len(geojson_list)}')

    # get all polygons on the image
    geoms = get_shapely_geom(geojson_list)

    lines = draw_dot_dashed_lines(blob_points, geoms)
    result_polys, dangles, cuts, invalids = polygonize_full(lines)
    result_polys = list(result_polys)

    print(f'Created final polygons: {len(result_polys)}')
    # for poly in result_polys:
    #     plt.plot(*poly.exterior.xy)

    # ----- Label Extraction
    ocr_result = recognize_texts(IMAGE_PATH, MODEL_PATH, TARGET_ALPHABETS)
    # plot_prediction_result(IMAGE_PATH, prediction_result)

    # affine transform for recognized labels
    labels = []
    transform = rasterio.open(IMAGE_PATH).transform
    for word, bbox in ocr_result:
        geo_box_coords = []
        for px_coord in bbox:
            geo_coord = affine_transform(transform, px_coord)
            geo_box_coords.append(geo_coord)
        labels.append((word, Polygon(geo_box_coords).centroid))

    # ----- Writing Shapefile
    write_shapefile(result_polys, labels, SHAPEFILE_PATH)


if __name__ == '__main__':
    main()


