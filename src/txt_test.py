from txt_proc import plot_prediction_result, recognize_texts
from settings import IMAGE_PATH, MODEL_PATH, TARGET_ALPHABETS
import rasterio
from shapely.geometry import Polygon
from line_proc import affine_transform

prediction_result = recognize_texts(IMAGE_PATH, MODEL_PATH, TARGET_ALPHABETS)
# plot_prediction_result(IMAGE_PATH, prediction_result)

# affine transform for recognized labels
transformed_result = []
transform = rasterio.open(IMAGE_PATH).transform
for word, bbox in prediction_result:
    geo_box_coords = []
    for px_coord in bbox:
        geo_coord = affine_transform(transform, px_coord)
        geo_box_coords.append(geo_coord)
    transformed_result.append((word, Polygon(geo_box_coords).centroid))

