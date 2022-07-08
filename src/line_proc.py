import cv2
import rasterio
from rasterio.features import shapes
from centerline.geometry import Centerline
from shapely.geometry import shape, MultiLineString, box, LineString
from shapely.ops import nearest_points, linemerge
from shapely.strtree import STRtree

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
    # simplify the polygon to prevent unnecessary branches
    geom = geom.buffer(0.00001, join_style=1)

    centerline_obj = Centerline(geom, 0.00006)                

    # merge the centerlines into a single centerline 
    multi_line = MultiLineString(centerline_obj)
    merged_line = linemerge(multi_line)        

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


# dots: dots as shapely Points 
# polygons: all the polygons on the map
# dash_range: the range where dashes will be detected around dots
# max_dot_length: maximum length of dot polygon; used to filter dashes from nearby polygons 
# max_dash_length: maximum length of dash polygon; used to filter dashes from nearby polygons 
# returns drawn lines (lines enclosing dot-dashed lines) merged as MultiLineString or LineString
# tentative default value :  dash_range=0.00035, max_dot_length=0.0004, max_dash_length=0.008
def draw_dot_dashed_lines (dots, polygons, dash_range=0.00030, max_dot_length=0.0004, max_dash_length=0.024):        
    # polygons that their centerlines are already drawn; key: id(polygon), value: drawn centerline
    centerlines = {} 

    # drawn lines including centerlines and connecting lines between dot and dashes
    # all lines should be LineString type since it is merged later    
    all_drawn_lines = []     

    tree_polygons = STRtree(polygons)

    # draw centerlines of dashes(lines near dots)        
    for dot in dots:        
        # find the polygons around each dot
        search_box = get_search_box(dot, dash_range)
        bbox_x, bbox_y = search_box.exterior.xy
        #plt.plot(bbox_x, bbox_y)        

        # check the polygons within range whether they are dashes or not
        polys_in_bbox = [poly for poly in tree_polygons.query(search_box) if search_box.intersects(poly)]
        dash_polys = [poly for poly in polys_in_bbox if max_dot_length < poly.length < max_dash_length]
        
        # filter out the polys that their centerline is already created 
        dash_polys_not_drawn = [poly for poly in dash_polys if not id(poly) in centerlines]
        dash_polys_drawn = [poly for poly in dash_polys if not poly in dash_polys_not_drawn]

        # lines already created by dashes that are near the current dot
        nearby_lines = [centerlines[id(poly)] for poly in dash_polys_drawn]

        # draw centerlines of dashes that are not created as lines yet        
        for poly in dash_polys_not_drawn:            
            dash_line = create_centerline(poly)                        
            centerlines[id(poly)] = dash_line
            nearby_lines.append(dash_line)
            
            if dash_line.type == "MultiLineString":
                for line in dash_line: all_drawn_lines.append(line)
            else:
                all_drawn_lines.append(dash_line)

        # draw lines connecting the dot to dashes                
        for line in nearby_lines:
            conn_line = LineString([dot, nearest_points(dot, line)[1]])             
            all_drawn_lines.append(conn_line)
                                
    return linemerge(all_drawn_lines)

def affine_transform (tif_file_path, pixel_coordinates):
    results = []
    with rasterio.open(tif_file_path) as src:        
        for p_x, p_y in pixel_coordinates:
            g_x, g_y = rasterio.transform.xy(transform = src.transform, 
                            rows = p_y, 
                            cols = p_x)            
            results.append((g_x,g_y))    
    return results