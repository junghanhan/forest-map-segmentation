import os

# --- Dot detection related parameters
# Minimum area of connected regions of black pixels that is regarded as dot
MIN_BLOB_AREA = 20
# Maximum area of connected regions of black pixels that is regarded as dot
MAX_BLOB_AREA = 140
# This parameter determines the distance from the endpoint of dashes and virtual dots.
# It is used when the algorithm determines the location of virtual dots to be created.
DASH_DOT_DIST = 0.00007

# --- Dot-dashed line extraction related parameters
# Maximum dot length (i.e., perimeter of a vector polygon that is regarded as a dot).
# This value is used to filter out polygons of dots when searching for dash polygons.
MAX_DOT_LEN = 0.00048
# The total width of the boxes for searching dash polygons
DASH_SEARCH_BOX_W = 0.000371 * 2
# The total height of the boxes for searching dash polygons
DASH_SEARCH_BOX_H = 0.00009 * 2
# The width of the small box for searching dash polygons.
# This small box is used for the additional search of dash polygons in the line extraction algorithm.
SML_DASH_SEARCH_BOX_W = 0.000466
# The height of the small box for searching dash polygons.
# This small box is used for the additional search of dash polygons in the line extraction algorithm.
SML_DASH_SEARCH_BOX_H = 0.00012
# The radius of the circle used to filter out redundant endpoints of centerline of dash polygons.
# It is common that the centerlines of dash polygon contain unnecessary branches when they are created,
# and due to these branches, there can be redundant endpoints of dashes.
# This value determines the size of a circle and if there are multiple endpoints in this circle,
# they are regarded as redundant except one.
ENDPOINT_FILTER_R = 0.00014
# The radius of the circle used to filter out redundant virtual dots.
# Normally, if a dot is not properly detected, the dashes around this dot try to create a virtual dot.
# In this case, multiple dashes end up creating multiple virtual dots for a single dot that is not detected.
# This value determines the size of a circle and if there are multiple virtual dots in this circle,
# they are regarded as redundant except one.
VDOT_FILTER_R = 0.0002
# The rotation step of dash search boxes in degree.
# After each iteration of searches, the boxes for searching dash polygons are rotated around a dot by this value.
SEARCH_STEP_DEGREE = 15
# Maximum swamp symbol length (i.e., perimeter of a vector polygon that is regarded as a swamp symbol).
# This value is used to filter out polygons of swamp symbols when searching for dash polygons.
MAX_SWAMP_SYMBOL_LEN = 0.0025
# minimum endpoints that a swamp symbol has
MIN_SWAMP_ENDPOINTS = 5
# Maximum length of the final line representing a dash and its connecting lines to nearby dots.
# This value is used to filter out solid lines.
# If this value is too small, more solid lines are regarded as dot-dashed lines.
# If this value is too large, some dot-dashed lines are regarded as solid lines and filtered out.
MAX_DASH_LINE_LEN = 0.0035
# Maximum dot-to-dot distance. If the distance between two dots with a dash in between is too far,
# it is likely a solid line crosses this dash and connects those two dots.
# This value is used to filter out this kind of cases.
MAX_D2D_DISTANCE = 0.003
# This value is used as an argument of a function called Centerline of centerline python package.
# The higher value, the less branches on the created centerline.
INTERPOLATION_DIST = 0.00006
# This value is to simplify a polygon before its centerline is created from the polygon.
# The reason to simplify a polygon is to have less branches on the created centerline.
# The higher the value, the bigger the polygon.
CENTERLINE_BUFFER = 0.00001
# This value is used to simplify the centerline created by a polygon. The higher value, the simpler centerline.
SIMPLIFY_TOLERANCE = 0.00001  # this is solely to reduce the execution time
# Image bounding box inside buffer value. The larger the value, the bigger buffer inside the image bounding box.
# This image bounding box is used to find intersection points between the edge of the image and the dot-dashed lines.
IMAGE_BBOX_BUFFER = -0.00005

# --- Other miscellaneous parameters
EPSG_CODE = 4326

# basic directory name constants
SRC_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SRC_DIR, os.pardir))

# AI related constants
MODEL_DIR = os.path.join(ROOT_DIR,'trained_model')
MODEL_FILE = 'map_labels_v3.h5'
MODEL_PATH = os.path.join(MODEL_DIR, MODEL_FILE)
# Target alphabet should match the trained model
TARGET_ALPHABETS = '()*+-.0123456789ABCDEFGHIKLMNOPRSTUVWXYabcdefghijklmnopqrstuvwxyz'
DEFAULT_ALPHABETS = '0123456789abcdefghijklmnopqrstuvwxyz'  # Keras-OCR pretrained model default alphabet set

# input image and input/output directories constants
INPUT_DIR = os.path.join(ROOT_DIR,'input')
OUTPUT_DIR = os.path.join(ROOT_DIR,'output')

DRIVER = 'ESRI Shapefile'
