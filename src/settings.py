import os

# algorithm related constants
# DASH_SEARCH_BOX_W = 0.0003 * 2
# DASH_SEARCH_BOX_H = 0.00005 * 2
DASH_SEARCH_BOX_W = 0.0004 * 2
DASH_SEARCH_BOX_H = 0.00007 * 2
# the distance between dash and dot; used to place virtual dots
# should be smaller than actual distance cause when the dots are not detected,
# the dots are normally merged into the near dash
#DASH_DOT_DIST = 0.0001
DASH_DOT_DIST = 0.00007
MAX_DOT_LENGTH = 0.0005
MIN_BRANCH_LEN = 0.00005 # minimum branch length that will be regarded as a valid branch of a dash
MAX_SWAMP_SYMBOL_LEN = 0.0025 # swamp symbol polygon's maximum length

# centerline creation related constants
INTERPOLATION_DIST = 0.00006
CENTERLINE_BUFFER = 0.00001

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

# input image and input/output directories constants
INPUT_DIR = os.path.join(ROOT_DIR,'input')
#IMAGE_FILE = '093C10f_1975_D_1_clipped_small.tif' # default test image
#IMAGE_FILE = '093C10e_1975_D_1_clipped1.tif'
IMAGE_FILE = '093C10e_1975_D_1_clipped1.tif' ########
#IMAGE_FILE = '093C10f_1975_D_1_clipped2.tif'
OUTPUT_DIR = os.path.join(ROOT_DIR,'output')

# Multithread version constants
IMAGE_FILES = ['092L3g_1969_D_1_clipped_1000x500.gtiff',
               '093C10e_1975_D_1_clipped1.tif',
               '093C10e_1975_D_1_clipped2.tif',
               '093C10f_1975_D_1_clipped.tif']

DRIVER = 'ESRI Shapefile'

MULTITHREAD = True
