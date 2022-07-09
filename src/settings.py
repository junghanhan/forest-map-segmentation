import os


# DASH_SEARCH_BOX_W = 0.0003 * 2
# DASH_SEARCH_BOX_H = 0.00005 * 2
DASH_SEARCH_BOX_W = 0.0004 * 2
DASH_SEARCH_BOX_H = 0.00006 * 2
MAX_DOT_LENGTH = 0.0005
MIN_BRANCH_LEN = 0.00001 # minimum branch length that will be regarded as a valid branch of a dash

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SRC_DIR, os.pardir))
RESOURCE_DIR = f'{ROOT_DIR}/res'
OUTPUT_DIR = f'{ROOT_DIR}/output'

DRIVER = 'ESRI Shapefile'
