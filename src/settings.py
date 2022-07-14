import os


# DASH_SEARCH_BOX_W = 0.0003 * 2
# DASH_SEARCH_BOX_H = 0.00005 * 2
DASH_SEARCH_BOX_W = 0.0004 * 2
DASH_SEARCH_BOX_H = 0.00007 * 2
MAX_DOT_LENGTH = 0.0005
MIN_BRANCH_LEN = 0.00001 # minimum branch length that will be regarded as a valid branch of a dash

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SRC_DIR, os.pardir))
RESOURCE_DIR = f'{ROOT_DIR}/res'
IMAGE_FILE = '093C10f_1975_D_1_clipped_small.tif'
IMAGE_PATH = os.path.join(RESOURCE_DIR, IMAGE_FILE)

MODEL_DIR = f'{ROOT_DIR}/trained_model'
MODEL_FILE = 'map_labels_v3.h5'
MODEL_PATH = os.path.join(MODEL_DIR, MODEL_FILE)
# Target alphabet should match the trained model
TARGET_ALPHABETS = '()*+-.0123456789ABCDEFGHIKLMNOPRSTUVWXYabcdefghijklmnopqrstuvwxyz'

OUTPUT_DIR = f'{ROOT_DIR}/output'


DRIVER = 'ESRI Shapefile'
