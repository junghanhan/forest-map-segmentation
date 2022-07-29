import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'  # suppress keras warning

import matplotlib.pyplot as plt
import keras_ocr
import cv2
import math
import numpy as np




def recognize_texts(image_path, model_path=None, target_alphabets='0123456789abcdefghijklmnopqrstuvwxyz'):
    """
    Recognize the designated alphabets from the image

    :param image_path: GeoTiff image path that OCR will be done
    :param model_path: keras-OCR trained model path
    :param target_alphabets: the characters that will be detected and recognized.
        They should match the trained model's target alphabet
    :return: (word, box) tuples recognized from the image.
    """

    if model_path is None:
        # initializing a new model
        recognizer = keras_ocr.recognition.Recognizer()
    else:
        # loading the previously trained model
        recognizer = keras_ocr.recognition.Recognizer(alphabet=target_alphabets)
        recognizer.model.load_weights(model_path)

    recognizer.compile()

    # loading the recognizer
    # keras-ocr will automatically download pretrained weights for the detector.
    pipeline = keras_ocr.pipeline.Pipeline(recognizer=recognizer)

    # Image that OCR will be done
    image = keras_ocr.tools.read(image_path)

    # predictions in prediction_group is a list of
    # (word, box) tuples.
    prediction_results = pipeline.recognize([image])
    # filter out empty("") prediction results
    result = [(word, box) for word, box in prediction_results[0]
              if word is not ""]

    return result


def plot_prediction_result(image_file_path, prediction_result):
    image = keras_ocr.tools.read(image_file_path)

    # plot the predictions
    _, ax = plt.subplots(nrows=1, figsize=(20, 20))
    keras_ocr.tools.drawAnnotations(image=image, predictions=prediction_result, ax=ax)
    plt.show()


# returns text removed image (cv2)
def get_text_removed_image(image_file_path, prediction_result):
    def midpoint(x1, y1, x2, y2):
        x_mid = int((x1 + x2) / 2)
        y_mid = int((y1 + y2) / 2)
        return (x_mid, y_mid)

    # read image
    img = keras_ocr.tools.read(image_file_path)

    mask = np.zeros(img.shape[:2], dtype="uint8")
    for box in prediction_result:
        x0, y0 = box[1][0]
        x1, y1 = box[1][1]
        x2, y2 = box[1][2]
        x3, y3 = box[1][3]

        x_mid0, y_mid0 = midpoint(x1, y1, x2, y2)
        x_mid1, y_mi1 = midpoint(x0, y0, x3, y3)

        thickness = int(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2))

        cv2.line(mask, (x_mid0, y_mid0), (x_mid1, y_mi1), 255,
                 thickness)
        text_removed_img = cv2.inpaint(img, mask, 7, cv2.INPAINT_NS)

    return text_removed_img