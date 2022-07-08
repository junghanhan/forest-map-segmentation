import matplotlib.pyplot as plt
import keras_ocr
import cv2
import math
import numpy as np

# returns (word, box) tuples recognized from the image
def recognize_texts(image_file_path, recognizer):
    # keras-ocr will automatically download pretrained
    # weights for the detector and recognizer.
    pipeline = keras_ocr.pipeline.Pipeline(recognizer=recognizer)

    # Image that OCR will be done
    image = keras_ocr.tools.read(image_file_path)

    # predictions in prediction_group is a list of
    # (word, box) tuples.
    prediction_results = pipeline.recognize([image])

    return prediction_results[0]


def plot_prediction_result(image_file_path, prediction_result):
    image = keras_ocr.tools.read(image_file_path)

    # plot the predictions
    _, ax = plt.subplots(nrows=1, figsize=(20, 20))
    keras_ocr.tools.drawAnnotations(image=image, predictions=prediction_result, ax=ax)


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