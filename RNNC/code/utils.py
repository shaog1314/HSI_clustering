import os

import numpy as np
import xlsxwriter
from skimage import img_as_ubyte
from skimage.color import label2rgb
from skimage.io import imsave
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score

from pipeline import MetricsEnum


def save_generated_map(gt_hat: np.ndarray, dest_path: str):
    """
    Save the predicted map.

    :param gt_hat: Predicted map.
    :param dest_path: Destination path.
    :return: None.
    """
    imsave(os.path.join(dest_path, 'generated-map.png'), img_as_ubyte(label2rgb(gt_hat.astype(int))))
    np.save(os.path.join(dest_path, 'generated-map.npy'), gt_hat.astype(int))


def save_gt_based_map(gt_hat: np.ndarray, gt: np.ndarray, dest_path: str):
    """
    Save predicted map based on ground-truth background indexes.
    :param gt_hat: Ground-truth map.
    :param gt: Predicted map.
    :param dest_path: Destination path.
    :return: None.
    """
    background = np.where(gt == 0)
    for row, column in zip(*background):
        gt_hat[row, column] = -1
    imsave(os.path.join(dest_path, 'gt-based-map.png'), img_as_ubyte(label2rgb(gt_hat.astype(int))))


def save_metrics(dest_path: str, artifacts: dict):
    """
    Save the metrics.

    :param dest_path: Destination path.
    :param artifacts: Artifacts.
    :return: None.
    """
    # Get the total time of execution:
    artifacts[MetricsEnum.EXECUTION_TIME].append(np.sum(artifacts[MetricsEnum.TRAIN_TIME] +
                                                        artifacts[MetricsEnum.EVAL_TIME]) / 60)
    workbook = xlsxwriter.Workbook(os.path.join(dest_path, 'metrics.xlsx'))
    worksheet = workbook.add_worksheet()
    row, col = 0, 0
    for key in artifacts.keys():
        row = 1
        col += 1
        worksheet.write(row, col, key.value)
        for item in artifacts[key]:
            row += 1
            worksheet.write(row, col, item)
    workbook.close()


def get_metrics(gt_hat: np.ndarray, gt: np.ndarray, metrics: dict):
    """
    Get the metrics from predicted map.
    :param gt_hat: Predicted map.
    :param gt: Ground-truth map.
    :param metrics: Metrics to save the results.
    :return: None.
    """
    gt_hat = gt_hat[gt != 0]
    gt = gt[gt != 0]
    ars_score = adjusted_rand_score(labels_true=gt, labels_pred=gt_hat)
    nmi_score = normalized_mutual_info_score(labels_true=gt, labels_pred=gt_hat)
    metrics[MetricsEnum.ARS_SCORE].append(ars_score)
    metrics[MetricsEnum.NMI_SCORE].append(nmi_score)
