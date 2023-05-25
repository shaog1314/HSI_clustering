import os
from enum import Enum

import torch
from torch import nn
from torch.utils.data import dataloader


class MetricsEnum(Enum):
    MSE_LOSS = 'mse_loss'
    NMI_SCORE = 'nmi_score'
    ARS_SCORE = 'ars_score'
    TRAIN_TIME = 'train_time'
    EVAL_TIME = 'eval_time'
    EXECUTION_TIME = 'execution_time'


class Pipeline(object):
    """
    Class for the pipeline of the entire algorithm.
    """

    def __init__(self, model: nn.Module, loader: dataloader.DataLoader, epochs: int, patience: int):
        assert isinstance(model, nn.Module)
        assert isinstance(loader, dataloader.DataLoader)
        assert hasattr(model, 'destination_path')
        assert hasattr(model, 'optimizer')
        self.metrics = {
            MetricsEnum.MSE_LOSS: [],
            MetricsEnum.NMI_SCORE: [],
            MetricsEnum.ARS_SCORE: [],
            MetricsEnum.TRAIN_TIME: [],
            MetricsEnum.EVAL_TIME: [],
            MetricsEnum.EXECUTION_TIME: [],
        }
        self.model = model
        self.loader = loader
        self.epochs = epochs
        self.patience = patience

    def stopping_condition(self) -> bool:
        """
        Return True if no improvement and should stop learning.

        :return: Boolean.
        """
        return min(self.metrics[MetricsEnum.MSE_LOSS][:-self.patience]) <= \
               min(self.metrics[MetricsEnum.MSE_LOSS][-self.patience:])

    def save_model(self):
        """
        Save the model.

        :return: None.
        """
        if len(self.metrics[MetricsEnum.MSE_LOSS]) < 2 or \
                self.metrics[MetricsEnum.MSE_LOSS][-1] < min(self.metrics[MetricsEnum.MSE_LOSS][:-1]):
            print('Saving improvement...')
            torch.save(dict(state_dict=self.model.state_dict(),
                            optimizer_state_dict=self.model.optimizer.state_dict()),
                       os.path.join(self.model.destination_path, "checkpoint"))

    def run_pipeline(self):
        """
        Run the entire pipeline.

        :return: None.
        """
        for epoch in range(self.epochs):
            print("Epoch -> {}".format(epoch + 1))
            self.model.train_epoch(self.loader, self.metrics)
            if epoch > self.patience and self.stopping_condition():
                break
            self.save_model()
        self.model.eval_epoch(self.loader, self.metrics)
        print('NMI: {}\tARS: {}'.format(self.metrics[MetricsEnum.NMI_SCORE],
                                        self.metrics[MetricsEnum.ARS_SCORE]))
