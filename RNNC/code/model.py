import os
from time import time

import numpy as np
import torch
from sklearn.mixture import GaussianMixture
from torch import nn
from torch.utils.data import dataloader
from tqdm import tqdm

from utils import MetricsEnum, get_metrics, save_generated_map, \
    save_gt_based_map


class RecurrentAutoencoder(nn.Module):

    def __init__(self,
                 n_clusters: int,
                 spatial_size: int,
                 seq_len: int,
                 img_height: int,
                 img_width: int,
                 hidden_size: int,
                 num_layers: int,
                 dropout: float,
                 destination_path: str,
                 type_rnn: str,
                 device: torch.device = 'cpu',
                 ground_truth: np.ndarray = None):
        super(RecurrentAutoencoder, self).__init__()
        if ground_truth is not None:
            assert isinstance(ground_truth, np.ndarray), 'Ground-truth should be a numpy array.'
        print('Num layers: {}\tHidden size: {}\t Type rnn: {}'.format(num_layers, hidden_size, type_rnn))
        self.img_height = img_height
        self.img_width = img_width
        self.ground_truth = ground_truth
        self.destination_path = destination_path
        if type_rnn == 'gru':
            self.encoder = nn.GRU(input_size=spatial_size ** 2,
                                  hidden_size=hidden_size,
                                  num_layers=num_layers,
                                  batch_first=True).to(device)
        if type_rnn == 'lstm':
            self.encoder = nn.LSTM(input_size=spatial_size ** 2,
                                   hidden_size=hidden_size,
                                   num_layers=num_layers,
                                   batch_first=True).to(device)
        self.type_rnn = type_rnn
        self.decoder = nn.Sequential(
            nn.Linear(in_features=hidden_size, out_features=128),
            nn.ReLU(),
            nn.Dropout(p=dropout),

            nn.Linear(in_features=128, out_features=256),
            nn.ReLU(),
            nn.Dropout(p=dropout),

            nn.Linear(in_features=256, out_features=spatial_size ** 2 * seq_len),
        ).to(device)
        self.cost_function = nn.MSELoss()
        self.optimizer = torch.optim.Adam(params=self.parameters())
        self.n_clusters = n_clusters
        self.device = device
        self.spatial_size = spatial_size
        self.seq_len = seq_len

    def forward(self, batch: torch.Tensor):
        isinstance(batch, torch.Tensor)
        batch = self._encode(batch)
        batch = self._decode(batch)
        return batch

    def eval_epoch(self, loader: dataloader.DataLoader, metrics: dict):
        assert isinstance(loader, dataloader.DataLoader)
        state = torch.load(os.path.join(self.destination_path, "checkpoint"))
        self.load_state_dict(state['state_dict'])
        self.optimizer.load_state_dict(state['optimizer_state_dict'])
        eval_time = time()
        self.eval()
        outputs, indexes = [], []
        with torch.no_grad():
            for batch, idx in tqdm(loader, total=len(loader)):
                if torch.cuda.is_available():
                    batch = batch.type(torch.cuda.DoubleTensor)
                batch = batch.view(batch.shape[0], batch.shape[1], -1)
                encoder_output = self._encode(batch)
                outputs.append(encoder_output.detach().cpu().numpy().squeeze())
                indexes.append(idx)
        outputs, indexes = np.vstack(outputs), np.vstack(indexes)

        print('Clustering {} components'.format(self.n_clusters))
        clusters = GaussianMixture(n_components=self.n_clusters, max_iter=500).fit_predict(outputs)
        prediction_map = self.get_new_prediciton_map()
        for y_predicted, idx in zip(clusters, indexes):
            prediction_map[idx[0], idx[1]] = y_predicted
        metrics[MetricsEnum.EVAL_TIME].append(time() - eval_time)
        save_generated_map(prediction_map, self.destination_path)

        if self.ground_truth is not None:
            save_gt_based_map(prediction_map.copy(), self.ground_truth.copy(), self.destination_path)
            get_metrics(gt_hat=prediction_map, gt=self.ground_truth.copy(), metrics=metrics)

    def train_epoch(self, loader: dataloader.DataLoader, metrics: dict):
        assert isinstance(loader, dataloader.DataLoader)
        train_time = time()
        self.train()
        epoch_loss = []
        for batch, _ in tqdm(loader, total=len(loader)):
            if torch.cuda.is_available():
                batch = batch.type(torch.cuda.DoubleTensor)
            batch = batch.view(batch.shape[0], batch.shape[1], -1)
            self.zero_grad()
            out = self.forward(batch)
            loss = self.cost_function(out, batch)
            loss.backward()
            self.optimizer.step()
            epoch_loss.append(loss.clone().detach().cpu().numpy())
        epoch_loss = np.mean(epoch_loss)
        print("Training MSE -> {}".format(epoch_loss))
        metrics[MetricsEnum.TRAIN_TIME].append(time() - train_time)
        metrics[MetricsEnum.MSE_LOSS].append(epoch_loss)

    def get_new_prediciton_map(self) -> np.ndarray:
        return np.zeros((self.img_height, self.img_width), dtype=np.int)

    def _encode(self, batch: torch.Tensor):
        hidden = None
        if self.type_rnn == 'lstm':
            out, (hidden, cell) = self.encoder(batch)
        if self.type_rnn == 'gru':
            out, hidden = self.encoder(batch)
        return hidden[-1]

    def _decode(self, batch: torch.Tensor):
        batch = self.decoder(batch)
        batch = batch.view(batch.shape[0], self.seq_len, self.spatial_size * self.spatial_size)
        return batch
