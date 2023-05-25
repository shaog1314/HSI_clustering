import math

import numpy as np
import torch
from scipy import io
from tifffile import imread
from torch.utils.data import dataset, dataloader
from sklearn import decomposition
from band_mapper import BandMapper
from sklearn.preprocessing import StandardScaler, MinMaxScaler


class SpectralDataset(dataset.Dataset):
    """
    Class representing the hyper/multi spectral dataset.
    """
    ROW_AXIS = 0
    COLUMN_AXIS = 1
    SPECTRAL_AXIS = -1
    SHUFFLE = True

    def __init__(self, spatial_size: int):
        """
        Constructor of the dataset class.

        :param spatial_size: Size of spatial patch, must be odd number.
        """
        super(SpectralDataset, self).__init__()
        assert spatial_size % 2 != 0, 'The \'spatial_size\' argument must be a odd value, passed {}.' \
            .format(spatial_size)
        self.spatial_size = spatial_size
        self.pad_size = int((self.spatial_size - 1) / 2)
        self.data_cube = None
        self.ground_truth = None
        self.spectral_size = None
        self.row_size = None
        self.column_size = None
        self.dataloader = None

    def __len__(self) -> int:
        return int(self.row_size * self.column_size)

    def __getitem__(self, index: int) -> list:
        row_index, col_index = int(index / self.column_size) + self.pad_size, \
                               int(index % self.column_size) + self.pad_size
        sample = None
        try:
            if self.spatial_size == 1:
                sample = torch.DoubleTensor(self.data_cube[row_index, col_index]).unsqueeze(dim=-1)
            elif self.spatial_size > 1:
                sample = torch.DoubleTensor(
                    self.data_cube[row_index - self.pad_size:row_index + self.pad_size + 1,
                    col_index - self.pad_size:col_index + self.pad_size + 1]
                ).permute(2, 0, 1)
                row_index, col_index = row_index - self.pad_size, col_index - self.pad_size
            else:
                raise ValueError('\'spatial_size\' is in incorrect form.')
        except ValueError as error:
            print(error)
        return [sample, torch.LongTensor(np.array([row_index, col_index]))]

    def get_dataloader(self, batch_size: int) -> dataloader.DataLoader:
        """
        Get the dataloader.

        :param batch_size: Size of the batch.
        :return: Dataloader.
        """
        return dataloader.DataLoader(dataset=self, batch_size=batch_size,
                                     shuffle=self.__class__.SHUFFLE, num_workers=4)

    def load_data(self, data_path: str, gt_path: str = None):
        """
        Load data method.

        :param data_path: Path to data file.
        :param gt_path: Path to ground truth if exists.
        :return: None.
        """
        try:
            if data_path.endswith('.npy'):
                self.data_cube = np.load(data_path).astype(np.float64)
            elif data_path.endswith('.tif'):
                self.data_cube = imread(data_path).astype(np.float64)
            elif data_path.endswith('.mat'):
                mat = io.loadmat(data_path)
                for key in mat.keys():
                    if "__" not in key:
                        self.data_cube = np.asarray(mat[key])
                        break
            else:
                raise ValueError('This type of file is not supported for loading.')
        except (IOError, ValueError) as error:
            print(error)
        if gt_path is not None:
            if gt_path.endswith('.npy'):
                self.ground_truth = np.load(gt_path).astype(np.uint8)
            elif gt_path.endswith('.mat'):
                mat = io.loadmat(gt_path)
                for key in mat.keys():
                    if "__" not in key:
                        self.ground_truth = np.asarray(mat[key])
                        break

    def pad_data_cube(self):
        """
        This method applies the zero padding on the data cube, if the spatial_size > 1.

        :return: None.
        """
        if self.spatial_size > 1:
            self.data_cube = np.pad(array=self.data_cube,
                                    pad_width=((self.pad_size, self.pad_size),
                                               (self.pad_size, self.pad_size),
                                               (0, 0)),
                                    mode='constant')
        # Set the shape of the data cube:
        self.spectral_size = self.data_cube.shape[self.SPECTRAL_AXIS]
        self.row_size = self.data_cube.shape[self.ROW_AXIS] - self.pad_size * 2
        self.column_size = self.data_cube.shape[self.COLUMN_AXIS] - self.pad_size * 2

    def preprocess_data_cube(self, method: str = 'normalize', spectral_size: int = None):
        """
        In this method, several preprocessing choices are included.

        :param method: Type of the preprocessing step.
        :param spectral_size: Number of bands.
        :return: None.
        """
        print(method)
        method = method.lower()
        input_shape = self.data_cube.shape[:-1]
        if method == 'normalize':
            print('Normalizing...')
            shape = self.data_cube.shape
            self.data_cube = MinMaxScaler() \
                .fit_transform(self.data_cube.reshape(-1, self.data_cube.shape[self.SPECTRAL_AXIS]))
            self.data_cube = self.data_cube.reshape(shape)
            return
        elif method == 'ica':
            print('ICA...')
            self.data_cube = self.data_cube.reshape(-1, self.data_cube.shape[self.SPECTRAL_AXIS])
            self.data_cube = decomposition.FastICA(n_components=spectral_size, max_iter=1000).fit_transform(
                self.data_cube)
            self.data_cube = StandardScaler().fit_transform(self.data_cube)
        elif method == 'pca':
            print('PCA...')
            self.data_cube = self.data_cube.reshape(-1, self.data_cube.shape[self.SPECTRAL_AXIS])
            self.data_cube = decomposition.PCA(n_components=spectral_size).fit_transform(self.data_cube)
            self.data_cube = StandardScaler().fit_transform(self.data_cube)
        elif method == 's-msi':
            print('S-MSI')
            mapper = BandMapper()
            self.data_cube = self.data_cube.reshape(-1, self.data_cube.shape[self.SPECTRAL_AXIS])
            self.data_cube = mapper.map(data=self.data_cube, resulting_bands_count=spectral_size)
            self.data_cube = StandardScaler().fit_transform(self.data_cube)
        self.data_cube = self.data_cube.reshape(*input_shape, spectral_size)
        print('Spectral size: {}'.format(self.data_cube.shape[self.SPECTRAL_AXIS]))

    @staticmethod
    def get_min(value: int) -> int:
        """
        Get the minimal size of batch.
        :param value: Number of samples.
        :return: Minimal size of batch.
        """
        divisors = []
        for divisor in range(2, math.ceil(math.sqrt(value))):
            if value % divisor == 0:
                divisors.append(int(value / divisor))
        min_ = 1
        if len(divisors) != 0:
            min_ = np.min(divisors)
        return int(min_)

    @staticmethod
    def get_median(value: int) -> int:
        """
        Get the median size of batch.
        :param value: Number of samples.
        :return: Median size of batch.
        """
        divisors = []
        for divisor in range(2, math.ceil(math.sqrt(value))):
            if value % divisor == 0:
                divisors.append(int(value / divisor))
        median = 1
        if len(divisors) != 0:
            if len(divisors) % 2 != 1:
                median = np.median(divisors)
            else:
                median = divisors[int(len(divisors) / 2)]
        return int(median)
