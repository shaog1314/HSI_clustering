import argparse
import os
import typing

import torch

import dataset
import model
import pipeline
import utils


class Arguments(typing.NamedTuple):
    data_path: str
    labels_path: str
    spatial: int
    destination_path: str
    n_clusters: int
    epochs: int
    patience: int
    num_layers: int
    hidden_size: int
    type_rnn: str
    dropout: float


def arguments() -> Arguments:
    """
    Parse arguments passed by the user.

    :return: Parsed arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', dest='data_path', type=str, help='Path to data file.')
    parser.add_argument('--labels_path', dest='labels_path', type=str, help='Path to labels, if exists.', default=None)
    parser.add_argument('--spatial', dest='spatial', type=int, help='Size of spatial patch.')
    parser.add_argument('--destination_path', dest='destination_path', type=str, help='Destination path.')
    parser.add_argument('--n_clusters', dest='n_clusters', type=int, help='Number of clusters.')
    parser.add_argument('--epochs', dest='epochs', type=int, help='Number of epochs.')
    parser.add_argument('--patience', dest='patience', type=int, help='Number of epochs without improvement.')
    parser.add_argument('--num_layers', dest='num_layers', type=int, help='Number of layers in RNN.')
    parser.add_argument('--hidden_size', dest='hidden_size', type=int, help='Hidden size of recurrent layer.')
    parser.add_argument('--type_rnn', dest='type_rnn', type=str,
                        help='Type of RNN gate mechanism: \"LSTM\" or \"GRU\".')
    parser.add_argument('--dropout', dest='dropout', type=float, help='Dropout used by each layer in decoder.')
    return Arguments(**vars(parser.parse_args()))


def main():
    """
    Main loop of the algorithm.

    :return: None.
    """
    # Get the device name:
    device = torch.device('cuda:0') if torch.cuda.is_available() else torch.device('cpu')
    print(device)
    # Get arguments:
    args = arguments()
    # Create destination directory:
    os.makedirs(args.destination_path, exist_ok=True)
    # Get the dataset:
    data = dataset.SpectralDataset(spatial_size=args.spatial)
    data.load_data(data_path=args.data_path, gt_path=args.labels_path)
    data.preprocess_data_cube(method='s-msi', spectral_size=5)
    data.pad_data_cube()
    # Get the dataloader:
    dataloader = data.get_dataloader(batch_size=data.row_size)
    # Instantiate the algorithm model:
    algorithm = model.RecurrentAutoencoder(n_clusters=args.n_clusters, spatial_size=data.spatial_size,
                                           seq_len=data.spectral_size, img_height=data.row_size,
                                           img_width=data.column_size, ground_truth=data.ground_truth,
                                           destination_path=args.destination_path, device=device,
                                           type_rnn=args.type_rnn, hidden_size=args.hidden_size,
                                           num_layers=args.num_layers, dropout=args.dropout).to(device).double()
    # Get the runner:
    runner = pipeline.Pipeline(algorithm, dataloader, args.epochs, args.patience)
    runner.run_pipeline()
    # Save metrics:
    utils.save_metrics(args.destination_path, runner.metrics)


if __name__ == '__main__':
    main()
