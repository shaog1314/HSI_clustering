import tensorflow as tf
from keras.losses import mse
from keras import backend as K

import numpy as np
import os


# reparameterization trick
# instead of sampling from Q(z|X), sample epsilon = N(0,I)
# z = z_mean + sqrt(var) * epsilon
import dataset


import numpy as np
def sampling(args):
    """Reparameterization trick by sampling from an isotropic unit Gaussian.
    # Arguments
        args (tensor): mean and log of variance of Q(z|X)
    # Returns
        z (tensor): sampled latent vector
    """

    z_mean, z_log_var = args
    batch = K.shape(z_mean)[0]
    dim = K.int_shape(z_mean)[1]
    # by default, random_normal has mean = 0 and std = 1.0
    epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
    return z_mean + K.sqrt(K.exp(z_log_var)) * epsilon


reducion = 'norm'
datasetname = 'houston'
original_dim = 50*5*5
epochs=10
data = dataset.SpectralDataset(spatial_size=5)
data.load_data(data_path=r'D:\Hypernet\data\{}.npy'.format(datasetname),
               gt_path=r'D:\Hypernet\data\{}_gt.npy'.format(datasetname))

data.reduce_dim(type_=reducion)
data.padd()

x_train = np.zeros((data.ground_truth.shape[0] * data.ground_truth.shape[1], original_dim))
counter = 0
for i in range(0, data.ground_truth.shape[0]):
    for j in range(0, data.ground_truth.shape[1]):
        x_train[counter] = data.data_cube[i:i+5, j:j+5].ravel()
        counter +=1
# x_train = samples.copy()
#
# x_train = np.stack(x_train, axis=0)
# x_train = tf.data.Dataset.from_tensor_slices(x_train).batch(batch_size=1024, drop_remainder=False)


# network parameters
batch_size = 1024
latent_dim = 25

# VAE model = encoder + decoder
# build encoder model
inputs = tf.keras.Input(shape=(original_dim, ), name='encoder_input')
x = tf.keras.layers.Dense(512, activation='relu')(inputs)

x = tf.keras.layers.Dense(latent_dim, activation='relu', name='before_last')(x)

z_mean = tf.keras.layers.Dense(latent_dim, name='z_mean')(x)
z_log_var = tf.keras.layers.Dense(latent_dim, name='z_log_var')(x)

# use reparameterization trick to push the sampling out as input
# note that "output_shape" isn't necessary with the TensorFlow backend
z = tf.keras.layers.Lambda(sampling, output_shape=(latent_dim,), name='z')([z_mean, z_log_var])


# instantiate encoder model
encoder = tf.keras.Model(inputs, [z_mean, z_log_var, z], name='encoder')
encoder.summary()


# build decoder model
latent_inputs = tf.keras.Input(shape=(latent_dim,), name='z_sampling')
x = tf.keras.layers.Dense(512, activation='relu')(latent_inputs)
outputs = tf.keras.layers.Dense(original_dim, activation='sigmoid')(x)

# instantiate decoder model
decoder = tf.keras.Model(latent_inputs, outputs, name='decoder')
decoder.summary()


# instantiate VAE model
outputs = decoder(encoder(inputs)[2])
vae = tf.keras.Model(inputs, outputs, name='vae_mlp')


import os
from sklearn import decomposition, mixture
from sklearn import cluster
from time import time
from utils import save_gt_based_map, save_gt, save_metrics, save_generated_map, get_metrics, MetricsEnum
if __name__ == '__main__':

    dest_path = r'D:\unsupervised-segmentation-rnn\vae'
    dest_path = os.path.join(dest_path, datasetname, reducion)
    os.makedirs(dest_path, exist_ok=True)

    reconstruction_loss = mse(inputs, outputs)
    reconstruction_loss *= original_dim

    kl_loss = 1 + z_log_var - K.square(z_mean) - K.exp(z_log_var) # distance from N(0, 1) of mu and sigma
    kl_loss = K.sum(kl_loss, axis=-1)
    kl_loss *= -0.5

    vae_loss = K.mean(reconstruction_loss + kl_loss)
    vae.add_loss(vae_loss)
    vae.compile(optimizer='adam')
    traintime = None
    artifacts = None
    if not os.path.exists(os.path.join(dest_path, 'vae-weights.h5')):
        traintime = time()
        artifacts = vae.fit(x_train,  epochs=epochs, batch_size=batch_size)
        vae.save_weights(os.path.join(dest_path, 'vae-weights.h5'))
        traintime = time()-traintime


    vae.load_weights(os.path.join(dest_path, 'vae-weights.h5'))

    model = tf.keras.Model(vae.inputs, vae.layers[1].layers[-4].output)
    predicted = model.predict(x_train,batch_size=batch_size)
    print('predicted shape: {}'.format(predicted.shape))
    metrics = {
        MetricsEnum.MSE_LOSS: [],
        MetricsEnum.NMI_SCORE: [],
        MetricsEnum.ARS_SCORE: [],
        MetricsEnum.TRAIN_TIME: [],
        MetricsEnum.EVAL_TIME: [],
        MetricsEnum.EXECUTION_TIME: []
    }
    if traintime is not None:
        metrics[MetricsEnum.TRAIN_TIME] = [traintime]
    if artifacts is not None:
        metrics[MetricsEnum.MSE_LOSS] = artifacts.history['loss']
    gt = data.ground_truth
    n_clusters = np.unique(gt.ravel()).size
    eval_time = time()

    predicted = mixture.GaussianMixture(n_components=n_clusters, random_state=0).fit_predict(predicted)

    preds = predicted.reshape(data.ground_truth.shape)
    path = os.path.join(dest_path, 'gm')
    os.makedirs(path, exist_ok=True)
    save_generated_map(preds, path)
    save_gt_based_map(preds, gt, path)
    get_metrics(gt_hat=preds,gt=gt, metrics=metrics)
    save_gt(gt, path)
    print("saved gt")
    metrics[MetricsEnum.EVAL_TIME].append(time() - eval_time)
    save_metrics(path, metrics)
    print("Done")

    # FOR PATCHES:
    # indiana_steps = [[(0, 0), (145, 73)],
    #                  [(0, 73), (145, 145)]]
    # pavia_steps = [[(0, 0), (152, 170)],
    #                [(0, 170), (152, 340)],
    #                [(152, 0), (304, 170)],
    #                [(152, 170), (304, 340)],
    #                [(304, 0), (456, 170)],
    #                [(304, 170), (456, 340)],
    #                [(456, 0), (610, 170)],
    #                [(456, 170), (610, 340)]]
    # # salinas_steps = [[(0, 0), (128, 217)],
    # #                  [(128, 0), (256, 217)],
    # #                  [(256, 0), (384, 217)],
    # #                  [(384, 0), (512, 217)]]
    # # houston_steps = []
    # # for i in range(0, 1201, 200):
    # #     for j in range(0, 4601, 200):
    # #         if i + 200 > 1200:
    # #             break
    # #         last = j + 200
    # #         if last == 4800:
    # #             last = 4768
    # #         xval = i + 200
    # #         if xval == 1200:
    # #             xval = 1202
    # #
    # from utils import MetricsEnum
    # for window_id, range_ in enumerate(pavia_steps):
    #     patch_indexes = []
    #     patch_outputes = []
    #
    #     for i in range(range_[0][0], range_[1][0]):
    #         for j in range(range_[0][1], range_[1][1]):
    #             sample = predicted[int(i*data.ground_truth.shape[1]) + j]
    #             patch_indexes.append([i, j])
    #             patch_outputes.append(np.squeeze(sample))
    #
    #     patch_outputes = np.asarray(patch_outputes)
    #
    #     metrics = {
    #         MetricsEnum.MSE_LOSS: [],
    #         MetricsEnum.NMI_SCORE: [],
    #         MetricsEnum.ARS_SCORE: [],
    #         MetricsEnum.TRAIN_TIME: [],
    #         MetricsEnum.EVAL_TIME: [],
    #         MetricsEnum.EXECUTION_TIME: [],
    #         MetricsEnum.DELTA: [],
    #         MetricsEnum.NPIXELS: [],
    #         MetricsEnum.NCLUSTERS: []
    #     }
    #
    #     eval_time = time()
    #     # patch_outputes = decomposition.PCA(2).fit_transform(patch_outputes)
    #     # patch_outputes = np.ascontiguousarray(patch_outputes)
    #     gt = data.ground_truth[range_[0][0]:range_[1][0], range_[0][1]:range_[1][1]]
    #     n_clusters = np.unique(gt.ravel()).size
    #     print('N clusters for this patch: {}'.format(n_clusters))
    #     #
    #     clu = cluster.AgglomerativeClustering(n_clusters=n_clusters,
    #                                             affinity='euclidean',
    #                                             linkage='ward').fit_predict(patch_outputes)
    #     # clu = Cluster(patch_outputes)
    #     # deltas = clu.delta
    #     # deltas = np.sort(deltas)[-(n_clusters+1)]
    #     # metrics[MetricsEnum.DELTA].append(deltas)
    #     # metrics[MetricsEnum.NCLUSTERS].append(n_clusters)
    #     # clu.assign(min_density=0, min_delta=deltas)
    #
    #     preds = clu#.membership.astype(int)
    #     gthatmap = np.zeros((range_[1][0]-range_[0][0], range_[1][1]-range_[0][1]))
    #     for y_predicted, idx in zip(preds, patch_indexes):
    #         gthatmap[idx[0]-range_[0][0], idx[1]-range_[0][1]] = y_predicted
    #     preds = gthatmap
    #     # preds = AffinityPropagation().fit_predict(patch_outputes)\
    #     #     .reshape(range_[1][0]-range_[0][0], range_[1][1]-range_[0][1]).astype(int)
    #     # pred_from_data = self.data_cube[range_[0][0]:range_[1][0], range_[0][1]:range_[1][1]]
    #     # preds = pred_map
    #     print('Clusters in predicted map:')
    #     print(np.unique(preds).size)
    #     # assert np.unique(preds).size == np.unique(gt.ravel()).size
    #     path = os.path.join(dest_path, 'hc')
    #     small_path = os.path.join(path, 'patch-{}'.format(window_id))
    #     os.makedirs(small_path, exist_ok=True)
    #     save_generated_map(preds, small_path)
    #     save_gt_based_map(preds, gt, small_path)
    #     get_metrics(gt_hat=preds,gt=gt, metrics=metrics)
    #     save_gt(gt, small_path)
    #     metrics[MetricsEnum.EVAL_TIME].append(time() - eval_time)
    #     save_metrics(small_path, metrics)

