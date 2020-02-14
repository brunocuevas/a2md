from a2mdnet.a2mdkeras import NaiveA2mdNN
import numpy as np

DATA_FOLDER = r'C:\Users\Bruno\ownCloud\main\a2mdnet\data\set500/'
MODEL_FOLDER = r'C:\Users\Bruno\ownCloud\main\a2mdnet\models/'
# DATA_FOLDER = r'F:\aamdlib\a2mdnet/'


if __name__ == '__main__':

    print("naive-a2mdnet - deep learning for prediction of electron density")

    feats_t = np.load(DATA_FOLDER + 'set_500_training.ani.npy')
    feats_v = np.load(DATA_FOLDER + 'set_500_validation.ani.npy')
    labels_t = np.load(DATA_FOLDER + 'set_500_training.lab.npy')
    labels_v = np.load(DATA_FOLDER + 'set_500_validation.lab.npy')
    targets_t = np.load(DATA_FOLDER + 'set_500_training.ppp.npy')
    targets_v = np.load(DATA_FOLDER + 'set_500_validation.ppp.npy')

    print("flattening tensors")

    feats_t = feats_t.reshape((-1, feats_t.shape[2]))
    feats_v = feats_v.reshape((-1, feats_v.shape[2]))
    labels_t = labels_t.reshape(-1)
    labels_v = labels_v.reshape(-1)
    targets_t = targets_t.reshape((-1, targets_t.shape[2]))
    targets_v = targets_v.reshape((-1, targets_v.shape[2]))

    na2mdnn = NaiveA2mdNN(
        training_feats=feats_t,
        training_labels=labels_t,
        training_targets=targets_t,
        validation_feats=feats_v,
        validation_targets=targets_v,
        validation_labels=labels_v,
        training_loss_function='mean_squared_error',
        validation_loss_function=None,
        hidden_layers=[25, 10, 2],
        natoms=[1, 6, 7, 8]
    )

    na2mdnn.build_model()

    print("model built")

    hist_b, hist = na2mdnn.train(epochs=5000, subnet_epochs=100)

    with open(MODEL_FOLDER + 'naive_a2mdnet_subnet_history_2.dat', 'w') as f:
        f.write('{:<8s}\t{:<12s}\t{:<12s}\t{:<12s}\t{:<12s}\n'.format('epoch', 'lossH', 'lossC', 'lossN', 'lossO'))
        for i, (lossH, lossC, lossN, lossO) in enumerate(
                zip(hist_b[0].history['loss'], hist_b[1].history['loss'], hist_b[2].history['loss'],
                    hist_b[3].history['loss'])
        ):
            f.write('{:<8d}\t{:<12.4e}\t{:<12.4e}\t{:<12.4e}\t{:<12.4e}\n'.format(i, lossH, lossC, lossN, lossO))

    with open(MODEL_FOLDER + 'naive_a2mdnet_history_2.dat', 'w') as f:
        f.write('{:<8s}\t{:<12s}'.format('epoch', 'loss'))
        for i, loss in enumerate(hist.history['loss']):
            f.write('{:<8d}\t{:<12.4e}'.format(i, loss))

    na2mdnn.validate(save_output_filename=MODEL_FOLDER + 'naive_a2mdnet_fe.csv')
    na2mdnn.save_model(save_model_filename=MODEL_FOLDER + 'naive_a2mdnet_fe')
