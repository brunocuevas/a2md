from a2mdnet.a2mdkeras import NaiveA2mdRF
import numpy as np

DATA_FOLDER = r'C:\Users\Bruno\ownCloud\main\a2mdnet\data/'
MODEL_FOLDER = r'C:\Users\Bruno\ownCloud\main\a2mdnet\models/'


if __name__ == '__main__':

    print("naive-a2mdrf- machine learning for prediction of electron density")

    feats_t = np.load(DATA_FOLDER + 'set_500_training.sym.npy')
    feats_v = np.load(DATA_FOLDER + 'set_500_validation.sym.npy')
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

    feats_t = feats_t[labels_t != 0, :]
    feats_v = feats_v[labels_v != 0, :]
    targets_t = targets_t[labels_t != 0, :]
    targets_v = targets_v[labels_v != 0, :]
    labels_t = labels_t[labels_t != 0]
    labels_v = labels_v[labels_v != 0]

    na2mdrf = NaiveA2mdRF(
        training_feats=feats_t,
        training_labels=labels_t,
        training_targets=targets_t,
        validation_feats=feats_v,
        validation_targets=targets_v,
        validation_labels=labels_v,
        natoms=[1, 6, 7, 8],
        max_depth=3,
        n_estimators=15
    )

    na2mdrf.train()
    na2mdrf.validate(save_output_filename=MODEL_FOLDER + 'naive_a2mdrf_validation.csv')
    na2mdrf.save_model(save_model_filename=MODEL_FOLDER + 'naive_a2mdrf')
