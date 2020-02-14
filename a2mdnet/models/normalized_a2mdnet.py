from a2mdnet.a2mdkeras import DirectA2mdNN
import numpy as np

# sys.path.append(r'C:\Users\Bruno\Dropbox\doctorado\aAMDlib\src')


if __name__ == '__main__':

    feats_np_t = np.load(r'C:\scratch\feats_training.npy')
    feats_np_v = np.load(r'C:\scratch\feats_validation.npy')

    labels_np_t = np.load(r'C:\scratch\labels_training.npy')
    labels_np_v = np.load(r'C:\scratch\labels_validation.npy')

    targets_np_t = np.load(r'C:\scratch\targets_training.npy')
    targets_np_v = np.load(r'C:\scratch\targets_validation.npy')

    charges_np_t = np.load(r'C:\scratch\charges_training.npy')
    charges_np_v = np.load(r'C:\scratch\charges_validation.npy')

    amd_charges_np_t = np.load(r'C:\scratch\amd_charges_training.npy')
    amd_charges_np_v = np.load(r'C:\scratch\amd_charges_validation.npy')

    a2mdnn = DirectA2mdNN(
        training_feats=feats_np_t,
        training_labels=labels_np_t,
        training_targets=targets_np_t,
        training_function_charge_tensor=amd_charges_np_t,
        training_charge_tensor=charges_np_t,

        validation_feats=feats_np_v,
        validation_labels=labels_np_v,
        validation_targets=targets_np_v,
        validation_function_charge_tensor=amd_charges_np_v,
        validation_charge_tensor=charges_np_v,
        hidden_layers=[15, 10, 5, 2]
    )

    hist = a2mdnn.train(epochs=100)

    zp = a2mdnn.validate('normalized_net.csv')
    a2mdnn.save_model('tmp_norm_nn')

    print(np.sum(zp[0, :, :] * amd_charges_np_v[0, :, :]))
    print(charges_np_v[0, :, :])

    with open('norm_train_history.csv', 'w') as f:
        f.write('{:<8s}\t{:<12s}\n'.format('epoch', 'lossH'))
        for i, h in enumerate(hist.history['loss']):
            f.write('{:<8d}\t{:<12.4f}\n'.format(i, h))
    #
    # na2mdnn.validate(save_output_filename='naive_net_validation.csv')
    #
    # na2mdnn.save_model(save_model_filename='tmp_naive_net')
    #
    # # na2mdnn.load_model(model_file_name='tmp_naive_net')
