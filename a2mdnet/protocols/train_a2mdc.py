import torch
import numpy as np
import json
from torch.utils import data
from a2mdnet.models.direct import DensityCoupled
from a2mdnet.data import MolecularElectronDensityDataset
from pathlib import Path
import pandas as pd
import time
import os

gdb9set = Path(os.environ['GDB9'])
fdaset = Path(os.environ['FDA'])

if gdb9set is None or fdaset is None:
    print("ERROR : please, define environment variables GDB9 and FDA")

configurations = dict(
    narrow_linear=dict(
        common_net=[384, 128, 112, 48], bond_net=[96, 2], atom_net=[48, 2],
        fe_net=None, fe_layer=None
    )
)


def model_generator(
    common_net, atom_net, bond_net, fe_net, fe_layer
):
    with open(gdb9set / 'gdb9set_coefficients_dist_grouped.json') as f:
        norm = json.load(f)
    DEFAULT_ARCHITECTURE = dict(
        common_net=common_net,
        atom_net=atom_net,
        bond_net=bond_net,
        coefficients_distribution=dict(isotropic=norm['isotropic'], anisotropic=norm['anisotropic']),
        fe_net=fe_net,
        fe_layer=fe_layer
    )
    model = DensityCoupled(
        architecture=DEFAULT_ARCHITECTURE,
        device=device
    ).to(device)
    return model


def coefficient_task(dataset, model, optimizer, number_epochs, save_name):
    """

    :param dataset:
    :param model:
    :param optimizer:
    :param number_epochs:
    :param save_name:
    :return:
    """
    history = np.zeros([number_epochs, 3], dtype='float64')
    for t in range(number_epochs):
        start = time.time()
        for batch in dataset:
            l_, c_, x_, q_, i_iso, i_aniso, iso, aniso = batch
            _, _, iso_pred, aniso_pred = model.forward_coefficients(l_, c_, x_, q_, i_iso, i_aniso)
            l2_iso = (iso - iso_pred).pow(2.0).sum()
            l2_aniso = (aniso - aniso_pred).pow(2.0).sum()
            l2 = l2_iso + l2_aniso
            n_iso = iso[iso != 0.0].size()[0]
            n_aniso = aniso[aniso != 0.0].size()[0]
            history[t, 0] += l2_iso.item() / n_iso
            history[t, 1] += l2_aniso.item() / n_aniso
            l2.backward()
            optimizer.step()
            model.zero_grad()
        history[t, 2] = time.time() - start

        history[t, 0] = np.sqrt(history[t, 0] / len(dataset))
        history[t, 1] = np.sqrt(history[t, 1] / len(dataset))

        if t % 10 == 0:
            torch.save(model, save_name.format(t))

    return model, history


def validate(dataset, model):
    """
    validate
    ---
    - validation_dl : object storing all the data used for validation
    - model: object to validate
    returns rmse, rmlse, rmse(iso), rmse(aniso), time
    """
    l2_iso_buffer = 0.0
    l2_aniso_buffer = 0.0
    l2_dens_buffer = 0.0
    l2_ldens_buffer = 0.0
    time_spent = 0.0
    start = time.time()
    for batch in dataset:
        with torch.no_grad():
            l_, c_, x_, q_, i_iso, i_aniso, iso, aniso, samx, samp = batch
            _, _, dens_pred, iso_pred, aniso_pred = model.forward(l_, c_, x_, q_, i_iso, i_aniso, samx)
            n_iso = iso[iso != 0.0].size()[0]
            n_aniso = aniso[aniso != 0.0].size()[0]
            l2_iso = (iso - iso_pred).pow(2.0).sum()
            l2_aniso = (aniso - aniso_pred).pow(2.0).sum()
            l2_dens = (samp - dens_pred).pow(2.0).mean()
            l2_ldens = (torch.log10(samp) - torch.log10(dens_pred)).pow(2.0).mean()
            l2_iso_buffer += l2_iso.item() / n_iso
            l2_aniso_buffer += l2_aniso.item() / n_aniso
            l2_dens_buffer += l2_dens.item()
            l2_ldens_buffer += l2_ldens.item()
            end = time.time()
            time_spent = end - start

    rmse_dens = np.sqrt(l2_dens_buffer)
    rmse_iso = np.sqrt(l2_iso_buffer)
    rmse_aniso = np.sqrt(l2_aniso_buffer)
    rmlse_dens = np.sqrt(l2_ldens_buffer)
    return rmse_dens, rmlse_dens, rmse_iso, rmse_aniso, time_spent


if __name__ == '__main__':

    print("-- training a2mdnet on coefficient reproduction")

    with open(gdb9set / 'gdb9_filtered.json') as f:
        training_idx = json.load(f)[:100]

    with open(fdaset / 'opt/rfdaset_opt.json') as f:
        validation_idx = json.load(f)[:100]

    training_data = MolecularElectronDensityDataset(
        ids=training_idx, device=torch.device('cuda:0'), dtype=torch.float, prefetch=True,
        density_data_path=None,
        molecular_data_path=gdb9set / 'mol2',
        model_parameters_path=gdb9set / 'pyp'
    )

    validation_data = MolecularElectronDensityDataset(
        ids=validation_idx, device=torch.device('cuda:0'), dtype=torch.float, prefetch=True,
        density_data_path=fdaset / 'opt/sampled/train/',
        molecular_data_path=fdaset / 'opt/mol2/',
        model_parameters_path=fdaset / 'opt/pyp/'
    )

    training_data_load_params = dict(
        batch_size=128, shuffle=True
    )
    validation_data_load_params = dict(
        batch_size=256, shuffle=False
    )

    training_dl = data.DataLoader(training_data, **training_data_load_params)
    validation_dl = data.DataLoader(validation_data, **validation_data_load_params)

    device = torch.device('cuda:0')
    # Optimization settings

    epochs = 200
    lr = 1e-4
    betas = (0.9, 0.999)
    weight_decay = 1e-2
    niterations = 5

    training_history_dict = dict()
    print("-- starting training")
    for name, confs in configurations.items():

        print("-- processing {:s} configuration".format(name))

        training_history_dict[name] = []
        validation_history = []
        for i in range(niterations):
            start = time.time()
            torch.manual_seed(42 + i)
            folder_name = './model_{:s}_{:03d}/'.format(name, i)
            try:
                os.mkdir(folder_name)
            except FileExistsError:
                pass
            model = model_generator(
                common_net=confs['common_net'],
                atom_net=confs['atom_net'],
                bond_net=confs['bond_net'],
                fe_net=None, fe_layer=None
            )

            opt = torch.optim.Adam(model.parameters(), betas=betas, lr=lr, weight_decay=weight_decay)

            model, history = coefficient_task(
                dataset=training_dl, model=model, number_epochs=101,
                save_name='{:s}/{:s}'.format(folder_name, name) + '_{:04d}.pt',
                optimizer=opt
            )

            training_history_dict[name].append(history)
            end = time.time()
            print("\t -- iteration {:d} completed, time elapsed : {:8.4f}".format(i + 1, end - start))

            for j in range(0, epochs, 10):

                try:
                    model = torch.load('{:s}/{:s}_{:04d}.pt'.format(folder_name, name, j))
                except FileNotFoundError:
                    continue
                rmse_dens, rmlse_dens, rmse_iso, rmse_aniso, ts = validate(validation_dl, model)
                validation_history.append([j, name, rmse_dens, rmlse_dens, rmse_iso, rmse_aniso, ts])

            history = pd.DataFrame(validation_history, columns=[
                'epoch', 'architecture', 'rmse_d', 'rmlse_d', 'rmse_i', 'rmse_a', 'ts'
            ])
            history['model'] = i
            history.to_csv(folder_name + '.training_history.csv')
            print("\t -- saving to {:s}".format(folder_name + '.training_history.csv'))