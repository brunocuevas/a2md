from a2mdnet.density_models import GenAMD
from a2mdnet.data import MonomerDataset
from a2mdnet.modules import QMDensityBatch
from a2mdnet.utils import CoordinatesSampler
from a2mdio.params import AMDParameters
from a2mdnet.models.monomernet import Monomernet
from torch.optim import Adam
from torch.utils import data
from pathlib import Path
from a2mdnet import LIBRARY_PATH
import torch
import time
import click
import json
import sys


@click.command()
@click.option('--output_path', default='./', help='path for output models')
@click.option('--learning_rate', default=1e-4, help='learning_rate')
@click.option('--device', default='cuda:0', help='device at which the network will be stored')
@click.option('--epochs', default=200, help='number of iterations around the dataset')
@click.option('--sampling', type=click.Choice(['spheres', 'random', 'box']), default='spheres',
              help='either spherical or random')
@click.option('--sampling_args', default=None, help='json containing a dict with args')
@click.option('--batch_size', default=16, help='number of molecules at each iteration')
@click.option('--aevs', default=None, help='atomic environment vector config')
@click.option('--seed', default=42, help='random seed number')
@click.option('--save_each', default=10, help='epoch intervals for saving nn model')
@click.argument('name')
@click.argument('dataset')
@click.argument('protomolecule')
@click.argument('molecule')
@click.argument('architecture')
def train(
        name: str, dataset: str, protomolecule: str, molecule: str, architecture: str,
        output_path: str, learning_rate: float, device: str,
        epochs: int, sampling: str, sampling_args: str, batch_size: int, aevs: str,
        seed: int, save_each: int
):
    print('-- monomernet training: {:s}'.format(name))
    print('-- opening dataset')
    print('-- setting random seed to'.format(seed))
    dataset = Path(dataset)
    output_path = Path(output_path)
    dataset_path = dataset.parent
    with dataset.open() as f:
        dataset_info = json.load(f)
    try:
        print('-- loading dataset {:s}'.format(dataset_info['name']))
        print('--    size : {:d}'.format(dataset_info['nmolecules']))
        print('--    natoms : {:d}'.format(dataset_info['natoms']))
        print('--    nbonds : {:d}'.format(dataset_info['nbonds']))
        print('--    molecules at {:s}'.format(dataset_info['mol']))
        print('--    wfn hdf5 at {:s}'.format(dataset_info['wfn']))
        index = dataset_info['index']
    except KeyError:
        print('>> ERROR : missing attribute at dataset information')
        sys.exit()

    mol_path = dataset_path / Path(dataset_info['mol'])
    wfn_path = dataset_path / Path(dataset_info['wfn'])

    print('-- using {:s} device for training'.format(device))
    device = torch.device(device)
    float_dtype = torch.float
    if aevs is None:
        print('-- using default polymer params for aev')
        aevs = LIBRARY_PATH / 'params/aev_polymer.params'

    print('-- processing architecture')
    print('--     architecture string : {:s}'.format(architecture))
    architecture = architecture.split('x')
    layers = [int(i) for i in architecture]

    print('-- checking sampler args')
    if sampling_args is not None:
        with open(sampling_args) as f:
            sampling_args = json.load(f)
    else:
        sampling_args = dict()

    print('-- setup over')
    print('-- using {:s} for protomolecular electron density'.format(protomolecule))
    print('-- using {:s} for molecular electron density'.format(molecule))
    proto_amd = AMDParameters.from_file(protomolecule)
    mol_amd = AMDParameters.from_file(molecule)
    gamd = GenAMD(proto_amd, device=torch.device('cuda:0'), dtype=torch.float)

    print('-- declaring sampler')
    cs = CoordinatesSampler(
        sampler=sampling,
        sampler_args=sampling_args,
        dtype=float_dtype, device=device
    )
    print("-- declaring monomernet")
    mn = Monomernet(
        layers=layers, density_model=mol_amd, symfeats_model=aevs,
        dtype=torch.float, device=device
    )

    print('-- loading dataset')
    ds_load_start = time.time()
    md = MonomerDataset(
        device=device, dtype=float_dtype, ids=index,
        molecular_data_path=mol_path,
        max_atoms=dataset_info['natoms'],
        max_bonds=dataset_info['nbonds']
    )
    ds_load_end = time.time()
    print('--     molecules loaded. time elapsed : {:12.4f}'.format(ds_load_end - ds_load_start))
    ds_load_start = time.time()
    qmbatch = QMDensityBatch(filename=wfn_path, index=index, device=device, dtype=float_dtype)
    ds_load_end = time.time()
    print('--     wavefunctions loaded. time elapsed : {:12.4f}'.format(ds_load_end - ds_load_start))
    mddl = data.DataLoader(md, batch_size=batch_size, shuffle=True)
    print('-- using batches of size : {:d}'.format(batch_size))

    print('-- load over')
    print('-- using ADAM optimizer with a lr : {:12.4e}'.format(learning_rate))
    opt = Adam(params=mn.parameters(), lr=learning_rate)

    print('--starting training for {:d} epochs'.format(epochs))
    print('@ {:12s} {:12s}'.format('loss', 'time'))
    j = 0
    training_start = time.time()
    for i in range(epochs):
        epoch_start = time.time()
        for index, labels, topo, coords, charge in mddl:
            itstart = time.time()
            sample = cs(coords)
            density = qmbatch.forward(index, sample)
            protodensity = gamd.protodensity(coordinates=sample, labels=labels, centers=coords)
            protointegrals = gamd.protointegrate(labels)
            moldensity = density - protodensity
            dcharge = charge - protointegrals
            preddensity = mn.forward(
                coordinates=sample, labels=labels, mol_coordinates=coords, charge=dcharge
            )
            loss = ((moldensity - preddensity).pow(2.0)).sum()
            loss.backward()
            opt.step()
            mn.zero_grad()
            itend = time.time()
            print("@ {:12.6e} {:12.6f}".format(loss.item(), itend - itstart))
            j += 1
        epoch_end = time.time()
        print('-- epoch {:d} over, time elapsed : {:12.4f}'.format(i, epoch_end - epoch_start))

        if i % save_each == 0:
            save_name = output_path / "{:s}_{:04d}.pt".format(name, i)
            print('--saving nn state at epoch {:d} to {:s}'.format(i, str(save_name)))
            torch.save(mn, save_name)
    training_end = time.time()
    print('-- training over. time elapsed: {:12.4f}'.format(training_end - training_start))
    print('-- done')

    save_name = output_path / "{:s}_{:04d}.pt".format(name, epochs - 1)
    print('--saving final nn state at epoch {:d} to {:s}'.format(epochs - 1, str(save_name)))
    torch.save(mn, save_name)


if __name__ == '__main__':
    train()
