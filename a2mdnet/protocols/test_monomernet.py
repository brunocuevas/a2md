from a2mdnet.density_models import GenAMD
from a2mdnet.data import MonomerDataset
from a2mdnet.modules import QMDensityBatch
from a2mdnet.utils import CoordinatesSampler
from a2mdio.params import AMDParameters
from torch.utils import data
from pathlib import Path
import torch
import time
import click
import json
import sys


@click.command()
@click.option('--path', default=None, help='location of the nns')
@click.option('--device', default='cuda:0', help='device at which the network will be stored')
@click.option('--sampling', type=click.Choice(['spheres', 'random']), default='spheres',
              help='either spherical or random')
@click.option('--batch_size', default=8, help='number of molecules per iteration')
@click.option('--sampling_args', default=None, help='json containing a dict with args')
@click.option('--replicas', default=1, help='repetitions of the dataset')
@click.argument('networks')
@click.argument('dataset')
@click.argument('protomolecule')
def test(
        networks: str, dataset: str, protomolecule: str, path: str, device: str,
        sampling: str, sampling_args: str, batch_size: int, replicas: int
):
    print('-- monomernet testing: {:s}'.format(networks))
    networks = Path(networks)
    if path is not None:
        print('-- networks located at {:s}'.format(path))
        path = Path(path)
    else:
        path = Path('./')
    print('-- opening network info')
    with networks.open() as f:
        nets = [i.strip() for i in f.readlines()]
    print('--     number of networks to validate {:d}'.format(len(nets)))

    print('-- opening dataset')
    dataset = Path(dataset)
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

    print('-- using {:s} device for testing'.format(device))
    device = torch.device(device)
    float_dtype = torch.float

    print('-- checking sampler args')
    if sampling_args is not None:
        with open(sampling_args) as f:
            sampling_args = json.load(f)
    else:
        sampling_args = dict()

    print('-- setup over')
    print('-- using {:s} for protomolecular electron density'.format(protomolecule))
    proto_amd = AMDParameters.from_file(protomolecule)
    gamd = GenAMD(proto_amd, device=torch.device('cuda:0'), dtype=torch.float)

    print('-- declaring sampler')
    cs = CoordinatesSampler(
        sampler=sampling,
        sampler_args=sampling_args,
        dtype=float_dtype, device=device
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

    print('-- MSE : mean squared error')
    print('-- relMSE : mean squared error divided by density')
    print('-- MLSE : mean logarithm squared error')
    print('-- neg : negative exponential score')
    print('@ {:24s} {:12s} {:12s} {:12s} {:12s} {:12s}'.format(
        'model_name', 'replica', 'MSE', 'rel_MSE', 'MLSE', 'time')
    )
    testing_start = time.time()
    for n in nets:
        net_name = path / n
        mn = torch.load(net_name)
        for i in range(replicas):
            mse = 0.0
            rel_mse = 0.0
            mlse = 0.0

            epoch_start = time.time()
            for index, labels, topo, coords, charge in mddl:

                sample = cs(coords)
                density = qmbatch.forward(index, sample)
                protodensity = gamd.protodensity(coordinates=sample, labels=labels, centers=coords)
                protointegrals = gamd.protointegrate(labels)
                dcharge = charge - protointegrals
                preddensity = mn.forward(
                    coordinates=sample, labels=labels, mol_coordinates=coords, charge=dcharge
                )

                molpreddensity = protodensity + preddensity
                del protodensity

                mse += (molpreddensity - density).pow(2.0).sum()
                rel_mse += ((molpreddensity - density).pow(2.0) / density).sum()
                mlse += (molpreddensity.log() - density.log()).sum()

            epoch_end = time.time()

            print(
                '@ {:24s} {:12d} {:12.4e} {:12.4e} {:12.4e} {:12.4e}'.format(
                    n, i, mse, rel_mse, mlse, epoch_end - epoch_start
                )
            )

    testing_end = time.time()
    print('-- testing over. time elapsed: {:12.4f}'.format(testing_end - testing_start))
    print('-- done')


if __name__ == '__main__':
    test()
