from a2mdnet.density_models import GenAMD
from a2mdnet.data import MonomerDataset
from a2mdnet.modules import QMDensityBatch
from a2mdnet.utils import CoordinatesSampler
from a2mdio.params import AMDParameters
from a2mdnet.functionals import mean_squared_error, kl_divergence
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
@click.option('--batch_size', default=1, help='number of molecules per iteration')
@click.option(
    '--grid', type=click.Choice(['extra_coarse', 'coarse', 'medium', 'tight']),
    default='coarse', help='type of lebdenev sphere'
)
@click.option('--resolution', type=int, default=20, help='gauss-chebysev radial resolution')
@click.argument('networks')
@click.argument('dataset')
@click.argument('protomolecule')
def test(
        networks: str, dataset: str, protomolecule: str, path: str, device: str,
        batch_size: int,  grid: str, resolution: int
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

    print('-- setup over')
    print('-- using {:s} for protomolecular electron density'.format(protomolecule))
    proto_amd = AMDParameters.from_file(protomolecule)
    gamd = GenAMD(proto_amd, device=torch.device('cuda:0'), dtype=torch.float)

    print('-- declaring sampler')
    cs = CoordinatesSampler(
        sampler='becke',
        sampler_args=dict(grid=grid, resolution=resolution),
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

    print('-- DKL : mean squared error')
    print('-- MSE : mean squared error divided by density')
    print('@ {:24s} {:12s} {:12s} {:12s}'.format(
        'model_name', 'MSE', 'DKL', 'time')
    )
    testing_start = time.time()
    for n in nets:
        net_name = path / n
        mn = torch.load(net_name)

        epoch_start = time.time()
        for index, labels, topo, coords, charge in mddl:

            z, w = cs(coords)
            q = qmbatch.forward(index, z)
            pt = gamd.protodensity(coordinates=z, labels=labels, centers=coords)
            pi = gamd.protointegrate(labels)
            dcharge = charge - pi
            p = mn.forward(
                coordinates=z, labels=labels, mol_coordinates=coords, charge=dcharge
            )
            p = p + pt

            Q = (q * w).sum()
            P = (p * w).sum()

            mse = mean_squared_error(p, q, w)
            dkl = kl_divergence(p / P, q/Q, w)

            epoch_end = time.time()

            print(
                '@ {:24s} {:12.4e} {:12.4e} {:12.4e}'.format(
                    n, mse, dkl, epoch_end - epoch_start
                )
            )

    testing_end = time.time()
    print('-- testing over. time elapsed: {:12.4f}'.format(testing_end - testing_start))
    print('-- done')


if __name__ == '__main__':
    test()
