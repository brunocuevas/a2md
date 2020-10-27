from a2mdnet.data import MonomerDataset
from a2mdnet.modules import QMDensityBatch
from a2mdnet.utils import DxGrid
from torch.utils import data
from pathlib import Path
import torch
import time
import click
import json
import sys


@click.command()
@click.option('--device', default='cuda:0', help='device at which the network will be stored')
@click.option('--resolution', type=float, default=0.5, help='gauss-chebysev radial resolution')
@click.option('--spacing', type=float, default=3.0, help='margins')
@click.argument('model')
@click.argument('dataset')
def test(
        model: str, dataset: str, device: str, resolution: float, spacing: float
):
    print('-- monomernet volume evaluation')
    print('-- model: {:s}'.format(model))
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

    print('-- loading dataset')
    ds_load_start = time.time()
    md = MonomerDataset(
        device=device, dtype=float_dtype, ids=index,
        molecular_data_path=mol_path,
        max_atoms=dataset_info['natoms'],
        max_bonds=dataset_info['nbonds'],
        batch_size=1, shuffle=False
    )
    ds_load_end = time.time()
    print('--     molecules loaded. time elapsed : {:12.4f}'.format(ds_load_end - ds_load_start))
    ds_load_start = time.time()
    qmbatch = QMDensityBatch(filename=wfn_path, index=index, device=device, dtype=float_dtype)
    ds_load_end = time.time()
    print('--     wavefunctions loaded. time elapsed : {:12.4f}'.format(ds_load_end - ds_load_start))
    mn = torch.load(model, map_location=torch.device(device)).to(torch.device(device))

    dxg = DxGrid(device=device, dtype=float_dtype, resolution=resolution, spacing=spacing)

    for i, labels, topo, coords, charge in md.epoch():
        print('-- \t processing {:s}'.format(index[i]))
        z, r0, basis, dims = dxg.generate_grid(coords)

        q = qmbatch.forward(i, z)
        p = mn.forward(
            coordinates=z, labels=labels, mol_coordinates=coords, charge=charge
        )
        d = p - q

        pdx = dxg.dx(grid=p, r0=r0, basis=basis, dims=dims)
        pdx.write('{:s}_prediction.dx'.format(index[i]))
        qdx = dxg.dx(grid=q, r0=r0, basis=basis, dims=dims)
        qdx.write('{:s}_reference.dx'.format(index[i]))
        ddx = dxg.dx(grid=d, r0=r0, basis=basis, dims=dims)
        ddx.write('{:s}_difference.dx'.format(index[i]))

    print('-- done')


if __name__ == '__main__':
    test()
    data.dataloader