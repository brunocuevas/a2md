from a2mdnet.models.monomernet import Monomernet
from a2mdnet import LIBRARY_PATH
from a2mdio.params import AMDParameters
from a2mdnet.data import MonomerDataset
from torch.utils import data
import torch
import json
from pathlib import Path
if __name__ == '__main__':
    device = torch.device('cpu')
    float_dtype = torch.float
    mol_path = Path('./')

    with open('.dataset.json') as f:
        dataset_info = json.load(f)

    md = MonomerDataset(
        device=device, dtype=float_dtype, ids=dataset_info['index'],
        molecular_data_path=mol_path,
        max_atoms=dataset_info['natoms'],
        max_bonds=dataset_info['nbonds']
    )
    dl = data.DataLoader(md, batch_size=16)
    print('-- testing persistance and device change')
    ppp = AMDParameters.from_file('model.json')
    aevs = LIBRARY_PATH / 'params/aev_polymer.params'
    mn = Monomernet(
        layers=[110, 55, 12],
        density_model=ppp,
        symfeats_model=aevs,
        device=torch.device('cuda:0'),
        dtype=torch.float
    )
    mn = mn.to(device=torch.device('cpu'))

    for index_, labels, topo, coords, charge in dl:
        r = torch.rand(16, 1000, 3, dtype=torch.float, device=torch.device('cpu'))
        p = mn.forward(r, labels, coords, charge)
        break

    print('DONE!')
