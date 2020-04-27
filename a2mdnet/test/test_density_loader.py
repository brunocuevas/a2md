from a2mdnet.data import MolecularElectronDensityDataset
from a2mdnet.test.datasets import fdaset
from a2md.utils import integrate_from_old_dict
import torch
from torch.utils import data
import time

import json

params = {'batch_size': 64,
          'shuffle': True,
          'num_workers': 0}

GLOBAL_START = time.time()
LOAD_START = time.time()

d = MolecularElectronDensityDataset(
    device=torch.device('cuda:0'),
    dtype=torch.float,
    ids=fdaset.idx,
    model_parameters_path=fdaset.parameters_path,
    molecular_data_path=fdaset.molecular_path,
    density_data_path=fdaset.density_path,
    mol_prop=fdaset.energy, prefetch=False,
    integration_method=integrate_from_old_dict
)

LOAD_END = time.time()
print("load time : {:8.4f}".format(LOAD_END - LOAD_START))
print("---")
data_generator = data.DataLoader(d, **params)

for i in range(50):
    ITERATION_START = time.time()
    for u in data_generator:
        pass
    ITERATION_END = time.time()
    print("iteration {:4d}, time: {:8.4f}".format(i, ITERATION_END - ITERATION_START))
print("---")
print("global time: {:8.4f}".format(time.time() - GLOBAL_START))