from a2mdnet.modules import TorchaniFeats
from a2mdnet.test.datasets import aniset
from a2mdnet.data import MolecularDataset
from a2md.utils import integrate_from_old_dict
from a2mdnet.data import match_old_fun_names
import torch
from torch.utils import data

device = torch.device('cpu')
params = {'batch_size': 64,
          'shuffle': False,
          'num_workers': 0}

fdaset_ds = MolecularDataset(
    device=torch.device('cuda:0'), dtype=torch.float32,
    ids=aniset.idx[:5], model_parameters_path=aniset.environ_path / "nppp",
    molecular_data_path=aniset.environ_path / "nppp", prefetch=False
)

data_generator = data.DataLoader(fdaset_ds, **params)

ta_feats = TorchaniFeats(device=torch.device('cuda:0')).to(torch.device('cuda:0'))

for labels, topology, coords, charge, iso_integral, aniso_integral, iso, aniso in data_generator:
    x = iso_integral[0, :, :]
    v = iso[0, :, :]
    y = aniso_integral[0, :, :]
    w = aniso[0, :, :]
    z = ((x * v).sum() + (y * w).sum())
    print(z)


