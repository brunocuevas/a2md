from a2mdtest.a2mdtests import benzene
from a2mdnet.modules import QMDensityFun
from a2mdnet.utils import CoordinatesSampler
from a2mdnet.data import mol2tensor
import h5py as h5
import torch

if __name__ == '__main__':

    print('-- testing density functionals')
    na, nb, r, l, t, q = mol2tensor(benzene.mol2, device=torch.device('cuda:0'), dtype=torch.float)

    with h5.File(benzene.path / 'benzene.wfn.hdf5', 'r') as f:
        benzene_qm = QMDensityFun(f['benzene.wfn'], dtype=torch.float, device=torch.device('cuda:0'))
    r = benzene_qm.coords.unsqueeze(0)
    z, w = CoordinatesSampler.integration_grid(r, torch.device('cuda:0'), torch.float, resolution=20)


