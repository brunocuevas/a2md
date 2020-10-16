from a2mdnet.density_models import GenAMD
from a2mdnet.data import MonomerDataset
from a2mdnet.modules import QMDensityBatch
from a2mdnet.utils import CoordinatesSampler
from a2mdio.params import AMDParameters
from torch.utils import data
import torch
import time
if __name__ == '__main__':

    print("--protodensity calculation")
    with open('F:/edip/extended/ata500/.ata.index') as f:
        index = [i.strip() for i in f.readlines()]

    device = torch.device('cuda:0')
    float_dtype = torch.float
    int_dtype = torch.long

    md = MonomerDataset(device=device, dtype=float_dtype, ids=index, molecular_data_path='F:/edip/extended/ata500/')
    mddl = data.DataLoader(md, batch_size=16, shuffle=True)

    qmbatch = QMDensityBatch(
        filename='F:/edip/extended/ata500/.ata.wfn.hdf5',
        index=index, device=device, dtype=float_dtype
    )
    amdp = AMDParameters.from_file('../a2mdt/test/spherical.json')
    gamd = GenAMD(amdp, device=torch.device('cuda:0'), dtype=torch.float)

    cs = CoordinatesSampler(
        sampler=CoordinatesSampler.random_box,
        sampler_args=dict(n_sample=10000, spacing=6.0),
        dtype=float_dtype, device=device
    )

    for index, labels, topo, coords, charge in mddl:
        itstart = time.time()
        sample = cs(coords)
        density = qmbatch.forward(index, sample)
        protodensity = gamd.protodensity(coordinates=sample, labels=labels, centers=coords)
        moldensity = density - protodensity
        itend = time.time()
        print("{:12.6f}".format(itend - itstart))

    print("--done")


