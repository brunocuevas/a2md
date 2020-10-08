from a2mdnet.modules import QPSegmentChargeNormalization, QPChargeNormalization
from a2mdnet.data import MolecularDataset, PolymerDataset
from a2mdnet.density_models import GenAMD
from a2mdio.params import AMDParameters
from torch.utils.data import DataLoader
import torch
import json

if __name__ == '__main__':

    print("-- charge normalization test")
    print("-- using two datasets")
    print("-- 1. aniset (first 6 molecules")
    print("-- 2. edip/ccc (first 6 molecules")

    qp = QPChargeNormalization(device=torch.device('cuda:0'))
    qps = QPSegmentChargeNormalization(device=torch.device('cuda:0'))

    with open('F:/aniset/aniset_all.json') as f:
        idx_monomers = json.load(f)[:5]

    with open('F:/edip/dynamic/ccc/.ccc.input') as f:
        idx_polymers = json.load(f)[:6]

    test_monomer_data = MolecularDataset(
        device=torch.device('cuda:0'), dtype=torch.float,
        ids=idx_monomers,
        model_parameters_path='F:/aniset/nppp/',
        molecular_data_path='F:/aniset/mol2/',
        prefetch=True
    )

    test_polymer_data = PolymerDataset(
        device=torch.device('cuda:0'), dtype=torch.float, ids=idx_polymers,
        molecular_data_path='F:/edip/dynamic/ccc/',
    )

    dl_monomer = DataLoader(test_monomer_data, batch_size=32, shuffle=False)
    dl_polymer = DataLoader(test_polymer_data, batch_size=32, shuffle=False)

    print("-- using spherical model of electron density")

    gamd = GenAMD(
        AMDParameters.from_file('../a2mdt/test/spherical.json'),
        device=torch.device('cuda:0'), dtype=torch.float
    )
    print("-- testing standard normalization")
    for _, _, _, q, i_iso, i_aniso, iso, aniso in dl_monomer:
        u, v = qp.forward(iso, aniso, i_iso, i_aniso, q)
        break

    qpred = (u * i_iso).sum((1, 2)) + (v * i_aniso).sum((1, 2))
    qmod = (iso * i_iso).sum((1, 2)) + (aniso * i_aniso).sum((1, 2))
    qreal = q.sum(1)
    print("{:12s} {:12s} {:12s}".format("prediction", "model", "real"))
    for i in range(u.size()[0]):
        print(
            "{:12.4f} {:12.4f} {:12.4f}".format(qpred[i].item(), qmod[i].item(), qreal[i].item())
        )

    print("-- charge check in qp finished")
    print("-- testing segment normalization")
    for i, l, t, x, q, s, sq in dl_polymer:
        i_iso = gamd.integrate(l)
        iso = torch.ones_like(i_iso)
        u = qps.forward(iso, i_iso, s, sq)
        qreal, qpred = qps.check(u, i_iso, s, sq)
        _, qmod = qps.check(iso, i_iso, s, sq)
    print("{:12s} {:12s} {:12s}".format("prediction", "model", "real"))
    for i in range(u.size()[0]):
        print(
            "{:12.4f} {:12.4f} {:12.4f}".format(qpred[i].item(), qmod[i].item(), qreal[i].item())
        )
    print("-- charge check in qps finished")

