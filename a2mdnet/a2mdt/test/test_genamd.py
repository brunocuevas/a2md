import torch
from a2mdnet.density_models import GenAMD, genamd_from_mol2
from a2mdio.params import AMDParameters
from a2mdio.molecules import Mol2
from a2mdnet import ALLOWED_SPECIES
from a2mdtest.a2mdtests import benzene
from a2md.models import a2md_from_mol
import json
import numpy as np
if __name__ == '__main__':

    print('-- gen amd')
    print('-- testing generalized amd to reproduce integrals')

    amdp = AMDParameters.from_file('spherical.json')
    gamd = GenAMD(amdp, device=torch.device('cuda:0'), dtype=torch.float)

    labels = torch.tensor([[0, 1, 2, 3, 4]], dtype=torch.long, device=torch.device('cuda:0'))

    ints = gamd.integrate(labels=labels)
    element_ints = ints.sum(2)

    labels = labels.data.cpu().numpy().reshape(-1)
    element_ints = element_ints.data.cpu().numpy().reshape(-1)

    for lab, int in zip(labels, element_ints):
        symbol = ALLOWED_SPECIES[lab]
        print("{:12d} {:8.4f}".format(symbol, int))

    print("-- testing exponential kernel functions")
    labels = torch.tensor([[0, 1, 2, 3, 4], [0, 1, 2, 3, 4]], dtype=torch.long, device=torch.device('cuda:0'))
    x = torch.rand((2, 1000, 3), dtype=torch.float, device=torch.device('cuda:0'))
    r = torch.rand((2, 5, 3), dtype=torch.float, device=torch.device('cuda:0'))
    c = torch.ones((2, 5, 4), dtype=torch.float, device=torch.device('cuda:0'))
    p = gamd.forward(x, c, labels, r)
    print("-- done")

    print("-- comparing amd vs genamd")
    benzene_mol2 = Mol2(benzene.mol2)
    benzene_amd = a2md_from_mol(benzene_mol2)
    with open('spherical.json') as f:
        spherical = json.load(f)
    benzene_amd.parametrize(param_dict=spherical)
    print('-- charge = {:6.4f}'.format(benzene_amd.integrate()))
    dx = benzene_amd.eval_volume(spacing=3.0, resolution=0.25, kind='density')
    dx.write('benzene.dx')
    print('-- writting benzene.dx')
    r = (np.random.rand(1000, 3) - 0.5) * 10
    l, x, c = genamd_from_mol2(benzene_mol2, device=torch.device('cuda:0'))
    r_torch = torch.tensor(r, device=torch.device('cuda:0'), dtype=torch.float).unsqueeze(0)
    p_ref = benzene_amd.eval(r)
    # p_val = gamd.forward(r_torch, c, l, x).data.cpu().numpy()
    p_val = gamd.protodensity(r_torch, l, x).data.cpu().numpy()

    l2 = np.power((p_ref - p_val), 2.0).sum()
    print('-- l2 = {:8.4f}'.format(l2))
    if l2 < 1e-4:
        print(" -- succesful test!")
