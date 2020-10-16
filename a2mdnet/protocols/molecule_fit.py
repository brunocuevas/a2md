from a2mdnet.modules import QMDensityBatch
from a2mdnet.utils import CoordinatesSampler
from a2mdnet.density_models import GenAMD

from a2mdio.params import AMDParameters
from a2mdio.molecules import Mol2
from a2mdtest.a2mdtests import benzene
from a2mdnet.data import convert_label2tensor
from a2mdnet.models.amd import AMDMolecule
from a2mdnet.utils import torch_eval_volume
import matplotlib.pyplot as plt
import torch
from torch.optim import SGD


if __name__ == '__main__':

    dtype = torch.float
    device = torch.device('cuda:0')
    benzene_mol2 = Mol2(benzene.mol2)
    molecule_coordinates = benzene_mol2.get_coordinates(units='au')
    molecule_charge = benzene_mol2.get_partial_charges()
    molecule_labels = convert_label2tensor(benzene_mol2.get_atomic_numbers(), device=device)
    proto_pars = AMDParameters.from_file(
        'C:/Users/Bruno/ownCloud/projects/a2md/a2mdnet/a2mdt/test/spherical.json'
    )
    mol_pars = AMDParameters.from_file(
        'C:/Users/Bruno/ownCloud/projects/a2md/a2mdnet/a2mdt/test/extension.json'
    )

    qm = QMDensityBatch(
        filename=benzene.path / "benzene.wfn.hdf5", index=['benzene.wfn'], dtype=dtype, device=device
    )

    cs = CoordinatesSampler(
        sampler='spheres',
        sampler_args=dict(resolution=10, max_radius=15.0),
        dtype=dtype, device=device
    )

    proto_amd = GenAMD(proto_pars, device=torch.device('cuda:0'), dtype=torch.float)

    benzene_t = AMDMolecule(
        coordinates=molecule_coordinates, labels=molecule_labels,
        parameters=mol_pars, charge=molecule_charge, device=device, dtype=dtype
    )

    # dx = torch_eval_volume(
    #     lambda x: proto_amd.protodensity(x, benzene_t.labels.unsqueeze(0), benzene_t.coordinates.unsqueeze(0)),
    #     resolution=0.25, steps=50, device=torch.device('cuda:0')
    # )
    # dx.write('C:/scratch/proto.dx')
    sgd = SGD(params=benzene_t.parameters(), lr=1e-4)
    for i in range(100):

        sample = cs(benzene_t.coordinates.unsqueeze(0))
        dens = qm.forward(torch.tensor([0]), sample)
        proto = proto_amd.protodensity(
            sample, benzene_t.labels.unsqueeze(0), benzene_t.coordinates.unsqueeze(0)
        )
        defpred = benzene_t.forward(sample)
        l2 = dens - proto - defpred
        l2 = l2.pow(2.0)
        l2 = (l2 / dens).sum()
        l2 = l2.sum() + (benzene_t.coefficients ** 2).sum() * 1e-2
        print(l2.item())
        l2.backward()
        sgd.step()
        benzene_t.zero_grad()

    dx = torch_eval_volume(
        lambda x: benzene_t.forward(x),
        resolution=0.25, steps=50, device=torch.device('cuda:0')
    )
    dx.write('C:/scratch/mol.dx')
