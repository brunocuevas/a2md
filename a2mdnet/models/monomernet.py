import torch
from torch import nn
from a2mdnet.modules import TorchElementSpecificA2MDNN, SymFeats
from a2mdnet.modules import QPChargeNormalization
from a2mdnet.density_models import GenAMD
from a2mdio.params import AMDParameters
from typing import List, Dict


class Monomernet(nn.Module):

    def __init__(
        self, layers: List[int], density_model: AMDParameters, symfeats_model: Dict,
        device: torch.device, dtype: torch.dtype
    ):
        """
        A2MDnet/monomer

        reproduces an electron density model from molecular coordinates using a generalized
        amd model.

        :param layers: list providing the architecture of the neural network
        :param density_model: AMDParameters instance containing AMD parameters
        :param symfeats_model: TorchANI constants for symfeats
        :param device: either cpu of cuda:0
        :param dtype: please, use torch.float

        Example:
            pars = AMDParameters('extended_valence.json')
            symfeat = Constants('monomer_aev')
            mnet = Monomernet(
                [110, 64, 64, 12], pars, symfeat, device=torch.device('cuda:0'), dtype=torch.float
            )

        """
        super(Monomernet, self).__init__()
        self.net = TorchElementSpecificA2MDNN(nodes=layers, elements=[1, 6, 7, 8, 16], device=device)
        self.feats = SymFeats(symfeats_model).to(device)
        self.norm = QPChargeNormalization(device=device)
        self.density = GenAMD(density_model, device=device, dtype=dtype)
        self.dtype = dtype
        self.device = device

    def forward(
            self, coordinates: torch.Tensor, labels: torch.Tensor,
            mol_coordinates: torch.Tensor, charge: torch.Tensor
    ) -> torch.Tensor:
        """

        calculates electron density at given coordinates (au)

        :param coordinates: Coordinates at which to calculate the electron density values. Use Bohrs
        :param labels: Element labels of the molecule. Use ELEMENT2NN
        :param mol_coordinates: Molecule oordinates. Use Bohrs
        :param charge: Atomic charges. Use electrons
        :returns density: in electron / (bohr ^ 3)
        """
        _, symfeats = self.feats(labels, mol_coordinates)
        _, unnormcoeff = self.net(labels, symfeats)
        integrals = self.density.integrate(labels)
        normcoeff = self.norm.forward_iso(unnormcoeff, integrals, charge)
        density = self.density.forward(coordinates, normcoeff, labels, mol_coordinates)
        return density
