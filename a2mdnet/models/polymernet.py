from torch import nn
from a2mdnet.modules import TorchElementSpecificA2MDNN, SymFeats
from a2mdnet.modules import QPSegmentChargeNormalization
from a2mdnet.density_models import GenAMD
from a2mdio.params import AMDParameters
from typing import List, Dict
import torch


class PolymerNet(nn.Module):

    def __init__(
        self, layers: List[int], density_model: AMDParameters, symfeats_model: Dict,
        device: torch.device, dtype: torch.dtype
    ):
        """
        A2MDnet/monomer

        reproduces an electron density model from molecular coordinates using a generalized
        amd model. Charges are restricted by segments, so it can be applied to polymers.
        Otherwise, we recommend to use Monomernet.

        :param layers: list providing the architecture of the neural network
        :param density_model: AMDParameters instance containing AMD parameters
        :param symfeats_model: TorchANI constants for symfeats
        :param device: either cpu of cuda:0
        :param dtype: please, use torch.float

        Example:
            pars = AMDParameters('extended_valence.json')
            symfeat = Constants('monomer_aev')
            mnet = Polymernet(
                [110, 64, 64, 12], pars, symfeat, device=torch.device('cuda:0'), dtype=torch.float
            )

        """
        super(PolymerNet, self).__init__()
        self.device = device
        self.feats = SymFeats(symfeats_model).to(device)
        self.net = TorchElementSpecificA2MDNN(nodes=layers, elements=[1, 6, 7, 8, 16], device=device)
        self.norm = QPSegmentChargeNormalization(device=device)
        self.density = GenAMD(density_model, device=device, dtype=dtype)
        self.dtype = dtype
        self.device = device

        for param in self.featurizer.parameters():
            param.requires_grad = False

    def forward(
            self, coordinates: torch.Tensor, labels: torch.Tensor,
            mol_coordinates: torch.Tensor, segment_charge: torch.Tensor,
            segment: torch.Tensor
    ):
        """

        calculates electron density at given coordinates (au)

        :param coordinates: Coordinates at which to calculate the electron density values. Use Bohrs
        :param labels: Element labels of the molecule. Use ELEMENT2NN
        :param mol_coordinates: Molecule oordinates. Use Bohrs
        :param segment_charge: Atomic charges. Use electrons
        :param segment: Molecule monomer index.
        :returns density: in electron / (bohr ^ 3)
        """
        _, symfeats = self.feats(labels, mol_coordinates)
        _, unnormcoeff = self.net(labels, symfeats)
        integrals = self.density.integrate(labels)
        normcoeff = self.norm.forward(unnormcoeff, integrals, segment, segment_charge)
        density = self.density.forward(coordinates, normcoeff, labels, mol_coordinates)
        return density

    def to(self, new_device, **kwargs):
        # Following https://stackoverflow.com/questions/54706146/moving-member-tensors-with-module-to-in-pytorch
        new_self = super(PolymerNet, self).to(new_device)
        new_self.device = new_device
        new_self.normalizer.device = new_device
        new_self.featurizer.device = new_device
        new_self.atom_net.device = new_device
        new_self.bond_net.device = new_device
        new_self.common_net.device = new_device
        new_self.density_model.device = new_device
        return new_self
