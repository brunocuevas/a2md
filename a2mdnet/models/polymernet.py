from torch import nn
from a2mdnet import modules


class PolymerNet(nn.Module):

    def __init__(self, architecture, device):
        super(PolymerNet, self).__init__()
        self.device = device

        self.net = modules.TorchElementSpecificA2MDNN(
            nodes=architecture, elements=[1, 6, 7, 8], device=self.device
        )

        self.featurizer = modules.SymFeats()

        for param in self.featurizer.parameters():
            param.requires_grad = False

        self.normalizer = modules.QPSegmentChargeNormalization(device=self.device)

    def forward(self, *x):

        l, t, x, s, sq, i_iso = x
        l, x = self.featurizer.forward(l, x)
        l, y = self.net.forward(l, x)
        y = self.normalizer.forward(
            y, i_iso,  s, sq
        )
        return y

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
