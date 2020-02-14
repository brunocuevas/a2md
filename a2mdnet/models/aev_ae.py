import torch
import torch.nn as nn
import torchani
import numpy as np
from a2mdnet.modules import TorchElementSpecificA2MDNN
from a2mdnet import ALLOWED_SPECIES

DEVICE = torch.device('cpu')
FLOAT_DTYPE = torch.float32
DEFAULT_NODES = [384, 96, 16, 6]


class AevAE(nn.Module):

    def __init__(self, nodes=DEFAULT_NODES):

        super(AevAE, self).__init__()

        reversed_nodes = nodes.copy()
        reversed_nodes.reverse()

        self.encoder = TorchElementSpecificA2MDNN(
            nodes=nodes, elements=ALLOWED_SPECIES, device=DEVICE
        )

        self.decoder = TorchElementSpecificA2MDNN(
            nodes=reversed_nodes, elements=ALLOWED_SPECIES, device=DEVICE
        )

    def forward(self, *x):
        """

        :param x:
        :return:
        """

        labels, features = x
        labels, encoded_feats = self.encoder(labels, features)
        labels, decoded_feats = self.decoder(labels, encoded_feats)
        return labels, decoded_feats


if __name__ == "__main__":

    from torchani.neurochem import Constants
    from a2mdnet.data import Set500Energy
    from torch.utils import data
    # import matplotlib.pyplot as plt

    device_cpu = torch.device('cpu')
    device_gpu = torch.device('cuda:0')
    set500 = Set500Energy(device=DEVICE, dtype=FLOAT_DTYPE, normalize=True)
    set500_dl = data.DataLoader(set500, batch_size=64, num_workers=0, shuffle=True)

    epochs = 100
    lr = 1e-5

    cx = Constants(r'C:/Users/Bruno/ownCloud/main/a2mdnet/data/aev_params_simp.params')
    aev = torchani.AEVComputer(**cx)
    aev_ae = AevAE()
    history = np.zeros(epochs)
    optim = torch.optim.SGD(params=aev_ae.parameters(), lr=lr, momentum=0.9)

    for t in range(epochs):
        buffer = []
        for labels, coordinates, targets in set500_dl:
            _, aev_feats = aev.forward((labels, coordinates))
            _, output = aev_ae.forward(labels, aev_feats)
            l2 = torch.pow(output - aev_feats, 2.0).sum()

            l2.backward()
            optim.step()
            aev_ae.zero_grad()
            buffer.append(l2.item())

        history[t] = np.mean(buffer)
        print("epoch {:4d} : {:12.4e}".format(t, history[t]))
