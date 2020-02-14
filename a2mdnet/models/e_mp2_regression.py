from a2mdnet.modules import TorchaniFeats, TorchElementSpecificA2MDNN
from a2mdnet.data import Set500Energy
import torch
from torch import nn
from torch.utils import data
import numpy as np
# -------------------------------------------------------------------------------------------------------------------- #
#   MP2 Energy Regression Task
#   Model: 1. TorchANI as feature extractor
#          2. Element specific MLP
# -------------------------------------------------------------------------------------------------------------------- #

LR = 1e-5
NODES = [96, 6, 1]
ELEMENTS = [1, 6, 7, 8]
EPOCHS = 10
BATCH_SIZE = 32
MOMENTUM = 0.9
DEVICE = torch.device('cpu')
MODEL = r'C:\Users\Bruno\ownCloud\main\a2mdnet\models\emp2r.pt'


class EMP2r(nn.Module):
    def __init__(self, architecture):
        """

        :param architecture:
        """
        super(EMP2r, self).__init__()
        self.fe = TorchaniFeats(device=DEVICE)
        for param in self.fe.parameters():
            param.requires_grad = False
        self.reg = TorchElementSpecificA2MDNN(nodes=architecture, elements=[1, 6, 7, 8], device=DEVICE)

    def forward(self, *x):
        """

        :param x:
        :return:
        """
        labels = x[0]
        coordinates = x[1]
        pred_buffer = []

        ta_feats = self.fe.forward(labels, coordinates)
        _, pred = self.reg.forward(labels, ta_feats)

        pred_buffer.append(pred.sum(dim=(1, 2)))

        return torch.stack(pred_buffer).reshape(-1, 1)


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    print("EMP2r")

    set500 = Set500Energy(device=DEVICE, dtype=torch.float, normalize=True)
    set500_size = len(set500)
    data_gen = data.DataLoader(set500,  batch_size=64, shuffle=True, num_workers=0)
    model = EMP2r(architecture=NODES)
    history = np.zeros(EPOCHS, dtype='float64')

    opt = torch.optim.SGD(
        params=model.parameters(),
        lr=LR, momentum=MOMENTUM
    )

    targets = None
    output = None

    for t in range(EPOCHS):
        buffer = []
        for labels, coords, targets in data_gen:
            output = model.forward(labels, coords)
            output = output.reshape(-1)
            l2 = torch.pow(output - targets, 2).sum()
            buffer.append(l2.item())

            l2.backward()
            opt.step()
            model.zero_grad()

        history[t] = np.mean(buffer)
        print("{:>4d}, {:<12.4e}".format(t, history[t]))

    fig, ax = plt.subplots(1, 2)
    ax[0].plot(np.arange(EPOCHS), history, marker='.')
    ax[1].scatter(targets.data, output.data)

    plt.show()

    # torch.save(model, MODEL)
