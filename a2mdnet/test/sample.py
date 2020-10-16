import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import torch
from a2mdnet.utils import CoordinatesSampler
from a2mdnet.data import MonomerDataset
from torch.utils import data


if __name__ == '__main__':

    print('-- coordinate sampler')

    cs = CoordinatesSampler(
        device=torch.device('cpu'), dtype=torch.float, sampler='spheres',
        sampler_args=dict(max_radius=2.0, resolution=2)
    )

    r = torch.rand(5, 10, 3) - 0.5
    r *= 10.0
    u = cs(r)

    u = u.reshape(-1, 3).data.numpy()
    r = r.reshape(-1, 3).data.numpy()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(r[:, 0], r[:, 1], r[:, 2], color='red')
    ax.scatter(u[:, 0], u[:, 1], u[:, 2], color='blue')
    plt.show()
    print('-- done!')
    print('-- starting qm sampling')
    with open('F:/edip/extended/ata500/.ata.index') as f:
        index = [i.strip() for i in f.readlines()]

    md = MonomerDataset(
        device=torch.device('cpu'), dtype=torch.float,
        ids=index, molecular_data_path='F:/edip/extended/ata500/'
    )
    mddl = data.DataLoader(md, batch_size=32, shuffle=True)

    for index, labels, topo, coords, charge in mddl:
        sample = cs(coords)
        break

    u = sample.reshape(-1, 3).data.numpy()
    r = coords.reshape(-1, 3).data.numpy()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(r[:, 0], r[:, 1], r[:, 2], color='red')
    ax.scatter(u[:, 0], u[:, 1], u[:, 2], color='blue', alpha=0.1)
    plt.show()
    print('-- done!')
