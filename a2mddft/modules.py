import torch
import numpy as np
from torch import nn
from torch.utils import data
from typing import List
from a2mdio.qm import WaveFunction
from a2mdio.qm import symetry_index
from a2md import LEBEDEV_DESIGN

class WaveFunctionReader:
    def __init__(self, idx : List[str]):

        self.idx = idx
        self.n = len(idx)


    def __len__(self):

        return len(self.idx)

    def __getitem__(self, item):

        wfn = WaveFunction(file=self.idx[item], verbose=False)

        centers = wfn.cent.tolist()
        sym = wfn.sym.tolist()
        exp = wfn.exp.tolist()
        r = torch.tensor(wfn.get_coordinates(), dtype=torch.float)

        dm = torch.tensor(wfn.density_matrix, dtype=torch.float)

        return r, centers, exp, sym, dm

    def generate_wfns(self, permute=False):
        order = np.arange(self.n)
        if permute:
            order = np.random.permutation(order)

        for i in order:
            yield self.__getitem__(i)


class DensityFunction(nn.Module):
    def __init__(self, device):

        super(DensityFunction, self).__init__()
        self.device = device
        self.r = None
        self.centers = None
        self.exp = None
        self.sym = None
        self.dm = None
        self.coords = None
        self.n = None

    @staticmethod
    def calculate_distance_vector(x, r):
        """
        calculates |x - r| distance
        :param x:
        :param r:
        :return:
        """

        x = x.clone()
        r = r.clone().unsqueeze(0).expand(x.size()[0], 3)
        return x - r

    @staticmethod
    def calculate_distance2(v):
        """

        :param v:
        :return:
        """
        return torch.pow(v, 2.0).sum(1)


    @staticmethod
    def gaussian(d2, alpha):
        """

        :param d2:
        :param alpha:
        :return:
        """

        return torch.exp(-alpha * d2)

    @staticmethod
    def symmetry(v, s):
        """

        :param v:
        :param s:
        :return:
        """
        lx = torch.pow(v[:, 0], symetry_index[s][0])
        ly = torch.pow(v[:, 1], symetry_index[s][1])
        lz = torch.pow(v[:, 2], symetry_index[s][2])
        return lx * ly * lz

    def calculate_density(self, coords):
        """

        :param coords:
        :return:
        """
        nprims = len(self.exp)
        ncoords = coords.size()[0]
        x_buffer = torch.zeros([nprims, ncoords], dtype=torch.float, device=self.device)
        p = torch.zeros(ncoords, dtype=torch.float, device=self.device)
        for i in range(nprims):
            c = self.centers[i]
            v = self.calculate_distance_vector(coords, self.coords[c, :])
            d = self.calculate_distance2(v)
            sh = self.symmetry(v, self.sym[i])
            g = self.gaussian(d, self.exp[i])
            x_buffer[i, :] = sh * g

        for i in range(nprims):
            for j in range(i):
                p += self.dm[i, j] * x_buffer[i, :] * x_buffer[j, :] * 2.0
            p += self.dm[i, i] * x_buffer[i, :] * x_buffer[i, :]
        return p

    def calculate_density_gradient(self, coords):
        """

        :param coords:
        :return:
        """
        dh = 1e-8
        grad = torch.zeros(coords.size()[0], dtype=torch.float, device=self.device)
        for i in range(3):
            coords_buffer = coords.clone()
            coords_buffer[:, i] = coords_buffer[:, i] + dh
            dens_up = self.calculate_density(coords_buffer)
            coords_buffer[:, i] = coords_buffer[:, i] - 2 * dh
            dens_down = self.calculate_density(coords_buffer)
            grad += torch.pow((dens_up - dens_down) / (2.0 * dh), 2.0)
        return torch.sqrt(grad)

    def parametrize(
            self,
            coords : torch.Tensor,
            centers : List[int], exp : torch.Tensor,
            sym : List[List[float]], dm : torch.Tensor
    ):
        """

        :param coords:
        :param centers:
        :param exp:
        :param sym:
        :param dm:
        :return:
        """
        self.coords = coords.to(device=self.device)
        self.centers = centers
        self.exp = exp
        self.sym = sym
        self.dm = dm
        self.n = len(centers)

    def forward(self, coords : torch.Tensor):
        """

        :param coords:
        :return:
        """
        return self.calculate_density(coords), self.calculate_density_gradient(coords)


class MLDensityFunctional(nn.Module):

    def __init__(self, architecture, device):
        super(MLDensityFunctional, self).__init__()
        layers = []
        self.device=device
        for i in range(len(architecture) - 1):
            layers.append(
                nn.Linear(architecture[i], architecture[i + 1])
            )
            layers.append(
                nn.CELU()
            )
        self.network = nn.ModuleList(layers).to(self.device)

    @staticmethod
    def tesellate(coords, ref):

        d = np.zeros((ref.shape[0], coords.shape[0]), dtype='float64')

        for i in range(ref.shape[0]):

            d[i, :] = np.linalg.norm(coords - ref[i, :], axis=1)
        d = np.argmin(d.T, axis=1)
        return d

    def generate_lebdenev_grid(self, ref):
        """

        :return:
        """
        ref = ref.data.numpy()
        n = ref.shape[0]
        ld = np.loadtxt(LEBEDEV_DESIGN['medium'])

        lebedev = np.zeros((ld.shape[0], 3), dtype='float64')
        lebedev[:, 0] = np.cos(np.deg2rad(ld[:, 0])) * np.sin(np.deg2rad(ld[:, 1]))
        lebedev[:, 1] = np.sin(np.deg2rad(ld[:, 0])) * np.sin(np.deg2rad(ld[:, 1]))
        lebedev[:, 2] = np.cos(np.deg2rad(ld[:, 1]))
        res = 150
        ld_w = ld[:, 2]
        u = np.log(10.0 + 1) / res
        r_grid = np.exp(np.arange(1, res+ 1) * u) - 1

        w = []
        coords = []
        for i, r in enumerate(r_grid):
            dr = r - (np.exp( u * i ) - 1 )
            dv = ((r ** 3) / 3) - (((r - dr) ** 3) / 3)
            coords.append(r * lebedev)
            w.append(ld_w * 4.0 * np.pi * dv)

        w = np.concatenate(w)
        coords = np.concatenate(coords)
        final_coords = []
        final_w = []
        for i in range(n):
            r0 = ref[i, :]
            coords_ = coords + r0
            vor = self.tesellate(coords_, ref)
            final_coords.append(coords_[vor == i, :])
            final_w.append(w[vor == i])

        final_coords = np.concatenate(final_coords)
        final_w = np.concatenate(final_w)
        final_coords = torch.tensor(final_coords, device=self.device, dtype=torch.float)
        final_w = torch.tensor(final_w, device=self.device, dtype=torch.float)
        return final_coords, final_w


    def forward(self, centers : torch.Tensor, density_function : nn.Module):
        """
        :param centers
        :param density_function:
        :return:
        """
        grid, w = self.generate_lebdenev_grid()
        p, dp, lp, zp = density_function.forward(grid)
        pp = torch.cat([p, dp, lp, zp], dim=1)
        df = self.network.forward(pp)
        dfi = (w * df).sum(1)
        return dfi