import torch
import torch.nn as nn
from torchani.aev import radial_terms

class APEV(nn.Module):
    def __init__(self, output_size, parameters, device):
        super(APEV, self).__init__()

        self.parameters = parameters
        self.EtaR = torch.tensor(parameters['EtaR'], dtype=torch.float, device=device)
        self.ShfR = torch.tensor(parameters['ShfR'], dtype=torch.float, device=device)
        self.cutoff = parameters['Rc']
        self.output_size = output_size
        self.device = device

    def forward(self, *x):
        connectivity, coords = x

        assert isinstance(connectivity, torch.Tensor)
        assert isinstance(coords, torch.Tensor)

        n_connectivity = connectivity.size(1)
        n_batch = connectivity.size(0)
        n_atoms = coords.size(1)
        n_atoms_arange = torch.arange(n_batch, device=self.device) * n_atoms
        connectivity_f = connectivity.flatten(0, 1)
        coords_f = coords.flatten(0, 1)
        con_mask = connectivity_f[:, 0] != -1
        con_mask_coords = con_mask.unsqueeze(1).expand(connectivity_f.size(0), 3)

        connectivity_f = connectivity_f + (
            n_atoms_arange
        ).reshape(-1, 1).repeat(1, n_connectivity).reshape(-1, 1)

        coords_aceptor = torch.zeros(
            connectivity_f.size(0), 3, dtype=torch.float, device=self.device
        )
        coords_donnor = torch.zeros(
            connectivity_f.size(0), 3, dtype=torch.float, device=self.device
        )

        coords_aceptor.masked_scatter_(
            con_mask_coords,
            coords_f.index_select(dim=0, index=connectivity_f[con_mask, 0])
        )
        coords_donnor.masked_scatter_(
            con_mask_coords,
            coords_f.index_select(dim=0, index=connectivity_f[con_mask, 1])
        )

        distances = (coords_aceptor - coords_donnor).pow(2.0).sum(dim=1).sqrt()

        y = torch.zeros(connectivity_f.size(0), self.output_size, dtype=torch.float, device=self.device)

        rad_terms = radial_terms(
            Rcr=self.cutoff,
            EtaR=self.EtaR,
            ShfR=self.ShfR,
            distances=distances
        )

        mask_output = con_mask.unsqueeze(1).expand(connectivity_f.size(0), self.output_size)
        y.masked_scatter_(mask_output, rad_terms)
        y = y.reshape(n_batch, n_connectivity, self.output_size)

        return connectivity, y


def sample_cell(fun, min_r, max_r, n_initial, n_step, n_max, runs=10, var_treshold=0.1):
    """
    sample cell
    ---
    parameters
    - fun : function to evaluate
    - min_r : starting point of the cell
    - max_r : maximum point of the cell
    - n_initial : number of initial sampling points
    - n_step : increase number of initial sampling points
    - n_max : maximum number of steps
    - runs : number of parallel evaluations
    - var_treshold : treshold for convergence
    """
    import numpy as np
    dims = max_r - min_r
    volume = dims.prod()
    d = np.zeros(runs, dtype='float64')
    n = n_initial
    for i in range(runs):
        r = (np.random.rand(n_initial, 3) * dims) + min_r
        d[i] = fun(r).sum() * volume / n

    while n < n_max and d.var() > var_treshold:
        for i in range(runs):
            r = (np.random.rand(n_step, 3) * dims) + min_r
            f = fun(r).sum() * volume / n_step
            d[i] = (d[i] * (n / (n + n_step))) + ((n_step / (n + n_step)) * f)
        n = n + n_step

    return d.var() < var_treshold, n, d.mean()


def generate_voxel(r_min, r_max, resolution):
    import numpy as np
    steps = ((r_max - r_min) / resolution).astype('int64')
    range_x = np.arange(r_min[0], r_max[0], resolution)
    range_y = np.arange(r_min[1], r_max[1], resolution)
    range_z = np.arange(r_min[2], r_max[2], resolution)
    voxel_r_min = np.zeros(3, dtype="float64")
    for ix in range(steps[0]):
        voxel_r_min[0] = (ix * resolution) + r_min[0]
        for iy in range(steps[1]):
            voxel_r_min[1] = (iy * resolution) + r_min[1]
            for iz in range(steps[2]):
                voxel_r_min[2] = (iz * resolution) + r_min[2]
                yield voxel_r_min, voxel_r_min + resolution

