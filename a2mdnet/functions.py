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