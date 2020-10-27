from a2mdio.molecules import Mol2
from a2mdnet.data import convert_label2tensor
from a2md.utils import integrate_from_dict
from a2mdnet.data import match_fun_names, Coordinates
from a2mdio.utils import eval_volume
from a2mdio.units import bohr, Unit
from a2md import LEBEDEV_DESIGN
import warnings
import torch
from typing import Callable, Dict, Tuple, List
import math
import numpy as np


def get_charges(l, t, i_iso, i_aniso, iso, aniso, device):
    """

    :param l:
    :param t:
    :param i_iso:
    :param i_aniso:
    :param iso:
    :param aniso:
    :param device:
    :return:
    """

    # Defining problems
    n_batch = l.size(0)
    n_atoms = l.size(1)
    n_bond = t.size(1)
    n_bond_funs = int(i_aniso.size(2) / 2)

    # Calculating function charges
    charge_iso = (i_iso * iso)
    charge_aniso = (i_aniso * aniso)
    q_iso = charge_iso.sum(2)
    # Calculating charges per atom
    # To do so, the topology tensor is expanded to match the shape of the anisotric charge tensor
    # Then a range tensor is added to allow to map index in flat vectors
    # Finally, two operations are employed to map back the anisotropic charge to each atom
    #       1. scatter add. Sums all the positions related to the same index and stores them in that given index
    #       2. masked_scatter. Gets all the non zero values and stores them in the output

    t_r = t.reshape(n_batch, n_bond_funs * n_bond, 1).expand(n_batch, n_bond_funs * n_bond, n_bond_funs).reshape(
        n_batch, n_bond, 2 * n_bond_funs)
    r = (torch.arange(0, n_batch, device=device) * n_atoms).unsqueeze(1).unsqueeze(2).expand(n_batch, n_bond,
                                                                                             2 * n_bond_funs)
    t_r.masked_scatter_(t_r != -1, t_r.masked_select(t_r != -1) + r.masked_select(t_r != -1))
    t_r[t_r == -1] = 0
    q_buffer = torch.zeros_like(charge_aniso).flatten()
    t_rf = t_r.flatten()
    ca_f = charge_aniso.flatten()
    q_buffer = q_buffer.scatter_add(0, t_rf, ca_f)
    q_aniso = torch.zeros_like(q_iso)
    q_aniso = q_aniso.masked_scatter(l != - 1, q_buffer.masked_select(q_buffer != 0.0))

    q_per_atom = q_aniso + q_iso
    return q_per_atom


class Parametrizer:

    def __init__(self, model, device):
        self.model = model
        self.device = device

    def parametrize(self, mol, params):
        """

        :param mol:
        :param params
        :return:
        """
        if not isinstance(mol, Mol2):
            try:
                mol = Mol2(file=mol)
            except ValueError:
                raise IOError("unknown format for mol. Please, use Mol2 or str")
        coords = torch.tensor(mol.get_coordinates(), device=self.device, dtype=torch.float).unsqueeze(0)
        labels = convert_label2tensor(mol.get_atomic_numbers(), device=self.device).unsqueeze(0)
        connectivity = torch.tensor(mol.get_bonds(), dtype=torch.long, device=self.device).unsqueeze(0)
        charge = torch.tensor(mol.charges + mol.atomic_numbers, dtype=torch.float, device=self.device).unsqueeze(0)
        natoms = mol.get_number_atoms()
        nbonds = mol.get_number_bonds()
        connectivity -= 1

        int_iso = torch.zeros(1, natoms, 2)
        int_aniso = torch.zeros(1, nbonds, 4)

        for fun in params:
            center = fun['center']
            funtype, pos = match_fun_names(fun)
            if funtype == 'core':
                charge[0, center] -= integrate_from_dict(fun)
            elif funtype in 'iso':
                int_iso[0, center, pos] = integrate_from_dict(fun)
            elif funtype in 'aniso':
                idx, col = self.match_bond(connectivity, fun)
                int_aniso[0, idx, col + pos] = integrate_from_dict(fun)

        int_iso = int_iso.to(self.device)
        int_aniso = int_aniso.to(self.device)

        _, _, iso_out, aniso_out = self.model.forward_coefficients(labels, connectivity, coords, charge, int_iso,
                                                                   int_aniso)

        if iso_out.is_cuda:
            iso_out = iso_out.squeeze(0).cpu().data.numpy()
            aniso_out = aniso_out.squeeze(0).cpu().data.numpy()
        else:
            iso_out = iso_out.squeeze(0).data.numpy()
            aniso_out = aniso_out.squeeze(0).data.numpy()

        for fun in params:
            center = fun['center']
            funtype, pos = match_fun_names(fun)
            if funtype == 'core':
                continue
            if funtype == 'iso':
                fun['coefficient'] = iso_out[center, pos].item()
            elif funtype == 'aniso':
                idx, col = self.match_bond(connectivity, fun)
                fun['coefficient'] = aniso_out[idx, col + pos].item()

        return params

    @staticmethod
    def match_bond(topology, fun):
        topology = topology.squeeze(0)
        for j in range(topology.shape[0]):
            c1 = (topology[j, 0] == fun['center']) and (
                    topology[j, 1] == fun['bond']
            )
            c2 = (topology[j, 1] == fun['center']) and (
                    topology[j, 0] == fun['bond']
            )
            if c1:
                return j, 0
            if c2:
                return j, 0 + 2

        warnings.warn("bond was not matched")


class CoordinatesSampler:
    def __init__(
            self, device: torch.device, dtype: torch.dtype,
            sampler: str, sampler_args: Dict, units: Unit = bohr
    ):
        """
        coordinates sampler
        ---
        performs a 3d coordinates sample given some initical coordinates
        """

        self.internal_methods = dict(
            random=CoordinatesSampler.random_box,
            spheres=CoordinatesSampler.spheres,
            box=CoordinatesSampler.regular_grid
        )
        self.sampler = self.internal_methods[sampler]
        self.sampler_args = sampler_args
        self.device = device
        self.dtype = dtype
        self.units = units

    @staticmethod
    def principal_components(coords: torch.Tensor):
        mean = coords.mean(1, keepdim=True)
        x = coords - mean
        n = coords.size()[1]
        c = x.transpose(1, 2) @ x / n
        eig, eiv = torch.symeig(c, eigenvectors=True)
        eivp = eiv.inverse()
        x = (eivp @ x.transpose(1, 2)).transpose(1, 2)
        # x += mean
        return x, eiv, mean

    @staticmethod
    def random_box(coords, device, dtype, n_sample=1000, spacing=6.0):
        r = torch.rand(coords.size()[0], n_sample, 3, device=device, dtype=dtype)
        rotcoords, eivp, mean = CoordinatesSampler.principal_components(coords)
        box_min = rotcoords.min(1, keepdim=True)[0] - spacing
        box_max = rotcoords.max(1, keepdim=True)[0] + spacing
        diff = box_max + (box_min * -1)
        r = (r * diff) + box_min
        r = (eivp @ r.transpose(1, 2)).transpose(1, 2)
        r += mean
        return r

    @staticmethod
    def spheres(coords, device, dtype, grid='coarse', resolution=10, max_radius=10.0):
        import numpy as np
        n = coords.size()[0]
        m = coords.size()[1]

        design = np.loadtxt(LEBEDEV_DESIGN[grid])
        design = torch.tensor(design, device=device, dtype=dtype)
        phi, psi, _ = design.split(1, dim=1)
        phi = (phi / 180.0) * math.pi
        psi = (psi / 180.0) * math.pi
        design_size = phi.size(0)
        sphere_x = (phi.cos() * psi.sin()).flatten()
        sphere_y = (phi.sin() * psi.sin()).flatten()
        sphere_z = (psi.cos()).flatten()

        sphere = torch.stack([sphere_x, sphere_y, sphere_z], dim=1)
        radius = torch.arange(resolution, dtype=torch.float, device=device) / resolution
        radius = 0.1 + (radius * (max_radius - 0.1))
        n_spheres = radius.size()[0]
        radius = radius.unsqueeze(1).unsqueeze(2)

        sphere = sphere.unsqueeze(0).expand(n_spheres, design_size, 3)
        concentric_spheres = radius * sphere
        concentric_spheres = concentric_spheres.reshape(-1, 3)
        concentric_spheres = concentric_spheres.unsqueeze(0)

        ms = concentric_spheres.size()[1]

        coords = coords.reshape(-1, 3)
        coords_x, coords_y, coords_z = torch.split(coords, 1, dim=1)
        mask_x = coords_x == 0
        mask_y = coords_y == 0
        mask_z = coords_z == 0
        mask = mask_x * mask_y * mask_z
        n_zeros = mask.sum()
        mean_coords = coords.mean(dim=0)
        std_coords = coords.std(dim=0)
        off_centers = mean_coords + (torch.randn((n_zeros, 3), dtype=dtype, device=device) * std_coords)

        coords = torch.masked_scatter(coords, mask, off_centers)
        coords = coords.unsqueeze(1)
        coords = coords.expand(n * m, ms, 3)

        concentric_spheres = concentric_spheres + coords
        concentric_spheres = concentric_spheres.reshape(n, m * ms, 3)
        return concentric_spheres

    @staticmethod
    def regular_grid(coords, device, dtype, grid='medium', spacing=2.0):
        sizes = dict(
            coarse=10, medium=20, fine=30
        )
        try:
            size = sizes[grid]
        except KeyError:
            raise IOError('use either coarse, medium or fine')

        n = coords.size()[0]
        x = torch.arange(size, dtype=dtype, device=device) / float(size)
        x, y, z = torch.meshgrid([x, x, x])
        r = torch.stack([x.flatten(), y.flatten(), z.flatten()], dim=1)
        r = r.unsqueeze(0).expand([n, size ** 3, 3])
        rotcoords, eivp, mean = CoordinatesSampler.principal_components(coords)
        box_min = rotcoords.min(1, keepdim=True)[0] - spacing
        box_max = rotcoords.max(1, keepdim=True)[0] + spacing
        diff = (box_max + (box_min * -1))
        r = (r * diff) + box_min
        r = (eivp @ r.transpose(1, 2)).transpose(1, 2)
        r += mean
        return r

    def __call__(self, coords: Coordinates, *args) -> torch.Tensor:

        coords = coords.get_cartessian(self.units)

        r = self.sampler(
            coords=coords, device=self.device, dtype=self.dtype, *args, **self.sampler_args
        )

        return Coordinates(r, units=self.units)


def torch_eval_volume(fun: Callable, resolution: float, steps: int, device: torch.device):
    modfun = lambda x: fun(
        torch.tensor(x, device=device, dtype=torch.float).unsqueeze(0)).data.cpu().numpy()
    dx = eval_volume(modfun, resolution, steps, shift=[0, 0, 0])
    return dx


class IntegrationGrid:

    def __init__(
            self, device: torch.device, dtype: torch.dtype, grid: str = 'coarse',
            radial_resolution: int = 15,
            units: Unit = bohr, softening: int = 3, rm: float = 5.0
    ):
        """
        Integration grid
        Generates a integration grid suitable for electron density problems. To do so,
        concentrical Lebdenev spheres are placed on top of the molecule coordinates, using
        radius from Gauss-Chebysev quadrature; then weights are calculated using the
        elliptic coordinates suggested by Becke (JCP, 1988).

        :param device:
        :param dtype:
        :param grid: use either coarse, medium or tight
        :param radial_resolution: number of bins. Values larger than 15 are ok
        :param units: use a unit with a bohr reference unit
        :param softening: number of softening passes
        :param rm: middle point for radiual integration. In coordinates units

        """
        self.device = device
        self.dtype = dtype
        self.sphere, self.sphere_weights, self.design_size = self.load_design(grid=grid)
        self.units = units
        self.radial_resolution = radial_resolution
        self.softening = softening
        self.rm = rm

    def load_design(self, grid) -> Tuple[torch.Tensor, torch.Tensor, int]:
        """
        Loads one of the precalculated designs for lebdenev spheres
        :param grid:
        :return:
        """
        import numpy as np
        design = np.loadtxt(LEBEDEV_DESIGN[grid])
        design = torch.tensor(design, device=self.device, dtype=self.dtype)
        phi, psi, weights = design.split(1, dim=1)
        phi = (phi / 180.0) * math.pi
        psi = (psi / 180.0) * math.pi
        design_size = phi.size(0)
        sphere_x = (phi.cos() * psi.sin()).flatten()
        sphere_y = (phi.sin() * psi.sin()).flatten()
        sphere_z = (psi.cos()).flatten()
        weights = weights.flatten()
        sphere = torch.stack([sphere_x, sphere_y, sphere_z], dim=1)
        return sphere, weights, design_size

    def integration_grid(self, coords: Coordinates) -> Tuple[Coordinates, torch.Tensor]:
        """
        Generates grid on top of the coordinates
        NOTE: There is no implementation right now to avoid the padding problem!
        :param coords: Coordinates
        :return:
        """

        dm = coords.get_distance_matrix(self.units)
        coords = coords.get_cartessian(self.units)

        n = coords.size()[0]
        m = coords.size()[1]

        # Defining spherical grids
        ms, concentric_spheres, weights = self.spheric_integration_grid()
        # Placing spheres at molecule coordinates
        coords_ = coords.reshape(-1, 3)
        coords_ = coords_.unsqueeze(1)
        coords_ = coords_.expand(n * m, ms, 3)
        concentric_spheres = concentric_spheres + coords_
        concentric_spheres = concentric_spheres.reshape(n, m * ms, 3)
        weights = weights.unsqueeze(1).expand(1, m, ms).reshape(1, m * ms)
        # Generating list of sampling centers
        sampling_centers = torch.arange(
            m, device=self.device, dtype=torch.long
        ).unsqueeze(0).unsqueeze(2).expand(n, m, ms).reshape(n, m * ms)
        # Tesellation
        v = self.becke_tesellation(
            x=concentric_spheres, r=coords, dm=dm, i=sampling_centers
        )
        w = v * weights

        concentric_spheres = Coordinates(values=concentric_spheres, units=self.units)

        return concentric_spheres, w

    def spheric_integration_grid(
            self
    ) -> Tuple[int, torch.Tensor, torch.Tensor]:
        """
        Generates concentric Lebdenev-Gauss Chebysev spheres
        NOTE: There is no implementation right now to avoid the padding problem!
        :return:
        """

        weights = self.sphere_weights.clone()
        sphere = self.sphere.clone()

        # Generating the radial component
        # Using Gauss-Chebysev with variable change
        #
        #       r = rm ( 1 + x ) / ( 1 - x )
        #

        i = torch.arange(1, self.radial_resolution + 1, dtype=self.dtype, device=self.device)
        z = - (math.pi * ((2.0 * i) - 1.0) / (2.0 * self.radial_resolution)).cos()
        dr = 2.0 * self.rm * torch.pow(1 - z, -2.0)
        r = self.rm * (1 + z) / (1 - z)
        w = torch.sqrt(1 - z.pow(2.0)) * dr * math.pi / self.radial_resolution
        w = r.pow(2.0) * 4.0 * math.pi * w

        # Stacking concentric spheres
        n_spheres = r.size()[0]
        r = r.unsqueeze(1).unsqueeze(2)
        sphere = sphere.unsqueeze(0).expand(n_spheres, self.design_size, 3)
        weights = weights.unsqueeze(0).expand(n_spheres, self.design_size)
        concentric_spheres = r * sphere
        concentric_spheres = concentric_spheres.reshape(-1, 3)
        concentric_spheres = concentric_spheres.unsqueeze(0)
        thinness = w.unsqueeze(1)
        weights = weights * thinness
        weights = weights.reshape(-1).unsqueeze(0)
        ms = concentric_spheres.size()[1]
        return ms, concentric_spheres, weights

    def becke_tesellation(
        self, x: torch.Tensor, r: torch.Tensor, dm: torch.Tensor,
        i: torch.Tensor
    ) -> torch.Tensor:
        """

        Generates a weight scheme to avoid overlap

        :param x: sample
        :param r: centers
        :param dm: centers distance matrix
        :param i: map sample -> center
        :return:
        """

        nx = x.size()[0]  # number of sample matches. It should match nx
        mx = x.size()[1]  # number of samples per molecule
        nr = r.size()[0]  # number of molecules
        mr = r.size()[1]  # number of atoms per molecule
        if nx != nr:
            raise IOError("batch dims don't not match: {:d} {:d}".format(nx, nr))

        # calculate r-r distance matrix

        # calculate x-r distance matrix
        # dim 1 will be r and dim 2 will be x
        # final dims should be [nr, mx, mr]

        x_exp = x.unsqueeze(2).expand(nr, mx, mr, 3)
        r_exp = r.unsqueeze(1).expand(nr, mx, mr, 3)
        rx_dm = torch.norm(x_exp - r_exp, dim=3)
        del x_exp, r_exp

        # calculate (x-r)-(x-r)T distance vectors

        rx_dm_expanded1 = rx_dm.unsqueeze(3).expand(nr, mx, mr, mr)
        rx_dm_expanded2 = rx_dm.unsqueeze(2).expand(nr, mx, mr, mr)
        xx_dm_expanded = dm.unsqueeze(0)
        mu = (rx_dm_expanded1 - rx_dm_expanded2) / xx_dm_expanded.clamp(min=1e-12)

        del rx_dm_expanded1, rx_dm_expanded2, xx_dm_expanded
        del dm, rx_dm
        # soft-cutoff function

        for _ in range(1, self.softening):
            mu = 0.5 * mu * (3.0 - (mu.pow(2.0)))

        mu = 0.5 * (1 - mu)

        # diagonal fill

        diag_fill = torch.eye(
            mr, mr, device=self.device, dtype=self.dtype
        ).unsqueeze(0).unsqueeze(1).expand(nr, mx, mr, mr)

        mu += (diag_fill * 0.5)

        # productory

        mu = mu.prod(dim=3)

        # mu should be now [nr, mx, mr] containing
        # the weights of each nuclei for each sampling point

        i = i.unsqueeze(2)

        v_i = torch.gather(mu, 2, i).squeeze(2)
        v_all = mu.sum(2)
        w = v_i / v_all

        # w should be now [nr, mx] containing the weight of each sample

        return w


class DxGrid:

    def __init__(self, device: torch.device, dtype: torch.dtype, resolution: float, spacing: float, units: Unit = bohr):
        """
        Dx grid
        Eases representation of electron density in a volumetric format that can be read by popular chemistries
        software as Chimera


        :param device:
        :param dtype:
        :param resolution: it shares units with the coordinates
        :param spacing: it shares units with the coordinates
        :param units: units of the coordinates
        """

        self.device = device
        self.dtype = dtype
        self.resolution = resolution
        self.spacing = spacing
        self.units = units
        from a2mdio.units import angstrom
        from a2mdio.volumes import Volume
        self.angstrom = angstrom
        self.volume = Volume

    def generate_grid(self, coords: Coordinates):
        """
        Generates a box that wraps coordinates using a voxel of size (resolution^3)
        :param coords:
        :return:
        """
        coords = coords.get_cartessian(unit=self.units).squeeze(0)

        # rotcoords, basis, mean = CoordinatesSampler.principal_components(coords)
        basis = torch.eye(3, 3, dtype=self.dtype, device=self.device) * (self.resolution * (self.units / self.angstrom))
        box_min = coords.min(0, keepdim=False)[0] - self.spacing
        box_max = coords.max(0, keepdim=False)[0] + self.spacing
        diff = (box_max + (box_min * -1))
        dims = ((diff / self.resolution).floor() + 1).to(torch.long)

        xx = torch.arange(dims[0], dtype=self.dtype, device=self.device) * self.resolution
        yy = torch.arange(dims[1], dtype=self.dtype, device=self.device) * self.resolution
        zz = torch.arange(dims[2], dtype=self.dtype, device=self.device) * self.resolution

        xg, yg, zg = torch.meshgrid(xx, yy, zz)
        xg = xg.flatten()
        yg = yg.flatten()
        zg = zg.flatten()

        rgrid = torch.stack([xg, yg, zg], dim=1)

        rgrid += box_min
        rgrid = rgrid.reshape(-1, 3).unsqueeze(0)

        box_min *= (self.units / self.angstrom)
        rgrid = Coordinates(rgrid, units=self.units)
        return rgrid, box_min.tolist(), basis.data.cpu().numpy(), dims.tolist()

    def dx(self, grid: torch.Tensor, r0: List[float], basis: np.ndarray, dims: List[int]):
        """
        converts a tensor into a grid that can be stored in dx format
        :param grid:
        :param r0:
        :param basis:
        :param dims:
        :return:
        """
        grid = grid.reshape(dims[0], dims[1], dims[2]).data.cpu().numpy()
        return self.volume(
            filename=None, dxvalues=grid, r0=r0, basis=basis
        )
