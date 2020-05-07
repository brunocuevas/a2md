import numpy as np
from a2mdnet import ELEMENT2NN
from a2mdio.molecules import Mol2
from a2md.utils import integrate_from_dict
from a2mdnet import FUNCTION_NAMES2POSITION
from torch.utils import data
import torch
import json
import warnings

def convert_label2tensor(label, device=torch.device('cpu')):
    u = [ELEMENT2NN[i] for i in label]
    return torch.tensor(u, device=device, dtype=torch.long)

def convert_targets2tensor(
        targets, device=torch.device('cpu'), dtype=torch.float
):
    return torch.tensor(targets, device=device, dtype=dtype).reshape(-1, 1)

def match_old_fun_names(fun):
    if fun['bond'] is not None and fun['support_type'] == "AG":
        if fun['params']['Psi'] == 0:
            return "aniso", 0
        elif fun['params']['Psi'] == 1:
            return "aniso", 1
        else:
            raise IOError("old parametrization does not allow Psi value {:d}".format(fun['parameters']['Psi']))
    elif fun['bond'] is None:
        if fun['support_type'] == 'ORC':
            return "core", None
        elif fun['support_type'] == 'ORCV':
            return "iso", 0
        elif fun['support_type'] == 'ORV':
            return "iso", 1
        elif fun['support_type'] == 'OR':
            return "iso", 0
        else:
            raise IOError("unknown support type {:s}".format(fun['support_type']))
    else:
        raise IOError("can not understand function {:s}".format(json.dumps(fun)))

def match_fun_names(fun):
    fn = fun['support_type']
    return FUNCTION_NAMES2POSITION[fn]


class MolecularDataset(data.Dataset):
    def __init__(
            self, device, dtype, ids=None,
            model_parameters_path=None,
            molecular_data_path=None,
            prefetch=False, integration_method=integrate_from_dict,
            match_method=match_fun_names
    ):
        """
        MolecularDataset
        ---
        Allows to load batches of coordinates + atom_types + topology + charge + function integrals +
        coefficients targets. It is used to train deep learning models.


        :param device: either cuda or cpu
        :param dtype: float32 or float64 to specify either single or double precission
        :param ids: names of the molecules that will be read by this data loader
        :param model_parameters_path: path where parameters files of specified molecules should be located
        :param molecular_data_path: path where Mol2 files of specified molecules should be located
        :param prefetch: wether to read all the inputs at once or not. It saves time in long term, but it takes time
        """

        self.ids = ids

        max_atoms = 0
        max_bonds = 0

        for i in ids:

            mol2_filename = molecular_data_path / '{:s}.mol2'.format(i)
            mm = Mol2(file=mol2_filename)

            natoms = mm.get_number_atoms()
            nbonds = mm.get_number_bonds()

            if natoms > max_atoms:
                max_atoms = natoms
            if nbonds > max_bonds:
                max_bonds = nbonds

        self.mol_path = molecular_data_path
        self.params_path = model_parameters_path
        self.device = device
        self.dtype = dtype
        self.prefetch = prefetch
        self.max_atoms = max_atoms
        self.max_bonds = max_bonds

        self.integration_method = integration_method
        self.match_method = match_method

        if prefetch:
            self.labels = []
            self.connectivity = []
            self.coordinates = []
            self.charge = []
            self.integrals_iso = []
            self.integrals_aniso = []
            self.iso = []
            self.aniso = []

            for idx in range(len(self.ids)):
                l, t, x, q, i_iso, i_aniso, iso, aniso = self.fetch(idx)
                self.labels.append(l)
                self.connectivity.append(t)
                self.coordinates.append(x)
                self.charge.append(q)
                self.integrals_iso.append(i_iso)
                self.integrals_aniso.append(i_aniso)
                self.iso.append(iso)
                self.aniso.append(aniso)
                if idx % 100 == 0:
                    if idx % 1000 != 0:
                        print("-", end="")
                    elif idx % 1000 == 0 and idx != 0:
                        print("-", end=" {:8d}\n".format(idx))


    def calculate_integral(self, fun):
        """

        :param fun:
        :return:
        """
        return self.integration_method(fun)

    @staticmethod
    def match_bond(connectivity_tensor, fun, nbonds):
        """

        :param connectivity_tensor:
        :param fun:
        :param nbonds:
        :return:
        """
        for j in range(nbonds):
            c1 = (connectivity_tensor[j, 0] == fun['center']) and (
                    connectivity_tensor[j, 1] == fun['bond']
            )
            c2 = (connectivity_tensor[j, 1] == fun['center']) and (
                    connectivity_tensor[j, 0] == fun['bond']
            )
            if c1:
                return j, 0
            elif c2:
                return j, 2

        warnings.warn("bond was not matched")
        return None, None

    def __len__(self):
        return len(self.ids)

    def fetch(self, item):
        """

        :param item:
        :return:
        """
        import math

        identifier = self.ids[item]
        mm = Mol2(self.mol_path / '{:s}.mol2'.format(identifier))
        with open(self.params_path / '{:s}.ppp'.format(identifier)) as f:
            params = json.load(f)

        na = mm.get_number_atoms()
        nb = mm.get_number_bonds()

        coords_tensor = torch.zeros(self.max_atoms, 3, dtype=self.dtype, device=torch.device('cpu'))
        labels_tensor = torch.ones(self.max_atoms, dtype=torch.long, device=torch.device('cpu')) * -1
        connectivity_tensor = torch.ones(self.max_bonds, 2, device=torch.device('cpu'), dtype=torch.long) * -1
        charges_tensor = torch.zeros(self.max_atoms, dtype=self.dtype, device=torch.device('cpu'))
        isotropic_params = torch.zeros(self.max_atoms, 2, dtype=self.dtype, device=torch.device('cpu'))
        anisotropic_params= torch.zeros(self.max_bonds, 4, dtype=self.dtype, device=torch.device('cpu'))
        isotropic_integrals= torch.zeros(self.max_atoms, 2, dtype=self.dtype, device=torch.device('cpu'))
        anisotropic_integrals = torch.zeros(self.max_bonds, 4, dtype=self.dtype, device=torch.device('cpu'))

        coords_tensor[:na, :] = torch.tensor(mm.get_coordinates(), device=torch.device('cpu'), dtype=self.dtype)
        labels_tensor[:na] = convert_label2tensor(mm.get_atomic_numbers(), device=self.device)
        connectivity_tensor[:nb, :] = torch.tensor(mm.get_bonds() - 1, device=torch.device('cpu'), dtype=torch.uint8)
        charges_tensor[:na] = torch.tensor(mm.get_absolute_charges(), device=self.device, dtype=self.dtype)

        for fun in params:

            if math.isnan(fun['coefficient']):
                raise IOError("issue with instance {:s}".format(identifier))

            funtype, pos = self.match_method(fun)
            if funtype == 'core':
                charges_tensor[fun['center']] -= self.calculate_integral(fun)
            elif funtype == 'iso':
                isotropic_params[fun['center'], pos] = fun['coefficient']
                isotropic_integrals[fun['center'], pos] = self.calculate_integral(fun)
            elif funtype == 'aniso':
                bond_idx, col = self.match_bond(connectivity_tensor, fun, nb)
                anisotropic_params[bond_idx, col + pos] = fun['coefficient']
                anisotropic_integrals[bond_idx, col + pos] = self.calculate_integral(fun)

        return [
            i.to(self.device) for i in [labels_tensor, connectivity_tensor, coords_tensor, charges_tensor,
                                        isotropic_integrals, anisotropic_integrals, isotropic_params,
                                        anisotropic_params]
        ]

    def __getitem__(self, item):

        if self.prefetch:

            return [
                self.labels[item], self.connectivity[item], self.coordinates[item],
                self.charge[item], self.integrals_iso[item], self.integrals_aniso[item],
                self.iso[item], self.aniso[item]
            ]

        else:

            return self.fetch(item)


class MolecularElectronDensityDataset(MolecularDataset):
    DENSITY_POINTS=1000
    def __init__(
        self, device, dtype, ids=None, model_parameters_path=None,
        molecular_data_path=None, density_data_path=None,
        mol_prop=None, prefetch=True, integration_method=integrate_from_dict,
        match_method=match_fun_names
    ):
        MolecularDataset.__init__(
            self, device, dtype, ids, model_parameters_path,
            molecular_data_path, prefetch, integration_method=integration_method,
            match_method=match_method
        )

        self.density_data_path = density_data_path
        self.density_buffer = []

        if self.density_data_path is not None and prefetch:
            for item in self.ids:
                self.density_buffer.append(
                    np.load(self.density_data_path / '{:s}.npy'.format(item))
                )

        if mol_prop is None:
            self.mol_prop = None
        else:
            with open(mol_prop) as f:
                self.mol_prop = json.load(f)

    def __len__(self):
        return len(self.ids)

    def __getitem__(self, item):
        if self.density_data_path is None and self.mol_prop is None:
            return super(MolecularElectronDensityDataset, self).__getitem__(item=item)
        else:
            output = super(MolecularElectronDensityDataset, self).__getitem__(item=item)
            if self.density_data_path is not None:

                if not self.prefetch:

                    density = np.load(self.density_data_path / '{:s}.npy'.format(self.ids[item]))
                    density = density[np.random.randint(0, density.shape[0], self.DENSITY_POINTS), :]
                    density = torch.tensor(density, dtype=torch.float, device=self.device)
                    output.append(density[:, :3])
                    output.append(density[:, 3])

                else:

                    density = self.density_buffer[item]
                    density = density[np.random.randint(0, density.shape[0], self.DENSITY_POINTS), :]
                    density = torch.tensor(density, dtype=torch.float, device=self.device)
                    output.append(density[:, :3])
                    output.append(density[:, 3])

            if self.mol_prop is not None:
                output.append(torch.tensor(self.mol_prop[self.ids[item]], dtype=self.dtype, device=self.device))

            return output
