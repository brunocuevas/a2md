import numpy as np
from a2mdnet import MAX_NUMBER_EVAL_POINTS, ELEMENT2NN
from a2mdio.molecules import Mol2
from torch.utils import data
import torch
import json
import warnings

# import rdkit.Chem.AllChem as Chem


def convert_label2tensor(label, device=torch.device('cpu')):
    u = [ELEMENT2NN[i] for i in label]
    return torch.tensor(u, device=device, dtype=torch.long)


def convert_targets2tensor(
        targets, device=torch.device('cpu'), dtype=torch.float
):
    return torch.tensor(targets, device=device, dtype=dtype).reshape(-1, 1)


class MolecularDataset(data.Dataset):
    def __init__(
            self, device, dtype, ids=None,
            model_parameters_path=None,
            molecular_data_path=None,
            prefetch=False
    ):
        """

        :param device:
        :param dtype:
        :param ids:
        :param path2targets:
        :param path2data:
        :param normalize_sample:
        :param charges
        :param integrals
        """
        import math
        from a2md.mathfunctions import expfun_integral, angular_gaussian_integral
        self.radial_integral = expfun_integral
        self.angular_integral = angular_gaussian_integral
        self.ids = ids

        max_atoms = 0
        max_bonds = 0

        for i in ids:

            if molecular_data_path is not None:
                mol2_filename = molecular_data_path + i + '.mol2'

            else:
                mol2_filename = i + '.mol2'

            mm = Mol2(file=mol2_filename, verbose=False)

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
        if fun['support_type'] in ['ORCV', 'OR', 'ORV', 'ORC']:
            a = fun['params']['A3']
            b = fun['params']['B3']
            return self.radial_integral(a, b) * 2
        elif fun['support_type'] == 'AG':
            g = fun['params']['G']
            u = fun['params']['U']
            alpha = fun['params']['Alpha']
            return self.angular_integral(g, u, alpha)

    @staticmethod
    def match_bond(anisotropic_params, connectivity_tensor, fun, nbonds):
        """

        :param anisotropic_params:
        :param connectivity_tensor:
        :param fun:
        :param nbonds:
        :return:
        """
        edge = int(fun['params']['Psi'] == 1)  # 0 if Psi = 0, 1 if Psi = 1
        flag = False
        j = 0
        for j in range(nbonds):
            c1 = (connectivity_tensor[j, 0] == fun['center']) and (
                    connectivity_tensor[j, 1] == fun['bond']
            )
            c2 = (connectivity_tensor[j, 1] == fun['center']) and (
                    connectivity_tensor[j, 0] == fun['bond']
            )
            if c1:
                flag = True
                anisotropic_params[j, edge] = fun['coefficient']
                break
            elif c2:
                flag = True
                edge = edge + 2
                anisotropic_params[j, edge] = fun['coefficient']
                break
        if not flag:
            warnings.warn("bond was not matched")
            return None
        else:
            return j, edge

    def __len__(self):
        return len(self.ids)

    def fetch(self, item):
        """

        :param item:
        :return:
        """
        import math

        identifier = self.ids[item]
        mm = Mol2(self.mol_path + identifier + '.mol2')
        with open(self.params_path + identifier + '.ppp') as f:
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
        charges_tensor[:na] = torch.tensor(mm.get_charge(kind='total'), device=self.device, dtype=self.dtype)

        for fun in params:

            if math.isnan(fun['coefficient']):
                raise IOError("issue with instance {:s}".format(identifier))

            if fun['support_type'] == 'ORC':
                charges_tensor[fun['center']] -= self.calculate_integral(fun)

            if fun['support_type'] == 'ORCV' or fun['support_type'] == 'OR':
                isotropic_params[fun['center'], 0] = fun['coefficient']
                isotropic_integrals[fun['center'], 0] = self.calculate_integral(fun)

            elif fun['support_type'] == 'ORV':
                isotropic_params[fun['center'], 1] = fun['coefficient']
                isotropic_integrals[fun['center'], 1] = self.calculate_integral(fun)

            elif fun['support_type'] == 'AG':
                current_center, current_type = self.match_bond(anisotropic_params, connectivity_tensor, fun, nb)
                anisotropic_integrals[current_center, current_type] = self.calculate_integral(fun)

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


class CompleteSetDensityParams(MolecularDataset):
    def __init__(
            self, device, dtype, kind='training', number=None,
            charges=False, integrals=False
    ):
        from a2mdnet import COMPLETE_SET_DENSITY_DATASET
        if kind not in ['curated_training', 'curated_validation', 'training', 'validation', 'all', 'test']:
            raise IOError("the requested complete set kind does not exist")

        idx_file = COMPLETE_SET_DENSITY_DATASET[kind + '_idx']

        with open(idx_file) as f:
            idx = json.load(f)

        if number is not None:
            idx = idx[:number]

        MolecularDataset.__init__(
            self,
            device=device,
            dtype=dtype,
            ids=idx,
            path2targets=COMPLETE_SET_DENSITY_DATASET['density_path'],
            path2data=COMPLETE_SET_DENSITY_DATASET['geometry_path'],
            charges=charges, integrals=integrals
        )


class MolecularElectronDensityDataset(MolecularDataset):
    DENSITY_POINTS=1000
    def __init__(
        self, device, dtype, ids=None, model_parameters_path=None,
        molecular_data_path=None, density_data_path=None,
        mol_prop=None, prefetch=True
    ):
        MolecularDataset.__init__(
            self, device, dtype, ids, model_parameters_path,
            molecular_data_path, prefetch
        )

        self.density_data_path = density_data_path
        self.density_buffer = []

        if self.density_data_path is not None and prefetch:
            for item in self.ids:
                self.density_buffer.append(
                    np.load(self.density_data_path + item + '.npy')
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

                    density = np.load(self.density_data_path + self.ids[item] + '.npy')
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

class Aniset(MolecularElectronDensityDataset):
    def __init__(
            self, device, dtype, kind='training', number=None,
            energy=False
    ):
        from a2mdnet import ANISET_DATASET
        if kind not in ['training', 'validation', 'all', 'test']:
            raise IOError("the requested complete set kind does not exist")

        idx_file = ANISET_DATASET[kind + '_idx']

        with open(idx_file) as f:
            idx = json.load(f)

        if number is not None:
            idx = idx[:number]

        mol_prop = None
        if energy:
            mol_prop = ANISET_DATASET['energy']


        MolecularElectronDensityDataset.__init__(
            self,
            device=device,
            dtype=dtype,
            ids=idx,
            model_parameters_path=ANISET_DATASET['density_path'],
            molecular_data_path=ANISET_DATASET['geometry_path'],
            density_data_path=ANISET_DATASET['density_values_path'],
            load_charges=load_charges, load_integrals=load_integrals, mol_prop=mol_prop
        )

class Fdaset(MolecularElectronDensityDataset):
    def __init__(
            self, device, dtype, kind='training', number=None,
            energy=False
    ):
        from a2mdnet import FDASET
        if kind not in ['all', 'test']:
            raise IOError("the requested complete set kind does not exist")

        idx_file = FDASET[kind + '_idx']

        with open(idx_file) as f:
            idx = json.load(f)

        if number is not None:
            idx = idx[:number]

        mol_prop = None
        if energy:
            mol_prop = FDASET['energy']

        MolecularElectronDensityDataset.__init__(
            self,
            device=device,
            dtype=dtype,
            ids=idx,
            model_parameters_path=FDASET['density_path'],
            molecular_data_path=FDASET['geometry_path'],
            density_data_path=FDASET['density_values_path'],
            load_charges=load_charges, load_integrals=load_integrals, mol_prop=mol_prop
        )
