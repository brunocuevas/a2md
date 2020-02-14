from a2md.baseclass import A2MD_basis
import copy
import numpy as np
import json
from importlib import import_module
import os
import sys

ANG2AU = 1.8897

ELEMENT2AN = dict(
    H=1,
    C=6,
    N=7,
    O=8
)

class a2md_general_parameters(A2MD_basis):
    def __init__(self, version, verbose = False, filename = None):
        A2MD_basis.__init__(self, verbose=verbose, name='a2md parameters store')
        self.__filename = filename
        self.__version = version
        self.__params = None

    def __check_fundamental_params(self, tmp_params):
        self.log("checking parameters")
        pars = tmp_params['_ATOMS']
        for k, atom in pars.items():
            for fun in  atom:
                name_check = '_NAME' in fun
                path_check = '_PATH' in fun
                kws_check = '_KWS' in fun
                specific_pars_check = '_PARS' in fun
                bonding_check = '_BOND' in fun
                global_check = name_check & path_check & kws_check & specific_pars_check & bonding_check
                if not global_check:
                    raise IOError("some fundamental keywords were missing. check documentation")
        self.log("no errors fundamental kws error detected")
        return True

    def get_atom(self, atom):
        all_atoms = self.__params['_ATOMS']
        return all_atoms[atom]

    def get_atoms(self):
        return list(self.__params['_ATOMS'].keys())

    def get_filename(self):
        return self.__filename

    def get_version(self):
        return self.__version

    def get_params(self):
        return self.__params

    def set_filename(self):
        return self.__filename

    def set_params_from_dict(self, pars):
        if pars['_VERSION'] == self.__version :
            self.log("version was correct")
        else:
            self.log("WARNING. Version is not correct. Errors may arise")
        if self.__check_fundamental_params(pars):
            self.__params = pars
        else:
            raise IOError

    def read_parameters(self):
        with open(self.__filename) as f:
            tmp_params = json.load(f)
        if tmp_params['_VERSION'] == self.__version :
            self.log("version was correct")
        else:
            self.log("WARNING. Version is not correct. Errors may arise")
        if self.__check_fundamental_params(tmp_params):
            self.__params = tmp_params
        else:
            raise IOError


    def set_support_library(self):
        from a2md.amd import SUPPORT_TYPE
        for k, atom in self.__params['_ATOMS'].items():
            for fun in atom:
                name = fun['_NAME']
                if fun['_NAME'] not in SUPPORT_TYPE.keys():
                    try :
                        module_name = os.path.basename(fun['_PATH']).replace('.py','')
                        current_module = import_module(module_name, fun['_PATH'])
                        SUPPORT_TYPE[name] = getattr(current_module, fun['_NAME'])
                    except AttributeError:
                        raise IOError("the provided function was not located in the provided file")
                    except ModuleNotFoundError:
                        sys.path.append(fun['_PATH'].replace(os.path.basename(fun['_PATH']), ''))
                        module_name = os.path.basename(fun['_PATH']).replace('.py', '')
                        current_module = import_module(module_name, fun['_PATH'])
                        SUPPORT_TYPE[name] = getattr(current_module, fun['_NAME'])


    def write_parameters(self, filename):
        if self.__params is None :
            raise RuntimeError("there were no parameters")
        with open(filename, 'w') as f:
            json.dump(self.__params, f, indent=4, sort_keys=True)


class preprocessor(A2MD_basis):
    def __init__(self, file=None, input_format=None, verbose=True):
        """
        preprocessor eases the parametrization of aAMD.
        - aAMD uses a list of individual atom features
        - parameters are specified as a dictionary in which entries appear as functions of element name
        - so aAMD makes a new list by reading the molecular coordinates and then it uses them to feed aAMD
        :param verbose:
        """
        from a2md import AMD_PARAMS_19
        A2MD_basis.__init__(self, name='aAMD preprocessor', verbose=verbose)
        self.__params = copy.copy(AMD_PARAMS_19)
        self.__molecule_coordinates = None
        self.__molecule_elements = None
        self.__molecule_topo = None
        self.__molecule_charge = None

        if file is not None:
            self.read_molecule(
                file=file,
                input_format=input_format
            )

    def __parametrize_frozen_core_gaussian_model(self):
        number_atoms = len(self.__molecule_elements)
        A = []
        Xi = []
        C = []
        PN = []
        F = []

        for atom_idx in range(number_atoms):

            elem = self.__molecule_elements[atom_idx]
            if elem == 'H':
                A.append(atom_idx)
                Xi.append(
                    dict(
                        A3=self.__params[elem]['AA'],
                        B3=self.__params[elem]['BA']
                    )
                )
                PN.append((atom_idx, 'OR'))
                F.append(False)
                C.append(1.0)
            else:
                # Adding the core function
                A.append(atom_idx)
                Xi.append(
                    dict(
                        A3=self.__params[elem]['A1'],
                        B3=self.__params[elem]['B1'],
                    )
                )
                PN.append((atom_idx, 'ORC'))
                F.append(True)
                C.append(1.0)
                # Adding the core-valence function
                A.append(atom_idx)
                Xi.append(
                    dict(
                        A3=self.__params[elem]['A2'],
                        B3=self.__params[elem]['B2'],
                    )
                )
                PN.append((atom_idx, 'ORCV'))
                F.append(False)
                C.append(1.0)
                # Adding the valence function
                A.append(atom_idx)
                Xi.append(
                    dict(
                        A3=self.__params[elem]['A3'],
                        B3=self.__params[elem]['B3'],
                    )
                )
                PN.append((atom_idx, 'ORV'))
                F.append(False)
                C.append(1.0)
                # Angular functions
            for bond in self.__molecule_topo[atom_idx]:
                A.append(atom_idx)
                Xi.append(
                    dict(
                        U=self.__params[elem]['U'],
                        G=self.__params[elem]['G'],
                        Alpha=self.__params[elem]['Alpha1'],
                        Psi=0
                    )
                )

                PN.append((atom_idx, 'AG', bond))
                F.append(False)
                C.append(1.0)

                A.append(atom_idx)
                Xi.append(
                    dict(
                        U=self.__params[elem]['U'],
                        G=self.__params[elem]['G'],
                        Alpha=self.__params[elem]['Alpha2'],
                        Psi=1
                    )
                )

                PN.append((atom_idx, 'AG', bond))
                F.append(False)
                C.append(1.0)

        number_functions = len(A)
        atom_parametrization_list = [dict()] * number_functions

        for xi_iter in range(number_functions):
            tmp_input = dict(
                center=A[xi_iter],
                support_type=PN[xi_iter][1],
                coefficient=C[xi_iter],
                frozen=F[xi_iter]
            )
            tmp_params = dict()
            for key, item in Xi[xi_iter].items():
                tmp_params[key] = item
            tmp_input['params'] = tmp_params
            if PN[xi_iter][1][0] == 'A':
                try:
                    tmp_input['bond'] = PN[xi_iter][2]
                except IndexError:
                    raise IOError("angular support functions should have a bonded atom")
            else:
                tmp_input['bond'] = None

            atom_parametrization_list[xi_iter] = tmp_input

        return atom_parametrization_list

    def __parametrize_frozen_core_model(self):
        number_atoms = len(self.__molecule_elements)
        A = []
        Xi = []
        C = []
        PN = []
        F = []

        for atom_idx in range(number_atoms):

            elem = self.__molecule_elements[atom_idx]
            if elem == 'H':
                A.append(atom_idx)
                Xi.append(
                    dict(
                        A3=self.__params[elem]['AA'],
                        B3=self.__params[elem]['BA']
                    )
                )
                PN.append((atom_idx, 'OR'))
                F.append(False)
                C.append(1.0)
            else:
                # Adding the core function
                A.append(atom_idx)
                Xi.append(
                    dict(
                        A3=self.__params[elem]['A1'],
                        B3=self.__params[elem]['B1'],
                    )
                )
                PN.append((atom_idx, 'ORC'))
                F.append(True)
                C.append(1.0)
                # Adding the core-valence function
                A.append(atom_idx)
                Xi.append(
                    dict(
                        A3=self.__params[elem]['A2'],
                        B3=self.__params[elem]['B2'],
                    )
                )
                PN.append((atom_idx, 'ORCV'))
                F.append(False)
                C.append(1.0)
                # Adding the valence function
                A.append(atom_idx)
                Xi.append(
                    dict(
                        A3=self.__params[elem]['A3'],
                        B3=self.__params[elem]['B3'],
                    )
                )
                PN.append((atom_idx, 'ORV'))
                F.append(False)
                C.append(1.0)
                # Angular functions
            for bond in self.__molecule_topo[atom_idx]:
                for psi in range(self.__params[elem]['Psi']):
                    A.append(atom_idx)
                    Xi.append(
                        dict(
                            U=self.__params[elem]['U'],
                            G=self.__params[elem]['G'],
                            Psi=psi
                        )
                    )
                    
                    PN.append((atom_idx, 'AC', bond))
                    F.append(False)
                    C.append(1.0)

        number_functions = len(A)
        atom_parametrization_list = [dict()] * number_functions

        for xi_iter in range(number_functions):
            tmp_input = dict(
                center=A[xi_iter],
                support_type=PN[xi_iter][1],
                coefficient=C[xi_iter],
                frozen=F[xi_iter]
            )
            tmp_params = dict()
            for key, item in Xi[xi_iter].items():
                tmp_params[key] = item
            tmp_input['params'] = tmp_params
            if PN[xi_iter][1][0] == 'A':
                try:
                    tmp_input['bond'] = PN[xi_iter][2]
                except IndexError:
                    raise IOError("angular support functions should have a bonded atom")
            else:
                tmp_input['bond'] = None

            atom_parametrization_list[xi_iter] = tmp_input

        return atom_parametrization_list

    @staticmethod
    def __read_mol2(file, xyz_array, elem_array, topo_array, charge_array):
        from a2mdlib.molecules import Mol2
        mol2_instance = Mol2(file, verbose=False)

        labels = mol2_instance.get_labels()
        charges = mol2_instance.get_charge(kind='total')
        coords = mol2_instance.get_coordinates() * ANG2AU
        topology_matrix = mol2_instance.get_bonds()

        for i in range(mol2_instance.get_number_atoms()):
            charge_array.append(float(charges[i]))
            elem_array.append(labels[i])
            xyz_array.append([float(coords[i, 0]), float(coords[i, 1]), float(coords[i, 2])])

        for i in range(mol2_instance.get_number_atoms()): topo_array.append([])
        for i in range(mol2_instance.get_number_bonds()):
            begin = int(topology_matrix[i, 0])
            end = int(topology_matrix[i, 1])
            topo_array[begin - 1].append(end - 1)
            topo_array[end - 1].append(begin - 1)

    def get_atomic_numbers(self):
        return [ELEMENT2AN[i] for i in self.__molecule_elements]

    def get_coordinates(self):
        return copy.copy(self.__molecule_coordinates)

    def get_labels(self):
        return copy.copy(self.__molecule_elements)

    def get_topology(self):
        return copy.copy(self.__molecule_topo)

    def get_charge(self):
        return copy.copy(self.__molecule_charge)

    def get_bond_median_points(self):
        mol_coords=self.__molecule_coordinates
        mol_topo=self.__molecule_topo
        topo_sorted_dict = []
        median_points_coords = []
        for a1, bonded in enumerate(mol_topo):
            for a2 in bonded:
                label = ''.join(sorted([str(a1), str(a2)]))
                if not label in topo_sorted_dict :
                    topo_sorted_dict.append(label)
                    median = (mol_coords[a1,:] + mol_coords[a2,:])/2.0
                    median_points_coords.append(median)
        return np.array(median_points_coords)

    def read_parameters(self, filename=None, dictionary=None):
        """
        Reads a collection of parameters
        :param filename: name of the parameters file. It will be read as .json
        :param dictionary:
        :type filename: str
        :type dictionary: dict
        :return:
        """
        if (filename is None) and (dictionary is None):
            raise IOError("nothing to read")
        elif dictionary is None:
            with open(filename) as f :
                self.log("file %s was opened" % filename)
                params = json.load(f)
                self.__params = params
        elif filename is None:
            self.__params = dictionary
        else:
            self.log("reading from default parameters")
            from . import AMD_PARAMS_18
            self.__params = copy.copy(AMD_PARAMS_18)
        return True

    def read_molecule(self, file, input_format ='mol2'):
        """
        Reads a molecule stored in a table
        :param file: name of the coordinates file
        :param input_format: type of file
        :type file: str
        :type input_format: str
        :return:
        """
        xyz_array = []
        elem_array = []
        topo_array = []
        charge_array = []


        if input_format == 'mol2':
            self.__read_mol2(file, xyz_array, elem_array, topo_array, charge_array)
        else:
            raise NotImplementedError("only mol2 files are implemented by now")

        self.__molecule_coordinates = np.array(xyz_array, dtype='float64')
        self.__molecule_elements = elem_array
        self.__molecule_topo = topo_array
        self.__molecule_charge = np.array(charge_array)

    def parametrize(self, kind = 'proto'):
        if not kind in [
            'proto', 'mol', 'bond', 'trigonometric', 'none', 'frozen_core', 'frozen_core_gaussian'
        ]:
            raise IOError("the specified kind of setup is not implemented")

        number_atoms = len(self.__molecule_elements)
        A  = []
        Xi = []
        C  = []
        PN = []

        if kind == 'none':
            return self.__molecule_coordinates, self.__molecule_charge, self.__molecule_topo, None
        elif kind == 'frozen_core':
            atom_parametrization_list = self.__parametrize_frozen_core_model()
            return self.__molecule_coordinates, self.__molecule_charge, self.__molecule_topo, atom_parametrization_list
        elif kind == 'frozen_core_gaussian':
            atom_parametrization_list = self.__parametrize_frozen_core_gaussian_model()
            return self.__molecule_coordinates, self.__molecule_charge, self.__molecule_topo, atom_parametrization_list
        else :
            for aim in range(number_atoms):
                elem = self.__molecule_elements[aim]
                if kind in ['proto', 'mol', 'trigonometric']:
                    A.append(aim)
                    Xi.append(
                        dict(
                            A1=self.__params[elem]['A1'],
                            A2=self.__params[elem]['A2'],
                            B1=self.__params[elem]['B1'],
                            B2=self.__params[elem]['B2'],
                        )
                    )
                    PN.append((aim, 'IR'))
                    C.append(1.0)

                    A.append(aim)
                    Xi.append(
                        dict(
                            A3=self.__params[elem]['A3'],
                            B3=self.__params[elem]['B3'],
                        )
                    )
                    PN.append((aim, 'OR'))
                    C.append(1.0)

                if kind == 'mol' or kind == 'bond':
                    for bond in self.__molecule_topo:

                        A.append(aim)
                        Xi.append(
                            dict(
                                F=self.__params[elem]['F'],
                                G=self.__params[elem]['G'],
                                U=self.__params[elem]['U'],
                                V=self.__params[elem]['V'],
                                Psi=self.__params[elem]['Psi'],
                            )
                        )
                        PN.append((aim, 'ACS', bond))
                        C.append(1.0)

                if kind == 'trigonometric':
                    for bond in self.__molecule_topo[aim]:
                        for psi in range(self.__params[elem]['Psi']):

                            A.append(aim)
                            Xi.append(
                                dict(
                                    U=self.__params[elem]['U'],
                                    G=self.__params[elem]['G'],
                                    Psi=psi
                                )
                            )
                            PN.append((aim, 'AC', bond))
                            C.append(1.0)

                            A.append(aim)
                            Xi.append(
                                dict(
                                    U=self.__params[elem]['V'],
                                    G=self.__params[elem]['G'],
                                    Psi=psi
                                )
                            )
                            PN.append((aim, 'AS', bond))
                            C.append(1.0)

            number_functions = len(A)

            atom_parametrization_list = [dict()] * number_functions

            for xi_iter in range(number_functions):
                tmp_input = dict(center=A[xi_iter], support_type=PN[xi_iter][1], coefficient=C[xi_iter])
                tmp_params = dict()
                for key, item in Xi[xi_iter].items():
                    tmp_params[key] = item
                tmp_input['params'] = tmp_params
                if PN[xi_iter][1][0] == 'A':
                    try:
                        tmp_input['bond'] = PN[xi_iter][2]
                    except IndexError:
                        raise IOError("angular support functions should have a bonded atom")
                else:
                    tmp_input['bond'] = None

                atom_parametrization_list[xi_iter] = tmp_input

            return self.__molecule_coordinates, self.__molecule_charge, self.__molecule_topo, atom_parametrization_list

    def parametrize_from_modelfile(self, modelfile, input_format='json'):
        """

        :param modelfile:
        :param input_format
        :return:
        """
        mp = model_parser(m=modelfile, input_format=input_format)
        A, Xi, PN, F, C = mp.build_model(
            coordinates=self.__molecule_coordinates,
            labels=self.__molecule_elements,
            topology=self.__molecule_topo,
            charges=self.__molecule_charge
        )

        number_functions = len(A)
        atom_parametrization_list = [dict()] * number_functions

        for xi_iter in range(number_functions):
            tmp_input = dict(
                center=A[xi_iter],
                support_type=PN[xi_iter][1],
                coefficient=C[xi_iter],
                frozen=F[xi_iter]
            )
            tmp_params = dict()
            for key, item in Xi[xi_iter].items():
                tmp_params[key] = item
            tmp_input['params'] = tmp_params
            if PN[xi_iter][1][0] == 'A':
                try:
                    tmp_input['bond'] = PN[xi_iter][2]
                except IndexError:
                    raise IOError("angular support functions should have a bonded atom")
            else:
                tmp_input['bond'] = None

            atom_parametrization_list[xi_iter] = tmp_input

        return self.__molecule_coordinates, self.__molecule_charge, self.__molecule_topo, atom_parametrization_list


    def set_topology(self, topo):
        self.__molecule_topo = topo


class generalizable_preprocessor(preprocessor):
    def __init__(self, verbose = False):
        preprocessor.__init__(self, verbose=verbose)
        self.custom_pars = None

    def read_custom_parameters(self, filename=None, dictionary=None):
        if filename is not None:
            a2md_pars = a2md_general_parameters(version=1.0, filename=filename)
        elif dictionary is not None:
            a2md_pars = a2md_general_parameters(version=1.0)
            a2md_pars.set_params_from_dict(dictionary)
        else:
            raise IOError
        a2md_pars.read_parameters()
        a2md_pars.set_support_library()
        self.custom_pars = a2md_pars


    def parametrize(self, kind = 'trigonometric'):
        """

        :param kind:
        :return:
        """
        if kind in ['proto', 'mol', 'bond', 'trigonometric', 'none']:
            return super(generalizable_preprocessor, self).parametrize(kind=kind)
        else:
            if self.custom_pars is None:
                raise IOError("custom parameters were not defined")
            cp = self.custom_pars # :type cp: a2md_general_parameters
            """:type cp: a2md_general_parameters"""
            atom_coordinates = self.get_coordinates()
            atom_labels = self.get_labels()
            topology = self.get_topology()
            charge = self.get_charge()
            A = []
            Xi = []
            C = []
            PN = []
            B = []
            number_atoms = atom_coordinates.shape[0]
            # ATOM LOOP
            for atom_index in range(number_atoms):
                atom_l = atom_labels[atom_index]
                atom_fun_set = cp.get_atom(atom_l)
                # SUPPORT FUNCTION LOOP
                for fun in atom_fun_set:
                    if fun['_BOND'] != "0":
                        # BONDED FUNCTIONS
                        for bond in topology[atom_index]:
                            if atom_labels[bond] == fun['_BOND'] or "A" == fun["_BOND"]:
                                A.append(atom_index)
                                Xi.append(fun['_PARS'])
                                PN.append((atom_index, fun['_NAME'], bond))
                                C.append(1.0)
                                B.append(bond)
                    else:
                        # NON-BONDED FUNCTIONS
                        A.append(atom_index)
                        Xi.append(fun['_PARS'])
                        PN.append((atom_index, fun['_NAME']))
                        C.append(1.0)
                        B.append(None)
            number_functions = len(A)

            atom_parametrization_list = [dict()] * number_functions
            # DICTIONARY OF PARAMETERS IN THE FORMAT THAT A2MD.AMD.AAMD LIKES
            for xi_iter in range(number_functions):
                tmp_input = dict(
                    center=A[xi_iter],
                    support_type=PN[xi_iter][1],
                    coefficient=C[xi_iter],
                )
                tmp_params = dict()
                for key, item in Xi[xi_iter].items():
                    tmp_params[key] = item
                tmp_input['params'] = tmp_params
                tmp_input['bond']   = B[xi_iter]
                atom_parametrization_list[xi_iter] = tmp_input
            return atom_coordinates, charge, topology, atom_parametrization_list


class model_parser(A2MD_basis):
    allowed_formats = ['json']
    def __init__(self, m, input_format='json'):
        """

        :param m:
        :param input_format:
        """
        A2MD_basis.__init__(self, name='model parser', verbose=False)
        if input_format == 'json':
            self.model_contents = self.read_json(m)
        elif input_format == 'dict':
            self.model_contents = m
        else:
            raise IOError("not allowed format. Please, use: {:s}".format(' '.join(self.allowed_formats)))


    @staticmethod
    def read_json(fname):
        import json
        with open(fname) as f:
            file_contents = json.load(f)
        return file_contents

    def build_model(
            self,
            coordinates,
            labels,
            topology,
            charges
    ):
        # Load model information
        model = self.model_contents['_MODEL']

        # Create the placeholders where the output data will be stored
        A = list()  # Function-atom map
        Xi = list()  # Parameters list
        PN = list()  # function-atom, function names, function-bond
        F = list()  # Frozen term (true/false)
        C = list()  # Coefficients of each function in the output

        # Loop along the atoms of the molecule
        for i, (element, atom_coords, atom_charge) in enumerate(zip(labels, coordinates, charges)):
            # Loop along the functions associated to each element
            for j, fun in enumerate(model[element]):

                if fun['_CONNECT'] == '_NONE':

                    A.append(i)
                    Xi.append(fun['_PARAMS'])
                    PN.append((i, fun['_NAME']))
                    F.append(fun['_FROZEN'])
                    C.append(fun['_COEFFICIENT'])

                elif fun['_CONNECT'] == '_TOPO':

                    for bond in topology[i]:

                        A.append(i)
                        Xi.append(fun['_PARAMS'])
                        PN.append((i, fun['_NAME'], bond))
                        F.append(fun['_FROZEN'])
                        C.append(fun['_COEFFICIENT'])

        return A, Xi, PN, F, C

class symmetry_features(A2MD_basis):
    """

    This class allows to define symmetry features between atoms and bonds, so their parameters
    can be collapsed, easing calculations.

    """
    def __init__(self, file=None):
        """

        :param file:
        """
        from a2mdlib.molecules import Mol2
        A2MD_basis.__init__(self, name='symfeats', verbose=False)
        if file is not None:
            self.__mol = Mol2(file=file, verbose=False)
            self.__mol.change_units('au')
        else:
            self.__mol = None
        self.__x = None
        self.__l = None

    @staticmethod
    def build_distance_matrix(coord_matrix):
        """

        :param coord_matrix:
        :return:
        """
        n = coord_matrix.shape[0]
        dm = np.zeros((n, n), dtype='float64')
        for i in range(n):
            d2 = coord_matrix - coord_matrix[i, :]
            dm[i, :] = np.sqrt(np.sum(d2 ** 2, axis=1))
        return dm

    @staticmethod
    def get_cosines(coord_matrix, i, j, k, t=None):
        """

        :param coord_matrix:
        :param i:
        :param j:
        :param k:
        :param t:
        :return:
        """
        r1 = coord_matrix[j, :] - coord_matrix[i, :]
        r2 = coord_matrix[k, :] - coord_matrix[i, :]
        cos = r1.dot(r2) / (np.linalg.norm(r1) * np.linalg.norm(r2))
        angle = np.arccos(cos)
        if t is None:
            return cos
        else:
            angle = angle - t
            return np.cos(angle)

    @staticmethod
    def cutoff(r, rc):
        """

        :param r:
        :param rc:
        :return:
        """
        res = np.zeros(r.shape, dtype='float64')
        res[r < rc] = 0.5 * (np.cos(np.pi * (r[r < rc] / rc)) + 1)
        return res

    def get_bonds(self):
        """

        :return:
        """
        from a2mdlib.molecules import Mol2
        if self.__mol is not None:
            assert isinstance(self.__mol, Mol2)
            topo = self.__mol.get_bonds()
            return topo

    def get_coordinates(self):
        """

        :return:
        """
        if self.__mol is not None:
            self.__mol.change_units('au')
            return self.__mol.get_coordinates()
        else:
            return self.__x

    def get_atomic_nums(self):
        """

        :return:
        """
        if self.__mol is not None:
            return self.__mol.get_atomic_numbers()
        else:
            return self.__l

    def get_total_charge(self):
        """

        :return:
        """
        if self.__mol is not None:
            return self.__mol.get_charge().sum()

    def feed_coordinates(self, x, l):
        """

        :param x:
        :param l:
        :return:
        """
        self.__x = x
        self.__l = l

    def get_radial_features(self, params):
        """

        :param params:
        :return:
        """
        coordinates = self.__mol.get_coordinates()
        dm = self.build_distance_matrix(coordinates)
        an = self.get_atomic_nums()
        feats = np.zeros((dm.shape[0], len(params)), dtype='float64')
        for i, (eta, rs, rc, a) in enumerate(params):
            co = self.cutoff(dm, rc) - np.identity(dm.shape[0])
            gr = np.sum(np.exp(- eta * ((dm - rs) * (dm - rs))) * co * (an == a), axis=1)
            feats[:, i] = gr
        return feats

    def get_angular_features(self, params):
        """

        :param params:
        :return:
        """
        coordinates = self.__mol.get_coordinates()
        dm = self.build_distance_matrix(coordinates)
        feats = np.zeros((dm.shape[0], len(params)), dtype='float64')
        for p, (theta, eta, rc, rs) in enumerate(params):
            co = self.cutoff(dm, rc)
            for i in range(dm.shape[0]):
                for j in range(dm.shape[0]):
                    for k in range(dm.shape[0]):
                        if i == j or j == k or i == k:
                            continue
                        else:
                            angular_dependent_term = (1 + self.get_cosines(coordinates, i, j, k, theta)) ** 6
                            radial_dependent_term = np.exp(-eta * (0.5 * (dm[i, j] + dm[i, k]) - rs) ** 2)
                            cutoff_term = co[i, j] * co[i, k]
                            prefactor = (2 ** -5)
                            feats[i, p] += prefactor * angular_dependent_term * radial_dependent_term * cutoff_term
        return feats
