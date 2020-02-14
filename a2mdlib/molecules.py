from a2mdlib import A2MDlib
import numpy as np

ALLOWED_ATOMS = ['C', 'N', 'O', 'H']
ATOM_LABELS = ['', 'h', '', '', '', '', 'c', 'n', 'o', 'f']
ATOMIC_NUMBERS = {"C": 6, "N": 7, "O": 8, "H": 1, "F": 9}

UNITS_TABLE = dict(
    angstrom=dict(
        au=1.8897261254535,
        nm=0.1
    ),
    au=dict(
        angstrom=0.52917721067121,
        nm=0.05291772106712
    ),
    nm=dict(
        au=18.897261254535,
        angstrom=10.0
    )
)


class Zinc15(A2MDlib):
    """
    Zinc 15

    the purpose of this class is to create a client that downloads information
    from the docking database ZINC15 by providing a list of ZINC15 identifiers.
    It works as.
        id_list = ['ZINC00000015','ZINC00000016',...]
        Zinc15(id_list)
        Zinc15.query(path='~/chemical_data/')
    See a2mdlib/test/test_zinc15.py
    """
    ZINC_url = 'http://zinc15.docking.org/substances/'

    def __init__(self, id_list, verbose=False):
        """

        :param id_list: a list of zinc15 codes
        :type id_list: list[str]
        :param verbose: display or not messages
        :type verbose: bool
        """
        import urllib3
        A2MDlib.__init__(self, name='zinc 15', verbose=verbose)
        self.__mol_input = id_list
        self.__http = urllib3.PoolManager()

    @staticmethod
    def __write(data, file):
        """

        :param data:
        :param file:
        :return:
        """
        with open(file, 'w') as f:
            f.write(data)

    def query(self, path="./"):
        """

        :param path: location where files should be saved
        :return: number of sucessful pulls
        """
        n_success = 0
        for i, mol in enumerate(self.__mol_input):
            self.log("downloading {0} from ZINC15".format(mol))
            r = self.__http.request('GET', self.ZINC_url + mol + '.mol2')
            if r.status == 200:
                self.log("\tsuccess downloading {0} - {1} len".format(mol, len(r.data)))
                self.__write(data=r.data.decode('UTF-8'), file=path + mol + '.mol2')
                n_success += 1
            else:
                self.log("\tcould not download {0} from ZINC15")
        self.log(
            "SUCCESSFUL PULLS : {:d} SUCESS RATIO : {:3.2f} %".format(
                n_success, (n_success / len(self.__mol_input))*100
            )
        )
        return n_success


class Mol2(A2MDlib):
    def __init__(self, file=None, basis='6-311G(d)', method='MP2', proc=8, units='angstrom', verbose=True):
        """

        :param file: file from which to read information
        :param basis: basis set used for quantum calculations. Only used if writting g09
        :param method: qm method used for quantum calculations. Only used if writting g09
        :param proc: number of processors to use for quantum calculations. Only used if writting g09
        :param units: either angstrom, au or nm. All other options will raise a IOError
        """
        A2MDlib.__init__(self, name='mol2 utils', verbose=verbose)
        self.__filename = file
        if units in UNITS_TABLE.keys():
            self.__units = units
        else:
            raise IOError(
                "please provide coordinates in format : {:s}".format(
                    ', '.join([i for i in UNITS_TABLE.keys()])
                )
            )
        if file is None:
            raise IOError
        else:
            coordinates, labs, top, charges, na, nb, st, bt = self.__read(fn=file)
        self.__coordinates = coordinates
        self.__labels = labs
        self.__topology = top
        self.__charges = charges
        self.__na = na
        self.__nb = nb
        self.__basis = basis
        self.__method = method
        self.__number_processors = proc
        self.__sybil = st
        self.__bond_type = bt

    def __calculate_total_charges(self):

        an = self.get_atomic_numbers()
        ch = self.__charges
        return an - ch

    @staticmethod
    def __coords2string(labels, matrix, sep=' '):
        """

        :param labels:
        :param matrix:
        :param sep:
        :return:
        """
        shape = matrix.shape
        buffer = [''] * shape[0]
        for i in range(shape[0]):
            tmp = ['{:12.6f}'.format(j) for j in matrix[i, :]]
            buffer[i] = labels[i] + sep + sep.join(tmp)
        return '\n'.join(buffer)

    @staticmethod
    def __read(fn):
        """

        :param fn:
        :return:
        """
        with open(fn) as f:
            raw_text = f.read()
            mol, atom, bond = raw_text.split(r'@')[1:4]

        if len(mol) == 0 or len(atom) == 0 or len(bond) == 0:
            raise IOError("the mol2 file could not be read")

        # PARSING MOLECULE PART
        head, name, numbers = mol.split('\n')[:3]
        numbers = numbers.split()
        number_atoms = int(numbers[0])
        number_bonds = int(numbers[1])

        # DEFINING MATRICES FOR COORDINATES AND LISTS FOR TOPOLOGY

        coordinate_matrix = np.zeros((number_atoms, 3), dtype='float64')
        labels_list = ['X'] * number_atoms
        bonding_list = np.zeros((number_bonds, 2), dtype='int32')
        charge_list = np.zeros(number_atoms, dtype='float64')
        sybil_type = [''] * number_atoms
        bond_type = [''] * number_bonds
        # PARSING ATOM PART

        for line in atom.split('\n'):
            try:
                split_line = line.split()
                atom_id, label, x, y, z = split_line[:5]
                charge = float(split_line[-1])
                sybil = split_line[5]
            except ValueError:
                continue
            atom_id = int(atom_id)
            coordinate_matrix[atom_id - 1, :] = float(x), float(y), float(z)
            labels_list[atom_id - 1] = str(label[0])
            charge_list[atom_id - 1] = charge
            sybil_type[atom_id - 1] = sybil

        # PARSING BONDING PART

        for line in bond.split('\n'):
            try:
                bond_id, atom1, atom2, bt = line.split()
            except ValueError:
                continue
            bond_id = int(bond_id)
            bonding_list[bond_id - 1, :] = int(atom1), int(atom2)
            bond_type[bond_id - 1] = bt

        return coordinate_matrix, labels_list, bonding_list, charge_list, number_atoms, number_bonds, \
            sybil_type, bond_type

    def __get_bound(self, atom_number):
        """

        :param atom_number:
        :return:
        """
        atom_number = atom_number + 1
        bound_list = []
        for bond_index in range(self.__nb):
            if self.__topology[bond_index, 0] == atom_number:
                bound_list.append("{:d}".format(self.__topology[bond_index, 1]-1))
            elif self.__topology[bond_index, 1] == atom_number:
                bound_list.append("{:d}".format(self.__topology[bond_index, 0]-1))
        return ','.join(bound_list)

    def __get_net_charge(self):
        """

        :return:
        """
        # expected_charge = np.sum([ATOMIC_CHARGE_DICT['elem_%1s' % i] for i in self.__labels])
        present_charge = np.sum(self.__charges)
        # dif = present_charge - expected_charge
        return int(np.round(present_charge))

    def __write_aamd(self):
        """

        :return:
        """
        if self.__units != 'au':
            raise IOError("WARNING : aamd default units are atomic units. this can give place to awfuld headaches!")
        amd_str_contents = '#id\tx({0})\ty({0})\tz({0})\tsymbol\ttopo\tq(d)\n'.format(self.__units)
        charge = self.__calculate_total_charges()
        for i in range(self.__na):
            amd_str_contents = amd_str_contents + '{:2d}\t{:8.6f}\t{:8.6f}\t{:8.6f}\t{:s}\t{:s}\t{:8.6f}\n'.format(
                (i + 1),  # index
                self.__coordinates[i, 0], self.__coordinates[i, 1], self.__coordinates[i, 2],  # coordinates x y z
                self.__labels[i], self.__get_bound(i), charge[i]  # element bounded elements charges
            )
        return amd_str_contents

    def __write_mol2(self):
        """

        :return:
        """
        import os
        if self.__units != 'angstrom':
            self.log("WARNING : Mol2 default units are angstroms. This can give place to awful headaches!")
        molecule_str = ""
        atoms_str = ""
        bonds_str = ""
        molecule_str = molecule_str + "@<TRIPOS>MOLECULE\n"
        molecule_str = molecule_str + os.path.basename(self.__filename).replace('.mol2', '') + "\n"
        molecule_str = molecule_str + "{:d} {:d} 0 0 0\n".format(self.__na, self.__nb)

        if self.__units != 'angstrom':
            comment_str = "@<TRIPOS>ATOM\n"
            comment_str = comment_str + "coordinates are in {:s}\n".format(self.__units)
        else:
            comment_str = None

        atoms_str = atoms_str + "@<TRIPOS>ATOM\n"
        for i in range(self.__na):
            atoms_str = atoms_str + "    {:d} {:<15s} {:>8.4f} {:>8.4f} {:>8.4f} {:<6s} 1 {:<6s} {:>8.4f}\n".format(
                i + 1, self.__labels[i], self.__coordinates[i, 0], self.__coordinates[i, 1], self.__coordinates[i, 2],
                self.__sybil[i], 'LIG1', self.__charges[i]
            )

        bonds_str = bonds_str + "@<TRIPOS>BOND\n"
        for i in range(self.__nb):
            a1, a2 = self.__topology[i, :]
            bt = self.__bond_type[i]
            bonds_str = bonds_str + "    {:d} {:d} {:d} {:s}\n".format(
                i + 1, a1, a2, bt
            )

        if self.__units != 'angstrom':
            mol2_contents_list = [molecule_str, comment_str, atoms_str, bonds_str]
        else:
            mol2_contents_list = [molecule_str, atoms_str, bonds_str]

        return '\n'.join(mol2_contents_list)

    def __write_gaussian09(self, filename, calculation_type='single'):
        """

        :param filename:
        :param calculation_type:
        :return:
        """
        import os
        filename = os.path.basename(filename)
        if calculation_type not in ['opt', 'single', 'opt-single']:
            raise NotImplementedError

        if self.__units == 'au':
            self.log("WARNING: since default units are not angstroms, a command indicating this unit will be included")
        elif self.__units not in ['au', 'angstrom']:
            raise IOError("gaussian does not accept these units")

        net_charge = self.__get_net_charge()

        if not calculation_type == 'opt-single':

            gaussian_input_contents = [''] * 6
            gaussian_input_contents[0] = '%NProcShared={:d}\n'.format(self.__number_processors)

            if self.__units == 'au':
                if calculation_type == 'single':
                    gaussian_input_contents[1] = '#{:s}/{:s} output=wfn density=current pop=npa units=(au)\n'.format(
                        self.__method, self.__basis
                    )
                elif calculation_type == 'opt':
                    gaussian_input_contents[1] = '#{:s}/{:s} Opt units=(au)\n'.format(
                        self.__method, self.__basis
                    )
            else:
                if calculation_type == 'single':
                    gaussian_input_contents[1] = '#{:s}/{:s} output=wfn density=current pop=npa\n'.format(
                        self.__method, self.__basis
                    )
                elif calculation_type == 'opt':
                    gaussian_input_contents[1] = '#{:s}/{:s} Opt\n'.format(
                        self.__method, self.__basis
                    )
            gaussian_input_contents[2] = filename.replace('.g09', '') + '\n'
            gaussian_input_contents[3] = '{:1d}, 1\n'.format(net_charge)
            gaussian_input_contents[4] = self.__coords2string(self.__labels, self.__coordinates) + '\n'
            if calculation_type == 'single':
                gaussian_input_contents[5] = filename.replace('.g09', '.wfn') + '\n'
            # ASSEMBLING FILE
            gaussian_str_contents = gaussian_input_contents[0] + gaussian_input_contents[1] + '\n'
            gaussian_str_contents = gaussian_str_contents + gaussian_input_contents[2] + '\n'
            gaussian_str_contents = gaussian_str_contents + gaussian_input_contents[3] + gaussian_input_contents[4]
            if calculation_type == 'single':
                gaussian_str_contents = gaussian_str_contents + '\n' + gaussian_input_contents[5] + '\n'
            else:
                gaussian_str_contents = gaussian_str_contents + '\n'
        else:

            check_file_name = filename.replace('.g09', '_chk')

            gaussian_input_contents = [''] * 11
            gaussian_input_contents[0] = '%NProcShared={:d}'.format(self.__number_processors)
            gaussian_input_contents[1] = '%Chk={:s}'.format(check_file_name)
            if self.__units == 'au':
                gaussian_input_contents[2] = '#T B3LYP/6-311++G opt units=(au)\n'.format(self.__method, self.__basis)
            elif self.__units == 'angstrom':
                gaussian_input_contents[2] = '#T B3LYP/6-311++G opt\n'.format(self.__method, self.__basis)
            else:
                raise IOError("gaussian does not accept these units")
            gaussian_input_contents[3] = '{:s}\n'.format(filename.replace('.g09', '_opt'))
            gaussian_input_contents[4] = '{:1d}, 1'.format(net_charge)
            gaussian_input_contents[5] = '{:s}\n'.format(self.__coords2string(self.__labels, self.__coordinates))
            gaussian_input_contents[6] = '--Link1--'
            gaussian_input_contents[7] = '{:s}\n{:s}\n#T {:s}/{:s} geom=check guess=read output=wfn density=current ' \
                                         'pop=npa\n'.format(
                gaussian_input_contents[0], gaussian_input_contents[1], self.__method, self.__basis
            )
            gaussian_input_contents[8] = '{:s}\n'.format(filename.replace('.g09', ''))
            gaussian_input_contents[9] = '{:1d}, 1\n'.format(net_charge)
            gaussian_input_contents[10] = '{:s}\n'.format(filename.replace('.g09', '.wfn'))

            gaussian_str_contents = '\n'.join(gaussian_input_contents)

        return gaussian_str_contents

    def calculate_partial_charges(self, charge):
        """

        The purpose of this function is to allow reparametrization of Mol2 files using charges
        obtained by other methods. For instance, it is usual to consider total charges when reading
        from Gaussian outputs. However, Mol2 always use partial charges. This function
        allows the conversion from total charges to partial charges

        :param charge:
        :return:
        """
        an = self.get_atomic_numbers()
        return an - charge

    def change_units(self, units):
        """

        this function allows to convert the units of the coordinates matrix. When a
        new unit is indicated, a change of units will be performed by consulting the UNITS_TABLE.
        New units and their conversions can be placed in that table to allow conversion.


        :param units:
        :return:
        """
        if units == self.__units:
            self.log("WARNING: old and new units are the same")
            return
        try:
            factor = UNITS_TABLE[self.__units][units]
        except KeyError:
            raise NotImplementedError("unit conversion not in table")
        self.log(
            "change of coordinates from {:s} to {:s}. Multiplying by {:8.4f}".format(
                self.__units, units, factor
            )
        )
        self.__coordinates *= factor
        self.__units = units
        return

    def get_atomic_numbers(self):
        labels = self.__labels
        return np.array([ATOMIC_NUMBERS[i] for i in labels])

    def get_bonds(self):
        return self.__topology

    def get_basis(self):
        return self.__basis

    def get_bond_type(self):
        return self.__bond_type

    def get_charge(self, kind='partial'):
        if kind == 'partial':
            return self.__charges
        elif kind == 'total':
            return self.__calculate_total_charges()

    def get_coordinates(self):
        return self.__coordinates

    def get_labels(self):
        return self.__labels

    def get_number_atoms(self):
        return self.__na

    def get_number_bonds(self):
        return self.__nb

    def get_sybil(self):
        return self.__sybil

    def get_units(self):
        return self.__units

    def sample_probes(self, region='box', resolution=0.1):
        """

        This function should allow to obtain samples of points around the atomic
        nuclei of the molecule. Notice that by now, only the "box" sampling
        is implemented

        :param region:
        :param resolution:
        :return:
        """
        if region not in ['box', 'bonding']:
            raise NotImplementedError
        if region == 'box':

            pass

        elif region == 'bonding':

            sampling = []
            for bond in self.__topology:
                bond_vector = self.__coordinates[bond[1] - 1] - self.__coordinates[bond[0] - 1]

                bond_distance = np.linalg.norm(bond_vector)
                bond_vector /= bond_distance

                sampled_distance = bond_distance / 0.75
                app_point = self.__coordinates[bond[0] - 1]
                # sampling.append(bond_vector.reshape(1, 3) + app_point)
                rot_matrix = np.zeros([3, 3], dtype='float64')
                rot_matrix[0, :] = bond_vector
                rot_matrix[1, :] = np.cross(
                    bond_vector, np.array([0.0, 0.0, 1.0], dtype='float64')
                )
                rot_matrix[2, :] = np.cross(
                    bond_vector, rot_matrix[:, 1]
                )

                x = np.arange(0, sampled_distance + resolution, resolution)

                buffer_ = np.zeros([x.size, 3], dtype='float64')
                buffer_[:, 0] = x
                rot_buffer = buffer_.copy().dot(rot_matrix) + app_point

                sampling.append(rot_buffer)
                for ir, rad in enumerate(np.arange(0.1, 1.3, 0.3)):
                    for ia, ang in enumerate(np.arange(0.0, 2*np.pi, np.pi/25)):
                        y = rad * np.sin(ang)
                        z = rad * np.cos(ang)
                        buffer_[:, 1] = y
                        buffer_[:, 2] = z
                        rot_buffer = buffer_.copy().dot(rot_matrix) + app_point
                        sampling.append(rot_buffer)

            stacked_sample = np.concatenate(sampling, axis=0)
            return stacked_sample

    def set_basis(self, basis):
        """

        :param basis:
        :return:
        """
        self.__basis = basis

    def set_bonds(self, bonds):
        """

        :param bonds:
        :return:
        """
        if isinstance(bonds, list):
            bonds = np.array(bonds, dtype='int64')
        self.__topology = bonds

    def set_coordinates(self, coordinates, units):
        """

        This method allows to change the coordinates matrix of the atom nuclei.
        It is mandatory to specify the coordinates of the new coordinates, to avoid
        confussions between files whose coordinates are in angstroms and files whose
        coordinates are in atomic units

        :param coordinates:
        :param units:
        :return:
        """
        if units in UNITS_TABLE.keys():
            self.__units = units
        else:
            raise IOError(
                "please provide coordinates in format : {:s}".format(
                    ', '.join([i for i in UNITS_TABLE.keys()])
                )
            )
        self.__coordinates = coordinates

    def set_charges(self, charges, kind='partial'):
        """

        :param charges:
        :param kind:
        :return:
        """
        if kind == 'partial':
            self.__charges = charges
        elif kind == 'total':
            self.__charges = self.calculate_partial_charges(charges)
        else:
            raise NotImplementedError("use partial or total")

    def set_method(self, method):
        """

        :param method:
        :return:
        """
        self.__method = method

    def set_number_processors(self, number_processors):
        """

        :param number_processors:
        :return:
        """
        self.__number_processors = number_processors

    def write(self, file, output_format='g09', calculation_type='single'):
        """

        :param file:
        :param output_format: Mol2 objects can write aamd files and gaussian inputs
        :param calculation_type: There are three calculation types predefined (single, opt, opt-single)
        :return:
        """
        if output_format == 'g09':
            str2write = self.__write_gaussian09(
                file,
                calculation_type=calculation_type
            )
        elif output_format == 'aamd':
            str2write = self.__write_aamd()
        elif output_format == 'mol2':
            str2write = self.__write_mol2()
        else:
            raise IOError
        with open(file, 'w') as f:
            f.write(str2write)


class Cleaner(A2MDlib):
    """
    The purpose of this class is to extract a subset from gdb9
    """

    def __init__(self, verbose=True):
        A2MDlib.__init__(self, name='cleaner', verbose=verbose)

    def clean(self, file, path, allowed_atoms=ALLOWED_ATOMS):
        """

        This object takes a large sdf file and extracts all molecules
        fulfilling some conditions:
        1. Not having any ring
        2. Not having any triple bond
        3. Not having any atom that is not included in the  allowed_atoms list
        4. That can be read by rdkit (some molecules can not be read!)

        :param file:
        :param path:
        :param allowed_atoms: list of atomic numbers
        :return:
        """
        from rdkit.Chem import AllChem as Chem
        supplier = Chem.SDMolSupplier(file, removeHs=False)
        self.log('sdf was read')
        for molecule in supplier:
            if molecule is None:
                continue
            name = molecule.GetProp('_Name')
            try:
                symbol_list = [i.GetSymbol() for i in molecule.GetAtoms()]
            except AttributeError:
                self.log("there was some issue with molecule {:s}".format(name))
                continue
            if all([i in allowed_atoms for i in symbol_list]):
                ri = molecule.GetRingInfo()
                if ri.NumRings() < 1:

                    triple_bond_control = False
                    for bond in molecule.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
                            triple_bond_control = True
                            break
                    if triple_bond_control:
                        continue
                    self.log("writting {:s} to {:s}".format(name, path))
                    Chem.MolToMolFile(
                        molecule,
                        "{:s}\{:s}.mol".format(path, name)
                    )
                else:
                    self.log("molecule {:s} was not discarded since it had rings".format(name))
            else:
                self.log("molecule {:s} was not discarded since it has not allowed elements".format(name))
