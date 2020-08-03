from a2mdio import A2MDlib
from a2mdio import PDB_PROTEIN_TOPOLOGY, PDB_PROTEIN_TYPES, PDB_PROTEIN_CHARGES, PDB_PROTEIN_TYPE_CHARGES
import numpy as np

ALLOWED_ATOMS = ['C', 'N', 'O', 'H']
ATOM_LABELS = [
    '',
    'h', '',
    '', '', '', 'c', 'n', 'o', 'f', '',
    '', '', '', 'si', 'p', 's', 'cl', ''
]
ATOMIC_NUMBERS = {
    "H": 1, "He":2,
    "Li": 3, "Be":4, "B":5, "C": 6, "N": 7, "O": 8,  "F": 9, "Ne":10,
    "Na":11, "Mg":12, "Al":13, "Si":14, "P":15, "S":16, "Cl":17, "Ar":18,
    "K":19, "Ca":20, "Ga":31, "Ge":32, "As":33, "Se":34, "Br":35, "Kr":36
}

UNITS_TABLE = dict(
    angstrom=dict(
        au=1.8897261254535,
        nm=0.1,
        angstrom=1.0
    ),
    au=dict(
        angstrom=0.52917721067121,
        nm=0.05291772106712,
        au=1.0
    ),
    nm=dict(
        au=18.897261254535,
        angstrom=10.0,
        nm=1.0
    )
)

class QmSetUp(A2MDlib):
    def __init__(
            self, basis, method, calculation_type='single', nprocs=1, disk=None, memory=None,
            additional_commands=None, verbose=None
    ):
        A2MDlib.__init__(self, name='G09 input writer', verbose=verbose)
        self.basis = basis
        self.method = method
        self.nprocs = nprocs
        self.disk = disk
        self.memory = memory
        self.calculation_type = calculation_type
        self.additional_commands = additional_commands

    @staticmethod
    def get_mol_info(mol):
        assert issubclass(type(mol), MolRepresentation)
        try:
            coordinates = mol.get_coordinates()
            labels = mol.get_symbols()
            size = mol.get_number_atoms()
            total_charge = mol.get_total_charge()
            units = mol.get_units()
            multiplicity = mol.get_multiplicity()
        except AttributeError:
            raise IOError("some of the requested fields was missing")
        return size, total_charge, units, multiplicity, coordinates, labels
    def write_g09(self, file, mol, wfn=False):
        import os

        size, total_charge, units, multiplicity, coordinates, labels = self.get_mol_info(mol)
        size: int

        file_basename = os.path.basename(file)
        if self.calculation_type not in ['opt', 'single']:
            raise NotImplementedError("use opt/single as calculation types")

        coordinates = coordinates * UNITS_TABLE[units]['angstrom']
        self.log("coordinates were transformed from {:s} to angstrom ".format(units))


        proc_line = '%NProcShared={:d}\n'.format(self.nprocs)

        if self.disk is not None :
            disk_line = '%MaxDisk={:s}\n'.format(self.disk)
        else:
            disk_line = ''

        if self.memory is not None :
            mem_line = '%MaxMem={:s}\n'.format(self.memory)
        else:
            mem_line = ''
        if self.calculation_type is "opt":
            calculation_type_str = 'opt'
        else:
            calculation_type_str = ''
        g09_command_line = '# {:s}/{:s} {:s}'.format(self.method, self.basis, calculation_type_str)
        if self.additional_commands is not None:
            g09_command_line = g09_command_line + ' {:s}\n\n'.format(' '.join(self.additional_commands))
        else:
            g09_command_line = g09_command_line + '\n\n'

        g09_run_name_line = '{:s}\n\n'.format(file_basename)
        g09_mult_and_charge_line = '{:d},{:d}\n'.format(total_charge, multiplicity)

        g09_coords_matrix = []
        for i in range(size):
            g09_coords_matrix.append('{:s} {:12.4f} {:12.4f} {:12.4f}'.format(
                labels[i], coordinates[i, 0], coordinates[i, 1], coordinates[i, 2])
            )
        g09_coords_str = '\n'.join(g09_coords_matrix) + '\n'
        if wfn:
            g09_output_line = '\n{:s}\n'.format(file.replace('.g09.output', '') + '.wfn')
        else:
            g09_output_line = '\n'

        g09_str = ''.join(
            [
                proc_line, disk_line, mem_line, g09_command_line, g09_run_name_line, g09_mult_and_charge_line,
                g09_coords_str, g09_output_line
            ]
        )

        with open(file, 'w') as f:
            f.write(g09_str)

        return g09_str
        
    def write_orca(self, filename, mol):
        import os
        size, total_charge, units, multiplicity, coordinates, labels = self.get_mol_info(mol)
        size: int
        filename = os.path.basename(filename)
        options = []
        if self.calculation_type == "opt":
            options.append('Opt')
        if self.nprocs > 1:
            options.append('PAL{:d}'.format(self.nprocs))
        header = "! {:s} {:s} {:s}".format(
            self.method, self.basis, ' '.join(options)
        )
        name = r"#" + "  {:s}".format(filename)
        xyz = ["*xyz {:d} {:d}".format(total_charge, multiplicity)]
        for i in range(size):
            xyz.append(
                "{:s} {:12.6f} {:12.6f} {:12.6f}".format(
                    labels[i], coordinates[i, 0],
                    coordinates[i, 1], coordinates[i, 2]
                )
            )
        xyz.append("*")
        xyz = '\n'.join(xyz)
        orca_input_contents = '\n'.join([header, name, xyz])
        with open(filename, 'w') as f:
            f.write(orca_input_contents)

        return orca_input_contents


class MolRepresentation(A2MDlib):
    def __init__(
            self, name, units, number_atoms, number_bonds, atomic_numbers, coordinates,
            charges, atom_types, atom_names, bonds, bond_types

    ):
        """

        :param name:
        :param units:
        :param number_atoms:
        :param number_bonds:
        :param coordinates:
        :param atomic_numbers:
        :param charges:
        :param atom_types:
        :param bond_types:
        """
        A2MDlib.__init__(self, name=name, verbose=False)
        self.coordinates = coordinates
        self.units = units
        self.natoms = number_atoms
        self.nbonds = number_bonds
        self.atomic_numbers = atomic_numbers
        self.charges = charges
        self.atom_types = atom_types
        self.atom_names = atom_names
        self.bonds = bonds
        self.bond_types = bond_types
        self.multiplicity = 1

    def get_bonds(self):
        return self.bonds

    def get_coordinates(self, units=None):
        if units is None: units = self.units
        return self.coordinates * UNITS_TABLE[self.units][units]

    def get_absolute_charges(self):
        return self.atomic_numbers - self.charges

    def get_partial_charges(self):
        return self.charges

    def get_atomic_numbers(self):
        return self.atomic_numbers

    def get_symbols(self):
        return [ATOM_LABELS[i].upper() for i in self.atomic_numbers]

    def get_atom(self, idx):
        return dict(
            atomic_number=self.atomic_numbers[idx],
            symbol=ATOM_LABELS[self.atomic_numbers[idx]].upper(),
            coordinates=self.coordinates[idx, :],
            atom_type=self.atom_types[idx],
            atom_name=self.atom_names[idx],
            charge=self.charges[idx]
        )
    def get_bond(self, idx):
        return dict(
            bonded_atoms=self.bonds[idx],
            bonded_elements=[
                self.atomic_numbers[self.bonds[idx][0]],
                self.atomic_numbers[self.bonds[idx][1]],
            ],
            bond_type=self.bond_types[idx],
            distance=np.linalg.norm(
                self.coordinates[self.bonds[idx][0]] - self.coordinates[self.bonds[idx][1]]
            )
        )

    def get_number_atoms(self):
        return self.natoms

    def get_number_bonds(self):
        return self.nbonds

    def get_total_charge(self):
        present_charge = np.sum(self.charges)
        total_charge =  int(np.round(present_charge, decimals=0))
        if total_charge - present_charge > 0.1:
            self.log("review charges; their sum is not accurate integer")
        return total_charge

    def get_multiplicity(self):
        return self.multiplicity

    def get_units(self):
        return self.units

class Mol2(MolRepresentation):
    def __init__(self, file=None):
        """

        :param file: file from which to read information
        """

        coordinates, atomic_numbers, bonds, charges, natoms, nbonds, \
        atom_types, atom_names, bond_types, segments, segment_idx = self.__read(fn=file)

        MolRepresentation.__init__(
            self,
            name='mol2 reading instance',
            units='angstrom',
            number_atoms=natoms, number_bonds=nbonds,
            coordinates=coordinates, atomic_numbers=atomic_numbers,
            charges=charges, atom_types=atom_types, atom_names=atom_names,
            bond_types=bond_types, bonds=bonds
        )
        self.file = file
        self.segments = segments
        self.segment_idx = segment_idx

    @staticmethod
    def __read(fn):
        """

        :param fn:
        :return:
        """
        with open(fn) as f:
            raw_text = f.read()
            mol, atom, bond = raw_text.split(r'@<TRIPOS>')[1:4]

        # PARSING MOLECULE PART
        head, name, numbers = mol.split('\n')[:3]
        numbers = numbers.split()
        number_atoms = int(numbers[0])
        number_bonds = int(numbers[1])

        # DEFINING MATRICES FOR COORDINATES AND LISTS FOR TOPOLOGY

        coordinate_matrix = np.zeros((number_atoms, 3), dtype='float64')
        atom_names = ['X'] * number_atoms
        atom_types = [''] * number_atoms
        atom_charges = np.zeros(number_atoms, dtype='float64')
        bonds = np.zeros((number_bonds, 2), dtype='int32')
        bond_type = [''] * number_bonds
        atom_segments = [''] * number_atoms
        atom_segments_idx = [0] * number_atoms
        segment_idx = 0

        for line in atom.split('\n'):
            try:
                split_line = line.split()
                atom_id, label, x, y, z = split_line[:5]
                charge = float(split_line[-1])
                tripos_type = split_line[5]
                segment = split_line[7]
            except ValueError:
                continue

            atom_id = int(atom_id)
            coordinate_matrix[atom_id - 1, :] = float(x), float(y), float(z)
            atom_names[atom_id - 1] = label
            atom_charges[atom_id - 1] = charge
            atom_types[atom_id - 1] = tripos_type

            try:
                segment_flag = segment == last_atom_segment
            except NameError:
                last_atom_segment = segment
                segment_flag = segment == last_atom_segment

            if not segment_flag:
                segment_idx += 1
                last_atom_segment = segment
            atom_segments[atom_id - 1] = segment
            atom_segments_idx[atom_id - 1] = segment_idx

        # PARSING BONDING PART

        for line in bond.split('\n'):
            try:
                bond_id, atom1, atom2, bt = line.split()
            except ValueError:
                continue
            bond_id = int(bond_id)
            bonds[bond_id - 1, :] = int(atom1), int(atom2)
            bond_type[bond_id - 1] = bt

        element_names = [i.split('.')[0] for i in atom_types]
        atomic_numbers = np.array([ATOMIC_NUMBERS[i] for i in element_names], dtype='int32')

        return coordinate_matrix, atomic_numbers, bonds, atom_charges, number_atoms, number_bonds, \
            atom_types, atom_names, bond_type, atom_segments, atom_segments_idx

    def write(self, file):
        """

        :return:
        """
        import os
        molecule_str = ""
        atoms_str = ""
        bonds_str = ""
        molecule_str = molecule_str + "@<TRIPOS>MOLECULE\n"
        molecule_str = molecule_str + os.path.basename(file).replace('.mol2', '') + "\n"
        molecule_str = molecule_str + "{:d} {:d} 0 0 0\n".format(self.natoms, self.nbonds)

        coordinates = self.coordinates * UNITS_TABLE[self.units]['angstrom']

        atoms_str = atoms_str + "@<TRIPOS>ATOM\n"
        for i in range(self.natoms):
            atoms_str = atoms_str + "    {:d} {:<15s} {:>8.4f} {:>8.4f} {:>8.4f} {:<6s} 1 {:<6s} {:>8.4f}\n".format(
                i + 1, self.atom_names[i], coordinates[i, 0], coordinates[i, 1], coordinates[i, 2],
                self.atom_types[i], self.segments[i], self.charges[i]
            )

        bonds_str = bonds_str + "@<TRIPOS>BOND\n"
        for i in range(self.nbonds):
            a1, a2 = self.bonds[i, :]
            bt = self.bond_types[i]
            bonds_str = bonds_str + "    {:d} {:d} {:d} {:s}\n".format(
                i + 1, a1, a2, bt
            )

        mol2_contents_list = [molecule_str, atoms_str, bonds_str]

        with open(file, 'w') as f:
            f.write('\n'.join(mol2_contents_list))

        return '\n'.join(mol2_contents_list)


class PDBLine:
    def __init__(self, line):
        self.atom_name = line[12:16].strip()
        # atom_alternate_location = line[16]
        self.atom_residue_name = line[17:20].strip()
        self.atom_chain = line[21].strip()
        self.atom_residue_idx = int(line[22:26])
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.atom_element = line[76:78].strip()
        self.atom_occupancy = float(line[54:60])
        self.atom_bfactor = float(line[60:66])
        self.charge = line[78:80].strip()


class PQRLine:
    def __init__(self, line):
        fields = line.split()
        self.atom_name = fields[2]
        # atom_alternate_location = line[16]
        self.atom_residue_name = fields[3]
        self.atom_chain = fields[4]
        self.atom_residue_idx = int(fields[5])
        self.x = float(fields[6])
        self.y = float(fields[7])
        self.z = float(fields[8])
        self.charge = float(fields[9])
        self.radius = float(fields[10])


class PDB(MolRepresentation):
    topology_source = PDB_PROTEIN_TOPOLOGY
    atom_types_source = PDB_PROTEIN_TYPES
    residue_charge_source = PDB_PROTEIN_CHARGES
    residue_charge_type_source = PDB_PROTEIN_TYPE_CHARGES
    def __init__(self, file, input_format='pdb'):
        if input_format in ['pdb', 'pqr']:
            self.input_format = input_format
        else:
            raise IOError("unrecognized input format, {:s}".format(input_format))
        na, nb, nr, nc, coordinates, topo, \
        atom_idx, atomic_numbers, atom_names, atom_residues_idx, atom_chains, atom_charges, \
        residues_idx, residues_names, residues_charges, residues_chains, residues_extent, \
        chain_names = self.read(file=file)

        MolRepresentation.__init__(
            self, name='pdb reading instance',
            units='angstrom', number_atoms=na, number_bonds=nb, atomic_numbers=atomic_numbers,
            coordinates=coordinates, charges=atom_charges, atom_types=atom_names,
            atom_names=atom_names,
            bonds=topo, bond_types=None
        )
        self.atom_idx = atom_idx
        self.atom_residues_idx = atom_residues_idx
        self.residue_idx = residues_idx
        self.residue_names = residues_names
        self.residue_charges = residues_charges
        self.residue_chains = residues_chains
        self.residue_extent = residues_extent
        self.chain_names = chain_names
        self.anotation = dict()
        self.sequence = self.make_sequences(self.residue_names, self.residue_chains)


    def read(self, file):


        with open(file) as f:
            contents = f.readlines()

        atom_lines = [i for i in contents if i[:4] == "ATOM"]

        # ATOM ATTRIBUTES
        atom_coordinates = []
        atom_elements = []
        atom_names = []
        atom_idx = []
        atom_residues_idx = []
        atom_chains = []
        atom_charges = []
        # RESIDUES ATTRIBUTES
        residues_idx = []
        residues_names = []
        residues_charges = []
        residues_chains = []
        residues_extent = []
        # TOPOLOGY
        topo = []

        # CHAIN ATTRIBUTES
        chain_names = []
        # ITERATION
        _residue_start_ = 0
        current_chain = ''
        current_residue = 0
        for i, line in enumerate(atom_lines):

            _atom_idx_ = i
            if self.input_format == 'pdb':
                line = PDBLine(line)
            elif self.input_format == 'pqr':
                line = PQRLine(line)
            else:
                raise IOError("unrecognized input format, {:s}".format(self.input_format))

            _atom_element_ = self.atom_types_source[line.atom_residue_name][line.atom_name]
            atom_coordinates.append([line.x, line.y, line.z])
            atom_elements.append(_atom_element_)
            atom_names.append(line.atom_name)
            atom_residues_idx.append(line.atom_residue_idx)
            atom_idx.append(_atom_idx_)
            atom_charges.append(self.assign_partial_charges(line.atom_name, line.atom_residue_name))
            atom_chains.append(line.atom_chain)

            if line.atom_residue_idx != current_residue:
                _residue_end_ = i

                if i > 1:
                    residues_extent.append((_residue_start_, _residue_end_))
                    topo = topo + self.set_residue_topology(
                        resname=residues_names[-1],
                        atom_idx=atom_idx[_residue_start_:_residue_end_],
                        names=atom_names[_residue_start_:_residue_end_]
                    )

                _residue_start_ = i
                residues_idx.append(line.atom_residue_idx)
                residues_names.append(line.atom_residue_name)
                residues_charges.append(self.residue_charge_source[line.atom_residue_name]['charge'])
                residues_chains.append(line.atom_chain)
                current_residue = line.atom_residue_idx

                if len(residues_extent) > 1 and residues_chains[-1] == residues_chains[-2]:
                    topo = topo + self.join_peptidic_bonds(
                        names_1=atom_names[residues_extent[-2][0]:residues_extent[-2][1]],
                        atom_idx_1=atom_idx[residues_extent[-2][0]:residues_extent[-2][1]],
                        names_2=atom_names[residues_extent[-1][0]:residues_extent[-1][1]],
                        atom_idx_2=atom_idx[residues_extent[-1][0]:residues_extent[-1][1]]
                    )

            if line.atom_chain != current_chain:
                chain_names.append(line.atom_chain)
                current_chain = line.atom_chain
        _residue_end_ = i
        residues_extent.append((_residue_start_, _residue_end_))
        topo = topo + self.set_residue_topology(
            resname=residues_names[-1],
            atom_idx=atom_idx[_residue_start_:_residue_end_],
            names=atom_names[_residue_start_:_residue_end_]
        )

        coordinates = np.array(atom_coordinates, dtype='float64')
        atomic_numbers = np.array(atom_elements, dtype='int64')
        na = coordinates.shape[0]
        nb = len(topo)
        nr = len(residues_names)
        nc = len(chain_names)

        return na, nb, nr, nc, coordinates, topo, \
               atom_idx, atomic_numbers, atom_names, atom_residues_idx, atom_chains, atom_charges, \
               residues_idx, residues_names, residues_charges, residues_chains, residues_extent, \
               chain_names

    def set_residue_topology(self, resname, atom_idx, names):
        topology_buffer = []
        res_topo = self.topology_source[resname]
        for i, (atom_id1, name1) in enumerate(zip(atom_idx, names)):
            for j, (atom_id2, name2) in enumerate(zip(atom_idx, names)):
                joint_name = "|".join([name1, name2])
                if joint_name in res_topo:
                    topology_buffer.append([atom_id1, atom_id2])
        return topology_buffer

    @staticmethod
    def join_peptidic_bonds(names_1, atom_idx_1, names_2, atom_idx_2):
        topology_buffer = []
        try:
            c_idx = names_1.index("C")
        except ValueError:
            raise RuntimeError("could not find the carboxyl carbon within selection")
        try:
            n_idx = names_2.index("N")
        except ValueError:
            raise RuntimeError("could not find the carboxyl carbon within selection")
        topology_buffer.append([atom_idx_1[c_idx], atom_idx_2[n_idx]])
        return topology_buffer

    @staticmethod
    def join_ss_bonds(names_1, atom_idx_1, names_2, atom_idx_2):
        topology_buffer = []
        try:
            ss1_idx = names_1.index("SG")
        except ValueError:
            raise  RuntimeError("could not find a sulphur group within selection")
        try:
            ss2_idx = names_2.index("SG")
        except ValueError:
            raise  RuntimeError("could not find a sulphur group within selection")
        topology_buffer.append([atom_idx_1[ss1_idx], atom_idx_2[ss2_idx]])
        return topology_buffer

    def add_ss_bonds(self):
        try:
            ssbonds = self.anotation['ssbonds']
        except KeyError:
            raise IOError("there are not defined ssbonds in the anotation")
        for i, (ss1, ss2) in enumerate(ssbonds):
            ss1_idx = self.residue_idx.index(ss1)
            ss2_idx = self.residue_idx.index(ss2)
            ss1_start, ss1_end = self.residue_extent[ss1_idx]
            ss2_start, ss2_end = self.residue_extent[ss2_idx]
            ssbond_topo = self.join_ss_bonds(
                names_1=self.atom_names[ss1_start:ss1_end],
                names_2=self.atom_names[ss2_start:ss2_end],
                atom_idx_1=self.atom_idx[ss1_start:ss1_end],
                atom_idx_2=self.atom_idx[ss2_start:ss2_end]
            )
            self.bonds = self.bonds + ssbond_topo

    def assign_partial_charges(self, atom_type, resname):
        if atom_type in self.residue_charge_type_source[resname].keys():
            return self.residue_charge_type_source[resname][atom_type]
        else:
            return 0.0

    def make_sequences(self, residue_names, residue_chains):

        chain_sequences = dict()
        for chain in residue_chains:
            if chain not in chain_sequences.keys():
                chain_sequences[chain] = ""


        for chain in chain_sequences.keys():
            chain_sequences[chain] = ""
            residues_in_chain = [i for i, j in zip(residue_names, residue_chains) if j == chain]
            for res in residues_in_chain:
                chain_sequences[chain] = chain_sequences[chain] + self.residue_charge_source[res]['symbol']

        return chain_sequences

    def read_anotation(self, file=None, dictionary=None):
        """
        PDB read anotation
        ---
        tries to solve the lack of annotations in PQR file by reading a custom dictionary in json format
        :param file:
        :param dictionary:
        :return:
        """
        import json
        if file is None and dictionary is None:
            raise IOError("lacking input")
        if file is None:
            anotation = dictionary
        else:
            with open(file) as f:
                anotation = json.load(f)

        if 'ssbonds' in anotation.keys():

            self.anotation['ssbonds'] = []
            for ss1, ss2 in anotation['ssbonds']:
                self.anotation['ssbonds'].append([ss1, ss2])
