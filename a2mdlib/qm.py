import numpy as np
import re
from a2mdlib.volumes import Volume

class A2MDlibQM:
    def __init__(self, name, verbose):
        self.__name = name
        if verbose in [True, False]:
            self.__verbose = verbose
        else:
            self.__verbose = True

    def log(self, mssg):
        if self.__verbose:
            print("[{:s}] {:s}".format(self.__name, mssg))


symetry_index = np.array(
    [
        [0, 0, 0],  # s
        [1, 0, 0],  # x
        [0, 1, 0],  # y
        [0, 0, 1],  # z
        [2, 0, 0],  # xx
        [0, 2, 0],  # yy
        [0, 0, 2],  # zz
        [1, 1, 0],  # xy
        [1, 0, 1],  # xz
        [0, 1, 1],  # yz
        [3, 0, 0],  # xxx
        [0, 3, 0],  # yyy
        [0, 0, 3],  # zzz
        [2, 1, 0],  # xxy
        [2, 0, 1],  # xxz
        [0, 2, 1],  # yyz
        [1, 2, 0],  # xyy
        [1, 0, 2],  # xzz
        [0, 1, 2],  # yzz
        [1, 1, 1]   # xyz
    ]
)

atom_names = list(
    ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']
)


def parse_fortran_scientific(fortran_number):
    """
    This functions turns a fortran scienfitic notation into a C scientific notation.
    For instance:
    8.45612D-05 -> 8.45612E-05
    so it can be converted to python float

    :param fortran_number: string containing a number in fortran scientific notations
    :return: floating number
    """
    c_number = fortran_number.replace('D', 'E')
    return float(c_number)


def set_nearest_atom(points, coordinates):
    """
    This function takes a set of sampling points and the coordinates of an
    object, creates a distance matrix, and then it takes the shortest distance
    for each point from each atom

    :param points: set of sampling points
    :type points: np.ndarray
    :param coordinates: coordinates of the nuclei of a molecule
    :type coordinates: np.ndarray
    :return: vector of shortest distances
    :rtype: np.ndarray
    """
    n_points = points.shape[0]
    n_coords = coordinates.shape[0]
    d = np.zeros((n_points, n_coords))
    for i in range(n_coords):
        d[:, i] = np.sqrt(
            ((points[:, 0] - coordinates[i, 0])**2) +
            ((points[:, 1] - coordinates[i, 1])**2) +
            ((points[:, 2] - coordinates[i, 2])**2)
        )
    return d.min(axis=1)


class Wavefunction(A2MDlibQM):
    def __init__(self, verbose=True, file=None, batch_size=1000, prefetch_dm=True):
        A2MDlibQM.__init__(self, verbose=verbose, name='wavefunction handler')
        self.__file = file
        self.__coeff = None
        self.__occ = None
        self.__exp = None
        self.__sym = None
        self.__cent = None
        self.__coord = None
        self.__types = None
        self.__charges = None
        self.__n_prim = None
        self.__n_orb = None
        self.__n_nuc = None
        self.__D = None
        self.__batch_size = batch_size
        if file is not None:
            self.read()
        if prefetch_dm:
            self.calculate_density_matrix()

    @staticmethod
    def __gaussian(r, s, exp):
        """

        Calculates a gaussian function (e^(-(x^2))) with a given
        harmonic, expressed in the function s. For instance, if
        s = [0,0,2], then the function calculates:

        f(x) = e^(-x^"2) * (z^2)

        :param r: radial distance to center
        :type r: np.ndarray
        :param s: symmetry type
        :type s: np.ndarray
        :param exp: exponent of the gassian
        :type exp: np.ndarray
        :return: gaussian function value
        :rtype: np.ndarray
        """
        r2 = np.sum(r * r, axis=1)
        smx = r[:, 0] ** s[0]
        smx *= r[:, 1] ** s[1]
        smx *= r[:, 2] ** s[2]
        return smx*np.exp(-r2*exp)

    @staticmethod
    def __parse_centers(wfn, primitives):
        """

        In charge of parsing the distribution of centers for the basis set

        :param wfn:
        :param primitives:
        :return:
        """
        centers = np.zeros(primitives, dtype='int64')
        line = wfn[0]
        i = 0
        j = 0
        while line[:6] == "CENTRE":
            sp = line.split()
            sp = sp[2:]
            for centre in sp:
                centers[i] = int(centre) - 1
                i += 1
            j += 1
            line = wfn[j]
        return centers, wfn[j:]

    @staticmethod
    def __parse_coordinates(wfn_instance, nuclei_number):
        """

        In charge of parsing the geometry of the molecule

        :param wfn_instance:
        :param nuclei_number:
        :return:
        """
        coordinates_lines = wfn_instance[:nuclei_number]
        atom_types = [None] * nuclei_number
        coordinates = np.zeros((nuclei_number, 3))
        charge = np.zeros(nuclei_number)
        i = 0
        for line in coordinates_lines:
            splitted_line = line.split()
            atom_types[i] = splitted_line[0]
            coordinates[i, 0] = float(splitted_line[4])
            coordinates[i, 1] = float(splitted_line[5])
            coordinates[i, 2] = float(splitted_line[6])
            charge[i] = float(splitted_line[-1])
            i += 1
        return coordinates, atom_types, charge, wfn_instance[nuclei_number:]

    @staticmethod
    def __parse_exponents(wfn, primitives):
        """

        In charge of parsing the exponents of the basisi set

        :param wfn:
        :param primitives:
        :return:
        """
        exponents = np.zeros(primitives, dtype='float64')
        line = wfn[0]
        i = 0
        j = 0
        while line[:9] == "EXPONENTS":
            sp = line.split()
            sp = sp[1:]
            for exp in sp:
                exponents[i] = parse_fortran_scientific(exp)
                i += 1
            j += 1
            line = wfn[j]
        return exponents, wfn[j:]

    @staticmethod
    def __parse_head(wfn_instance):
        """

        In charge of parsing the head

        :param wfn_instance:
        :return:
        """
        head_line = wfn_instance[0]
        head_terms = head_line.split()
        head = dict(
            molecular_orbitals=int(head_terms[1]),
            primitives=int(head_terms[4]),
            nuclei=int(head_terms[6])
        )
        return head, wfn_instance[1:]

    @staticmethod
    def __parse_orbitals(wfn, primitives, orbitals):
        """

        In charge of parsing the coefficients of the wave function for each molecular orbital

        :param wfn:
        :param primitives:
        :param orbitals:
        :return:
        """
        orbital_text = ''.join(wfn)
        orbital_text = re.split('END DATA', orbital_text)[0]
        orbital_text = re.split('MO\s*\d{1,3}\s*MO\s\d\.\d\s*OCC\sNO\s=\s*', orbital_text)
        orbital_text = orbital_text[1:]
        coeff = np.zeros((orbitals, primitives))
        occ = np.zeros(orbitals)
        i = 0
        j = 0
        for orb in orbital_text:
            orb_lines = orb.split('\n')
            header = orb_lines[0].split()
            occ[i] = float(header[0])
            orb_lines = orb_lines[1:]
            for line in orb_lines:
                sp = line.split()
                for item in sp:
                    coeff[i, j] = parse_fortran_scientific(item)
                    j += 1
            j = 0
            i += 1
        return coeff, occ

    @staticmethod
    def __parse_symmetry(wfn, primitives):
        """

        In charge of parsing the symmetry coefficients of the basis set

        :param wfn:
        :param primitives:
        :return:
        """
        symetry = np.zeros(primitives, dtype='int64')
        line = wfn[0]
        i = 0
        j = 0
        while line[:4] == "TYPE":
            sp = line.split()
            sp = sp[2:]
            for sym in sp:
                symetry[i] = int(sym) - 1
                i += 1
            j += 1
            line = wfn[j]
        return symetry, wfn[j:]

    # PUBLIC METHODS

    def calculate_density(self, coordinates):
        """

        Performs a calculation of electron density for the given coordinates by applying:

        rho = \sum_i \sum_j D_ij Xi_i Xi_j

        where D is the value of the density matrix (obtained by multiplying the coefficents i and j
        and the occupation numbers), and Xi_i and Xi_j are gaussian functions

        :param coordinates: coordinates upon which to provide the electron density value
        :type coordinates: np.ndarray
        :return: electron density
        :rtype: np.ndarray
        """
        if self.__D is None:
            raise IOError("density matrix has not been set")
        else:
            rho = np.zeros(coordinates.shape[0])
            batch_size = self.__batch_size
            n_batches = int(coordinates.shape[0]/batch_size)
            b = 0
            for b in range(n_batches):
                x = np.zeros((self.__n_prim, batch_size))
                for p in range(self.__n_prim):
                    cntr = self.__coord[self.__cent[p], :]
                    d = coordinates[b*batch_size:(b+1)*batch_size, :] - cntr
                    x[p, :] = self.__gaussian(
                        d,
                        symetry_index[self.__sym[p], :],
                        self.__exp[p]
                    )
                for p in range(self.__n_prim):
                    for q in range(self.__n_prim):
                        rho[b*batch_size:(b+1)*batch_size] += self.__D[p, q] * x[p, :] * x[q, :]

            if coordinates.shape[0] % batch_size != 0:
                x = np.zeros((self.__n_prim, coordinates.shape[0] % batch_size))
                for p in range(self.__n_prim):
                    cntr = self.__coord[self.__cent[p], :]
                    d = coordinates[b*batch_size:, :] - cntr
                    x[p, :] = self.__gaussian(
                        d,
                        symetry_index[self.__sym[p], :],
                        self.__exp[p]
                    )
                for p in range(self.__n_prim):
                    for q in range(self.__n_prim):
                        rho[b*batch_size:] += self.__D[p, q] * x[p, :] * x[q, :]
            return rho

    def calculate_density_by_atom(self, coordinates, center_i):
        """

        Performs a calculation of the density by only multiplying basis functions belonging
        only to the specific center

        :param coordinates:
        :type coordinates: np.ndarray
        :param center_i:
        :type center_i: int
        :return:
        """
        selected_centers = [i for i in range(len(self.__cent)) if self.__cent[i] == center_i]
        cntr = self.__coord[self.__cent[selected_centers[0]], :]
        n_sel = len(selected_centers)
        d = coordinates - cntr
        rho = np.zeros(coordinates.shape[0])

        x = np.zeros((len(selected_centers), coordinates.shape[0]))
        for i, p in zip(range(n_sel), selected_centers):
            trm = self.__sym[p]
            x[i, :] = self.__gaussian(d, symetry_index[trm, :], self.__exp[p])
        for i, p in zip(range(n_sel), selected_centers):
            for j, q in zip(range(n_sel), selected_centers):
                rho += self.__D[p, q]*x[i, :]*x[j, :]
        return rho

    def calculate_density_by_pairs(self, coordinates, center_i, center_j):
        """

        Performs a calculation of the density by only multiplying basis functions belonging
        to the centers i and j.

        :param coordinates:
        :type coordinates: np.ndarray
        :param center_i:
        :type center_i: int
        :param center_j:
        :type center_j: int
        :return:
        """
        selected_centers_i = [i for i in range(len(self.__cent)) if self.__cent[i] == center_i]
        selected_centers_j = [i for i in range(len(self.__cent)) if self.__cent[i] == center_j]
        coordinates_ori_i = self.__coord[center_i]
        coordinates_ori_j = self.__coord[center_j]

        di = coordinates - coordinates_ori_i
        dj = coordinates - coordinates_ori_j
        rho = np.zeros(coordinates.shape[0])
        x = np.zeros((len(selected_centers_i), coordinates.shape[0]))
        y = np.zeros((len(selected_centers_j), coordinates.shape[0]))

        for i, p in enumerate(selected_centers_i):
            trm = self.__sym[p]
            x[i, :] = self.__gaussian(di, symetry_index[trm, :], self.__exp[p])
        for j, q in enumerate(selected_centers_j):
            trm = self.__sym[q]
            y[j, :] = self.__gaussian(dj, symetry_index[trm, :], self.__exp[q])
        for i, p in enumerate(selected_centers_i):
            for j, q in enumerate(selected_centers_j):
                rho += self.__D[p, q]*x[i, :]*y[j, :]
        return rho

    def calculate_density_matrix(self):
        """

        Calculates the density matrix

        :return:
        """
        dm = np.zeros((self.__n_prim, self.__n_prim))
        for i in range(self.__n_orb):
            for p in range(self.__n_prim):
                dm[p, :] += self.__occ[i]*self.__coeff[i, p]*self.__coeff[i, :]
        self.__D = dm

    def get_atom_labels(self):
        """

        returns a list of the labels of the atoms ['C','H','O','C',...]

        :return: list of atomic labels
        :rtype: list
        """
        atom_labels = [atom_names[int(i) - 1] for i in self.__charges]
        return atom_labels

    def get_atomic_numbers(self):
        """

        returns a list of atomic numbers in format np.ndarry

        :return: atomic numbers
        :rtype: np.ndarray
        """
        return np.array([int(i) for i in self.__charges])

    def get_coordinates(self):
        """

        returns coordinates in Atomic Units

        1AU = 0.5292 Angstroms

        :return:
        """
        self.log("default units are AU")
        return self.__coord.copy()

    def get_number_molecular_orbitals(self):
        """

        returns the number of molecular orbitals

        :return: number of molecular orbitals
        """
        return self.__n_orb

    def get_number_primitives(self):
        """

        returns the number of primitive functions

        :return: number of primitive functions
        """
        return self.__n_prim

    def read(self):
        """

        reads the wave function

        """
        if self.__file is None:
            raise IOError("file is not defined")
        else:
            with open(self.__file) as fh:
                wfn = fh.readlines()
            wfn = wfn[1:]
            head, wfn = self.__parse_head(wfn)
            coords, types, charges, wfn = self.__parse_coordinates(wfn, head['nuclei'])
            centers, wfn = self.__parse_centers(wfn, head['primitives'])
            syms, wfn = self.__parse_symmetry(wfn, head['primitives'])
            exponents, wfn = self.__parse_exponents(wfn, head['primitives'])
            coeff, occ = self.__parse_orbitals(wfn, head['primitives'], head['molecular_orbitals'])

            self.__coeff = coeff
            self.__occ = occ
            self.__exp = exponents
            self.__sym = syms
            self.__cent = centers
            self.__coord = coords
            self.__types = types
            self.__charges = charges
            self.__n_prim = head['primitives']
            self.__n_orb = head['molecular_orbitals']
            self.__n_nuc = head['nuclei']

    def set_file(self, file):
        """

        Allows to change the file

        :param file:
        :return:
        """
        self.__file = file


class ElectronDensity(Volume):
    def __init__(self, verbose=False, file=None):
        Volume.__init__(self, verbose=verbose, filename=file)

    def sample(self, number_points):
        # SAMPLING THROUGH PERMUTATION TO AVOID REPETITIONS
        size = self.shape[0] * self.shape[1] * self.shape[2]
        if number_points > size:
            raise ArithmeticError("can't sample more data than there is")
        sampled = np.arange(size)
        sampled = np.random.permutation(sampled)[:number_points]
        tnsor = self.get_volume()
        ##
        cx, cy, cz = np.unravel_index(sampled, self.shape)
        itr = np.arange(number_points)
        sampled_density = np.zeros(number_points)
        for i, x, y, z in zip(itr, cx, cy, cz):
            sampled_density[i] = tnsor[x, y, z]
        r0 = self.get_r0()
        basis = self.get_basis()
        x = np.array([cx, cy, cz], dtype='float64').T
        x += r0
        x = x.dot(basis)
        return x, sampled_density


class GaussianLog(A2MDlibQM):

    def __init__(self, file, verbose=True):
        """

        :param file:
        :param verbose:
        """
        A2MDlibQM.__init__(self, verbose=verbose, name='gaussianLog')
        self.__basis = None
        self.__coordinates = None
        self.__atom_numbers = None
        self.__dipolar = None
        self.__E = None
        self.__E2 = None
        self.__energy_terms = None
        self.__fname = file
        self.__NPA = None
        self.__NPT = None
        self.read()

    def calculate_dipole(self):
        """

        Calculates dipole from NPA charges

        :return:
        """
        return np.sum(self.__coordinates * self.__NPA.reshape(
            self.__coordinates.shape[0], 1
        ), axis=0)

    def read(self):
        """

        :return:
        """
        coordinates_flag = False
        coordinates_barrier = 0
        charges_flag = False
        charges_lock = True

        coordinates = []
        atom_numbers = []

        natural_charges = []
        natural_total_charges = []
        with open(self.__fname) as f:
            self.log('file opened')
            for line in f:
                line = line.strip()
                if re.match('Standard\sorientation', line):
                    coordinates = []
                    atom_numbers = []
                    coordinates_flag = True
                    self.log('Reading coordinates')
                    continue
                if coordinates_flag:
                    if re.match('-{5}', line):
                        if coordinates_barrier == 2:
                            coordinates_barrier = False
                            coordinates_flag = False

                            self.__coordinates = np.array(coordinates)
                            self.__atom_numbers = np.array(atom_numbers)
                            self.log('Reading coordinates - FINISHED - %d ATOMS' % len(atom_numbers))
                        else:
                            coordinates_barrier += 1
                    elif re.match('Center\s*Atomic\s*Atomic\s*Coordinates\s\(Angstroms\)', line):
                        pass
                    elif re.match('Number\s*Number\s*Type\s*X\s*Y\s*Z', line):
                        pass
                    else:
                        cnter, atmnmbr, atmtype, x, y, z = line.split()
                        atmnmbr = int(atmnmbr)
                        coordinates.append([
                            float(x),
                            float(y),
                            float(z)
                        ])
                        atom_numbers.append(atmnmbr)
                if re.match('Standard\sbasis:\s(.*)', line):
                    m = re.match('Standard\sbasis:\s(.*)', line)
                    self.__basis = m.group(1)
                    self.log('Basis : %s' % m.group(1))
                    continue
                if re.match('SCF\sDone:\s*E\(RHF\)\s=(.*)\s*A\.U\.\safter\s*\d*\scycles', line):
                    m = re.match('SCF\sDone:\s*E\(RHF\)\s=(.*)\s*A\.U\.\safter\s*\d*\scycles', line)
                    self.__E = float(m.group(1).replace('D', 'E'))
                    self.log('Hartree-Fock Energy = %8.4f' % float(m.group(1).replace('D', 'E')))
                    continue
                if re.match('E2\s=\s*(.*)\sEUMP2\s=\s*(.*)', line):
                    m = re.match('E2\s=\s*(.*)\sEUMP2\s=\s*(.*)', line)
                    self.__E2 = float(m.group(2).replace('D', 'E'))
                    self.log('MP2 = %8.4f' % float(m.group(2).replace('D', 'E')))
                    continue
                if re.match('Natural Population', line) and charges_lock:
                    charges_flag = True
                    self.log('Reading NPA')
                    continue
                if charges_flag and charges_lock:
                    if re.match('\s*Atom\s*No\s*Charge\s*Core\s*Valence\s*Rydberg\s*Total', line):
                        pass
                    elif re.match('-{5}', line):
                        pass
                    elif re.match('Natural -*', line):
                        pass
                    elif re.match('={5}', line):
                        charges_flag = False
                        charges_lock = False
                        self.__NPA = np.array(natural_charges)
                        self.__NPT = np.array(natural_total_charges)
                        self.log('Reading NPA  - FINISHED')
                    else:

                        atmsymbl, no, charge, core, valence, rydberg, total = line.split()
                        try:
                            natural_charges.append(float(charge))
                        except ValueError:
                            pass
                        natural_total_charges.append(float(total))
                    continue
                if re.match(r'^\s*X=\s*(-?\d*.\d*)\s*Y=\s*(-?\d*.\d*)\s*Z=\s*(-?\d*.\d*)', line):
                    m = re.match(r'^\s*X=\s*(-?\d*.\d*)\s*Y=\s*(-?\d*.\d*)\s*Z=\s*(-?\d*.\d*)', line)
                    dx = float(m.group(1))
                    dy = float(m.group(2))
                    dz = float(m.group(3))
                    dipole = np.array([dx, dy, dz])
                    self.__dipolar = dipole
                    self.log('Dipole : %4.3f %4.3f %4.3f' % (dx, dy, dz))
                    continue
                if re.match(r'N-N=\s*(-?\d*.\d*D[+,-]\d{2})\s*E-N=\s*(-?\d*.\d*D[+,-]\d{2})\s*KE=\s*(-?\d*.\d*D[+,-]\d{2})\s*', line):
                    m = re.match(r'N-N=\s*(-?\d*.\d*D[+,-]\d{2})\s*E-N=\s*(-?\d*.\d*D[+,-]\d{2})\s*KE=\s*(-?\d*.\d*D[+,-]\d{2})\s*', line)

                    nn = float(m.group(1).replace('D', 'E'))
                    ne = float(m.group(2).replace('D', 'E'))
                    kin = float(m.group(3).replace('D', 'E'))

                    self.__energy_terms = dict(
                        nuclei_nuclei_potential=nn,
                        nuclei_electron_potential=ne,
                        kinetic=kin
                    )

                    self.log('Energy terms found :')
                    self.log("\t{:30s}: {:12.4f} Ha".format(
                            "Nuclei nuclei repulsion", self.__energy_terms['nuclei_nuclei_potential']
                    ))
                    self.log("\t{:30s}: {:12.4f} Ha".format(
                        "Nuclei electron attraction", self.__energy_terms['nuclei_electron_potential']
                    ))
                    self.log("\t{:30s}: {:12.4f} Ha".format(
                        "Kinetic energy", self.__energy_terms['kinetic']
                    ))
                    continue
        return True

    def seek_ep(self):
        """

        Iterates over the file until it founds the electrostatic potential block. It does not read
        the lines belonging to atoms. Output units are generally hartrees.

        :return: electostatic potential
        :rtype: np.ndarray
        """
        flag = False
        ep = []
        with open(self.__fname) as f:
            i = 0
            contents = f.readlines()
            for i, line in enumerate(contents):
                line = line.strip()
                if re.match('Electrostatic Properties \(Atomic Units\)', line):
                    self.log("found ep block")
                    flag = True
                    break
            if flag:
                for line in contents[i + 6:]:
                    line = line.strip()
                    if re.search(r'Atom', line):
                        continue
                    elif re.search(r'------', line):
                        break
                    else:
                        ep.append(float(line.split()[1]))

        return np.array(ep, dtype='float64')

    def get_atomic_numbers(self):
        """

        :return: np.ndarray
        """
        return self.__atom_numbers.copy()

    def get_energy_terms(self):
        """

        :return:
        """
        return self.__energy_terms

    def get_basis(self):
        """

        :return: basis employed in the QM calculation
        """
        return self.__basis

    def get_charges(self, kind='partial'):
        """

        :param kind: either partial or total
        :return:
        """
        if kind == 'partial':
            return self.__NPA.copy()
        elif kind == 'total':
            return self.__NPT.copy()
        else:
            raise NotImplementedError("unknown charge type {:s}".format(kind))

    def get_coordinates(self):
        """

        :return:
        """
        return self.__coordinates.copy()

    def get_dipole(self):
        """

        :return:
        """
        return self.__dipolar.copy()

    def get_energy(self, calctype='MP2'):
        """

        :param calctype: type of energy (HF or post-HF)
        :return:
        """
        if calctype == 'MP2':
            return self.__E2
        elif calctype == 'RHF':
            return self.__E

    def get_total_charges(self):
        """

        :return:
        """
        return self.__NPT

class CubeFile(A2MDlibQM):
    def __init__(self, file, verbose=False):
        """

        CubeFile allows to read the cube representation of electron density and store
        it in a numpy tensor.

        :param file:
        :param verbose:
        """
        A2MDlibQM.__init__(self, name='CubeFile', verbose=verbose)
        self.file = file
        self.cube_tensor = None
        self.origin = None
        self.basis = None
        self.read()

    def read(self):
        """
        Reads a cube file
        :return:
        """

        x_ori = None
        y_ori = None
        z_ori = None
        n_atoms = 1000
        x_n_values = None
        y_n_values = None
        z_n_values = None

        unit_cell = np.zeros((3, 3), dtype='float64')
        row_counter = 0
        row_buffer = None

        cube_tensor = None
        i = 0
        j = 0

        with open(self.file) as f:
            line = f.readline()
            cnt = 1
            while line:
                line = f.readline()
                cnt += 1

                if cnt == 3:

                    n_atoms, x_ori, y_ori, z_ori = line.split()
                    n_atoms = int(n_atoms)
                    x_ori = float(x_ori)
                    y_ori = float(y_ori)
                    z_ori = float(z_ori)

                    self.log(
                        "origin at : {:12.6f} {:12.6f} {:12.6f} (Angstroms)".format(
                            x_ori, y_ori, z_ori
                        )
                    )

                if cnt == 4:
                    x_n_values, xi, yi, zi = line.split()
                    x_n_values = int(x_n_values)
                    unit_cell[0, :] = float(xi), float(yi), float(zi)
                if cnt == 5:
                    y_n_values, xj, yj, zj = line.split()
                    y_n_values = int(y_n_values)
                    unit_cell[1, :] = float(xj), float(yj), float(zj)
                if cnt == 6:
                    z_n_values, xk, yk, zk = line.split()
                    z_n_values = int(z_n_values)
                    unit_cell[2, :] = float(xk), float(yk), float(zk)
                    n_values = z_n_values * x_n_values * y_n_values
                    row_buffer = np.zeros(z_n_values, dtype='float64')
                    cube_tensor = np.zeros((x_n_values, y_n_values, z_n_values), dtype='float32')
                    self.log("total points : {:12d}".format(n_values))
                if cnt > (6 + n_atoms):
                    row_length = len(line.split())
                    row_buffer[row_counter:row_counter + row_length] = list(map(float, line.split()))
                    row_counter += row_length
                    if row_counter == z_n_values:
                        cube_tensor[i, j, :] = row_buffer
                        row_counter -= row_counter
                        row_buffer -= row_buffer
                        j += 1
                        if j == y_n_values:
                            i += 1
                            j -= j

        self.cube_tensor = cube_tensor
        self.origin = np.array((x_ori, y_ori, z_ori), dtype='float64')
        self.basis = unit_cell

    def get_plane(self, value, axis = 'z'):
        """

        gives back the nearest layer of the tensor to the given value. Useful
        for representation of contours.

        :param value:
        :param axis:
        :return:
        """

        if axis not in ['x', 'y', 'z']:
            raise IOError("not found axis {:s}".format(axis))

        dim_idx = ['x', 'y', 'z'].index(axis)

        tensor_idx = int((value - self.origin[dim_idx]) / self.basis[dim_idx, dim_idx])

        if dim_idx == 0:
            return self.cube_tensor[tensor_idx, :, :]
        elif dim_idx == 1:
            return self.cube_tensor[:, tensor_idx, :]
        elif dim_idx == 2:
            return self.cube_tensor[:, :, tensor_idx]




if __name__ == '__main__':
    pass
