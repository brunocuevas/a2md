import numpy as np
import re
from a2mdio.volumes import Volume

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


class WaveFunction(A2MDlibQM):
    def __init__(self, file=None, verbose=True, batch_size=1000, prefetch_dm=True):
        A2MDlibQM.__init__(self, verbose=verbose, name='wavefunction handler')

        self.__file = file
        coeff, occ, exponents, syms, centers, coords, types, charges, primitives, molecular_orbitals, nuclei = self.read()
        self.coeff = coeff
        self.occ = occ
        self.exp = exponents
        self.sym = syms
        self.cent = centers
        self.coords = coords
        self.types = types
        self.charges = charges
        self.nprims = primitives
        self.norbs = molecular_orbitals
        self.ncenters = nuclei
        self.density_matrix = None
        self.batchsize = batch_size
        if prefetch_dm:
            self.density_matrix= self.calculate_density_matrix()

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

    def eval(self, coordinates):
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
        if self.density_matrix is None:
            raise IOError("density matrix has not been set")

        rho = np.zeros(coordinates.shape[0])
        batch_size = self.batchsize
        n_batches = int(coordinates.shape[0]/batch_size)
        b = 0
        for b in range(n_batches):
            x = np.zeros((self.nprims, batch_size))
            for p in range(self.nprims):
                cntr = self.coords[self.cent[p], :]
                d = coordinates[b*batch_size:(b+1)*batch_size, :] - cntr
                x[p, :] = self.__gaussian(
                    d,
                    symetry_index[self.sym[p], :],
                    self.exp[p]
                )
            for p in range(self.nprims):
                for q in range(self.nprims):
                    rho[b*batch_size:(b+1)*batch_size] += self.density_matrix[p, q] * x[p, :] * x[q, :]

        if coordinates.shape[0] % batch_size != 0:
            x = np.zeros((self.nprims, coordinates.shape[0] % batch_size))
            for p in range(self.nprims):
                cntr = self.coords[self.cent[p], :]
                d = coordinates[b*batch_size:, :] - cntr
                x[p, :] = self.__gaussian(
                    d,
                    symetry_index[self.sym[p], :],
                    self.exp[p]
                )
            for p in range(self.nprims):
                for q in range(self.nprims):
                    rho[b*batch_size:] += self.density_matrix[p, q] * x[p, :] * x[q, :]
        return rho


    def calculate_density_matrix(self):
        """

        Calculates the density matrix

        :return:
        """
        dm = np.zeros((self.nprims, self.nprims))
        for i in range(self.norbs):
            for p in range(self.nprims):
                dm[p, :] += self.occ[i] * self.coeff[i, p] * self.coeff[i, :]
        return dm

    def get_atom_labels(self):
        """

        returns a list of the labels of the atoms ['C','H','O','C',...]

        :return: list of atomic labels
        :rtype: list
        """
        atom_labels = [atom_names[int(i) - 1] for i in self.charges]
        return atom_labels

    def get_atomic_numbers(self):
        """

        returns a list of atomic numbers in format np.ndarry

        :return: atomic numbers
        :rtype: np.ndarray
        """
        return np.array([int(i) for i in self.charges])

    def get_coordinates(self):
        """

        returns coordinates in Atomic Units

        1AU = 0.5292 Angstroms

        :return:
        """
        self.log("default units are AU")
        return self.coords.copy()

    def get_number_molecular_orbitals(self):
        """

        returns the number of molecular orbitals

        :return: number of molecular orbitals
        """
        return self.norbs

    def get_number_primitives(self):
        """

        returns the number of primitive functions

        :return: number of primitive functions
        """
        return self.nprims

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

            return [
                coeff, occ, exponents, syms, centers, coords, types, charges, head['primitives'],
                head['molecular_orbitals'], head['nuclei']
            ]


class WaveFunctionGPU(WaveFunction):
    def __init__(self, file, dtype, device='cuda:0'):
        import torch
        self.torch = torch
        self.device = torch.device(device)
        self.dtype = dtype
        WaveFunction.__init__(self, file=file, verbose=True, batch_size=100000, prefetch_dm=True)
        self.density_matrix = self.calculate_density_matrix()
        self.coords = torch.tensor(self.coords, dtype=dtype, device=self.device)


    def convert_symmetry_to_tensor(self, sym, dims):

        sym_tensor = self.torch.tensor(symetry_index[sym], dtype=self.dtype, device=self.device)
        sym_tensor = sym_tensor.repeat(dims).reshape(dims, 3)
        return sym_tensor

    def basis_function(self, x, i):

        center = self.cent[i]
        coords = self.coords[center, :]
        sym_vector = symetry_index[self.sym[i], :]
        exp = self.exp[i]
        rv = (x - coords)
        r = rv.pow(2.0).sum(1)
        g = self.torch.exp(-exp * r)
        l = rv[:,0].pow(sym_vector[0]) * rv[:,1].pow(sym_vector[1]) * rv[:,2].pow(sym_vector[2])
        return g * l

    def distance_vector(self, x, i):
        dims = x.size(0)
        coords = self.coords[i, :].repeat(dims).reshape(dims, 3)
        rv = (x - coords)
        r = rv.pow(2.0).sum(1)
        return rv, r


    def eval(self, x):
        p = self.torch.zeros(x.size(0), device=self.device, dtype=self.dtype)
        buffer = self.torch.zeros(self.nprims, x.size(0), dtype=self.dtype, device=self.device)

        for i in range(self.nprims):
            buffer[i, :] = self.basis_function(x, i)

        for i in range(self.nprims):

            for j in range(self.nprims):

                p += self.density_matrix[i, j] * buffer[i, :] * buffer[j, :]
        return p


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

    def __init__(self, file, method, charges, ep=False, verbose=True):
        """

        :param file:
        :param verbose:
        """
        from a2mdio.parsers import std_coordinates
        from a2mdio.parsers import npa_charges, mk_charges
        from a2mdio.parsers import energy_decomposition, dipole
        from a2mdio.parsers import hf_energy, mp2_energy, dft_energy
        from a2mdio.parsers import electrostatic_potential
        A2MDlibQM.__init__(self, verbose=verbose, name='gaussianLog')
        self.fname = file
        tokkens = [
            std_coordinates, dipole, energy_decomposition
        ]
        tokkens_labels = ['coordinates', 'dipole', 'energy_decomposition']
        if method == 'MP2':
            tokkens.append(mp2_energy)
            tokkens_labels.append('energy')
        elif method == 'HF' :
            tokkens.append(hf_energy)
            tokkens_labels.append('energy')
        elif method[:3] == 'dft' :
            tokkens.append(lambda x : dft_energy(x, functional=method[4:]))
            tokkens_labels.append('energy')

        if charges == 'NPA':
            tokkens.append(npa_charges)
            tokkens_labels.append('charges')
        elif charges == 'MK' :
            tokkens.append(mk_charges)
            tokkens_labels.append('charges')

        if ep:
            tokkens.append(electrostatic_potential)
            tokkens_labels.append('ep')

        self.tokkens = tokkens
        self.tokken_labels = tokkens_labels

    def read(self):
        """

        :return:
        """
        output_dict = dict()
        for label, tk in zip(self.tokken_labels, self.tokkens):
            output_dict[label] = self.seek(tk)
        return output_dict

    def seek(self, fun):
        """

        :param fun:
        :return:
        """
        with open(self.fname) as f:
            return fun(f.readlines())



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
