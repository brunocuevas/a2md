import numpy as np
import re
from a2md.baseclass import A2MD_basis

symetry_index = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 2],
        [1, 1, 0],
        [1, 0, 1],
        [0, 1, 1],
    ]
)

atom_names = list(
    ['H', 'He', 'Li', 'Be', 'B' ,'C','N','O','F','Ne', 'Na','Mg','Al', 'Si','P','S', 'Cl','Ar']
)


def element2an(an):
    return atom_names.index(an) + 1

def an2element(elem):
    return atom_names[elem-1]

def parseFortranSciNot(fortran_number) :
    fortran_number = fortran_number.replace('D', 'E')
    return float(fortran_number)

def setNearestAtom(points, coordinates):
    n_points = points.shape[0]
    n_coords = coordinates.shape[0]
    D = np.zeros((n_points, n_coords))
    for i in range(n_coords):
        D[:,i] = np.sqrt(
            ((points[:,0] - coordinates[i,0])**2) +
            ((points[:,1] - coordinates[i,1])**2) +
            ((points[:,2] - coordinates[i,2])**2)
        )
    return D.min(axis=1)

def convert_connectivity_tree_to_pairs(connectivity_tree):
    pairs = []
    for i, item1 in enumerate(connectivity_tree):
        for j, item2 in enumerate(item1):
            if sorted([i, item2]) not in pairs:
                pairs.append(sorted([i, item2]))
    return pairs

def create_all2all_topology(n_items):
    topo_array = []
    for i in range(n_items):
        topo_array.append([])
    for i in range(n_items):
        for j in range(n_items):
            if i != j :
                topo_array[i].append(j)
    return topo_array

class Integrator(A2MD_basis):
    def __init__(self, verbose=False):
        """

        :param verbose:
        """
        A2MD_basis.__init__(self, verbose=verbose, name='Integrator')

    def integrate(self, params):
        from a2md.amd import SUPPORT_TYPE
        integrals = []
        for i, fun in enumerate(params):
            sup_type = fun['support_type']
            if sup_type == 'ORC':
                continue

            if sup_type in SUPPORT_TYPE.keys():
                input_dict = fun['params']
                input_dict['coordinates'] = [0, 0, 0]
                f = SUPPORT_TYPE[sup_type](input_dict)
                integral = f.integral()
                self.log("function {:d} : {:12.6e}".format(i, integral))
                integrals.append(f.integral())
        return integrals

    @staticmethod
    def integrate_fun(fun):
        from a2md.amd import SUPPORT_TYPE
        input_dict = fun['params']
        input_dict['coordinates'] = [0, 0, 0]
        f = SUPPORT_TYPE[fun['support_type']](input_dict)
        return f.integral()

class Wavefunction(A2MD_basis) :
    def __init__(self, verbose = True, file = None, batch_size = 1000):
        A2MD_basis.__init__(self, verbose= verbose, name='wavefunction handler')
        self.__file    = file
        self.__coeff   = None
        self.__occ     = None
        self.__exp     = None
        self.__sym     = None
        self.__cent    = None
        self.__coord   = None
        self.__types   = None
        self.__charges = None
        self.__n_prim  = None
        self.__n_orb   = None
        self.__n_nuc   = None
        self.__D       = None
        self.__batch_size = batch_size
    ## STATIC METHODS
    @staticmethod
    def __gaussian(r,s,exp):
        r2 = np.sum(r * r, axis = 1)
        smx  = r[:, 0] ** s[0]
        smx *= r[:, 1] ** s[1]
        smx *= r[:, 2] ** s[2]
        return smx*np.exp(-r2*exp)
    @staticmethod
    def __parseCenters(wfn, primitives):
        centers = np.zeros(primitives, dtype='int64')
        line = wfn[0]
        i = 0
        j = 0
        while line[:6] == "CENTRE" :
            sp = line.split()
            sp = sp[2:]
            for centre in sp :
                centers[i] = int(centre) - 1
                i += 1
            j += 1
            line = wfn[j]
        return centers, wfn[j:]
    @staticmethod
    def __parseCoordinates(wfn_instance, nuclei_number):
        coordinates_lines = wfn_instance[:nuclei_number]
        atom_types = [None] * nuclei_number
        coordinates = np.zeros((nuclei_number, 3))
        charge      = np.zeros(nuclei_number)
        i = 0
        for line in coordinates_lines :
            splitted_line = line.split()
            atom_types[i] = splitted_line[0]
            coordinates[i, 0] = float(splitted_line[4])
            coordinates[i, 1] = float(splitted_line[5])
            coordinates[i, 2] = float(splitted_line[6])
            charge[i] = float(splitted_line[-1])
            i += 1
        return coordinates, atom_types, charge, wfn_instance[nuclei_number:]
    @staticmethod
    def __parseExponents(wfn, primitives):
        exponents = np.zeros(primitives, dtype='float64')
        line = wfn[0]
        i = 0
        j = 0
        while line[:9] == "EXPONENTS":
            sp = line.split()
            sp = sp[1:]
            for exp in sp:
                exponents[i] = parseFortranSciNot(exp)
                i += 1
            j += 1
            line = wfn[j]
        return exponents, wfn[j:]
    @staticmethod
    def __parseHead(wfn_instance):
        head_line = wfn_instance[0]
        head_terms = head_line.split()
        head = dict(
            molecular_orbitals = int(head_terms[1]),
            primitives         = int(head_terms[4]),
            nuclei             = int(head_terms[6])
        )
        return head, wfn_instance[1:]
    @staticmethod
    def __parseOrbitals(wfn, primitives, orbitals):
        orbital_text = ''.join(wfn)
        orbital_text = re.split('END DATA', orbital_text)[0]
        orbital_text = re.split('MO\s*\d{1,3}\s*MO\s\d\.\d\s*OCC\sNO\s=\s*', orbital_text)
        orbital_text = orbital_text[1:]
        coeff        = np.zeros((orbitals, primitives))
        occ          = np.zeros(orbitals)
        i = 0 ; j = 0
        for orb in orbital_text :
            orb_lines = orb.split('\n')
            header = orb_lines[0].split()
            occ[i] = float(header[0])
            orb_lines = orb_lines[1:]
            for line in orb_lines :
                sp = line.split()
                for item in sp:
                    coeff[i,j] = parseFortranSciNot(item)
                    j += 1
            j  = 0
            i += 1
        return coeff, occ
    @staticmethod
    def __parseSimetry(wfn, primitives):
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
    ## PUBLIC METHODS
    def setDensityMatrix(self):
        D = np.zeros((self.__n_prim, self.__n_prim))
        for i in range(self.__n_orb) :
            for p in range(self.__n_prim):
                D[p,:] += self.__occ[i]*self.__coeff[i,p]*self.__coeff[i,:]
        self.__D = D
    def calculateDensity(self, coordinates):
        if self.__D is None :
            raise IOError("density matrix has not been set")
        else :
            Q = np.zeros(coordinates.shape[0])
            batch_size = self.__batch_size
            n_batches = int(coordinates.shape[0]/batch_size)
            b = 0

            # BATCHES

            for b in range(n_batches) :
                X = np.zeros((self.__n_prim, batch_size))
                for p in range(self.__n_prim):
                    cntr = self.__coord[self.__cent[p],:]
                    d = coordinates[b*batch_size:(b+1)*batch_size,:] - cntr
                    X[p,:] = self.__gaussian(
                        d,
                        symetry_index[self.__sym[p],:],
                        self.__exp[p]
                    )
                for p in range(self.__n_prim):
                    for q in range(self.__n_prim):
                         Q[b*batch_size:(b+1)*batch_size] += self.__D[p,q] * X[p,:] * X[q,:]
            ## LEFTOVERS

            if coordinates.shape[0] % batch_size != 0:
                X = np.zeros((self.__n_prim, coordinates.shape[0] % batch_size))
                for p in range(self.__n_prim):
                    cntr = self.__coord[self.__cent[p],:]
                    d = coordinates[b*batch_size:,:] - cntr
                    X[p,:] = self.__gaussian(
                        d,
                        symetry_index[self.__sym[p],:],
                        self.__exp[p]
                    )
                for p in range(self.__n_prim):
                    for q in range(self.__n_prim):
                        Q[b*batch_size:] += self.__D[p,q] * X[p,:] * X[q,:]


            return Q

    def calculateDensityByPairs(self, coordinates, center_i, center_j):
        selected_centers_i = [i for i in range(len(self.__cent)) if self.__cent[i] == center_i]
        selected_centers_j = [i for i in range(len(self.__cent)) if self.__cent[i] == center_j]
        coordinates_ori_i  = self.__coord[center_i]
        coordinates_ori_j  = self.__coord[center_j]


        di = coordinates - coordinates_ori_i
        dj = coordinates - coordinates_ori_j
        Q = np.zeros(coordinates.shape[0])
        X = np.zeros((len(selected_centers_i), coordinates.shape[0]))
        Y = np.zeros((len(selected_centers_j), coordinates.shape[0]))

        for i,p in enumerate(selected_centers_i):
            trm    = self.__sym[p]
            X[i,:] = self.__gaussian(di, symetry_index[trm, :], self.__exp[p])
        for j,q in enumerate(selected_centers_j):
            trm    = self.__sym[q]
            Y[j,:] = self.__gaussian(dj, symetry_index[trm, :], self.__exp[q])
        for i, p in enumerate(selected_centers_i) :
            for j, q in enumerate(selected_centers_j):
                Q += self.__D[p,q]*X[i,:]*Y[j,:]
        return Q

    def calculateDensityByAtom(self, coordinates, atom):
        selected_centers = [i for i in range(len(self.__cent)) if self.__cent[i] == atom]
        cntr = self.__coord[self.__cent[selected_centers[0]], :]
        n_sel = len(selected_centers)
        d = coordinates - cntr
        Q = np.zeros(coordinates.shape[0])

        X = np.zeros((len(selected_centers), coordinates.shape[0]))
        for i, p in zip(range(n_sel), selected_centers) :
            trm = self.__sym[p]
            X[i,:] = self.__gaussian(d, symetry_index[trm,:], self.__exp[p])
        for i, p in zip(range(n_sel), selected_centers):
            for j, q in zip(range(n_sel), selected_centers):
                Q += self.__D[p,q]*X[i,:]*X[j,:]
        return Q
    ## GETTERS

    def getCoordinates(self):
        atom_labels = [atom_names[int(i) -1] for i in self.__charges]
        return self.__coord.copy(), atom_labels

    def getMolecularOrbitals(self):
        return self.__n_orb

    def getPrimitives(self):
        return self.__n_prim


    def read(self):
        if self.__file is None :
            raise IOError("file is not defined")
        else:
            with open(self.__file) as fh :
                wfn = fh.readlines()
            wfn = wfn[1:]
            head, wfn = self.__parseHead(wfn)
            coords, types, charges, wfn = self.__parseCoordinates(wfn, head['nuclei'])
            centers, wfn = self.__parseCenters(wfn, head['primitives'])
            syms, wfn = self.__parseSimetry(wfn, head['primitives'])
            exponents, wfn = self.__parseExponents(wfn, head['primitives'])
            coeff, occ = self.__parseOrbitals(wfn, head['primitives'], head['molecular_orbitals'])

            self.__coeff = coeff
            self.__occ   = occ
            self.__exp   = exponents
            self.__sym   = syms
            self.__cent  = centers
            self.__coord = coords
            self.__types = types
            self.__charges = charges
            self.__n_prim = head['primitives']
            self.__n_orb  = head['molecular_orbitals']
            self.__n_nuc  = head['nuclei']
    def setFile(self, file):
        self.__file = file


# class ElectronDensity() :
#     def __init__(self, verbose=False, filename = None):
#         pbl.volumetric.__init__(self, verbose=verbose, filename = filename)
#     def sample(self, number_points):
#         ## SAMPLING THROUGH PERMUTATION TO AVOID REPETITIONS
#         size = self.shape[0] * self.shape[1] * self.shape[2]
#         if number_points > size :
#             raise ArithmeticError("can't sample more data than there is")
#         sampled = np.arange(size)
#         sampled = np.random.permutation(sampled)[:number_points]
#         tnsor = self.get_volume()
#         ##
#         cx, cy , cz = np.unravel_index(sampled, self.shape)
#         itr = np.arange(number_points)
#         sampled_density = np.zeros(number_points)
#         for i,x,y,z in zip(itr, cx, cy, cz) :
#             sampled_density[i] = tnsor[x, y, z]
#         r0 = self.get_r0()
#         basis = self.get_basis()
#         X = np.array([cx, cy, cz], dtype='float64').T
#         X += r0
#         X = X.dot(basis)
#         return X, sampled_density

class GaussianLog(A2MD_basis):
    def __init__(self, fname, verbose = True):
        A2MD_basis.__init__(self, verbose=verbose, name='gaussianLog')
        self.__basis = None
        self.__coordinates = None
        self.__atom_numbers = None
        self.__dipolar = None
        self.__E = None
        self.__E2 = None
        self.__fname = fname
        self.__NPA   = None
        self.__NPT   = None

    def calculate_dipole(self):
        return np.sum(self.__coordinates * self.__NPA.reshape(
            self.__coordinates.shape[0], 1
        ), axis = 0)

    def read(self):
        coordinates_flag = False
        coordinates_barrier = 0
        charges_flag = False
        charges_lock = True

        coordinates = []
        atom_numbers     = []

        natural_charges = []
        natural_total_charges = []
        with open(self.__fname) as f :
            self.log('file opened')
            for line in f:
                line = line.strip()
                if re.match('Standard\sorientation', line):
                    coordinates = []
                    atom_numbers = []
                    coordinates_flag = True
                    self.log('Reading coordinates')
                    continue
                if coordinates_flag :
                    if re.match('-{5}', line):
                        if coordinates_barrier  == 2:
                            coordinates_barrier = False
                            coordinates_flag    = False

                            self.__coordinates  = np.array(coordinates)
                            self.__atom_numbers = np.array(atom_numbers)
                            self.log('Reading coordinates - FINISHED - %d ATOMS' % len(atom_numbers))
                        else:
                            coordinates_barrier += 1
                    elif re.match('Center\s*Atomic\s*Atomic\s*Coordinates\s\(Angstroms\)', line):
                        pass
                    elif re.match('Number\s*Number\s*Type\s*X\s*Y\s*Z', line):
                        pass
                    else :

                        cnter, atmnmbr, atmtype, x,y,z = line.split()

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
                    m = re.match('E2\s=\s*(.*)\sEUMP2\s=\s*(.*)',line)
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
                    else :

                        atmsymbl, no, charge, core, valence, rydberg, total = line.split()
                        try:
                            natural_charges.append(float(charge))
                        except ValueError :
                            pass
                        natural_total_charges.append(float(total))
                    continue
                if re.match(r'^\s*X=\s*(-?\d*.\d*)\s*Y=\s*(-?\d*.\d*)\s*Z=\s*(-?\d*.\d*)', line):
                    m = re.match(r'^\s*X=\s*(-?\d*.\d*)\s*Y=\s*(-?\d*.\d*)\s*Z=\s*(-?\d*.\d*)', line)
                    dx = float(m.group(1))
                    dy = float(m.group(2))
                    dz = float(m.group(3))
                    dipole = np.array([dx,dy,dz])
                    self.__dipolar = dipole
                    self.log('Dipole : %4.3f %4.3f %4.3f' % (dx, dy, dz))
                    continue
        return True
    def getAtomicNumbers(self):
        return self.__atom_numbers.copy()
    def getBasis(self):
        return self.__basis
    def getCharges(self):
        return self.__NPA.copy()
    def getCoordinates(self):
        return self.__coordinates.copy()
    def getDipole(self):
        return self.__dipolar.copy()
    def getEnergy(self, calctype ='MP2'):
        if calctype == 'MP2' :
            return self.__E2
        elif calctype == 'RHF' :
            return self.__E
    def getTotalCharges(self):
        return self.__NPT


class SamplingBox:
    def __init__(self, distribution='uniform'):
        self.__dist_type = distribution
        self.__X = None

    def read_coordinates(self, X):
        if isinstance(X, np.ndarray):
            if X.shape[1] == 3:
                self.__X = X
            else:
                raise IOError("can not understand coordinates")
        else:
            raise IOError("use np.ndarray for coordinates")

    def create_box(self, spacing, rotate=True, npoints = 1000):
        X = np.copy(self.__X)
        P = None
        if rotate:
            C = X.T.dot(X)
            eig, P = np.linalg.eig(C)
            X = X.dot(np.linalg.inv(P))

        min_x = X[:,0].min() - spacing
        min_y = X[:,1].min() - spacing
        min_z = X[:,2].min() - spacing

        max_x = X[:,0].max() + spacing
        max_y = X[:,1].max() + spacing
        max_z = X[:,2].max() + spacing

        u = np.random.rand(npoints, 3)
        u[:, 0] = (u[:, 0] * (max_x - min_x)) + min_x
        u[:, 1] = (u[:, 1] * (max_y - min_y)) + min_y
        u[:, 2] = (u[:, 2] * (max_z - min_z)) + min_z

        if rotate:
            u = u.dot(P)
        return u

if __name__ == '__main__':
    pass