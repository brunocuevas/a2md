import numpy as np
import copy
from a2md.support import support_angular
from a2md.support import support_inner_radial
from a2md.support import support_outer_radial
from a2md.support import support_fund_trigonometric
from a2md.support import support_fund_trigonometric_exp
from a2md.support import support_fund_trigonometric_x2exp
from a2md.support import support_angular_gaussian
from a2md.support import SupportEnsemble
from a2md.support import support
from a2md.mathfunctions import trapezoidal_3d_integral
from a2md.baseclass import A2MD_basis
from a2md.utils import convert_connectivity_tree_to_pairs
CLUSTERING_TRESHOLD_VALUE = 0.02

# ----------------------------------------------------------------------------------#
# Classes:
# -aAMD_basis class : includes some methods to display information
# -support class : includes the calls to the mathematical functions involved
#   -inner_radial
#   -outer_radial
#   -angular
# -aAMD class : calls the support functions and optimizes the corresponding coeffs.
# ----------------------------------------------------------------------------------#


# To ease saving the parameters of the simulation
#

SUPPORT_TYPE = dict(
    IR = lambda args : support_inner_radial(**args),
    OR = lambda args : support_outer_radial(**args),
    ORC = lambda args : support_outer_radial(**args),
    ORCV = lambda args : support_outer_radial(**args),
    ORV = lambda args : support_outer_radial(**args),
    ACS = lambda args : support_angular(**args),
    AS = lambda args : support_fund_trigonometric(trigo = np.sin, **args),
    AC = lambda args : support_fund_trigonometric(trigo = np.cos, **args),
    AES = lambda args : support_fund_trigonometric_exp(trigo=np.sin, **args),
    AEC = lambda args : support_fund_trigonometric_exp(trigo=np.cos, **args),
    A2S = lambda args : support_fund_trigonometric_x2exp(trigo=np.sin, **args),
    A2C = lambda args : support_fund_trigonometric_x2exp(trigo=np.cos, **args),
    AG = lambda args : support_angular_gaussian(**args)
)



class A2MD(A2MD_basis) :
    def __init__(
            self, coordinates, atomic_numbers, charge, topology, parameters=None, verbose=False, predictor=None,
    ):
        """

        :param coordinates: nuclei coordinates
        :param charge: atoms in molecule charges
        :param topology: list of linked nuclei
        :param parameters: parameters for bonding
        :param verbose:
        :param predictor:
        :param atomic_numbers
        :type coordinates: np.ndarray
        :type charge: np.ndarray
        :type topology: list
        :type parameters: list
        :type verbose: bool
        :type atomic_numbers: np.ndarray

        """
        A2MD_basis.__init__(self, name ='A2MD', verbose=verbose)
        # ------------------------------------------------#
        # SETUP PARAMETERS
        # a -> number of atoms
        # N -> array of nuclei coordinates
        # Q -> array of aim charges
        # A -> number of nuclei / aims
        # con -> connectivity
        # pars -> function parameters
        # ------------------------------------------------#
        self.__a = coordinates.shape[0]
        self.__N = coordinates.copy()
        self.__Q = charge
        self.__con = topology
        self.__pars = parameters
        self.__an = atomic_numbers
        # ------------------------------------------------#
        # OPTIMIZATION PARAMETERS
        # Xi -> list support functions
        # A  -> map chi function -> aim index
        # C -> array of weights
        # o -> number of chi functions
        # ------------------------------------------------#
        self.__A = None
        self.__F = None
        self.__Xi = None
        self.__C = None
        self.__o = None
        self.__par_name = []
        # self.log("aAMD instance declared")
        self.__gamma = 0.0001
        self.__internal_methods = ['restricted', 'unrestricted']

        self.__predictor = predictor

        if parameters is not None:
            self.read(parameters)

    @staticmethod
    def __find_in_set(index, setlist):
        i = 0
        for setlist_instance in setlist:
            if index in setlist_instance:
                return i
            i += 1
        raise RuntimeError("could not find it among setlist")

    def __gather_functions(self, atom_features):

        sa, sb = self.gather_parameters(atom_features)

        frozen_ensemble = SupportEnsemble([], 'frozen ensamble')
        isotropic_ensemble = [dict() for i in range(len(sa))]
        anisotropic_ensemble = [dict() for i in range(len(sa))]
        iso_ens2sym = []
        aniso_ens2sym = []
        retrieve_list = []
        map_ens2sym_atoms = []
        # Functions gathering

        for i, (a, xi) in enumerate(zip(self.__A, self.__Xi)):
            assert issubclass(type(xi), support)

            if self.__F[i]:
                frozen_ensemble.append(fun=xi)
                retrieve_list.append('frozen')
                continue

            fun_name = self.__par_name[i][1]
            if xi.is_anisotropic():
                b =  self.__par_name[i][2]
                ab = sorted([a, b])
                sym_group_bond = self.__find_in_set(ab, sb)
                sym_group = self.__find_in_set(a, sa)
                psi = xi.get_params()['Psi']
                label = 'aniso_{:03d}_{:03d}_{:s}_{:03d}'.format(
                    sym_group, sym_group_bond, fun_name, psi
                )
            else:
                sym_group = self.__find_in_set(a, sa)
                label = 'iso_{:03d}_{:s}'.format(
                    sym_group, fun_name
                )

            if not xi.is_anisotropic():

                if label in isotropic_ensemble[sym_group].keys():
                    assert isinstance(
                        isotropic_ensemble[sym_group][label],
                        SupportEnsemble
                    )
                    isotropic_ensemble[sym_group][label].append(xi)
                    retrieve_list.append(label)
                else:
                    isotropic_ensemble[sym_group][label] = SupportEnsemble(
                        functions=[xi], name=label
                    )
                    retrieve_list.append(label)
                    iso_ens2sym.append(sym_group)

            else:
                if label in anisotropic_ensemble[sym_group].keys():
                    assert isinstance(
                        anisotropic_ensemble[sym_group][label],
                        SupportEnsemble
                    )
                    anisotropic_ensemble[sym_group][label].append(xi)
                    retrieve_list.append(label)
                else:
                    anisotropic_ensemble[sym_group][label] = SupportEnsemble(
                        functions=[xi], name=label
                    )
                    retrieve_list.append(label)
                    aniso_ens2sym.append(self.__find_in_set(a, sa))

        # Flattening
        map_index2label = dict()
        map_index2label['frozen'] = -1
        gathered_support_functions = []
        i = 0
        i_iso = 0
        i_aniso = 0
        for sym_group in isotropic_ensemble:
            for k, ens in sym_group.items():
                assert isinstance(ens, SupportEnsemble)
                gathered_support_functions.append(ens)
                map_index2label[ens.name] = i
                i += 1
                map_ens2sym_atoms.append(iso_ens2sym[i_iso])
                i_iso += 1
            #
        for sym_group in anisotropic_ensemble:
            for k, ens in sym_group.items():
                assert isinstance(ens, SupportEnsemble)
                gathered_support_functions.append(ens)
                map_index2label[ens.name] = i
                i += 1
                map_ens2sym_atoms.append(aniso_ens2sym[i_aniso])
                i_aniso += 1

        map_fun2ensamble = [map_index2label[i] for i in retrieve_list]

        # Charges gathering

        charges = self.__Q.copy()
        gathered_charges = []

        for i, sym_atom_list in enumerate(sa):
            charge_buffer = 0
            for atom in sym_atom_list:
                charge_buffer += charges[atom]

            gathered_charges.append(
                charge_buffer
            )

        for i, (a, xi) in enumerate(zip(self.__A, self.__Xi)):
            if self.__F[i]:
                sym_group = self.__find_in_set(a, sa)
                tmp = xi.integral()
                gathered_charges[sym_group] -= tmp

        return frozen_ensemble, gathered_support_functions, map_fun2ensamble, map_ens2sym_atoms, gathered_charges

    def __not_gather_functions(self):
        """

        :return:
        """
        frozen_ensemble = SupportEnsemble([], name='frozen')
        opt_ensembles = []
        charges = self.__Q.copy()
        retrieve_list = []
        map2atoms = []
        j = 0
        for i, (a, xi) in enumerate(zip(self.__A, self.__Xi)):
            if self.__F[i]:
                frozen_ensemble.append(xi)
                charges[a] -= xi.integral()
                retrieve_list.append(-1)
                continue
            opt_ensembles.append(xi)
            retrieve_list.append(j)
            map2atoms.append(a)
            j += 1

        return frozen_ensemble, opt_ensembles, retrieve_list, map2atoms, charges

    def eval(self, x):
        """
        Evaluation of aAMD at defined coordinates

        :param x: cartessian coordinates
        :type x: np.ndarray
        :return: density values
        :type: np.ndarray
        """
        d = np.zeros(x.shape[0])
        for sup, c in zip(self.__Xi, self.__C) :
            d += c * sup.eval(x)
        return d

    def eval_core(self, x):
        d = np.zeros(x.shape[0])
        for sup, c, f in zip(self.__Xi, self.__C, self.__F):
            if f: d += c * sup.eval(x)
        return d

    def eval_ep(self, x, units='au'): # Future change: au -> HA
        """
        Evaluation of aAMD at defined coordinates

        :param x: cartessian coordinates
        :param units: au
        :type x: np.ndarray
        :return: electrostatic potential
        :type: np.ndarray
        """
        ep_electronic = np.zeros(x.shape[0])
        ep_nuclei = np.zeros(x.shape[0])
        for sup, c, a in zip(self.__Xi, self.__C, self.__A) :
            buffer = sup.eval_ep(x) * c
            ep_electronic = ep_electronic - buffer

        for i in range(self.__a):
            d = np.linalg.norm(x - self.__N[i,:], axis=1)
            mask = d < 1e-3
            d[mask] = 1e-3
            ep_nuclei += self.__an.astype('float')[i] / d

        if units == 'au':
            self.log("warning : atomic units -> Hartrees")
            return ep_electronic + ep_nuclei
        elif units == 'hartrees':
            return ep_electronic + ep_nuclei
        elif units =='kcal/mol':
            return (ep_electronic + ep_nuclei) * 627.5
        else:
            Warning("the unit was not understood. The result is in atomic units")
            return ep_electronic + ep_nuclei

    def eval_potential_ne(self, units='hartrees'):
        """

        :param units:
        :return:
        """
        v = np.zeros(self.__N.shape[0], dtype='float64')

        for sup, c, a in zip(self.__Xi, self.__C, self.__A):
            v -= sup.eval_ep(self.__N) * c * self.__an.astype('float')

        if units == 'hartrees':
            return v.sum()
        elif units == 'kcal/mol':
            return v.sum() * 627.509391

    def eval_volume(self, extend, resolution, field='density'):
        """

        This method allows to create .dx volume files so density can be represented
        in software as chimera pymol or vmd

        :param extend: size that is summed to maximum and subtracted from minimum to each coordinate
        to setup the size of the box
        :param field
        :type extend: float
        :param resolution: size of the side of each cubic voxel
        :type resolution: float
        :return: 3D tensor with density values
        """
        from a2mdlib.qm import ElectronDensity
        if field == 'density':
            fun = self.eval
        elif field == 'ep':
            fun = lambda x : self.eval_ep(x, units='au')
        else:
            raise IOError("not found")

        minx = self.__N[:, 0].min() - extend
        miny = self.__N[:, 1].min() - extend
        minz = self.__N[:, 2].min() - extend
        maxx = self.__N[:, 0].max() + extend
        maxy = self.__N[:, 1].max() + extend
        maxz = self.__N[:, 2].max() + extend
        xx = np.arange(minx, maxx, resolution)
        yy = np.arange(miny, maxy, resolution)
        zz = np.arange(minz, maxz, resolution)
        DX = np.zeros((xx.size, yy.size, zz.size))
        for ix in range(xx.size) :
            r = np.zeros((zz.size, 3), dtype='float64')
            r[:,0] = xx[ix]
            for iy in range(yy.size) :
                r[:,1] = yy[iy]
                r[:,2] = zz[:]
                DX[ix,iy,:] = fun(r)

        vol_density = ElectronDensity(
            verbose=False
        )
        vol_density.set_r0(np.array([minx, miny, minz]) * 0.5292)
        vol_density.set_basis(np.identity(3) * resolution * 0.5292)
        vol_density.set_volume(DX)
        return vol_density

    def gather_parameters(self, atom_features):
        """

        This function's purpose is to create a list of equivalent atoms and equivalent bonds based in an arbitrary
        matrix of features. When inside self.optimize it allows to start an optimization in which coefficents
        are fitted the same way for all atoms belonging

        :param atom_features:
        :type atom_features: np.array
        :return:
        """
        labels = self.__an
        bonds = np.array(convert_connectivity_tree_to_pairs(self.__con))

        symmetric_atoms = []
        symmetric_bonds = []
        # Distance matrices for clustering of bonds and atoms
        atom_distance_matrix = np.zeros((atom_features.shape[0], atom_features.shape[0]), dtype='float64')
        bond_distance_matrix = np.zeros((bonds.shape[0], bonds.shape[0]), dtype='float64')
        n_feats = atom_features.shape[1]
        comp_vector_1 = np.zeros(2 * n_feats, dtype='float64')
        comp_vector_2 = np.zeros(2 * n_feats, dtype='float64')

        for i in range(atom_features.shape[0]):
            for j in range(atom_features.shape[0]):
                if labels[i] == labels[j]:
                    atom_distance_matrix[i, j] = np.linalg.norm(atom_features[i, :] - atom_features[j, :])
                else:
                    atom_distance_matrix[i, j] = 1000.0

        for i, b1 in enumerate(bonds):
            comp_vector_2[:n_feats] = atom_features[b1[0], :]
            comp_vector_2[n_feats:] = atom_features[b1[1], :]
            bond1_labels = [labels[bonds[i][0]], labels[bonds[i][1]]]
            for j, b2 in enumerate(bonds):
                bond2_labels = [labels[bonds[j][0]], labels[bonds[j][1]]]
                if sorted(bond1_labels) != sorted(bond2_labels):
                    bond_distance_matrix[i, j] = 1000.0
                else:
                    # aceptor-end
                    comp_vector_1[:n_feats] = atom_features[b2[0], :]
                    comp_vector_1[n_feats:] = atom_features[b2[1], :]
                    ae = np.linalg.norm(comp_vector_1 - comp_vector_2)
                    # end-aceptor
                    comp_vector_1[n_feats:] = atom_features[b2[0], :]
                    comp_vector_1[:n_feats] = atom_features[b2[1], :]
                    ea = np.linalg.norm(comp_vector_1 - comp_vector_2)

                    bond_distance_matrix[i, j] = np.min([ae, ea])


        atom_distance_matrix[atom_distance_matrix > CLUSTERING_TRESHOLD_VALUE] = 1.0
        atom_distance_matrix[atom_distance_matrix < CLUSTERING_TRESHOLD_VALUE] = 0.0
        atom_distance_matrix = 1 - atom_distance_matrix

        bond_distance_matrix[bond_distance_matrix > CLUSTERING_TRESHOLD_VALUE] = 1.0
        bond_distance_matrix[bond_distance_matrix < CLUSTERING_TRESHOLD_VALUE] = 0.0
        bond_distance_matrix = 1 - bond_distance_matrix

        placed_atoms = []
        placed_bonds = []

        for i in range(atom_features.shape[0]):
            if i in placed_atoms:
                continue
            ngroup = np.sum(atom_distance_matrix[:, i])
            if ngroup == 1:
                symmetric_atoms.append([i])
            else:
                tmp = []
                for j in range(i, atom_features.shape[0]):
                    if j not in placed_atoms:
                        if atom_distance_matrix[i, j] == 1.0:
                            tmp.append(j)
                            placed_atoms.append(j)
                symmetric_atoms.append(tmp)
        #bonds = bonds - 1
        for i in range(bonds.shape[0]):
            if i in placed_bonds:
                continue
            ngroup = np.sum(bond_distance_matrix[:, i])
            if ngroup == 1:

                symmetric_bonds.append([bonds[i].tolist()])
            else:
                tmp = []
                for j in range(i, bonds.shape[0]):
                    if j not in placed_bonds:
                        if bond_distance_matrix[i, j] == 1.0:
                            tmp.append(bonds[j].tolist())
                            placed_bonds.append(j)
                symmetric_bonds.append(tmp)
        return symmetric_atoms, symmetric_bonds

    def get_a2md_atomic_charge(self, atom):
        """

        :param atom:
        :return:
        """
        if type(atom) is int:
            if atom < self.__a:
                charge = 0
                for i, (xi, c) in enumerate(zip(self.__Xi, self.__C)):
                    if self.__A[i] == atom:
                        charge += c * xi.integral()
                    else:
                        pass
            else:
                raise IOError(f"there is no {atom} in this a2md system")
        else:
            raise IOError("only integer atom index are acepted")
        return charge

    def get_a2md_molecule_charge(self):
        """

        :return:
        """

        charge = np.zeros(self.__a, dtype='float64')
        for i, (xi, c) in enumerate(zip(self.__Xi, self.__C)):
            charge[self.__A[i]] += c * xi.integral()
        return charge

    def get_atomic_numbers(self):
        """

        :return:
        """
        return self.__an.copy()

    def get_coordinates(self):
        """

        :return:
        """
        return self.__N.copy()

    def get_dipole(self):
        """

        :return:
        """
        if self.__an is None:
            raise IOError("atomic numbers are not defined")

        atomic_dipole  = np.zeros((self.__a, 3), dtype='float64') # \mu
        atomic_charges_neg = np.zeros((self.__a, 3), dtype='float64') # q
        atomic_charges_pos = np.zeros((self.__a, 3), dtype='float64') # n
        for i in range(self.__a):
            atomic_charges_pos[i, :] = self.__N[i, :] * self.__an[i]

        for i, (c, xi, fn) in enumerate(zip(self.__C, self.__Xi, self.__par_name)):
            j = fn[0]
            assert issubclass(type(xi), support)

            if isinstance(xi, support_fund_trigonometric):
                if xi.get_params()['Psi'] == 0:
                    k = fn[2]
                    bonding_vector = self.__N[k, :] - self.__N[j, :]
                    bonding_vector /= np.linalg.norm(bonding_vector)
                    atomic_dipole[j, :] -= xi.first_moment()[0] * c * bonding_vector
                else:
                    atomic_charges_neg[j, :] -= (xi.integral() * c) * self.__N[j, :]
            elif isinstance(xi, support_angular_gaussian):
                k = fn[2]
                bonding_vector = self.__N[k, :] - self.__N[j, :]
                bonding_vector /= np.linalg.norm(bonding_vector)
                atomic_dipole[j, :] -= xi.first_moment() * c * bonding_vector
            else:
                atomic_charges_neg[j, :] -= (xi.integral() * c) * self.__N[j, :]

        atomic_dipole = atomic_dipole + atomic_charges_pos + atomic_charges_neg
        return  atomic_dipole.sum(axis=0) / 0.4

    def get_function_names(self):
        """

        :return:
        """
        return self.__par_name

    def get_gamma(self):
        """

        :return:
        """
        return self.__gamma

    def get_parametrization(self):
        """
        This method allows to create a list of function parameters

        :return:
        :rtype: list
        """
        atom_parametrization_list = [None] * self.__o

        xi_iter = 0
        for xi, fn, cc in zip(self.__Xi, self.__par_name, self.__C) :
            xi_params = xi.get_params()
            tmp_input = dict(center = fn[0], support_type = fn[1], coefficient=cc)
            tmp_params = dict()
            for key, item in xi_params.items() :
                tmp_params[key] = item
            tmp_input['params'] = tmp_params
            if fn[1][0] == 'A' :
                try :
                    tmp_input['bond'] = fn[2]
                except IndexError:
                    raise IOError("angular support functions should have a bonded atom")
            else :
                tmp_input['bond'] = None

            atom_parametrization_list[xi_iter] = tmp_input
            xi_iter += 1

        return atom_parametrization_list

    def get_topology(self):
        """

        :return:
        """
        return self.__con.copy()

    def load_coefficients(self, params):
        """
        This function is a workaround to the problem of serialization of generalized support
        functions. It works:
        1. Parametrizing using the generalizable proprocessor
        2. Once parametrized, open a json file of a previous run using this method
        :return:
        """
        import json
        with open(params) as f:
            params_list = json.load(f)

        for i, (lparams, fn) in enumerate(zip(params_list, self.__par_name)):
            if fn[0] == lparams['center'] and fn[1] == lparams['support_type']:
                self.__C[i] = params_list[i]['coefficient']

    def optimize(
            self, target_coordinates, target_density,
            weigths=None, method='restricted',
            atom_features=None
    ):
        """
        sets the coefficients associated to each support function by
        minimizing the lagrangian function:
        L = (\rho - \rho')^2 - \sum \lambda_{i} (\int \rho_{i}' dV - q_{i})
        where :
            - rho : ab-initio density
            - rho' : aAMD density
            - q_{i} : charge of atom i
        :param target_coordinates: coordinates of nuclei atoms
        :param target_density: density values
        :param weigths
        :param method
        :param atom_features: matrix of any atom feature from which symmetry between atoms of the system can be found
        :type target_coordinates: np.ndarray
        :type target_density: np.ndarray
        :type method: str
        :type atom_features: np.ndarray
        :return: coefficients
        :rtype: np.ndarray
        """
        if method not in ['restricted', 'unrestricted']:
            raise IOError("unidentified option. use restricted or unrestricted")

        # Gathering protocol

        # The main idea here is that some atoms are symmetric and, therefore, have the same properties.
        # We don't get into heavy group theory; we just observe the atomic environment vectors of
        # the different atoms, and we cluster those that have very similar environments.

        # Once clustered, the functions are grouped in ensembles of functions, and operations upon them
        # are done as if they were an only function. The symmetric atoms charges are also clustered.

        if atom_features is not None:

            F, Xi, MF, A, Q = self.__gather_functions(atom_features)

        else:

            F, Xi, MF, A, Q = self.__not_gather_functions()


        m = target_coordinates.shape[0] # number of sampling points
        X = target_coordinates.copy() # coordinates of those sampling points
        if weigths is not None:
            W = weigths.copy() # weights
        else:
            W = np.ones(m, dtype='float64') / m
        R = target_density.copy() # density of the sampling points

        n_coefficients = len(Xi)
        n_restrictions = len(Q)

        # Substract the charge of the frozen part of the density

        R -= F.eval(X)

        # Creating arrays to store intermediate values during optimization

        D = np.zeros((n_coefficients, m), dtype='float64')
        if method == 'restricted':
            P = np.zeros((n_coefficients + n_restrictions), dtype='float64')
            B = np.zeros(
                (n_coefficients + n_restrictions, n_coefficients + n_restrictions), dtype='float64'
            )
            Q = np.array(Q)
        elif method == 'unrestricted':
            P = np.zeros((n_coefficients + 1), dtype='float64')
            B = np.zeros(
                (n_coefficients + 1, n_coefficients + 1), dtype='float64'
            )
            Q = np.array(Q).sum()
        else:
            raise NotImplementedError("only restricted or unrestricted")

        # Calculating denisty upon the sampling coordinates

        for i, xi in enumerate(Xi):
            D[i, :] = xi.eval(X)

        # Filling array with L2 norm derivatives and regularization (upon pred density)
        effective_gamma = self.__gamma / n_coefficients
        B[:n_coefficients, :n_coefficients] = 2*(D*W).dot(D.T)
        B[:n_coefficients, :n_coefficients] = B[:n_coefficients, :n_coefficients] + \
            2 * effective_gamma * np.identity(n_coefficients, dtype='float64')

        # Filling array with charges

        if method == 'restricted':
            for i, xi in enumerate(Xi):
                j = A[i]
                B[n_coefficients + j, i] += xi.integral()
                B[i, n_coefficients + j] += xi.integral()
        elif method == 'unrestricted':
            for i, xi in enumerate(Xi):
                B[n_coefficients, i] += xi.integral()
                B[i, n_coefficients] += xi.integral()

        # Filling array with L2 norm derivatives (upon sampled density)

        P[:n_coefficients] = 2*(D*W).dot(R)
        P[n_coefficients:] = Q

        # Solving problem:
        #   BC = P
        # where:
        #   B is the matrix containing the L2 derivatives that depend upon the
        #       coefficients
        #   P is the matrix containing the L2 derivatives that do not depend upon
        #       the coefficients
        #   C is the coefficients vector that we want to solve

        C = np.linalg.solve(B, P)[:n_coefficients]

        C = self.unpack(C, MF)

        self.__C = C

        return self.__C.copy()


    def predict_parameters(self, predictor=None, normalized=False, integrate=False):
        """
        This method allows to use a regression function to predict
        the value of the internal coefficients of the model

        :param predictor: an object with an eval method that takes the coordinates and the atomic numbers,
        :param normalized: whether to fit or not the output to keep the charge constant
        and returns a matrix with two coefficents per atom, the inner and the outter one
        :param integrate:

        :return:
        """
        if predictor is None and self.__predictor is None:
            raise IOError("there is no predictor set")
        else:
            if predictor is None:
                current_pred = self.__predictor
            else:
                current_pred = predictor

        for par_name in self.__par_name:
            if par_name[1] not in ['IR', 'OR']:
                raise NotImplementedError("only prediction of isotropic component is available")

        assert isinstance(self.__C, np.ndarray)
        C_matrix = current_pred.eval(self.__N, self.__an)
        C_vector = np.zeros(self.__C.shape, dtype='float64')

        for i, par_name in enumerate(self.__par_name):
            if par_name[1] == 'IR':
                C_vector[i] = C_matrix[par_name[0], 0]
            else:
                C_vector[i] = C_matrix[par_name[0], 1]

        if not normalized:
            q = np.sum(self.__Q)
            n_c = C_vector.shape[0]
            A = np.zeros((n_c + 1, n_c + 1), dtype='float64')

            for i, fun in enumerate(self.__Xi):
                A[i, i] = 2.0
                A[-1, i] = fun.integral()
                A[i, -1] = fun.integral()

            b = np.zeros((n_c + 1), dtype='float64')
            b[:n_c] = 2*C_vector
            b[-1] = q
            x = np.linalg.solve(A, b)
            self.__C = x[:-1]
        else:
            self.__C = C_vector

        if integrate:
            extend = 3.0
            resolution = 5e-2
            domain = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            domain[0] = self.__N[:, 0].min() - extend
            domain[2] = self.__N[:, 1].min() - extend
            domain[4] = self.__N[:, 2].min() - extend
            domain[1] = self.__N[:, 0].max() + extend
            domain[3] = self.__N[:, 1].max() + extend
            domain[5] = self.__N[:, 2].max() + extend
            I, N = trapezoidal_3d_integral(
                domain=domain,
                fun=self.eval,
                resolution=resolution,
                sign_operator=True
            )
            print(I)
            return I

    def read(self, params):
        """

        This method allows to reparametrize from previous runs

        :param params: a list of parameters for support functions
        :type params: list
        :return:
        """
        if self.__Xi is not None:
            if not len(self.__Xi) == 0 :
                self.log("existing support functions will be erased")

        Xi = []
        C = []
        A = []
        F = []
        for xi_iter in range(len(params)) :
            ppp = params[xi_iter]
            A.append(ppp['center'])
            C.append(ppp['coefficient'])
            try:
                F.append(ppp['frozen'])
            except KeyError:
                F.append(False)
            input_dict = copy.copy(ppp['params'])
            input_dict['coordinates'] = self.__N[ppp['center'],:]
            Xi.append(SUPPORT_TYPE[ppp['support_type']](input_dict))
            if not ppp['bond'] is None :

                bonding_axis = self.__N[ppp['bond'],:] - self.__N[ppp['center'],:]
                Xi[-1].set_reference_frame(bonding_axis)
                self.__par_name.append((ppp['center'], ppp['support_type'], ppp['bond']))
            else:
                self.__par_name.append((ppp['center'], ppp['support_type']))
        C = np.array(C, dtype='float64')
        self.__C = np.array(C)
        self.__Xi = Xi
        self.__A = A
        self.__F = F
        self.__o = len(C)

    def save_model(self, file):
        """

        :param file:
        :type file: str
        :return:
        """

        # def ca2md_line(center, kind, bonded2, A, B, k, n, C):
        #     return
        #     )
        if self.__an is None:
            raise IOError("missing atomic numbers")
        n_radials = 0
        n_sines = 0
        n_cosines = 0
        for i, fn in enumerate(self.__par_name):
            if fn[1] in ['OR', 'IR']:
                n_radials += 1
            elif fn[1] in ['AS', 'AES']:
                n_sines += 1
            elif fn[1] in ['AC', 'AEC']:
                n_cosines +=1

        file_list = [
            "<NUMBER_ITEMS>",
            '\t{:d},{:d},{:d},{:d},{:d}'.format(self.__a, self.__o, n_radials, n_sines, n_cosines),
            "</NUMBER_ITEMS>",
            '<COORDS>'
        ]

        for i in range(self.__a):
            file_list.append(
                '\t{:>12.4f},{:>12.4f},{:>12.4f},{:>3d}'.format(
                    self.__N[i, 0], self.__N[i, 1], self.__N[i, 2],
                    self.__an[i]
                )
            )
        file_list.append('</COORDS>')

        file_list.append('<FUNCTIONS>')

        for i, (fn, xi, c) in enumerate(zip(self.__par_name, self.__Xi, self.__C)):

            pars = xi.get_params()

            if fn[1] == 'IR':
                file_list.append(
                    '\t{:<d},{:>3s},{:>3d},{:>12.4e},{:>12.4e},{:>6.3f},{:>6d},{:>12.4e}'.format(
                        fn[0], 'R', -1, pars['A1'], pars['B1'], 0, 0, c) + \
                    '\t{:<d},{:>3s},{:>3d},{:>12.4e},{:>12.4e},{:>6.3f},{:>6d},{:>12.4e}'.format(
                        fn[0], 'R', -1, pars['A2'], pars['B2'], 0, 0, c)
                    )
            elif fn[1] in ['OR', 'ORC', 'ORCV', 'ORV']:
                file_list.append(
                    '\t{:<d},{:>3s},{:>3d},{:>12.4e},{:>12.4e},{:>6.3f},{:>6d},{:>12.4e}'.format(
                        fn[0], 'R', -1, pars['A3'], pars['B3'], 0, 0, c)
                )
            elif fn[1] == 'AS':
                file_list.append(
                    '\t{:<d},{:>3s},{:>3d},{:>12.4e},{:>12.4e},{:>6.3f},{:>6d},{:>12.4e}'.format(
                        fn[0], 'S', fn[2], pars['U'], pars['G'], pars['Psi'], 1, c)
                )
            elif fn[1] == 'AC':
                file_list.append(
                    '\t{:<d},{:>3s},{:>3d},{:>12.4e},{:>12.4e},{:>6.3f},{:>6d},{:>12.4e}'.format(
                        fn[0], 'C', fn[2], pars['U'], pars['G'], pars['Psi'], 1, c)
                )
            elif fn[1] == 'AG':
                if pars['Psi'] == 0:
                    file_list.append(
                        '\t{:<d},{:>3s},{:>3d},{:>12.4e},{:>12.4e},{:>6.3f},{:>6d},{:>12.4e}'.format(
                            fn[0], 'G', fn[2], pars['U'], pars['G'], pars['Alpha'], 1, c)
                    )
                else:
                    file_list.append(
                        '\t{:<d},{:>3s},{:>3d},{:>12.4e},{:>12.4e},{:>6.3f},{:>6d},{:>12.4e}'.format(
                            fn[0], 'G', fn[2], pars['U'], pars['G'], pars['Alpha'], 1, c)
                    )

        file_list.append('</FUNCTIONS>')
        file_list.append('<META>')
        file_list.append('\tBy a2md')
        file_list.append('</META>')

        with open(file, 'w', newline='\n') as f:
            f.write('\n'.join(file_list))

    def set_gamma(self, gamma):
        """

        :param gamma:
        :return:
        """
        self.__gamma = gamma

    def unpack(self, C, map_fun2ens):
        C_f = np.ones(self.__o, dtype='float64')
        for i, f in enumerate(map_fun2ens):
            if f == -1:
                C_f[i] = 1.0
                continue
            C_f[i] = C[f]
        return C_f