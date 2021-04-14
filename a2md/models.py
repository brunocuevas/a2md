import numpy as np
import copy
from a2md.support import SupportRadial
from a2md.support import SupportAngular
from a2md.support import SupportEnsemble
from a2md.support import SupportHarmonic
from a2md import TOPO_RESTRICTED_PARAMS, HARMONIC_TOPO_RESTRICTED_PARAMS, EXTENDED_TOPO_RESTRICTED_PARAMS
from a2md import SPHERICAL_PARAMS
from a2md.baseclass import A2MDBaseClass
from a2mdio.molecules import Mol2, PDB
from a2md.utils import convert_connectivity_tree_to_pairs
from a2mdio import PDB_PROTEIN_TYPE_CHARGES, PDB_PROTEIN_CHARGES, PDB_PROTEIN_TYPES, PDB_PROTEIN_TOPOLOGY
from  typing import List, Union, Dict, Callable
CLUSTERING_TRESHOLD_VALUE = 0.02


ELEMENT2AN = dict(
    H=1,
    C=6,
    N=7,
    O=8
)
AN2ELEMENT = [
    '',
    'H', 'He',
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    '',  '',  '',  '', 'P', 'S',  '',   ''
]

SUPPORT_TYPE = {
    "_SPHERIC"  : lambda args : SupportRadial(**args),
    "_GAUSSIAN" : lambda args : SupportAngular(**args),
    "_HARMONIC" : lambda  args : SupportHarmonic(**args)
}

def a2md_from_mol(mol : Mol2):
    """
    returns a a2md.models.Molecule associated to a Mol2
    ---
    :param mol:
    :type mol: Mol2
    :return:
    :rtype: Molecule
    """
    an = mol.get_atomic_numbers()
    coords = mol.get_coordinates(units="au")
    charge = mol.get_absolute_charges()
    topology = mol.get_bonds()
    segments = mol.segment_idx
    topo_array = []
    for i in range(mol.get_number_atoms()): topo_array.append([])
    for i in range(mol.get_number_bonds()):
        begin = int(topology[i, 0])
        end = int(topology[i, 1])
        topo_array[begin - 1].append(end - 1)
        topo_array[end - 1].append(begin - 1)
    return Molecule(
        coordinates=coords, atomic_numbers=an, charge=charge, topology=topo_array,
        segments=segments
    )

def polymer_from_pdb(mol : PDB, chain : str ):
    """
    returns the polymer density model associated to a Mol2
    Still in development
    :param mol:
    :param chain:
    :return:
    """
    an = mol.get_atomic_numbers()
    coords = mol.get_coordinates(units='au')
    charge = mol.get_absolute_charges()
    topology = mol.get_bonds()
    topo_array = []
    for i in range(mol.get_number_atoms()): topo_array.append([])
    for i in range(mol.get_number_bonds()):
        begin = int(topology[i][0])
        end = int(topology[i][1])
        topo_array[begin].append(end)
        topo_array[end].append(begin)

    return Polymer(
        coordinates=coords, atomic_numbers=an, atom_labels=mol.atom_names,
        topology=topo_array, charge=charge, atom_residx=mol.atom_residues_idx,
        sequence=mol.sequence[chain], residue_names=mol.residue_names,
        residue_idx=mol.residue_idx, parameters=None, verbose=False
    )



class Molecule(A2MDBaseClass) :
    parametrization_default = TOPO_RESTRICTED_PARAMS
    parametrization_harmonic = HARMONIC_TOPO_RESTRICTED_PARAMS
    parametrization_extended = EXTENDED_TOPO_RESTRICTED_PARAMS
    parametrization_spherical = SPHERICAL_PARAMS
    def __init__(
            self, coordinates : np.ndarray, atomic_numbers : Union[List[int], np.ndarray],
            charge : Union[List[int], np.ndarray], topology : List[List[int]],
            parameters=None, verbose=False,
            atom_labels=None, segments=None
    ):
        """
        A2MD,models.molecule

        This object stores the information and the methods to generate models of electron density
        based on the linear combination of exponential-family functions.

        :param coordinates: nuclei coordinates
        :param charge: charge of each of the atoms
        :param topology: list of atoms joined to each atom
        :param parameters: electron density model parameters
        :param verbose:
        :param atomic_numbers: atomic number of each atom (integer)
        :param segments: to use in polymers
        :type coordinates: np.ndarray
        :type charge: np.ndarray
        :type topology: list
        :type parameters: list
        :type verbose: bool
        :type atomic_numbers: np.ndarray

        """
        A2MDBaseClass.__init__(self, name ='A2MD', verbose=verbose)
        self.natoms = coordinates.shape[0]
        self.coordinates = coordinates.copy()
        self.atom_charges = charge
        self.topology = topology
        self.parameters = parameters
        self.atomic_numbers = atomic_numbers
        self.map_function2center = None
        self.map_frozenfunctions = None
        self.functions = None
        self.opt_params = None
        self.nfunctions = None
        self.atom_labels = atom_labels
        self.function_names = []
        self.function_types = []
        self.regularization = 0.0001
        self.is_clusterized = False
        self.is_optimized = False
        self.segments = segments  # segments allows to perform semi-restricted optimizations, in which
        # a part of the molecule corresponding to a given segment must match the sum of charges of that
        # fragment. It's an intermediate step between fully atom-charge restriction and molecule-wide single
        # restriction.
        if segments is not  None:
            self.nsegments = np.unique(np.array(self.segments)).size

        if parameters is not None:
            self.read(parameters)

    def clustering(self, clusterizer : Callable):
        """
        a2md.models.Molecule.clusterize
        ---
        This function converts a problem of N atoms and F functions to a smaller problem
        by impossing relationships between the atoms and bonds of the molecule. These relationships
        are impossed by the clusterizer, which returns a list of clustered atoms and bonds according
        to some policy (mostly, symmetry relationships). It is important to notice: the clustering of bonds must take
        place over the clustering of atoms. Otherwise, charge normalization becomes difficult to track.

        Once done, clusterize converts the parameterization into a smaller one by joining charges and by joining
        functions into ensssemble functions. This function is transparent to optimize.

        :param clusterizer: callable. Takes as input (labels, topology, coordinates)

        """
        import copy

        x = self.get_coordinates()
        l = self.get_symbols()
        t = convert_connectivity_tree_to_pairs(self.get_topology())
        sa, sb = clusterizer(l, t, x)

        atom_charges = copy.copy(self.atom_charges)
        map_frozenfunction = copy.copy(self.map_frozenfunctions)
        functions = copy.copy(self.functions)
        function_names = copy.copy(self.function_names)
        function_types = copy.copy(self.function_types)
        expanded_fns = []

        for i in range(len(function_names)):

            bond_ = function_names[i][2]
            name_ = function_names[i][1]
            type_ = function_types[i]

            expanded_fns.append(
                dict(
                    type=type_,
                    center=function_names[i][0],
                    bond=bond_,
                    name=name_
                )
            )
            del bond_, name_, type_

        # Clusterize charges

        new_atom_charges = []
        for cluster in sa:
            current_charge = 0.0
            for atom_idx in cluster:
                current_charge += atom_charges[atom_idx]
            new_atom_charges.append(current_charge)
            del current_charge

        # Clusterize functions
        new_functions = []
        new_map2atom = []
        new_names = []
        new_frozenfunction = []
        new_functiontypes = []

        # Part 1, cluster isotropic functions

        for cluster_idx, cluster in enumerate(sa):

            # Part 1.1. Find types of functions

            current_cluster_function_names = []
            for atom_idx in cluster:
                for i, fun_info in enumerate(expanded_fns):
                    c1 = fun_info['center'] == atom_idx
                    c2 = fun_info['name'] not in current_cluster_function_names
                    c3 = fun_info['bond'] is None
                    if all([c1, c2, c3]):
                        current_cluster_function_names.append(fun_info['name'])
                    del c1, c2, c3

            # Part 1.1. Cluster isotropic functions

            for fun_name in current_cluster_function_names:
                is_frozen = []
                current_cluster_functions = []
                current_map2atoms = []
                current_map2bonds = []
                for atom_idx in cluster:
                    for i, fun_info in enumerate(expanded_fns):
                        c1 = fun_info['center'] == atom_idx
                        c2 = fun_info['name'] == fun_name
                        c3 = fun_info['bond'] is None
                        if all([c1, c2, c3]):
                            current_cluster_functions.append(functions[i])
                            is_frozen.append(map_frozenfunction[i])
                            current_map2atoms.append(atom_idx)
                            current_map2bonds.append(None)
                        del c1, c2, c3


                new_functions.append(
                    SupportEnsemble(
                        functions=current_cluster_functions,
                        name="ensemble_c{:d}_f{:s}".format(cluster_idx, fun_name),
                        map2atoms=current_map2atoms, map2bonds=current_map2bonds
                    )
                )

                if not all(is_frozen) and any(is_frozen):
                    self.log("there are hetereogeneous functions in the current group")

                new_frozenfunction.append(all(is_frozen))
                new_names.append([cluster_idx, "ENS_{:s}".format(fun_name)])
                new_map2atom.append(cluster_idx)
                new_functiontypes.append("_SPHERICAL")  # Harcoded


        # Part 2, cluster  anisotropic functions
        # bond cluster must contain the idx of the clusters of sa
        # so we can map atoms involved in this bond cluster.
        # NOTE: this operation is assymetric!!!
        for bond_cluster in sb:

            cluster_0 = sa[bond_cluster[0]]
            cluster_1 = sa[bond_cluster[1]]
            # only functions that come from cluster 0 are required
            current_cluster_function_names = []
            current_cluster_function_types = []
            for atom_idx in cluster_0:

                for i, fun_info in enumerate(expanded_fns):
                    c1 = fun_info['center'] == atom_idx
                    c2 = fun_info['name'] not in current_cluster_function_names
                    c3 = fun_info['bond'] in cluster_1
                    if all([c1, c2, c3]):
                        current_cluster_function_names.append(fun_info['name'])
                        current_cluster_function_types.append(fun_info['type'])
                    del c1, c2, c3


            for fun_name, fun_type in zip(current_cluster_function_names, current_cluster_function_types):
                current_cluster_functions = []
                current_map2atoms = []
                current_map2bonds = []
                for atom_idx in cluster_0:
                    for i, fun_info in enumerate(expanded_fns):
                        c1 = fun_info['center'] == atom_idx
                        c2 = fun_info['name'] == fun_name
                        c3 = fun_info['bond'] in cluster_1
                        if all([c1, c2, c3]):
                            current_cluster_functions.append(functions[i])
                            current_map2atoms.append(atom_idx)
                            current_map2bonds.append(fun_info['bond'])
                        del c1, c2, c3

                new_functions.append(
                    SupportEnsemble(
                        functions=current_cluster_functions,
                        name="ensemble_c{:d}_f{:s}_b{:d}".format(bond_cluster[0], fun_name, bond_cluster[1]),
                        map2atoms=current_map2atoms, map2bonds=current_map2bonds
                    )
                )
                new_names.append([bond_cluster[0], "ENS_{:s}".format(fun_name), bond_cluster[1]])
                new_map2atom.append(bond_cluster[0])
                new_frozenfunction.append(False)
                new_functiontypes.append(fun_type)

        self.natoms = len(sa)
        self.atom_charges = new_atom_charges
        self.map_function2center = new_map2atom
        self.map_frozenfunctions = new_frozenfunction
        self.function_names = new_names
        self.function_types = new_functiontypes
        self.functions = new_functions
        self.opt_params = np.zeros(len(sa), dtype="float64")
        self.nfunctions = len(new_functions)
        self.is_clusterized = True

    def eval(self, x : np.ndarray, kind='density'):
        """
        Evaluation of A2MD at defined coordinates.
        Electron density (-density-) or electrostatic potential (-ep-)

        :param x: cartessian coordinates
        :param kind: either density or ep
        :type x: np.ndarray
        :return: density values  (e-/Bohr^3)
        :type: np.ndarray
        """
        d = np.zeros(x.shape[0])
        for i, (sup, c) in enumerate(zip(self.functions, self.opt_params)):
            if kind == 'density':
                d += c * sup.eval(x)
            elif kind == 'ep':
                d -= c * sup.eval_ep(x)
            else:
                raise NotImplementedError("only density or ep")

        if kind == 'ep':
            d += self.eval_nuclear_potential(x)

        return d

    def eval_by_fun(self, x : np.ndarray, i : int):
        """
        Evaluation of a specific A2MD function at defined coordinates.
        Function must be defined by its index

        :param x: cartessian coordinates
        :param i: function index
        :type x: np.ndarray
        :type i: int
        :return: density values (e-/Bohr^3)
        :type: np.ndarray
        """
        sup = self.functions[i]
        c = self.opt_params[i]
        d = c * sup.eval(x)
        return d

    def eval_core(self, x : np.ndarray):
        """
        Evaluation of Frozen function (representing core function) at defined coordinates.

        :param x: cartessian coordinates
        :type x: np.ndarray
        :return: density values
        :type: np.ndarray
        """
        d = np.zeros(x.shape[0])
        for sup, c, f in zip(self.functions, self.opt_params, self.map_frozenfunctions):
            if f: d += c * sup.eval(x)
        return d

    def eval_by_name(self, x, fun_name):
        """
        Evaluation of a specific A2MD function at defined coordinates.
        Function must be defined by its name

        :param x: cartessian coordinates
        :param fun_name: function name
        :type x: np.ndarray
        :type fun_name: str
        :return: density values (e-/Bohr^3)
        :type: np.ndarray
        """
        d = np.zeros(x.shape[0])
        for sup, c, name in zip(self.functions, self.opt_params, self.function_names):
            if name[1] in fun_name:
                d += c * sup.eval(x)
        return d

    def eval_nuclear_potential(self, x):
        """
        Evaluation of nuclear electrostatic potential

        :param x: cartessian coordinates

        :type x: np.ndarray
        :return: nuclear potential (Ha)
        :type: np.ndarray
        """
        r = self.coordinates.copy()
        v = np.zeros(x.shape[0], dtype='float64')
        for i, z in enumerate(self.atomic_numbers):
            d = x - r[i, :].reshape(1, -1)
            d = np.linalg.norm(d, axis=1) + 1e-18
            v += z/d
        return v

    def eval_volume(self, spacing, resolution, kind='density', cutoff=None):
        """
        Evaluation of electron density in a grid. Useful to visualize electron density.
        Resulting object can be saved as a .dx file and visualized in Chimera and other
        tools.

        :param spacing: amount of space around the min and max of the molecule.
        :param resolution: use Bohr as unit of length
        :param kind: either electron density (-density-) or electrostatic potential (-ep-)
        :param cutoff: truncates values exceding some value.
        :return:
        """
        from a2mdio.qm import ElectronDensity
        minx, miny, minz = (self.coordinates - spacing).min(axis=0)
        maxx, maxy, maxz = (self.coordinates + spacing).max(axis=0)

        xx = np.arange(minx, maxx, resolution)
        yy = np.arange(miny, maxy, resolution)
        zz = np.arange(minz, maxz, resolution)

        dx = np.zeros((xx.size, yy.size, zz.size))
        if cutoff is None:
            f = lambda x : self.eval(x, kind=kind)
        else:
            def f(x):
                p = self.eval(x, kind=kind)
                p[p>cutoff] = cutoff
                return p

        for ix in range(xx.size):
            r = np.zeros((zz.size, 3), dtype='float64')
            r[:, 0] = xx[ix]
            for iy in range(yy.size):
                r[:, 1] = yy[iy]
                r[:, 2] = zz[:]
                dx[ix, iy, :] = f(r)

        vol_density = ElectronDensity(verbose=False)
        vol_density.set_r0(np.array([minx, miny, minz]) * 0.5292)
        vol_density.set_basis(np.identity(3) * resolution * 0.5292)
        vol_density.set_volume(dx)
        return vol_density

    def inflate(self):
        """
        This operation aims to revert the clustering. This is specially useful to write parametrizations. While in
        the clustering process, the coefficient value is lost, here it is preserved.
        Atom charge can not be recovered, so don't perform restricted optimizatons on inflated models.
        :return:
        """
        if not self.is_clusterized:
            return True

        new_functions = []
        new_map2atom = []
        new_names = []
        new_frozenfunction = []
        new_coefficients = []
        new_types = []

        for i, (ensemble, name, isfrozen, function_type) in enumerate(
                zip(self.functions, self.function_names, self.map_frozenfunctions, self.function_types)
        ):
            function_name = name[1].split("_")[1]
            for j, (fun, fun_atom, fun_bond) in enumerate(zip(ensemble.fun, ensemble.map2atoms, ensemble.map2bonds)):
                new_functions.append(fun)
                new_map2atom.append(fun_atom)
                new_frozenfunction.append(isfrozen)
                new_names.append([fun_atom, function_name, fun_bond])
                new_types.append(function_type)
                if self.is_optimized:
                    new_coefficients.append(self.opt_params[i])
                else:
                    new_coefficients = None

        self.natoms = self.coordinates.shape[0]
        self.atom_charges = sum(self.atom_charges)
        self.map_function2center = new_map2atom
        self.map_frozenfunctions = new_frozenfunction
        self.function_names = new_names
        self.function_types = new_types
        self.functions = new_functions
        self.nfunctions = len(new_functions)
        self.opt_params = new_coefficients
        self.is_clusterized = False

    def integrate(self):
        """
        Provides the sum of the integrals of the density functions.
        Useful for debugging.

        :return:
        """
        z = 0.0
        for i, (c, fun) in enumerate(zip(self.opt_params, self.functions)):
            z += c * fun.integral()

        return z

    def get_atomic_numbers(self):
        """

        :return: atomic numbers
        """
        return self.atomic_numbers.copy()

    def get_a2md_charges(self):
        """
        Uses a2md fitting procedure as a method to describe electronic populations.
        Its performance remains untested, this function has been developped with
        debugging purpose.

        :return:
        """
        q = np.zeros(self.natoms, dtype='float64')
        for an, fun, c in zip(self.map_function2center, self.functions, self.opt_params):
            q[an] += fun.integral() * c
        q = self.atomic_numbers - q
        return q

    def get_coordinates(self):
        """

        :return:
        """
        return self.coordinates.copy()

    def get_function_names(self):
        """

        :return:
        """
        return self.function_names

    def get_gamma(self):
        """
        Gamma is the regularization value. Useful for debugging.
        :return:
        """
        return self.regularization

    def get_integrals(self):
        """
        Calculates the integrals of each of the functions and returns
        a vector with such integral.s
        :return:
        """
        integrals = np.zeros(len(self.functions), dtype='float64')
        for i, fun in enumerate(self.functions):
            integrals[i] = fun.integral()
        return integrals

    def get_frozen_integrals(self):
        """

        :return:
        """
        integrals = []
        for i, fun in enumerate(self.functions):
            if self.map_frozenfunctions[i]:
                integrals.append(fun.integral())
        return np.array(integrals, dtype='float64')

    def get_unfrozen_integrals(self):
        """

        :return:
        """
        integrals = []
        for i, fun in enumerate(self.functions):
            if not self.map_frozenfunctions[i]:
                integrals.append(fun.integral())
        return np.array(integrals, dtype='float64')

    def get_number_functions(self):
        return len(self.functions)

    def get_number_optimizable_functions(self):
        """
        returns the number of functions that have not got frozen
        :return:
        """
        j = 0
        for i in self.map_frozenfunctions:
            if not i:
                j += 1
        return j

    def get_opt_coefficients(self):
        f = np.array(self.map_frozenfunctions)
        return self.opt_params[f == False]

    def get_parametrization(self):
        """
        Creates a list with the parameters of the function. Useful for
        model persistence.

        :return:
        :rtype: list
        """
        atom_parametrization_list = [None] * self.nfunctions

        xi_iter = 0
        for xi, fn, cc, tp in zip(self.functions, self.function_names, self.opt_params, self.function_types):

            xi_params = xi.get_params()
            tmp_input = dict(center = fn[0], support_type = fn[1], coefficient=cc, function_type=tp)
            tmp_params = dict()
            for key, item in xi_params.items() :
                tmp_params[key] = item
            tmp_input['params'] = tmp_params
            tmp_input['bond'] = fn[2]


            atom_parametrization_list[xi_iter] = tmp_input
            xi_iter += 1

        return atom_parametrization_list

    def get_symbols(self):
        """
        returns the atomic element symbols associated to each atomic number
        :return:
        """
        if self.atom_labels is None:
            return [AN2ELEMENT[i] for i in self.atomic_numbers]
        else:
            return self.atom_labels

    def get_topology(self):
        """

        :return:
        """
        return self.topology.copy()

    def prepare_restricted(self):
        n_restrictions = len(self.atom_charges)
        unfrozen_map2center = [self.map_function2center[i] for i in range(self.nfunctions) if
                               not self.map_frozenfunctions[i]]
        frozen_map2center = [self.map_function2center[i] for i in range(self.nfunctions) if
                             self.map_frozenfunctions[i]]
        q = self.atom_charges.copy()
        return n_restrictions, unfrozen_map2center, frozen_map2center, q

    def prepare_unrestricted(self):
        n_restrictions = 1
        unfrozen_map2center = [0 for i in range(self.nfunctions) if not self.map_frozenfunctions[i]]
        frozen_map2center = [0 for i in range(self.nfunctions) if self.map_frozenfunctions[i]]
        q = np.array([self.atom_charges.copy().sum()])
        return n_restrictions, unfrozen_map2center, frozen_map2center, q

    def prepare_semirestricted(self):
        try:
            n_restrictions = self.nsegments
        except ValueError:
            raise RuntimeError("there are no defined segments")
        tmp_unfrozen_map2center = [self.map_function2center[i] for i in range(self.nfunctions) if
                                   not self.map_frozenfunctions[i]]
        tmp_frozen_map2center = [self.map_function2center[i] for i in range(self.nfunctions) if
                                 self.map_frozenfunctions[i]]
        unfrozen_map2center = [self.segments[i] for i in tmp_unfrozen_map2center]
        frozen_map2center = [self.segments[i] for i in tmp_frozen_map2center]
        q = [0.0] * self.nsegments
        for i, atom_q in zip(self.segments, self.atom_charges):
            q[i] += atom_q
        return n_restrictions, unfrozen_map2center, frozen_map2center, q

    def optimize(
            self, training_coordinates, training_density, optimization_mode='restricted'
    ):
        """
        sets the coefficients associated to each support function by
        minimizing the lagrangian function:
        L = (\rho - \rho')^2 - \sum \lambda_{i} (\int \rho_{i}' dV - q_{i})
        where :
            - rho : ab-initio density
            - rho' : aAMD density
            - q_{i} : charge of atom i
        :param training_coordinates: coordinates of nuclei atoms
        :param training_density: density values
        :param optimization_mode: either restricted, unrestricted or semirestricted
        :type training_coordinates: np.ndarray
        :type training_density: np.ndarray
        :type optimization_mode: str
        :return: coefficients
        :rtype: np.ndarray
        """

        x = training_coordinates.copy()
        r = training_density.copy()

        frozen_ensemble = [self.functions[i] for i in range(self.nfunctions) if self.map_frozenfunctions[i]]
        unfrozen_ensemble = [self.functions[i] for i in range(self.nfunctions) if not self.map_frozenfunctions[i]]

        n_coefficients = len(unfrozen_ensemble)

        if optimization_mode == 'restricted':
            n_restrictions, unfrozen_map2center, frozen_map2center, q = self.prepare_restricted()
        elif optimization_mode == 'unrestricted':
            n_restrictions, unfrozen_map2center, frozen_map2center, q = self.prepare_unrestricted()
        elif optimization_mode == 'semirestricted':
            n_restrictions, unfrozen_map2center, frozen_map2center, q = self.prepare_semirestricted()
        else:
            raise IOError("optimization mode must be either restricted, unrestricted or semirestricted")

        n_training = training_coordinates.shape[0]
        w = np.ones(n_training, dtype='float64') / n_training

        for i, frozen_fun in enumerate(frozen_ensemble):
            j = frozen_map2center[i]
            q[j] -= frozen_fun.integral()
        frozen_ensemble = SupportEnsemble(functions=frozen_ensemble, name='frozen_set')
        r = r - frozen_ensemble.eval(x)

        d = np.zeros((n_coefficients, n_training), dtype='float64')
        p = np.zeros((n_coefficients + n_restrictions), dtype='float64')
        b = np.zeros((n_coefficients + n_restrictions, n_coefficients + n_restrictions), dtype='float64')

        for i, xi in enumerate(unfrozen_ensemble):
            d[i, :] = xi.eval(x)

        effective_gamma = self.regularization / n_coefficients
        b[:n_coefficients, :n_coefficients] = 2*(d * w).dot(d.T)
        b[:n_coefficients, :n_coefficients] = b[:n_coefficients, :n_coefficients] + \
            2 * effective_gamma * np.identity(n_coefficients, dtype='float64')

        for i, xi in enumerate(unfrozen_ensemble):
            j = unfrozen_map2center[i]
            b[n_coefficients + j, i] += xi.integral()
            b[i, n_coefficients + j] += xi.integral()

        p[:n_coefficients] = 2*(d * w).dot(r)
        p[n_coefficients:] = q

        c = np.linalg.solve(b, p)[:n_coefficients]
        mask = np.array(self.map_frozenfunctions)
        self.opt_params = np.ones(self.nfunctions, dtype='float64')
        self.opt_params[mask == False] = c
        integral = 0.0
        for i, (c, xi) in enumerate(zip(self.opt_params, self.functions)):
            integral += c * xi.integral()
        self.is_optimized = True
        return self.opt_params.copy()

    def parametrize(self, param_dict=None):
        """
        parametrize
        ---
        Assigns the different kinds of functions to each atom and bond by reading the parameters
        dictionary (param_dict). See guidelines to create a parameters dictionary.

        It is possible to use the dictionaries stored in the a2md instance. For
        instance:

            model.parametrize(param_dict=model.parametrization_spherical)

        Default is parametrization_default.

        :param param_dict: Dictionary including parameters in the field "_MODEL"
        :return:
        """

        if param_dict is None:
            param_dict = self.parametrization_default

        try:
            model = param_dict['_MODEL']
        except KeyError:
            raise IOError('missing the _MODEL entry on the parametrization dict')

        self.map_function2center = []
        self.functions = []
        self.function_names = []
        self.map_frozenfunctions = []
        self.function_types = []
        self.opt_params = []

        symbols = [i for i in self.get_symbols()]

        for i, (element, atom_coords, atom_charge) in enumerate(
                zip(symbols, self.coordinates, self.atom_charges)
        ):
            try:
                model_element = model[element]
            except KeyError:
                raise RuntimeError("the requested element is not present in the parameters")
            for j, fun in enumerate(model_element):

                if fun['_CONNECT'] == '_NONE':
                    self.map_function2center.append(i)
                    self.function_names.append((i, fun['_NAME'], None))
                    self.map_frozenfunctions.append(fun['_FROZEN'])
                    self.opt_params.append(fun['_COEFFICIENT'])
                    self.function_types.append(fun['_TYPE'])
                    ppp = copy.copy(fun['_PARAMS'])
                    ppp['coordinates'] = self.coordinates[i, :]
                    funtype = SUPPORT_TYPE[fun['_TYPE']]
                    self.functions.append(funtype(ppp))

                elif fun['_CONNECT'] == '_TOPO':

                    for bond in self.topology[i]:
                        self.map_function2center.append(i)
                        self.function_names.append((i, fun['_NAME'], bond))
                        self.map_frozenfunctions.append(fun['_FROZEN'])
                        self.opt_params.append(fun['_COEFFICIENT'])
                        self.function_types.append(fun['_TYPE'])
                        bonding_axis = self.coordinates[bond, :] - self.coordinates[i, :]
                        ppp = copy.copy(fun['_PARAMS'])
                        ppp['coordinates'] = self.coordinates[i, :]
                        funtype = SUPPORT_TYPE[fun['_TYPE']]
                        self.functions.append(funtype(ppp))
                        self.functions[-1].set_reference_frame(bonding_axis)


        self.nfunctions = len(self.functions)
        self.opt_params = np.array(self.opt_params, dtype='float64')
        return True

    def read(self, params):
        """

        This method allows to reparametrize from previous runs

        :param params: a list of parameters for support functions
        :type params: list
        :return:
        """
        if self.functions is not None:
            if not len(self.functions) == 0 :
                self.log('existing support functions will be erased')

        function_list = []
        opt_params = []
        map_fun2center = []
        map_funfrozen = []
        function_names = []
        function_types = []
        for xi_iter in range(len(params)) :
            ppp = params[xi_iter]
            map_fun2center.append(ppp['center'])
            opt_params.append(ppp['coefficient'])
            try:
                map_funfrozen.append(ppp['frozen'])
            except KeyError:
                map_funfrozen.append(False)
            input_dict = copy.copy(ppp['params'])
            input_dict['coordinates'] = self.coordinates[ppp['center'], :]
            if ppp['bond'] is None:
                support_function = SUPPORT_TYPE['_SPHERIC']
                function_types.append('_SPHERIC')  # Hardcoded
            else:
                try:
                    support_function = SUPPORT_TYPE[ppp['function_type']]
                    function_types.append(ppp['function_type'])
                except KeyError:
                    support_function = SUPPORT_TYPE['_GAUSSIAN']  # for the sake of retrocompatibility
                    function_types.append('_GAUSSIAN')
            function_list.append(support_function(input_dict))
            if not ppp['bond'] is None :

                bonding_axis = self.coordinates[ppp['bond'], :] - self.coordinates[ppp['center'], :]
                function_list[-1].set_reference_frame(bonding_axis)
                function_names.append((ppp['center'], ppp['support_type'], ppp['bond']))
            else:
                function_names.append((ppp['center'], ppp['support_type'], None))
        opt_params = np.array(opt_params, dtype='float64')
        self.opt_params = np.array(opt_params)
        self.functions = function_list
        self.map_function2center = map_fun2center
        self.map_frozenfunctions = map_funfrozen
        self.function_names = function_names
        self.function_types = function_types
        self.nfunctions = len(opt_params)

    def set_regularization_constant(self, gamma):
        """

        :param gamma:
        :return:
        """
        self.regularization = gamma

    def set_opt_coefficients(self, c):
        f = np.array(self.map_frozenfunctions)
        self.opt_params[f == False] = c

    def use_atomic_number_as_charge(self):
        """

        :return:
        """
        self.atom_charges = self.atomic_numbers.astype('float64')

    def modify_charge_by_segment(self, charge_str):
        """
        Useful to parametrize polymers

        :param charge_str:
        :return:
        """
        if len(charge_str) > self.nsegments:
            raise IOError("charge str must be equal to the number of segments")

        charges = self.atom_charges.copy()
        segment_idx = self.segments.copy()
        for i, character in enumerate(charge_str):
            if character == 'n': segment_charge = 0.0
            elif character == '+' : segment_charge = -1.0
            elif character == '-' : segment_charge = 1.0
            else:
                raise IOError("unrecognized character!")
            len_segment = len([j for j in segment_idx if j == i])
            segment_charge /= len_segment
            for j, segment in enumerate(self.segments):
                if segment == i:
                    charges[j] += segment_charge
        self.atom_charges = charges


class Polymer(Molecule):
    topology_source = PDB_PROTEIN_TOPOLOGY
    atom_types_source = PDB_PROTEIN_TYPES
    residue_charge_source = PDB_PROTEIN_CHARGES
    residue_charge_type_source = PDB_PROTEIN_TYPE_CHARGES

    def __init__(
        self,
        coordinates, atomic_numbers, atom_labels, topology, charge, atom_residx,
        sequence, residue_names, residue_idx,  parameters=None, verbose=False
    ):
        """
        A2MD.models.Polymer
        ---
        Similar to Molecule, it contains some extra tricks to allow parametrization of polymer molecules
        as proteins.
        :param coordinates:
        :param atomic_numbers:
        :param atom_labels:
        :param topology:
        :param charge:
        :param atom_residx:
        :param sequence:
        :param residue_names:
        :param residue_idx:
        :param parameters:
        :param verbose:
        """

        Molecule.__init__(
            self,
            coordinates=coordinates, atomic_numbers=atomic_numbers, topology=topology, charge=charge,
            parameters=parameters, verbose=verbose, atom_labels=atom_labels
        )

        self.sequence = sequence
        self.residue_names = residue_names
        self.residue_idx = residue_idx
        self.atom_residx = atom_residx


    def parametrize_as_polymer(self, param_dict=None):
        """
        parametrize as polymer
        ---
        this function takes a polymer contituted by the jonction of different residues, and parametrizes
        each of these units using a separated collection of parameters.


        :param param_dict:
        :return:
        """

        parameters = []
        #
        missing_keys = []
        for i in range(self.natoms):

            atom_name = self.atom_labels[i]
            atomic_number = self.atomic_numbers[i]
            atom_resname = self.residue_names[self.atom_residx[i] - 1]

            try:

                iso_funs = param_dict[atom_resname][atom_name]['_ISO']
                aniso_funs = param_dict[atom_resname][atom_name]['_ANISO']

            except KeyError:

                missing_keys.append((atom_resname, atom_name))
                current_label = AN2ELEMENT[atomic_number]
                for fun in self.parametrization_default['_MODEL'][current_label]:
                    if fun['_CONNECT'] == '_NONE':
                        parameters.append(
                            dict(
                                coefficient=1.0,
                                params=fun['_PARAMS'].copy(),
                                frozen=fun['_FROZEN'],
                                bond=None,
                                support_type=fun['_NAME'],
                                center=i
                            )
                        )
                continue


            for j, fun in enumerate(iso_funs):
                current_fun = fun.copy()
                current_fun['center'] = i
                parameters.append(current_fun)

            for k in self.topology[i]:
                bond_name = self.atom_labels[k]
                try:
                    bond_fun_list = aniso_funs[bond_name]
                except KeyError:
                    print("WARNING: could not find a parameter for {:s}|{:s}|{:s}".format(
                        atom_resname, atom_name, bond_name
                    ))
                    continue

                for fun in bond_fun_list:
                    current_fun = fun.copy()
                    current_fun['center'] = i
                    current_fun['bond'] = k
                    parameters.append(current_fun)

        self.read(parameters)
        return missing_keys

    def get_residue_charge(self, idx, kind='total'):
        """

        :return:
        """
        if kind not in ['partial', 'total']:
            raise IOError("kind must be either partial or total")
        charge = 0.0
        for i, (c, fun, center) in enumerate(zip(self.opt_params, self.functions, self.map_function2center)):
            res_idx = self.atom_residx[center]
            if res_idx == idx:
                charge += fun.integral() * c
        if kind == 'total':
            return charge
        elif kind == 'partial':
            charge = - charge
            for i, (an, res_idx) in enumerate(zip(self.atomic_numbers, self.atom_residx)):
                if res_idx == idx:
                    charge += an
            return charge

    def get_residue_functions(self, idx):
        for i, (c, fun, center, name) in enumerate(zip(self.opt_params, self.functions, self.map_function2center, self.function_names)):
            res_idx = self.atom_residx[center]
            if res_idx == idx:
                yield c, fun, center, name

class ConformerCollection(Molecule):
    """
    ConformerCollection
    ---
    Molecule instance that stores a group of conformations in its variable --conformers--.
    The method parametrize_conformer allows to change the positions of the centers of the molecule
    to one of the centers.

    The method eval conformers allows to use the different conformers at the same time for different sets of coordinates

    The method conformer optimize optimizes model coefficients considering all the possible conformations. Requires
    a set of coordinates and densities that is the same size than the number of conformations.

    """
    def __init__(
            self,
            coordinates, atomic_numbers, charge, topology,
            parameters=None, verbose=False, atom_labels=None, segments=None
    ):
        Molecule.__init__(
            self,
            coordinates=coordinates[0], atomic_numbers=atomic_numbers, topology=topology, charge=charge,
            parameters=parameters, verbose=verbose, atom_labels=atom_labels, segments=segments
        )
        self.nconformers = len(coordinates)
        self.conformers = coordinates
        self.current_conformer = 0

    def parametrize_conformer(self, i):

        if i > self.nconformers:
            raise RuntimeError("there is not conformer {:d}".format(i))
        current_coordinates = self.conformers[i]
        self.current_conformer = i

        for i, ((atom_idx, _, bond_idx), fun) in enumerate(zip(self.function_names, self.functions)):
            fun.coordinates = current_coordinates[atom_idx, :]
            if bond_idx is not None:
                bv = current_coordinates[bond_idx, :] - current_coordinates[atom_idx, :]
                bv /= np.linalg.norm(bv)
                fun.set_reference_frame(bv)

    def eval_conformers(self, x):
        p = []
        for i, x_c in enumerate(x):
            self.parametrize_conformer(i)
            p.append(np.zeros(x_c.shape[0], dtype="float64"))
            for j, (c, fun) in enumerate(zip(self.opt_params, self.functions)):
                p[-1] += c * fun.eval(x_c)

        return p

    def conformer_optimize(
            self, training_coordinates, training_densities, optimization_mode='restricted'
    ):
        if len(training_coordinates) != len(training_densities):
            raise RuntimeError(
                "the number of conformer densities does not match the number of training coordinates")

        if len(training_densities) != self.nconformers:
            raise RuntimeError("the number of conformer does not match the number of training coordinates")

        unfrozen_ensemble = [self.functions[i] for i in range(self.nfunctions) if not self.map_frozenfunctions[i]]
        frozen_ensemble = [self.functions[i] for i in range(self.nfunctions) if self.map_frozenfunctions[i]]
        n_coefficients = len(unfrozen_ensemble)

        if optimization_mode == 'restricted':
            n_restrictions, unfrozen_map2center, frozen_map2center, q = self.prepare_restricted()
        elif optimization_mode == 'unrestricted':
            n_restrictions, unfrozen_map2center, frozen_map2center, q = self.prepare_unrestricted()
        elif optimization_mode == 'semirestricted':
            n_restrictions, unfrozen_map2center, frozen_map2center, q = self.prepare_semirestricted()
        else:
            raise IOError("optimization mode must be either restricted, unrestricted or semirestricted")

        p = np.zeros((n_coefficients + n_restrictions), dtype='float64')
        b = np.zeros((n_coefficients + n_restrictions, n_coefficients + n_restrictions), dtype='float64')

        for i, frozen_fun in enumerate(frozen_ensemble):
            j = frozen_map2center[i]
            q[j] -= frozen_fun.integral()

        d = []
        r = []
        for i, (ctc, ctd) in enumerate(zip(training_coordinates, training_densities)):

            self.parametrize_conformer(i)

            frozen_ensemble = [self.functions[i] for i in range(self.nfunctions) if self.map_frozenfunctions[i]]
            unfrozen_ensemble = [self.functions[i] for i in range(self.nfunctions) if
                                 not self.map_frozenfunctions[i]]

            x = ctc.copy()
            current_r = ctd.copy()

            n_training = current_r.shape[0]

            frozen_ensemble = SupportEnsemble(functions=frozen_ensemble, name='frozen_set')
            current_r = current_r - frozen_ensemble.eval(x)

            current_d = np.zeros((n_coefficients, n_training), dtype='float64')

            for j, xi in enumerate(unfrozen_ensemble):
                current_d[j, :] = xi.eval(x)

            d.append(current_d)
            r.append(current_r)

        d = np.concatenate(d, axis=1)
        r = np.concatenate(r, axis=0)

        w = np.ones(d.shape[1], dtype='float64') / d.shape[1]
        effective_gamma = self.regularization / n_coefficients
        b[:n_coefficients, :n_coefficients] = 2 * (d * w).dot(d.T)
        b[:n_coefficients, :n_coefficients] = b[:n_coefficients, :n_coefficients] + \
                                              2 * effective_gamma * np.identity(n_coefficients, dtype='float64')


        for i, xi in enumerate(unfrozen_ensemble):
            j = unfrozen_map2center[i]
            b[n_coefficients + j, i] += xi.integral()
            b[i, n_coefficients + j] += xi.integral()

        p[:n_coefficients] = 2 * (d * w).dot(r)
        p[n_coefficients:] = q

        c = np.linalg.solve(b, p)[:n_coefficients]
        mask = np.array(self.map_frozenfunctions)
        self.opt_params = np.ones(self.nfunctions, dtype='float64')
        self.opt_params[mask == False] = c
        integral = 0.0
        for i, (c, xi) in enumerate(zip(self.opt_params, self.functions)):
            integral += c * xi.integral()
        self.is_optimized = True
        return self.opt_params.copy()
