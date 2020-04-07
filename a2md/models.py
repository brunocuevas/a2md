import numpy as np
import copy
from a2md.support import SupportRadial
from a2md.support import SupportAngular
from a2md.support import SupportEnsemble
from a2md import TOPO_RESTRICTED_PARAMS
from a2md.baseclass import A2MDBaseClass
from a2mdio.molecules import MolRepresentation, PDB
from a2md.utils import convert_connectivity_tree_to_pairs
from a2mdio import PDB_PROTEIN_TYPE_CHARGES, PDB_PROTEIN_CHARGES, PDB_PROTEIN_TYPES, PDB_PROTEIN_TOPOLOGY
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

SUPPORT_TYPE = dict(
    radial  = lambda args : SupportRadial(**args),
    angular = lambda args : SupportAngular(**args)
)

def a2md_from_mol(mol : MolRepresentation):
    """
    returns the density model associated to a Mol2
    :param mol:
    :type mol: Mol2
    :return:
    :rtype: Molecule
    """
    an = mol.get_atomic_numbers()
    coords = mol.get_coordinates(units="au")
    charge = mol.get_absolute_charges()
    topology = mol.get_bonds()
    topo_array = []
    for i in range(mol.get_number_atoms()): topo_array.append([])
    for i in range(mol.get_number_bonds()):
        begin = int(topology[i, 0])
        end = int(topology[i, 1])
        topo_array[begin - 1].append(end - 1)
        topo_array[end - 1].append(begin - 1)
    return Molecule(
        coordinates=coords, atomic_numbers=an, charge=charge, topology=topo_array
    )

def polymer_from_pdb(mol : PDB, chain : str ):
    """
    returns the polymer density model associated to a Mol2
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
    def __init__(
            self, coordinates, atomic_numbers, charge, topology, parameters=None, verbose=False,
            atom_labels=None
    ):
        """

        :param coordinates: nuclei coordinates
        :param charge: atoms in molecule charges
        :param topology: list of linked nuclei
        :param parameters: parameters for bonding
        :param verbose:
        :param atomic_numbers
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
        self.regularization = 0.0001
        self.is_clusterized = False
        self.is_optimized = False

        if parameters is not None:
            self.read(parameters)

    def clustering(self, clusterizer):
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

        expanded_fns = []

        for i in range(len(function_names)):

            bond_ = function_names[i][2]
            name_ = function_names[i][1]

            expanded_fns.append(
                dict(
                    type=name_,
                    center=function_names[i][0],
                    bond=bond_,
                    name=name_
                )
            )
            del bond_, name_

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

        # Part 1, cluster isotropic functions

        for cluster_idx, cluster in enumerate(sa):

            # Part 1.1. Find types of functions

            current_cluster_function_types = []
            for atom_idx in cluster:
                for i, fun_name in enumerate(expanded_fns):
                    c1 = fun_name['center'] == atom_idx
                    c2 = fun_name['type'] not in current_cluster_function_types
                    c3 = fun_name['bond'] is None
                    if all([c1, c2, c3]):
                        current_cluster_function_types.append(fun_name['type'])
                    del c1, c2, c3

            # Part 1.1. Cluster isotropic functions

            for fun_type in current_cluster_function_types:
                is_frozen = []
                current_cluster_functions = []
                current_map2atoms = []
                current_map2bonds = []
                for atom_idx in cluster:
                    for i, fun_name in enumerate(expanded_fns):
                        c1 = fun_name['center'] == atom_idx
                        c2 = fun_name['type'] == fun_type
                        c3 = fun_name['bond'] is None
                        if all([c1, c2, c3]):
                            current_cluster_functions.append(functions[i])
                            is_frozen.append(map_frozenfunction[i])
                            current_map2atoms.append(atom_idx)
                            current_map2bonds.append(None)
                        del c1, c2, c3


                new_functions.append(
                    SupportEnsemble(
                        functions=current_cluster_functions,
                        name="ensemble_c{:d}_f{:s}".format(cluster_idx, fun_type),
                        map2atoms=current_map2atoms, map2bonds=current_map2bonds
                    )
                )

                if not all(is_frozen) and any(is_frozen):
                    self.log("there are hetereogeneous functions in the current group")

                new_frozenfunction.append(all(is_frozen))
                new_names.append([cluster_idx, "ENS_{:s}".format(fun_type)])
                new_map2atom.append(cluster_idx)


        # Part 2, cluster  anisotropic functions
        # bond cluster must contain the idx of the clusters of sa
        # so we can map atoms involved in this bond cluster.
        # NOTE: this operation is assymetric!!!
        for bond_cluster in sb:

            cluster_0 = sa[bond_cluster[0]]
            cluster_1 = sa[bond_cluster[1]]
            # only functions that come from cluster 0 are required
            current_cluster_function_types = []
            for atom_idx in cluster_0:

                for i, fun_name in enumerate(expanded_fns):
                    c1 = fun_name['center'] == atom_idx
                    c2 = fun_name['type'] not in current_cluster_function_types
                    c3 = fun_name['bond'] in cluster_1
                    if all([c1, c2, c3]):
                        current_cluster_function_types.append(fun_name['type'])
                    del c1, c2, c3


            for fun_type in current_cluster_function_types:
                current_cluster_functions = []
                current_map2atoms = []
                current_map2bonds = []
                for atom_idx in cluster_0:
                    for i, fun_name in enumerate(expanded_fns):
                        c1 = fun_name['center'] == atom_idx
                        c2 = fun_name['type'] == fun_type
                        c3 = fun_name['bond'] in cluster_1
                        if all([c1, c2, c3]):
                            current_cluster_functions.append(functions[i])
                            current_map2atoms.append(atom_idx)
                            current_map2bonds.append(fun_name['bond'])
                        del c1, c2, c3

                new_functions.append(
                    SupportEnsemble(
                        functions=current_cluster_functions,
                        name="ensemble_c{:d}_f{:s}_b{:d}".format(bond_cluster[0], fun_type, bond_cluster[1]),
                        map2atoms=current_map2atoms, map2bonds=current_map2bonds
                    )
                )
                new_names.append([bond_cluster[0], "ENS_{:s}".format(fun_type), bond_cluster[1]])
                new_map2atom.append(bond_cluster[0])
                new_frozenfunction.append(False)

        self.natoms = len(sa)
        self.atom_charges = new_atom_charges
        self.map_function2center = new_map2atom
        self.map_frozenfunctions = new_frozenfunction
        self.function_names = new_names
        self.functions = new_functions
        self.opt_params = np.zeros(len(sa), dtype="float64")
        self.nfunctions = len(new_functions)
        self.is_clusterized = True

    def eval(self, x):
        """
        Evaluation of aAMD at defined coordinates

        :param x: cartessian coordinates
        :type x: np.ndarray
        :return: density values
        :type: np.ndarray
        """
        d = np.zeros(x.shape[0])
        for sup, c in zip(self.functions, self.opt_params) :
            d += c * sup.eval(x)
        return d

    def eval_core(self, x):
        d = np.zeros(x.shape[0])
        for sup, c, f in zip(self.functions, self.opt_params, self.map_frozenfunctions):
            if f: d += c * sup.eval(x)
        return d

    def eval_volume(self, spacing, resolution):
        """

        :param spacing:
        :param resolution:
        :return:
        """
        from a2mdio.qm import ElectronDensity
        minx, miny, minz = (self.coordinates - spacing).min(axis=0)
        maxx, maxy, maxz = (self.coordinates + spacing).max(axis=0)

        xx = np.arange(minx, maxx, resolution)
        yy = np.arange(miny, maxy, resolution)
        zz = np.arange(minz, maxz, resolution)

        dx = np.zeros((xx.size, yy.size, zz.size))

        for ix in range(xx.size):
            r = np.zeros((zz.size, 3), dtype='float64')
            r[:, 0] = xx[ix]
            for iy in range(yy.size):
                r[:, 1] = yy[iy]
                r[:, 2] = zz[:]
                dx[ix, iy, :] = self.eval(r)

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

        for i, (ensemble, name, isfrozen) in enumerate(
                zip(self.functions, self.function_names, self.map_frozenfunctions)
        ):
            function_type = name[1].split("_")[1]
            for j, (fun, fun_atom, fun_bond) in enumerate(zip(ensemble.fun, ensemble.map2atoms, ensemble.map2bonds)):
                new_functions.append(fun)
                new_map2atom.append(fun_atom)
                new_frozenfunction.append(isfrozen)
                new_names.append([fun_atom, function_type, fun_bond])
                if self.is_optimized:
                    new_coefficients.append(self.opt_params[i])
                else:
                    new_coefficients = None

        self.natoms = self.coordinates.shape[0]
        self.atom_charges = sum(self.atom_charges)
        self.map_function2center = new_map2atom
        self.map_frozenfunctions = new_frozenfunction
        self.function_names = new_names
        self.functions = new_functions
        self.nfunctions = len(new_functions)
        self.opt_params = new_coefficients
        self.is_clusterized = False

    def integrate(self):
        """

        :return:
        """
        z = 0.0
        for i, (c, fun) in enumerate(zip(self.opt_params, self.functions)):

            z += c * fun.integral()

        return z

    def get_atomic_numbers(self):
        """

        :return:
        """
        return self.atomic_numbers.copy()

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

        :return:
        """
        return self.regularization

    def get_parametrization(self):
        """
        This method allows to create a list of function parameters

        :return:
        :rtype: list
        """
        atom_parametrization_list = [None] * self.nfunctions

        xi_iter = 0
        for xi, fn, cc in zip(self.functions, self.function_names, self.opt_params) :

            xi_params = xi.get_params()
            tmp_input = dict(center = fn[0], support_type = fn[1], coefficient=cc)
            tmp_params = dict()
            for key, item in xi_params.items() :
                tmp_params[key] = item
            tmp_input['params'] = tmp_params
            tmp_input['bond'] = fn[2]


            atom_parametrization_list[xi_iter] = tmp_input
            xi_iter += 1

        return atom_parametrization_list

    def get_symbols(self):
        if self.atom_labels is None:
            return [AN2ELEMENT[i] for i in self.atomic_numbers]
        else:
            return self.atom_labels

    def get_topology(self):
        """

        :return:
        """
        return self.topology.copy()

    def optimize(
            self, training_coordinates, training_density
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
        :type training_coordinates: np.ndarray
        :type training_density: np.ndarray
        :return: coefficients
        :rtype: np.ndarray
        """

        x = training_coordinates.copy()
        r = training_density.copy()

        frozen_ensemble = [self.functions[i] for i in range(self.nfunctions) if self.map_frozenfunctions[i]]
        frozen_map2center = [self.map_function2center[i] for i in range(self.nfunctions) if self.map_frozenfunctions[i]]
        unfrozen_ensemble = [self.functions[i] for i in range(self.nfunctions) if not self.map_frozenfunctions[i]]
        unfrozen_map2center = [self.map_function2center[i] for i in range(self.nfunctions) if
                               not self.map_frozenfunctions[i]]

        n_coefficients = len(unfrozen_ensemble)
        n_restrictions = len(self.atom_charges)
        n_training = training_coordinates.shape[0]
        w = np.ones(n_training, dtype='float64') / n_training

        q = self.atom_charges.copy()
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
        This function aims to replace the preprocessor by allowing A2MD.Molecule
        to parametrize itself. Just need to provide a set of parameters.

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
        self.opt_params = []

        symbols = [i for i in self.get_symbols()]

        for i, (element, atom_coords, atom_charge) in enumerate(
                zip(symbols, self.coordinates, self.atom_charges)
        ):
            for j, fun in enumerate(model[element]):

                if fun['_CONNECT'] == '_NONE':
                    self.map_function2center.append(i)
                    self.function_names.append((i, fun['_NAME'], None))
                    self.map_frozenfunctions.append(fun['_FROZEN'])
                    self.opt_params.append(fun['_COEFFICIENT'])
                    ppp = copy.copy(fun['_PARAMS'])
                    ppp['coordinates'] = self.coordinates[i, :]
                    self.functions.append(SUPPORT_TYPE['radial'](ppp))

                elif fun['_CONNECT'] == '_TOPO':

                    for bond in self.topology[i]:
                        self.map_function2center.append(i)
                        self.function_names.append((i, fun['_NAME'], bond))
                        self.map_frozenfunctions.append(fun['_FROZEN'])
                        self.opt_params.append(fun['_COEFFICIENT'])
                        bonding_axis = self.coordinates[bond, :] - self.coordinates[i, :]
                        ppp = copy.copy(fun['_PARAMS'])
                        ppp['coordinates'] = self.coordinates[i, :]
                        self.functions.append(SUPPORT_TYPE['angular'](ppp))
                        self.functions[-1].set_reference_frame(bonding_axis)

        self.nfunctions = len(self.functions)
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
                support_function = SUPPORT_TYPE['radial']
            else:
                support_function = SUPPORT_TYPE['angular']
            function_list.append(support_function(input_dict))
            if not ppp['bond'] is None :

                bonding_axis = self.coordinates[ppp['bond'], :] - self.coordinates[ppp['center'], :]
                function_list[-1].set_reference_frame(bonding_axis)
                self.function_names.append((ppp['center'], ppp['support_type'], ppp['bond']))
            else:
                self.function_names.append((ppp['center'], ppp['support_type']))
        opt_params = np.array(opt_params, dtype='float64')
        self.opt_params = np.array(opt_params)
        self.functions = function_list
        self.map_function2center = map_fun2center
        self.map_frozenfunctions = map_funfrozen
        self.nfunctions = len(opt_params)

    def set_regularization_constant(self, gamma):
        """

        :param gamma:
        :return:
        """
        self.regularization = gamma


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

        for i in range(self.natoms):

            atom_name = self.atom_labels[i]
            atomic_number = self.atomic_numbers[i]
            atom_resname = self.residue_names[self.atom_residx[i] - 1]

            try:

                iso_funs = param_dict[atom_resname][atom_name]['_ISO']
                aniso_funs = param_dict[atom_resname][atom_name]['_ANISO']

            except KeyError:

                print("WARNING: could not find the combination {:s} {:s} in our parameters".format(
                    atom_resname, atom_name)
                )
                print("\ta standard isotropic electron distribution will be given to this atom")
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
                    print("WARNING: could not find a parameter for ...")
                    continue

                for fun in bond_fun_list:
                    current_fun = fun.copy()
                    current_fun['center'] = i
                    current_fun['bond'] = k
                    parameters.append(current_fun)



        self.read(parameters)
