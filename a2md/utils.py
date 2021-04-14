import numpy as np
from a2md.baseclass import A2MDBaseClass
from a2md import mathfunctions
import math
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


def set_nearest_atom(points, coordinates):
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

def topology_from_bonds(bonds, natoms, nbonds):
    topology = []
    for i in range(natoms): topology.append([])
    for i in range(nbonds):
        begin = int(bonds[i, 0])
        end = int(bonds[i, 1])
        topology[begin - 1].append(end - 1)
        topology[end - 1].append(begin - 1)
    return topology

def integrate_from_dict(fun_dict):
    """
    integrate from dict
    ---
    this function allows to integrate a function without having to declare it as a support function,
    directly from the ppp.

    Important: it does not multiply by the coefficient!!

    :param fun_dict: contains all the function parameters
    :return: float
    """
    try:
        a = fun_dict['params']['A']
        b = fun_dict['params']['B']
    except KeyError:
        raise IOError("could not find the parameters of the radial part")
    try:
        p = fun_dict['params']['P']
    except KeyError:
        if fun_dict['bond'] is None: p = 0
        else: p = 1
    if fun_dict['bond'] is None:
        return 4 * np.pi * mathfunctions.generalized_exponential_integral(a, b, p)
    else:
        alpha = fun_dict['params']['alpha']
        return 2 * np.pi * mathfunctions.generalized_exponential_integral(
            a, b, p) * mathfunctions.angular_gaussian_integral(alpha)

def integrate_from_old_dict(fun):
    """

    :return:
    """
    if fun['support_type'] in ['ORCV', 'OR', 'ORV', 'ORC']:
        a = fun['params']['A3']
        b = fun['params']['B3']
        return mathfunctions.generalized_exponential_integral(a, b, 0) * 4 * math.pi
    elif fun['support_type'] == 'AG':
        g = fun['params']['G']
        u = fun['params']['U']
        alpha = fun['params']['Alpha']
        v = mathfunctions.angular_gaussian_integral(alpha)
        v *= mathfunctions.generalized_exponential_integral(u, g, P=1)
        v *= 2 * math.pi
        return v


class ClusterTool(A2MDBaseClass):
    def __init__(self, name, verbose=False):
        A2MDBaseClass.__init__(self, name="cluter tool / {:s}".format(name), verbose=verbose)
        self.cluster_method = None
    def cluster(self, labels, topology, coordinates):
        return self.cluster_method(labels, topology, coordinates)

class RBFSymmetryCluster(ClusterTool):
    centers = [1.0, 2.7, 4.4, 6.1]
    vars = [1.0, 1.0, 1.0, 1.0]
    heights = [1.0, 0.5, 0.25, 0.1]
    n = 4
    def __init__(self, verbose=False):

        ClusterTool.__init__(self, name="Radial Gaussian Basis Function cluster", verbose=verbose)
        self.cluster_method = self.rbf_clustering

    def rbf(self, coordinates1, coordinates2):
        """

        :param coordinates1:
        :param coordinates2:
        :return:
        """
        d = np.zeros((coordinates1.shape[0], coordinates2.shape[0]), dtype="float64")
        for i in range(coordinates1.shape[0]):

            r = coordinates2 - coordinates1[i, :]
            d[i, :] = np.linalg.norm(r, axis=1)

        z = np.zeros((coordinates1.shape[0], self.n), dtype="float64")
        for i in range(self.n):
            z[:, i] = (np.exp(-((d - self.centers[i])**2)/self.vars[i]) * self.heights[i]).sum(1)

        return z


    def rbf_clustering(self, labels, topology, coordinates):
        """

        :param labels:
        :param topology:
        :param coordinates:
        :return:
        """
        from sklearn.cluster import AgglomerativeClustering
        # order by elements
        natoms = len(labels)
        present_elements = np.unique(labels)
        sa = []
        sb = []
        for i, element in enumerate(present_elements):

            index = []
            current_coordinates = []

            for j in range(natoms):
                if labels[j] == element:

                    current_coordinates.append(coordinates[j, :])
                    index.append(j)

            if len(current_coordinates) == 1:
                sa.append([index[0]])
                continue

            current_coordinates = np.array(current_coordinates, dtype="float64")
            rbf_values = self.rbf(current_coordinates, coordinates)


            atom_cluster = AgglomerativeClustering(
                n_clusters=None, affinity="euclidean", linkage="ward", distance_threshold=1e-3
            ).fit(rbf_values)

            cluster_labels = atom_cluster.labels_
            all_clusters = np.unique(cluster_labels)
            for l in all_clusters:
                sa.append([])
                filtered_labels = cluster_labels == l
                for j, f in enumerate(filtered_labels):
                    if f:
                        sa[-1].append(index[j])

        for cidx1, cluster_1 in enumerate(sa):
            for cidx2, cluster_2 in enumerate(sa):

                for i, b in enumerate(topology):
                    b0 = b[0]
                    b1 = b[1]
                    if b0 in cluster_1 and b1 in cluster_2:
                        joint_id = "|".join([str(cidx1), str(cidx2)])
                        if joint_id not in sb:
                            sb.append(joint_id)
                    if b0 in cluster_2 and b1 in cluster_1:
                        joint_id = "|".join([str(cidx1), str(cidx2)])
                        if joint_id not in sb:
                            sb.append(joint_id)
        sb = [i.split("|") for i in sb]
        sb = [[int(i[0]), int(i[1])] for i in sb]
        return sa, sb

def maptoconstraints(cp, x, q):

    a = np.identity(cp.size + 1)
    b = np.zeros(cp.size + 1)
    b[:-1] = 2*cp
    b[-1] = q
    a[-1, -1] = 0.0
    a = a * 2
    a[-1, :-1] = x
    a[:-1, -1] = x
    c = np.linalg.solve(a, b)
    return c[:-1]

def project(g, x):
    t1 = g.dot(x)
    mod = np.linalg.norm(x)
    return g - (t1 * x / (mod ** 2))