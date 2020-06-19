import numpy as np
import networkx as nx
from a2mdio.molecules import Mol2


def angle(x1, x2, x3):

    v1 = x1 - x2
    v2 = x3 - x2
    mv1 = np.linalg.norm(v1)
    mv2 = np.linalg.norm(v2)

    return np.arccos(v1.dot(v2)/(mv1 * mv2))

def dihedral(x1, x2, x3, x4):

    v1 = x2 - x1
    v2 = x3 - x2
    v3 = x4 - x3

    u1 = np.cross(v1, v2)
    u2 = np.cross(v2, v3)
    mu1 = np.linalg.norm(u1)
    mu2 = np.linalg.norm(u2)
    cosphi = u1.dot(u2) /(mu1 * mu2)
    sinphi = v2.dot(np.cross(u1, u2)) / (np.linalg.norm(u2) * mu1 * mu2)
    return np.arccos(cosphi) * np.sign(sinphi)

def distance(x1, x2):
    return np.linalg.norm(x2 - x1)


def decompose_molecule_graph(G):
    bonds = []
    angles = []
    dihedrals = []

    for i in range(G.number_of_nodes()):
        neighs = G[i]
        for j in neighs:
            pair = "|".join(sorted(["%d" % i, "%d" % j]))
            if pair not in bonds:
                bonds.append(pair)
            for k in G[j]:
                if k == i: continue
                triple = "|".join(sorted(["%d" % i, "%d" % k]))
                triple = "%s,%d" % (triple, j)
                if triple not in angles:
                    angles.append(triple)

                if G.nodes[i]["symbol"] == "H":
                    continue
                for l in G[k]:
                    if G.nodes[l]["symbol"] == "H":
                        continue

                    if l == i or l == j: continue

                    quad1 = "|".join(sorted(["%d" % i, "%d" % l]))
                    quad2 = "|".join(sorted(["%d" % j, "%d" % k]))
                    quad = "%s,%s" % (quad1, quad2)

                    if quad not in dihedrals:
                        dihedrals.append(quad)

    return sorted(bonds), sorted(angles), sorted(dihedrals)

def convert(mol : Mol2):
    G = nx.Graph()
    for i in range(mol.get_number_atoms()):
        G.add_node(i, **mol.get_atom(i))
    bonds = mol.get_bonds()
    bonds = bonds - 1
    for i in range(mol.get_number_bonds()):
        G.add_edge(bonds[i, 0], bonds[i, 1])

    bonds, angles, dihedrals = decompose_molecule_graph(G)
    dbonds = np.zeros(len(bonds), dtype='float64')
    dangles = np.zeros(len(angles), dtype='float64')
    ddihedrals = np.zeros(len(dihedrals), dtype='float64')

    for i in range(len(bonds)):
        atom1, atom2 = bonds[i].split("|")
        atom1 = int(atom1)
        atom2 = int(atom2)
        x1 = G.nodes[atom1]['coordinates']
        x2 = G.nodes[atom2]['coordinates']
        dbonds[i] = distance(x1, x2)

    for i in range(len(angles)):
        ends, center = angles[i].split(",")
        atom1, atom3 = ends.split("|")
        atom1 = int(atom1)
        atom3 = int(atom3)
        atom2 = int(center)
        x1 = G.nodes[atom1]['coordinates']
        x2 = G.nodes[atom2]['coordinates']
        x3 = G.nodes[atom3]['coordinates']
        dangles[i] = angle(x1, x2, x3)

    for i in range(len(dihedrals)):
        ends, center = dihedrals[i].split(",")
        atom1, atom4 = ends.split("|")
        atom2, atom3 = center.split("|")
        atom1 = int(atom1)
        atom3 = int(atom3)
        atom2 = int(atom2)
        atom4 = int(atom4)
        x1 = G.nodes[atom1]['coordinates']
        x2 = G.nodes[atom2]['coordinates']
        x3 = G.nodes[atom3]['coordinates']
        x4 = G.nodes[atom4]['coordinates']
        ddihedrals[i] = dihedral(x1, x2, x3, x4)

    return dbonds, dangles, ddihedrals

def convert_from_reference(coordinates, ref):
    t = np.zeros(len(ref), dtype='float64')
    for i, dhatoms in enumerate(ref):
        i1, i2, i3, i4 = dhatoms
        x1 = coordinates[i1 - 1, :]
        x2 = coordinates[i2 - 1, :]
        x3 = coordinates[i3 - 1, :]
        x4 = coordinates[i4 - 1, :]
        t[i] = dihedral(x1, x2, x3, x4)

    return t
