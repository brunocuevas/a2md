import numpy as np
from a2mdio.molecules import Mol2
from a2mdio.qm import WaveFunction, GaussianLog


def sample_box(
        coordinates,
        spacing,
        rotate=True,
        npoints=1000
):
    """

    Sample random points in a given box. The axis of the box
    can be rotated to the maximum covariance axis of the system
    to provide a more efficient sampling

    :param coordinates:
    :param spacing:
    :param rotate:
    :param npoints:
    :return:
    """
    p = None

    if not isinstance(coordinates, np.ndarray):
        raise TypeError("use np.ndarray for coordinates")
    x = np.copy(coordinates)

    if rotate:
        c = x.T.dot(x)
        eig, p = np.linalg.eig(c)
        x = x.dot(np.linalg.inv(p))
    min_x = x[:, 0].min() - spacing
    min_y = x[:, 1].min() - spacing
    min_z = x[:, 2].min() - spacing

    max_x = x[:, 0].max() + spacing
    max_y = x[:, 1].max() + spacing
    max_z = x[:, 2].max() + spacing

    u = np.random.rand(npoints, 3)
    u[:, 0] = (u[:, 0] * (max_x - min_x)) + min_x
    u[:, 1] = (u[:, 1] * (max_y - min_y)) + min_y
    u[:, 2] = (u[:, 2] * (max_z - min_z)) + min_z

    if rotate:
        u = u.dot(p)
    return u


def merge_sources(mol2, gaussian_output, wavefunction):
    """

    Merges data coming from different sources to ease the work.

    :param mol2:
    :type mol2: Mol2
    :param gaussian_output:
    :type gaussian_output: GaussianLog
    :param wavefunction:
    :type wavefunction: WaveFunction
    :return: Mol2 object with charges obtained from gaussian output and
    coordinates obtained from wavefunction file
    """
    mol2.set_coordinates(
        wavefunction.get_coordinates(), units='au'
    )
    mol2.set_charges(gaussian_output.get_charges(kind='partial'))
    mol2.change_units(units='angstrom')
    return mol2


def get_normal_modes(g09_output, natoms):
    """
    Reads a string containing the normal modes of a frequency calculation, and returns
    a numpy tensor contaning all the modes, and their frequencies and force constants.

    :param g09_output: string containing the normal modes (g09 calculation)
    :type g09_output: str
    :param natoms: number of atoms involved
    :type natoms: int
    :return: frequencies, forces and normal modes vectors
    """

    nmodes = (natoms * 3) - 6
    modes = np.zeros((nmodes, natoms, 3), dtype='float64')
    frequencies = np.zeros(nmodes, dtype='float64')
    forces = np.zeros(nmodes, dtype='float64')
    nrows = nmodes // 3 + (nmodes % 3 != 0)
    nlines = 7 + natoms

    for mode_row in range(nrows):
        row = g09_output[mode_row*nlines:(mode_row+1)*nlines]
        index = [int(i) - 1 for i in row[0].split()]
        freq_row = row[2].split()
        force_row = row[4].split()
        frequencies[index] = [float(i) for i in freq_row[2:]]
        forces[index] = [float(i) for i in  force_row[3 :]]
        vector_rows = row[7:]
        for i in range(natoms):
            vector = vector_rows[i].split()
            modes[index[0], i, :] = [float(j) for j in vector[2:5]]
            try:
                modes[index[1], i, :] = [float(j) for j in vector[5:8]]
            except IndexError:
                pass
            try:
                modes[index[2], i, :] = [float(j) for j in vector[8:]]
            except IndexError:
                pass

    return frequencies, forces, modes

def write_xyz(filename, coords, labels):

    with open(filename, 'w') as f:
        f.write('{:d}\n'.format(len(labels)))
        f.write('\n')
        for i in range(len(labels)):
            f.write(
                "{:s} {:12.4f} {:12.4f} {:12.4f}\n".format(
                    labels[i], coords[i, 0], coords[i, 1], coords[i, 2]
                )
            )