from a2mdio.molecules import Mol2
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from a2mdtests import gdb_test
print(gdb_test)
if __name__ == '__main__':
    methane = Mol2(file=gdb_test['gdb_000001']['path'] / gdb_test['gdb_000001']['mol2'])
    coords = methane.get_coordinates()
    res = methane.sample_probes(region='bonding', resolution=0.1)
    points = np.array([
        [-5.0, 0.0, 0.0],
        [5.0, 0.0, 0.0],
        [0.0, 5.0, 0.0],
        [0.0, -5.0, 0.0],
        [0.0, 0.0, 5.0],
        [0.0, 0.0, -5.0]
    ], dtype='float64') /3
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(res[:, 0], res[:, 1], res[:, 2])
    ax.scatter(coords[:,0], coords[:, 1], coords[:, 2], c = 'green', alpha=0.1)
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], c = 'red')
    plt.show()
