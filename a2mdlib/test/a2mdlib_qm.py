from a2mdlib.qm import Wavefunction, GaussianLog
from a2mdtests import gdb_test
import numpy as np

if __name__ == '__main__':
    print("testing to read a wfn")
    wfn = Wavefunction(
        verbose=True,
        file=gdb_test['gdb_000001']['path'] / gdb_test['gdb_000001']['wfn'],
        batch_size=1000
    )
    print("\tdone!")
    print("testing to calculate density upon wave function")
    random_coords = np.random.rand(1000, 3)
    wfn.calculate_density(random_coords)
    print("\tdone!")
    print("testing result from Cdens against Wavefunction")
    sampled = np.loadtxt(
        gdb_test['gdb_000001']['path'] / list(gdb_test['gdb_000001']['surfaces'].keys())[0],
        delimiter=',', skiprows=1
    )
    ro = wfn.calculate_density(sampled[:, :3])
    rmsd = np.sqrt((sampled[:, 3] - ro).dot((sampled[:, 3] - ro)))/sampled.shape[0]

    print("\tRMSD = {:8.4e}".format(rmsd))
    print("testing to obtain atomic labels and atomic numbers")
    for i, (al, an) in enumerate(zip(wfn.get_atom_labels(), wfn.get_atomic_numbers())):
        print("{:>8d},{:>3s},{:>3d}".format(i, al, an))
    print("testing to obtain coordinates")
    coords = wfn.get_coordinates()
    for i in range(coords.shape[0]):
        print(
            "{:>d},{:<8.4f},{:<8.4f},{:<8.4f}".format(
                i, coords[i, 0], coords[i, 1], coords[i, 2]
            )
        )
    print("-----------------------------------------------------")
    print("Wave function test - DONE")
    print("testing to read a gaussian log")
    gl = GaussianLog(
        verbose=False,
        file=gdb_test['gdb_000001']['path'] / gdb_test['gdb_000001']['out']
    )
    print("done")
    for i, (pc, tc) in enumerate(
            zip(gl.get_charges(kind='partial'), gl.get_charges(kind='total'))
    ):
        print("{:>3d}, {:<8.4f}, {:<8.4f}".format(i, pc, tc))


    from a2mdlib.utils import write_xyz

    coords = coords / 1.8897
    ans = wfn.get_atom_labels()

    write_xyz('test_07.xyz', coords, ans)