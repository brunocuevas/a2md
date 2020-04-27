from a2md.models import ConformerCollection
from a2mdio.molecules import Mol2
from a2md.utils import topology_from_bonds
import numpy as np

if __name__ == '__main__':

    conformers_training_coords = []
    conformers_training_density = []
    conformers_coords = []
    np.random.seed(42)
    for i in range(20):  # using the first 20 conformations
        mm = Mol2("F:/edip/asa/ASA_{:06d}_upd.mol2".format(i))
        conformers_coords.append(mm.get_coordinates(units='au'))
        tmp = np.concatenate(
            [
                np.loadtxt("F:/edip/asa/ASA_{:06d}_s1.csv".format(i)),
                np.loadtxt("F:/edip/asa/ASA_{:06d}_s2.csv".format(i)),
                np.loadtxt("F:/edip/asa/ASA_{:06d}_s3.csv".format(i))
            ]
        )

        conformers_training_coords.append(tmp[:, :3])
        conformers_training_density.append(tmp[:, 3])


    mm = Mol2("F:/edip/asa/ASA_{:06d}.mol2".format(0))

    charge = mm.get_absolute_charges()
    an = mm.get_atomic_numbers()

    topology = topology_from_bonds(
        bonds=mm.get_bonds(), natoms=mm.get_number_atoms(),
        nbonds=mm.get_number_bonds()
    )

    asa_conf20 = ConformerCollection(
        coordinates=conformers_coords, atomic_numbers=an, charge=charge,
        topology=topology, segments=mm.segments
    )
    asa_conf20.parametrize()

    asa_conf20.set_regularization_constant(1e-4)
    c = asa_conf20.conformer_optimize(
        training_coordinates=conformers_training_coords,
        training_densities=conformers_training_density,
        optimization_mode='semirestricted'
    )
    print(c)

    asa_conf20.parametrize_conformer(0)
    dx = asa_conf20.eval_volume(spacing=2.0, resolution=0.5)
    dx.write("asa.dx")