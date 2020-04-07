from a2md.models import a2md_from_mol
from a2mdio.molecules import Mol2
from a2mdtest.a2mdtests import aca
import json
import numpy as np
if __name__ == "__main__":

    print("a2md/molecules")
    print("---")
    m = Mol2(file=aca.mol2)
    dm = a2md_from_mol(m)

    with open(r"C:\Users\Bruno\ownCloud\projects\a2md\a2md\parameters\a2md_topo_bonded_model_extended.json") as f:
        dm.parametrize(json.load(f))

    training_coords  = np.loadtxt(aca.surfaces[1], skiprows=1, delimiter=',')
    training_density = training_coords[:,  3]
    training_coords  = training_coords[:, :3]
    #
    c = dm.optimize(training_coordinates=training_coords, training_density=training_density)
    #
    # dm.inflate()
    #
    # z = dm.integrate()
    # print("integral : {:8.4f}".format(z))
    #
    with open("aca.ppp", "w") as f:
        json.dump(dm.get_parametrization(), f, indent=4, sort_keys=True)
    #
    dm.read(dm.get_parametrization())
    print("DONE!")
