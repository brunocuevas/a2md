from a2md.models import polymer_from_pdb
from a2mdio.molecules import PDB
from a2mdtest.a2mdtests import prup3
import json
if __name__ == "__main__":

    acapdb = PDB(prup3.path / '2alg_b.pqr', input_format='pqr')
    acadm = polymer_from_pdb(acapdb, "B")

    with open("../parameters/all20.gp") as f:
        all20 = json.load(f)

    with open("../parameters/a2md_topo_bonded_model_extended.json") as f:
        acadm.parametrization_default = json.load(f)

    acadm.parametrize_as_polymer(all20)

    for i, res_name in enumerate(acadm.residue_names):
        q = acadm.get_residue_charge(i + 1, kind='partial')
        print("{:3d} {:3s} {:8.4f}".format(i, res_name, q))

    # dx = acadm.eval_volume(spacing=2.0, resolution=0.5)
    # dx.write("prup3.dx")

    print("DONE")
