from a2md.models import polymer_from_pdb
from a2mdio.molecules import PDB
from a2mdtest.a2mdtests import all20peptide
import json
if __name__ == "__main__":

    aca_mol = PDB(all20peptide.path / "all20peptide.pqr", input_format='pqr')
    aca_polymer = polymer_from_pdb(aca_mol, "A")

    with open("../parameters/all20.gp") as f:
        all20 = json.load(f)

    with open("../parameters/a2md_topo_bonded_model_extended.json") as f:
        aca_polymer.parametrization_default = json.load(f)

    aca_polymer.parametrize_as_polymer(all20)

    aca_dx = aca_polymer.eval_volume(spacing=2.0, resolution=0.5)
    aca_dx.write("all20.dx")

    print("DONE!")
