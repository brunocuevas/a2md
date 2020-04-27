from a2mdio.molecules import PDB
from a2mdtest.a2mdtests import prup3

if __name__ == "__main__":

    p = PDB("C:/scratch/foopeptide.pqr", input_format='pqr')

    # p.read_anotation(prup3.path / '2alg_b.anotation.json')
    # p.add_ss_bonds()
    # q = p.get_absolute_charges()[0:12]
    # xyz = p.get_coordinates()[0:12]
    topo = p.get_bonds()
    # seq = p.sequence['B']

    print("done")

