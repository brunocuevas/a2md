from a2mdio.molecules import PDB
from a2mdtest.a2mdtests import all20peptide

if __name__ == "__main__":

    p = PDB(all20peptide.path / 'all20peptide.pqr', input_format='pqr')


    q = p.get_absolute_charges()[0:12]
    xyz = p.get_coordinates()[0:12]
    topo = p.get_bonds()
    seq = p.sequence['A']

    print("done")

