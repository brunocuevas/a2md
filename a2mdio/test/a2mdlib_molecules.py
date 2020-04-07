from a2mdio.molecules import Mol2
from a2mdtest.a2mdtests import aca

# In this test we expect that:
# - Mol2 files can be read
# - Mol2 files can be written
# - Mol2 files can be read back
# - G09 can be written by using Mol2 coordinates
# - AAMD can be written by using Mol2 coordinates after using a change of units
if __name__ == '__main__':
    print("test a2md library")

    # Testing on ANIset

    aca = Mol2(
        file=aca.mol2
    )
    aca.write("aca.mol2")

    print("DONE!")
