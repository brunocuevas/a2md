from a2mdio.molecules import Mol2
from a2mdio.geometry import convert
from a2mdtest.a2mdtests import aca

if __name__ == "__main__":

    mm = Mol2(aca.mol2)
    b, a, d = convert(mm)

    print("DONE!")

