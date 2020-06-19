from a2mdio.molecules import Mol2, UNITS_TABLE
from a2mdtest.a2mdtests import benzene
from a2md.models import a2md_from_mol
from a2mdio.volumes import Volume
if __name__ == '__main__':

    mm = Mol2(benzene.mol2)
    dm = a2md_from_mol(mm)
    dm.parametrize()
    w = Volume('C:/Users/Bruno/ownCloud/main/research/edipff/apbs/apbs.charge.dx.dx')
    w.read()
    w.eval(lambda x : dm.eval(x * UNITS_TABLE['angstrom']['au'], kind='density'))
    w.write('test.dx')