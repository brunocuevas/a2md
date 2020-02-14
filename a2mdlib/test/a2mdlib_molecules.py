from a2mdlib.molecules import Mol2
from a2md.utils import Wavefunction
from a2mdtests import gdb_test
import glob

# In this test we expect that:
# - Mol2 files can be read
# - Mol2 files can be written
# - Mol2 files can be read back
# - G09 can be written by using Mol2 coordinates
# - AAMD can be written by using Mol2 coordinates after using a change of units
if __name__ == '__main__':
    print("test a2md library")
    # TESTING READING OF DIFFERENT MOL2 FILES OBTAINED FROM BABEL
    print("reading Babel Mol2 files")
    mol2_list = glob.glob('../../inp/qm9/set500/*.mol2')
    success_flag = True
    for i, mol2_fn in enumerate(mol2_list):
        try:
            tmp = Mol2(file=mol2_fn)
        except IOError:
            print("\treading test of Babel Mol2 files was not passed. Issues with {:s}".format(mol2_fn))
            success_flag = False
    if success_flag:
        print("success! all files could be read")
    # TESTING TO WRITE A FILE AND READ IT BACK
    print("writting and reading")
    methane = Mol2(file=gdb_test['gdb_000001']['path'] / gdb_test['gdb_000001']['mol2'])
    methane.write(file='test_01.mol2', output_format='mol2')
    try:
        methane = Mol2(file='test_01.mol2')
    except IOError:
        print("\tcould not read a Mol2 file written by the same object")
    # TESTING TO WRITE A G09 FILE
    # NOTE: You have to manually run g09 to test if it works
    print("writting gaussian 09 files")
    # methane.write(file='test_02.g09', output_format='g09', calculation_type='opt-single')
    # methane.write(file='test_03.g09', output_format='g09', calculation_type='single')
    # methane.write(file='test_04.g09', output_format='g09', calculation_type='opt')
    # TESTING TO WRITE AN AAMD FILE
    methane.change_units(units='au')
    wfn = Wavefunction(file=gdb_test['gdb_000001']['path'] / gdb_test['gdb_000001']['wfn'])
    wfn.read()
    coords, labs = wfn.getCoordinates()
    methane.set_coordinates(coordinates=coords, units='au')
    methane.change_units('angstrom')
    methane.write('test_05.mol2', output_format='mol2')


    # Testing on ANIset

    ani_479 = Mol2(
        file='C:/Users/Bruno/ownCloud/main/a2mdtests/ani_000479/ani_000479_0000.mol2'
    )
    ani_479.write('test_06.g09', output_format='g09')
