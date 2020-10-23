import click
from a2mdio.molecules import Mol2
from a2mdnet.data import DatasetCard
import h5py
from pathlib import Path

ALLOWED_SYMBOLS = ['H', 'C', 'N', 'O', 'S']


@click.command()
@click.option('remove_fs', type=bool)
@click.argument('name', type=str)
@click.argument('wfn_path', type=str)
@click.argument('mol_path', type=str)
@click.argument('index_file', type=str)
def build_dataset(
    name, wfn_path, mol_path, index_file
):
    mol_path = Path(mol_path)
    print('-- building dataset: {:s}'.format(name))
    with open(index_file) as f:
        index = [i.strip() for i in f.readlines()]
    print('-- index : {:d} files '.format(len(index)))

    print('-- reading hdf5 keys')
    f = h5py.File(wfn_path, 'r')
    wfn_keys = list(f.keys())
    f.close()

    print('-- starting check')
    cleared_ = []
    checked = []
    max_atoms = 0
    max_bonds = 0
    for i in index:
        if i not in wfn_keys:
            cleared_.append(i)
            continue

        m = Mol2(mol_path / i)
        symbols = m.get_symbols()
        if any([s not in ALLOWED_SYMBOLS for s in symbols]):
            cleared_.append(i)
            continue

        if m.get_number_atoms() > max_atoms:
            max_atoms = m.get_number_atoms()
        if m.get_number_bonds() > max_bonds:
            max_bonds = m.get_number_bonds()

        checked.append(i)

    n = len(checked)

    for cl in cleared_:
        print('-- discarded entry : {:s}'.format(cl))

    ds = DatasetCard(
        name, description="", n_molecules=n, max_atoms=max_atoms, max_bonds=max_bonds,
        path_dict=dict(wfn=wfn_path, mol=mol_path), properties=None, index=checked
    )
    print('-- saving to: {:s}'.format(name + '.dataset'))
    ds.to_json(name + '.dataset')
