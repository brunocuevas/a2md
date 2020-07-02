import click
from a2md.models import a2md_from_mol
from a2mdio.molecules import Mol2, UNITS_TABLE
from a2mdio.qm import WaveFunction
from a2mdio.volumes import Volume
import json

@click.group()
def cli():
    pass

@click.command()
@click.option('--charge', default=0)
@click.argument('name')
@click.argument('dx')
@click.argument('output')
def dx_add_wfn_charge(name, dx, charge, output):
    wfn = WaveFunction(file=name, verbose=False, batch_size=1000000)
    convert2au = lambda x: x*UNITS_TABLE['angstrom']['au']
    r3 = UNITS_TABLE['angstrom']['au'] ** 3
    fun = lambda x : -wfn.eval(convert2au(x))*r3

    dx1 = Volume(filename=dx)
    dx1.read()
    dx1.eval(fun)
    dx2 = Volume(filename=dx)
    dx2.read()

    q = dx2._Volume__dx.sum() + charge

    qt = dx1._Volume__dx.sum() * (dx1.get_basis()[0, 0] ** 3)
    dx1._Volume__dx = dx1._Volume__dx * (q/ qt)
    dx1._Volume__dx = -dx1._Volume__dx + dx2._Volume__dx
    dx1.write(output)
    return

@click.command()
@click.option('--charge', default=0)
@click.argument('name')
@click.argument('params')
@click.argument('dx')
@click.argument('output')
def dx_add_a2md_charge(name, params, dx, charge, output):
    mm = Mol2(name)
    dm = a2md_from_mol(mm)
    with open(params) as f:
        dm.read(json.load(f))


    convert2au = lambda x: x*UNITS_TABLE['angstrom']['au']
    r3 = UNITS_TABLE['angstrom']['au'] ** 3
    fun = lambda x : -dm.eval(convert2au(x))*r3

    dx1 = Volume(filename=dx)
    dx1.read()
    dx1.eval(fun)
    dx2 = Volume(filename=dx)
    dx2.read()

    q = dx2._Volume__dx.sum() * (dx1.get_basis()[0, 0] ** 3) + charge

    qt = dx1._Volume__dx.sum() * (dx1.get_basis()[0, 0] ** 3)
    dx1._Volume__dx = dx1._Volume__dx * (q/ qt)
    dx1._Volume__dx = -dx1._Volume__dx + dx2._Volume__dx
    dx1.write(output)
    return

@click.command()
@click.argument('name')
@click.argument('output')
def zmol2(name, output):
    mm = Mol2(file=name)
    mm.charges = mm.get_atomic_numbers()
    mm.write(output)
    return

cli.add_command(dx_add_wfn_charge)
cli.add_command(dx_add_a2md_charge)
cli.add_command(zmol2)

if __name__ == '__main__':
    cli()

