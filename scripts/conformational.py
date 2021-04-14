from ase.md.langevin import Langevin
from ase.optimize import BFGS
from ase import units
from ase import Atoms
import torchani
import torch
import time
import click
import sys
from a2mdio.molecules import Mol2

@click.group()
def cli():
    pass

@click.command()
@click.option('--model', default='ani1x', help='predictor used')
@click.option('--device', default='cpu', help='where data is stored')
@click.argument('molecule')
def energy(molecule, model, device):
    mm = Mol2(molecule)
    pos = mm.get_coordinates()
    labels = ''.join(mm.get_labels())
    device = torch.device(device)
    if model.lower() == 'ani1x':
        model = torchani.models.ANI1x()
    elif model.lower() == 'ani1ccx':
        model = torchani.models.ANI1ccx()
    else:
        print("unknwon model {:s}".format(model))
        sys.exit(1)
    assert isinstance(model, torchani.models.ANI1x) or isinstance(model, torchani.models.ANI1ccx)
    atom_group = Atoms(symbols=labels, positions=pos, calculator=model.ase())
    print("{:26s} {:18.6e}".format(molecule, atom_group.get_potential_energy() / units.Hartree))
    sys.exit(0)

@click.command()
@click.option('--model', default='ani1x', help='predictor used')
@click.option('--device', default='cpu', help='where data is stored')
@click.option('--output', default=None, help='output file')
@click.option('--tolerance', default=1e-2, help='tolerance for optimization')
@click.option('--steps', default=1000, help='maximum number of steps')
@click.argument('molecule')
def optimize(molecule, model, device, output, tolerance, steps):
    start = time.time()
    mm = Mol2(molecule)
    pos = mm.get_coordinates()
    labels = ''.join(mm.get_symbols())
    device=torch.device(device)
    if model.lower() == 'ani1x':
        model = torchani.models.ANI1x()
    elif model.lower() == 'ani1ccx':
        model = torchani.models.ANI1ccx()
    else:
        print("unknwon model {:s}".format(model))
        sys.exit(1)
    assert isinstance(model, torchani.models.ANI1x) or isinstance(model, torchani.models.ANI1ccx)
    atom_group = Atoms(symbols=labels, positions=pos, calculator=model.ase())

    print("\tenergy {:12.4f}".format(atom_group.get_total_energy()))
    opt = BFGS(atom_group, logfile=None)
    res = opt.run(fmax=tolerance, steps=steps)

    if res:
        print("optimization completed, time {:12.4f}".format(time.time() - start))
        print("initial energy {:12.4f}".format(atom_group.get_total_energy()))
        opt_pos = atom_group.get_positions()
        mm.coordinates = opt_pos
        if output is not None:
            mm.write(output)
        sys.exit(0)
    else:
        print("maximum number of steps reached, time {:12.4f}")
        print("final energy {:12.4f}".format(atom_group.get_total_energy()))
        sys.exit(1)

@click.command()
@click.option('--model', default='ani1x', help='predictor used')
@click.option('--device', default='cpu', help='where data is stored')
@click.option('--output', default=None, help='output file')
@click.option('--nsteps', default=10, help='number of iterations on the md')
@click.option('--temperature', default=300, help='temperature (K)')
@click.option('--step', default=20000, help='time per step (fs)')
@click.option('--friction', default=0.2, help='friction')
@click.argument('molecule')
def md(molecule, model, device, output, nsteps, temperature, step, friction):
    start = time.time()
    mm = Mol2(molecule)
    pos = mm.get_coordinates()
    labels = ''.join(mm.get_symbols())
    device = torch.device(device)
    if model.lower() == 'ani1x':
        model = torchani.models.ANI1x()
    elif model.lower() == 'ani1ccx':
        model = torchani.models.ANI1ccx()
    else:
        print("unknwon model {:s}".format(model))
        sys.exit(1)
    assert isinstance(model, torchani.models.ANI1x) or isinstance(model, torchani.models.ANI1ccx)
    atom_group = Atoms(symbols=labels, positions=pos, calculator=model.ase())
    print("step {:12d}\tenergy {:12.4f}".format(0, atom_group.get_total_energy()))
    dyn = Langevin(atom_group, 1*units.fs, temperature*units.kB, friction)

    for i in range(nsteps):
        dyn.run(step)
        print("step {:12d}\tenergy {:12.4f}".format((i + 1)*step, atom_group.get_total_energy()))
        mm.coordinates = atom_group.get_positions()
        if output is not None:
            mm.write(output + '_{:06d}.mol2'.format(i))
    print("time {:12.4f}".format(time.time() - start))
    sys.exit(0)

cli.add_command(optimize)
cli.add_command(md)
cli.add_command(energy)

if __name__ == '__main__':

    cli()
