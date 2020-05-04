from a2md.models import a2md_from_mol
from a2md.integrate import integrate_density_functional
from a2md.integrate import mse_functional, kullback_leibler_functional, vdwvolume_functional
from a2mdio.molecules import Mol2
from a2mdio.qm import WaveFunction
import json
import click
import time
import sys
import logging
logger = logging.getLogger('')

def admin_sources(mm, reference, reference_type):
    if reference_type == 'wfn':
        reference_d = WaveFunction(verbose=False, file=reference, batch_size=10000)
    elif reference_type == 'a2md':
        reference_d = a2md_from_mol(mm)
        with open(reference) as f:
            reference_d.read(json.load(f))
    else:
        logger.error("use either wfn or a2md as reference format")
        sys.exit()
    return reference_d


@click.group()
def cli():
    """
    A2MDintegrate
    Comparison of electron densities
    """

    pass

@click.command()
@click.option('--reference_type', default='wfn', help='type of reference')
@click.option('--candidate_type', default='wfn', help='type of candidate')
@click.option('--grid', default='coarse', help='coarse, medium or tight')
@click.option('--resolution', default=100, help="radial resolution")
@click.argument("name")
@click.argument("reference")
@click.argument("candidate")
def mse(name, reference, candidate, reference_type, candidate_type, grid, resolution):
    start = time.time()
    mm = Mol2(name)
    reference_d = admin_sources(mm, reference, reference_type)
    candidate_d = admin_sources(mm, candidate, candidate_type)

    msef = mse_functional(ref=reference_d.eval, fun=candidate_d.eval)
    msev = integrate_density_functional(msef, mm, res=resolution, grid=grid)
    print("{:24s} {:24s} MSE {:18.8e} {:8.4f}".format(reference, candidate, msev, time.time() - start))

@click.command()
@click.option('--reference_type', default='wfn', help='type of reference')
@click.option('--candidate_type', default='wfn', help='type of candidate')
@click.option('--grid', default='coarse', help='coarse, medium or tight')
@click.option('--resolution', default=100, help="radial resolution")
@click.argument("name")
@click.argument("reference")
@click.argument("candidate")
def dkl(name, reference, candidate, reference_type, candidate_type, grid, resolution):
    start = time.time()
    mm = Mol2(name)
    reference_d = admin_sources(mm, reference, reference_type)
    candidate_d = admin_sources(mm, candidate, candidate_type)

    dklf = kullback_leibler_functional(ref=reference_d.eval, fun=candidate_d.eval)
    dklv = integrate_density_functional(dklf, mm, res=resolution, grid=grid)
    print("{:24s} {:24s} DKL {:18.8e} {:8.4f}".format(reference, candidate, dklv, time.time() - start))

@click.command()
@click.option('--reference_type', default='wfn', help='wfn or a2md')
@click.option('--epsilon', default=1e-3, help='volume surface density value')
@click.option('--grid', default='coarse', help='coarse, medium, tight')
@click.option('--resolution', default=100, help='radial resolution')
@click.argument('name')
@click.argument('reference')
def volume(name, reference, reference_type, epsilon, grid, resolution):
    start = time.time()
    mm = Mol2(name)
    reference_d = admin_sources(mm, reference, reference_type)
    volf = vdwvolume_functional(reference_d.eval, eps=epsilon)

    vol = integrate_density_functional(volf, mm, res=resolution, grid=grid)
    print("{:24s} VOL(au) {:18.8e} {:8.4f}".format(reference, vol, time.time() - start))


cli.add_command(mse)
cli.add_command(dkl)
cli.add_command(volume)

if __name__ == '__main__':
    cli()