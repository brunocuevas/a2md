from a2md.models import a2md_from_mol
from a2md import utils
from a2mdio.molecules import QmSetUp
from a2mdio.molecules import Mol2
import numpy as np
import json
import click
import time
import sys
import logging

logger = logging.getLogger('')

@click.group()
def cli():
    """
    A2MDrun
    Optimization/prediction/evaluation of a2md models of electron density.

    """

    pass

@click.command()
@click.option('--output', default=None, help='file to save the info')
@click.option('--expand', default=2.0, help='file to save the info')
@click.option('--res', default=0.25, help='file to save the info')
@click.argument('name')
@click.argument('param_file')
def write_dx(name, param_file, output, expand, res):
    """

    writes a dx volume file with density from an a2md fit

    """

    start = time.time()

    mm = Mol2(name)
    dm = a2md_from_mol(mm)
    dm.parametrize()
    with open(param_file) as f:
        dm.read(json.load(f))

    dx = dm.eval_volume(spacing=expand, resolution=res)
    if output is None:
        dx.write(name.replace('.mol2', '') + '.dx')
    else:
        dx.write(output)
    print("writting to : {:s}".format(output))
    print("TE : {:12.4f}".format(time.time() - start))

@click.command()
@click.option('--output', default=None, help='file to save the info')
@click.argument('name')
@click.argument('param_file')
@click.argument('coordinates')
def evaluate(name, param_file, coordinates, output):
    """
    reads a model and runs an evaluation upon the specified coordinates
    """

    start = time.time()

    mm = Mol2(name)
    dm = a2md_from_mol(mm)
    dm.parametrize()
    with open(param_file) as f:
        dm.read(json.load(f))

    try:
        coordinates = np.loadtxt(coordinates, dtype='float64')
    except ValueError:
        print("could not read file. please, use space separated values")
        print("sed -i \"s/,/  /g\" COORDINATES")
        sys.exit()
    except FileNotFoundError:
        print("file not found")
        sys.exit()
    assert (type(coordinates) is np.ndarray)
    prediction = dm.eval(coordinates)
    assert (type(prediction) is np.ndarray)
    if output is None:
        for i in range(prediction.size):
            print(
                "{:12:4e} {:12:4e} {:12:4e} {:12:4e}",
                coordinates[i, 0], coordinates[i, 1],
                coordinates[i, 2], prediction[i]
            )
    else:
        # noinspection PyTypeChecker
        np.savetxt(
            output,
            np.stack([coordinates, prediction], axis=1)
        )

    print("TE : {:12.4f}".format(time.time() - start))

@click.command()
@click.option('--charge', default=0, help='charge for qm simulation')
@click.option('--multiplicity', default=1, help='spin multiplicity for qm simulation')
@click.option('--wfn', default=True, help='return a wfn file', type=bool)
@click.option('--population', default='npa', help='population analysis')
@click.option('--basis', default='6-311++G(d,p)', help='basis set')
@click.option('--method', default='MP2', help='qm method')
@click.option('--nprocs', default=4)
@click.argument('name')
def prepare_qm(name, charge, multiplicity, wfn, population, basis, method, nprocs):
    """
    sets a custom qm calculation for a given molecule
    """

    start = time.time()
    mm = Mol2(name)
    adcs = []
    if population != "none":
        adcs.append('pop={:s}'.format(population))
    if wfn:
        adcs.append('output=wfn')
        if method == 'MP2':
            adcs.append('density=current')
    mm.charges = (charge / mm.get_number_atoms())  * np.ones(mm.get_number_atoms(), dtype='float64')
    mm.multiplicity = multiplicity
    qmstp = QmSetUp(
        basis=basis, method=method, calculation_type='single', nprocs=nprocs,
        output='{:s}'.format(name.replace('.mol2', '.wfn')),
        additional_commands=adcs, verbose=False
    )
    qmstp.write_g09(name.replace('.mol2', '.g09.input'), mm)

    print("TE : {:12.4f}".format(time.time() - start))

@click.command()
@click.option('--opt_mode', default='restricted', help='either restricted, unrestricted or semirestricted')
@click.option('--regularization_constant', default=None, help='defines penalty on coefficient norm')
@click.option('--output', default=None, help="file where to store the output parameters")
@click.option('--cluster', default=None, help="use rbf to clusterize by distance signature") # to modify in the future
@click.option('--verbose', default=0, help="0 for no output, 1 for error, 2 for info")
@click.argument('name')
@click.argument('sample')
def fit(name, sample, opt_mode, regularization_constant, output, cluster, verbose):
    """
    ajusts the parameters of a density model to a sample of electron density
    """
    if verbose == 0: logging.basicConfig(level=logging.CRITICAL)
    elif verbose == 1 : logging.basicConfig(level=logging.ERROR)
    else: logging.basicConfig(level=logging.INFO)

    start = time.time()
    logger.info("reading inputs {:s} {:s} ".format(name, sample))
    mm = Mol2(name)
    sample_file = sample
    try:
        logger.info("reading sample file as npy")
        sample = np.load(sample)
    except FileNotFoundError:
        logger.error("could not find the {:s} file".format(sample))
    except ValueError:
        try:
            logger.info("reading csv was unsuccesful. Trying csv")
            sample = np.loadtxt(sample)
        except ValueError:
            logger.error("could not read sample file neither as npy nor csv")
            sys.exit()

    logger.info("reading of mol2 and sample file was succesful")
    logger.info("defining model")
    dm = a2md_from_mol(mm)
    logger.info("parametrizing")
    try:
        dm.parametrize()
    except RuntimeError:
        logger.error("there was some element which is not present in the input parameters")
        sys.exit()
    if cluster is not None:
        logger.info("using clustering of atoms and bonds")
        if cluster == 'rbf':
            logger.info("radial basis functions are used as symmetry function")
            rbf = utils.RBFSymmetryCluster(verbose=False)
            dm.clustering(rbf.cluster)
        else:
            logger.error("no found cluster method {:s}. Aborting".format(cluster))
            sys.exit(1)

    if regularization_constant is not None:
        logger.info("regularization constat is changed to {:12.4e}".format(regularization_constant))
        dm.set_regularization_constant(regularization_constant)

    logger.info("starting optimization, using a opt_mode={:s}".format(opt_mode))
    dm.optimize(sample[:, :3], sample[:, 3], optimization_mode=opt_mode)
    logger.info("finished optimization")

    if output is None:
        output = name.replace('.mol2', '.ppp')

    if cluster is not None:
        logger.info("inflating")
        dm.inflate()

    logger.info("saving to {:s}".format(output))
    with open(output, 'w') as f:
        json.dump(dm.get_parametrization(), f, indent=4)

    print("FIT NAME:{:s} SAMPLE:{:s} MODE:{:s}, TE:{:12.4f}".format(name, sample_file, opt_mode, time.time() - start))

@click.command()
@click.option('--opt_mode', default='restricted', help='either restricted, unrestricted or semirestricted')
@click.option('--regularization_constant', default=None, help='defines penalty on coefficient norm')
@click.option('--segment_charges', default=None, help="allows to define position specific charges")
@click.argument('conformation_data')
@click.argument('fitting_data')
@click.argument('output_name')
def fit_collection(
        conformation_data, fitting_data, output_name, opt_mode, regularization_constant,
        segment_charges
):
    """
    adjusts the parameters of a density model using a wide range of conformations and densities
    """
    start = time.time()

    from a2md.models import ConformerCollection
    from a2md.utils import topology_from_bonds
    from a2md import PARAMETERS_PATH

    with open(conformation_data) as f:
        conformation_data = [i.strip() for i in f.readlines()]
    with open(fitting_data) as f:
        fitting_data = [i.strip() for i in f.readlines()]

    coords_list = []
    fit_coords_list = []
    fit_dens_list = []

    for i, (name, sample) in enumerate(zip(conformation_data, fitting_data)):

        mm = Mol2(name)
        coords_list.append(mm.get_coordinates(units='au'))

        try:
            sample = np.load(sample)
        except FileNotFoundError:
            print("could not find the {:s} file".format(sample))
        except ValueError:
            sample = np.loadtxt(sample)

        sample_contents = np.loadtxt(sample)
        fit_coords_list.append(sample_contents[:, :3])
        fit_dens_list.append(sample_contents[:, 3])

    mm = Mol2(conformation_data[0])

    topo = topology_from_bonds(
        bonds=mm.get_bonds(), natoms=mm.get_number_atoms(),
        nbonds=mm.get_number_bonds()
    )

    dm = ConformerCollection(
        coordinates=coords_list, atomic_numbers=mm.get_atomic_numbers(),
        charge=mm.get_absolute_charges(), topology=topo,
        segments=mm.segment_idx
    )

    if regularization_constant is not None:

        dm.set_regularization_constant(regularization_constant)

    with open(PARAMETERS_PATH + "a2md_topo_bonded_model_extended.json") as f:
        params = json.load(f)

    dm.parametrize(param_dict=params)
    if opt_mode == 'semirestricted':
        dm.use_atomic_number_as_charge()
    if segment_charges is not None:
        dm.modify_charge_by_segment(segment_charges)

    dm.conformer_optimize(
        fit_coords_list, fit_dens_list, optimization_mode=opt_mode
    )

    with open(output_name, 'w') as f:
        json.dump(dm.get_parametrization(), f, indent=4)

    print("TE : {:12.4f}".format(time.time() - start))

@click.command()
@click.argument('name')
@click.argument('wfn')
@click.argument('gaussian_log')
@click.argument('output')
def update_mol2(name, wfn, gaussian_log, output):
    """
    updates a mol2 file using npa charges from gaussian output and coordinates of a wfn
    """
    from a2mdio.qm import WaveFunction, GaussianLog
    from a2mdio.molecules import UNITS_TABLE
    mm = Mol2(name)
    wfn_instance = WaveFunction(file=wfn, prefetch_dm=False, verbose=False)
    glog = GaussianLog(file=gaussian_log, verbose=False)
    mm.coordinates = wfn_instance.get_coordinates() * UNITS_TABLE['au']['angstrom']
    mm.charges = glog.get_charges(kind='partial')
    mm.write(output)

@click.command()
@click.option('--separator', default=',', help="symbol used to delimitate columns")
@click.option('--output_type', default='npy', help="either csv or npy")
@click.option('--random', default=False, help="permutates sample")
@click.option('--npoints', default=None, help="chooses the first N points")
@click.argument('name')
def convert_sample(name, random, separator, output_type, npoints):
    """
    simple and fast conversion from csv to npy and viceversa
    """
    try:
        sample = np.load(name)
    except ValueError:
        sample = np.loadtxt(name)

    if random:
        n = sample.shape[0]
        perm = np.random.permutation(np.arange(n))
        sample = sample[perm]

    if npoints is not None:
        sample = sample[:npoints]

    if output_type == 'csv':
        np.savetxt(name.replace('.npy', '.csv'), sample)
    elif output_type == 'npy':
        np.save(name.replace('.csv', '.npy'), sample)

cli.add_command(evaluate)
cli.add_command(write_dx)
cli.add_command(prepare_qm)
cli.add_command(fit)
cli.add_command(fit_collection)
cli.add_command(update_mol2)

if __name__ == "__main__":

    cli()