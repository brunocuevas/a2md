from a2md.models import a2md_from_mol
from a2md import utils
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
@click.option('--kind', default='density', help='either density or ep')
@click.argument('name')
@click.argument('param_file')
def write_dx(name, param_file, output, expand, res, kind):
    """

    writes a dx volume file with density from an a2md fit

    """

    start = time.time()

    mm = Mol2(name)
    dm = a2md_from_mol(mm)
    dm.parametrize()
    with open(param_file) as f:
        dm.read(json.load(f))

    dx = dm.eval_volume(spacing=expand, resolution=res, kind=kind)
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
                "{:12.4e} {:12.4e} {:12.4e} {:12.4e}".format(
                    coordinates[i, 0], coordinates[i, 1],
                    coordinates[i, 2], prediction[i]
                )
            )
    else:
        # noinspection PyTypeChecker
        np.savetxt(
            output,
            np.stack([coordinates, prediction], axis=1)
        )

    print("TE : {:12.4f}".format(time.time() - start))


@click.command()
@click.option('--opt_mode', default='restricted', help='either restricted, unrestricted or semirestricted')
@click.option('--scheme', default='default', help='either default, harmonic, extended, spheric')
@click.option('--regularization_constant', default=None, help='defines penalty on coefficient norm', type=float)
@click.option('--output', default=None, help="file where to store the output parameters")
@click.option('--cluster', default=None, help="use rbf to clusterize by distance signature") # to modify in the future
@click.option('--verbose', default=0, help="0 for no output, 1 for error, 2 for info")
@click.argument('name')
@click.argument('sample')
def fit(name, sample, opt_mode, scheme, regularization_constant, output, cluster, verbose):
    __fit_call(name, sample, opt_mode, scheme, regularization_constant, output, cluster, verbose)

def __fit_call(name, sample, opt_mode, scheme, regularization_constant, output, cluster, verbose):
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
    except OSError:
        try:
            logger.info("reading csv was unsuccesful. Trying csv")
            sample = np.loadtxt(sample)
        except ValueError:
            logger.error("could not read sample file neither as npy nor csv")
            sys.exit()
    except ValueError:
        sample = np.loadtxt(sample)

    logger.info("reading of mol2 and sample file was succesful")
    logger.info("defining model")
    dm = a2md_from_mol(mm)
    logger.info("parametrizing")
    try:
        if scheme == 'default': dm.parametrize()
        elif scheme == 'extended' : dm.parametrize(dm.parametrization_extended)
        elif scheme == 'harmonic' : dm.parametrize(dm.parametrization_harmonic)
        elif scheme == 'spheric': dm.parametrize(dm.parametrization_spherical)
        else:
            print("use a default, extended, harmonic or spheric scheme")
            sys.exit(1)

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
@click.option('--regularization_constant', default=None, help='defines penalty on coefficient norm', type=float)
@click.option('--scheme', default='default', help='either default, harmonic, extended, spheric')
@click.option('--cluster', default=None, help="use rbf to clusterize by distance signature") # to modify in the future
@click.option('--verbose', default=0, help="0 for no output, 1 for error, 2 for info")
@click.argument('names_file')
def fit_many(names_file, opt_mode, regularization_constant, scheme, cluster, verbose):
    """
    fits many compounds
    """
    with open(names_file) as f:
        names = json.load(f)
    mol2_ = [i['mol2'] for i in names]
    out_ = [i['output'] for i in names]
    sample_ = [i['sample'] for i in names]
    n = len(mol2_)
    global_start = time.time()
    for i, (m, s, o) in enumerate(zip(mol2_, sample_, out_)):
        __fit_call(
            m, s, opt_mode=opt_mode, regularization_constant=regularization_constant, output=o,
            cluster=cluster, verbose=verbose, scheme=scheme
        )
    global_end = time.time()
    time_elapsed = global_end - global_start
    print("FIT_MANY, {:8d} molecules, TE={:8.4f} s, TEPM={:8.4f} s/mol".format(n, time_elapsed, time_elapsed/n))

@click.command()
@click.option('--mol2_path', default=None, help="path for mol2 files")
@click.option('--sample_path', default=None, help="path for npy/csv files")
@click.option('--output_path', default=None, help="path for ppp files")
@click.option('--filter_names', default=None, help="list of names to avoid")
@click.option('--output_format', default='json', help="either json or txt")
@click.option('--split', default=1, help="divide the output in n files")
@click.argument('names')
@click.argument('output_file')
def prepare_fit_many(names, output_file, mol2_path, sample_path, output_path, filter_names, output_format, split):
    """
    prepares a set of files to run fit many
    """
    start = time.time()
    with open(names) as f:
        names = json.load(f)
    if filter_names is not None:
        with open(filter_names) as f:
            filter_names = json.load(f)
    else:
        filter_names=[]
    output = []
    if mol2_path is not None : mol2_str = mol2_path + "/{:s}.mol2"
    else: mol2_str = "{:s}.mol2"
    if sample_path is not None : sample_str = sample_path + "/{:s}.npy"
    else: sample_str = "{:s}.npy"
    if output_path is not None : output_str = output_path + "/{:s}.ppp"
    else: output_str = "{:s}.ppp"

    for i, n in enumerate(names):
        if n in filter_names:
            continue
        output.append(dict(
            mol2=mol2_str.format(n), sample=sample_str.format(n), output=output_str.format(n)
        ))

    len_div = len(output) // split
    lo = len(output) % split
    if lo != 0:
        split = split + 1

    if split != 1:
        output_file_format = output_file + '.{:02d}'
    else:
        output_file_format = output_file

    for i in range(split):
        current_chunk = output[i*len_div : (i+1)*len_div]
        output_file = output_file_format.format(i)

        if output_format == 'json':
            with open(output_file, "w") as f:
                json.dump(current_chunk, f, indent=4, sort_keys=True)
        elif output_format == 'txt':
            with open(output_file, "w") as f:
                for line in current_chunk:
                    f.write('{:26s} {:26s} {:26s}\n'.format(line['mol2'], line['sample'], line['output']))
    print("TE = {:8.4f}".format(time.time() - start))

@click.command()
@click.option('--opt_mode', default='restricted', help='either restricted, unrestricted or semirestricted')
@click.option('--regularization_constant', default=None, help='defines penalty on coefficient norm', type=float)
@click.option('--segment_charges', default=None, help="allows to define position specific charges")
@click.option('--scheme', default='default', help='either default, harmonic, extended, spheric')
@click.argument('conformation_data')
@click.argument('fitting_data')
@click.argument('output_name')
def fit_collection(
        conformation_data, fitting_data, output_name, opt_mode, regularization_constant,
        segment_charges, scheme
):
    """
    adjusts the parameters of a density model using a wide range of conformations and densities
    """
    start = time.time()

    from a2md.models import ConformerCollection
    from a2md.utils import topology_from_bonds

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


    if scheme == 'default':
        dm.parametrize()
    elif scheme == 'extended':
        dm.parametrize(dm.parametrization_extended)
    elif scheme == 'harmonic':
        dm.parametrize(dm.parametrization_harmonic)
    elif scheme == 'iso':
        dm.parametrize(dm.parametrization_spherical)
    else:
        print("use a default, extended, harmonic or spheric scheme")
        sys.exit(1)
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


cli.add_command(evaluate)
cli.add_command(write_dx)
cli.add_command(fit)
cli.add_command(fit_many)
cli.add_command(fit_collection)
cli.add_command(prepare_fit_many)

if __name__ == "__main__":

    cli()