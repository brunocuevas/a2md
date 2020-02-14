from a2md.preprocessor import preprocessor
from a2md.amd import A2MD
import numpy as np
import json
import click
import time
import sys

@click.group()
def cli():
    pass

@click.command()
@click.option('--output', default=None, help='file to save the info')
@click.option('--expand', default=2.0, help='file to save the info')
@click.option('--res', default=0.25, help='file to save the info')
@click.argument('name')
@click.argument('param_file')
def write_dx(name, param_file, output, expand, res):

    start = time.time()
    ppp = preprocessor(file=name, input_format='mol2')
    wx, wc, top, pars = ppp.parametrize(kind='frozen_core_gaussian')
    an = ppp.get_atomic_numbers()
    with open(param_file) as f:
        params = json.load(f)
    density_model = A2MD(
        coordinates=wx, atomic_numbers=an, charge=wc, topology=top,
        parameters=params, verbose=False
    )
    dx = density_model.eval_volume(extend=expand, resolution=res, field='density')
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

    start = time.time()
    ppp = preprocessor(file=name, input_format='mol2')
    wx, wc, top, pars = ppp.parametrize(kind='frozen_core_gaussian')
    an = ppp.get_atomic_numbers()
    with open(param_file) as f:
        params = json.load(f)
    density_model = A2MD(
        coordinates=wx, atomic_numbers=an, charge=wc, topology=top,
        parameters=params, verbose=False
    )
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
    prediction = density_model.eval(coordinates)
    assert (type(prediction) is np.ndarray)
    if output is None:
        for i in range(prediction.size):
            print(
                "{:8:4f} {:8:4f} {:8:4f} {:8:4f}",
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

cli.add_command(evaluate)
cli.add_command(write_dx)

if __name__ == "__main__":

    cli()