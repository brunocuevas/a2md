import torch
from a2mdnet.utils import Parametrizer
from a2md.preprocessor import preprocessor
from a2md.amd import A2MD
from a2mdlib.molecules import Mol2
import click
import json
import time

@click.group()
def cli():
    pass

@click.command()
@click.option('--predictor', default='x96x24_0200.pt', help='name of the network')
@click.option('--device', default='cpu', help='device to perform the inference')
@click.option('--output', default=None, help='file to save the info')
@click.argument('name')
def predict_model(name, predictor, device, output):


    start = time.time()
    dev = torch.device(device=device)
    ppp = preprocessor(file=name, input_format='mol2')
    wx, wc, top, pars = ppp.parametrize(kind='frozen_core_gaussian')
    model = torch.load(predictor, map_location=dev).to(dev)
    param = Parametrizer(model, device=dev)
    pars = param.parametrize(name, pars)
    if output is None:
        print(json.dumps(pars, indent=4))
    else:
        with open(output, 'w') as f:
            json.dump(pars, f, indent=4, sort_keys=True)

    print("TE : {:12.4f}".format(time.time() - start))


@click.command()
@click.option('--predictor', default='x96x24_0200.pt', help='name of the network')
@click.option('--device', default='cpu', help='device to perform the inference')
@click.option('--output', default=None, help='file to save the info')
@click.argument('name')
def predict_charges(name, predictor, device, output):


    start = time.time()
    dev = torch.device(device=device)
    ppp = preprocessor(file=name, input_format='mol2')
    wx, wc, top, pars = ppp.parametrize(kind='frozen_core_gaussian')
    an = ppp.get_atomic_numbers()
    model = torch.load(predictor, map_location=dev).to(dev)
    param = Parametrizer(model, device=dev)
    pars = param.parametrize(name, pars)
    density_model = A2MD(
        coordinates=wx, atomic_numbers=an, charge=wc, topology=top,
        parameters=pars
    )
    charges = density_model.get_a2md_molecule_charge()
    symbols = ppp.get_labels()
    if output is None:

        for i in range(len(charges)):
            print(
                "{:8d} {:8s} {:8.4e} {:8.4e} {:8.4e} {:8.4e}".format(
                    i, symbols[i], wx[i, 0], wx[i, 1], wx[i, 2], charges[i]
                )
            )
    else:
        mm = Mol2(file=name)
        mm.set_charges(charges, kind="total")
        mm.write(file=output, output_format='mol2')

    print("TE : {:12.4f}".format(time.time() - start))



cli.add_command(predict_model)
cli.add_command(predict_charges)

if __name__ == '__main__':

    cli()
    print("finished!")