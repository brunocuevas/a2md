import torch
from a2mdnet.utils import Parametrizer
from a2mdnet import MODELS
from a2md.models import a2md_from_mol
from a2mdio.molecules import Mol2
import click
import json
import time

@click.group()
def cli():
    pass

@click.command()
@click.option('--predictor', default='a2mdc', help='use a2mdc')
@click.option('--device', default='cpu', help='device to perform the inference')
@click.option('--output', default=None, help='file to save the info')
@click.argument('name')
def predict_model(name, predictor, device, output):


    start = time.time()
    dev = torch.device(device=device)
    mm = Mol2(name)
    dm = a2md_from_mol(mm)
    dm.parametrize()
    pars = dm.get_parametrization()
    model = torch.load(MODELS[predictor], map_location=dev).to(dev)
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
    mm = Mol2(name)
    dm = a2md_from_mol(mm)
    dm.parametrize()
    pars = dm.get_parametrization()
    model = torch.load(predictor, map_location=dev).to(dev)
    param = Parametrizer(model, device=dev)
    pars = param.parametrize(name, pars)
    dm.read(pars)
    charges = dm.get_a2md_charges()
    if output is None:
        for i, (sym, coords, q) in enumerate(zip(dm.get_symbols(), dm.get_coordinates(), charges)):
            print(
                "{:8d} {:8s} {:8.4e} {:8.4e} {:8.4e} {:8.4e}".format(
                    i, sym, coords[0], coords[1], coords[2], q
                )
            )
    else:
        mm = Mol2(file=name)
        mm.charges = charges
        mm.write(file=output)

    print("TE : {:12.4f}".format(time.time() - start))



cli.add_command(predict_model)
cli.add_command(predict_charges)

if __name__ == '__main__':

    cli()
    print("finished!")