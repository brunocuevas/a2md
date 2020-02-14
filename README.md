# A2MD

Anisotropic Analitically Modelled Density

<img src="a2mdnet_brochure.png" width=500px alt="brochure">

## Summary
A2MD is an ensamble of modules that aim to allow to work with electron density, 
model it, and predict it. It is based on the use of atom centered exponential
functions with anisotropic modifications to model bonding density.

The A2MD repository contains:

- **a2mdlib**: read/write chemical formats.
- **a2md**: fit models of electron electron density
- **a2mdnet**: prediction of model parameters using a deep learning approach
- **cdens**: a fast C executable to sample density
- **a2mdtests**: small ensemble of molecules to test the methodology
- **scripts**: the methodology turned into a cli set of tools

## Dependencies

- numpy
- torch
- torchani (1.1)

## Why separated modules?

The libraries were developped separately, at different times. There is not
a real reason to keep them separated, so it is likely that in the future
they get merged.

## How to

We suggest ti check the test folders within each module. Here we suggest how
to use some of the useful scripts:

Prediction of a density model can be made at the CLI easily using:

	python3 ./scripts/a2mdnet predict --output foo.ppp foo.mol2

Prediction of the electron density of an ensemble of models can be made using:

	python3 ./scripts/a2md evaluate --coordinates coords.xyz foo.mol2 foo.ppp

To visualize electron density, you can use the dx functionality

	python3 ./scripts/a2md write_dx --output volume.dx foo.mol2 foo.ppp

To obtain an electron density isosurface from an electron structure density, 
you can use the Cdens executable:

	./scripts/Cdens -i foo.wfn -u surface.csv -p 0.1 -r 10

To fit an electron density model to a reference calculation, you can use
the a2md script again.

	python3 ./scripts/a2md fit --output foo.ppp --gaussian_output foo.out --reference surface.csv --symmetrize foo.mol2


## Contact

There is a lot of things to polish, and we'd love feedback. Don't hesitate to contact us:

	bruno.czvuria@upm.es

## Licence

This project is distributed under a GNU GPL v3.0, so the only restriction
is to distribute closed source versions.

