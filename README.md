# A2MD

Anisotropic Analitically Modelled Density

**NOTE: We are retraining the networks, to ensure reproducibility**

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
- **examples**

Read the docs folder to learn to use the library and the CLI tools. 

## Dependencies

Some of the most important dependencies are:
- numpy
- torch
- torchani (1.2)

We provide a "requirements.txt" to ease installation through pip. 

## Contact

There is a lot of things to polish, and we'd love feedback. Don't hesitate to contact us:

	bruno.czvuria@upm.es

## Licence

This project is distributed under a GNU GPL v3.0, so the only restriction
is to distribute closed source versions.

