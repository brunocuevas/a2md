# A2MD manual


## Summary

A2MD is a new methodology based on linear models of electron density plus 
machine learning methods to parametrize it. We aim to develop a new 
methodology in the machine-learning for chemistry field that considers
electron structure as a way to unify the different descriptions (MM, QM)
and to derive properties without forcing the method to learn them.

### What can I do with a2md?

You can:
- prepare QM calculations to run in either Gaussian or ORCA
- sample density isosurfaces from AIM wave-function files.
- fit a linear model to the sampled densities
- obtain properties from the density model (shape, electrostatic potential)
- represent the molecular isosurface in a grid
- evaluate electron density at given points
- symmetrize the molecule using a custom function
- use a neural network to parametrize a model
- obtain atomic charges from an ML prediction

and we are still working to bring the calculation of energies to the method.

### How can I use the library?

We have developped two main methods to use A2MD:
- through the scripts located at a2md/scripts/.
- calling the library.

You can simply add the library location to your PYTHONPATH.

## A2MD

### Fundaments

A2MD is a method based in exponential radial functions model [AMD]() that
includes anisotropic functions to address the bonding-derived concentration
of density in molecules. This model can be fitted by using a restricted linear
model that accounts for charge normalization and regularization.

The A2MD module is divided in different modules:
- *models*. The core of the library. Here you can find the linear models we have
developped.
- *support*. Contains the functions that are used in models.
- *mathfunctions*. Mathematical functions used in the support.
- *integrate*. Numerical integration of electron density.
- *utils*. Clustering and other analysis tools.

## A2MDio

A2MDio is a library to handle input and output files. It has some custom
functionality related to A2MD applications.

Contains:
- *molecules*: module to read Mol2 and PDB files.
- *qm*: module to read wavefunction and gaussian09 log files.
- *volumes*: grid representation (dx files)
- *utils*


## A2MDnet

A2MDio is a library to make predictions of electrond density models. It uses
PyTorch as backend.

