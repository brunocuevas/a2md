# a2md-torch

This library is an implementation of the a2md models of electron density, using the pytorch 
framework as vectorial reference library, to allow the training of neural networks

## Summary

A2MD is a model of electron density based in the decomposition of electron density as
isotropic and anisotropic contributions, each of them given by specific functions that
share a distance-exponential decay. The reason to implement A2MD again, but in Torch, is
to allow:

- Fast computation in GPU
- Implicit gradient computation

This two features are very important to train neural networks.

## Contents:

- **functions.py**: contains elementary operations as distance and angles calculation.
- **modules.py**: organize the elementary functions into higher order objects that perform
the density calculation when coefficients are provided.
- **test/**: scripts to test this objects. Contain a test of the functions, a test of the modules,
and a test of the autograd/gpu Support.
- **a2md_topo_bonded_model.json**: a modified version of the same file in a2md. The reason for this
modification is to ease some issues regarding the labelling of the functions. Besides that, 
parameters are equivalent.    


