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


## Contents

### A2MD

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

### A2MDio

A2MDio is a library to handle input and output files. It has some custom
functionality related to A2MD applications.

Contains:
- *molecules*: module to read Mol2 and PDB files.
- *qm*: module to read wavefunction and gaussian09 log files.
- *volumes*: grid representation (dx files)
- *utils*

### A2MDnet

A2MDio is a library to make predictions of electron density models. It uses
PyTorch as backend.

### Scripts

If you want to avoid learning how the library works and scripting, don't worry:
we have written some useful scripts with a nice interface so you can run most of the
standard operations. 

Contact us if you think that somehting is missing here!

## Usage

### Installation

To ease installation, you can install a Python environment containing all the
modules that A2MD requires. The file can be found at the root folder of this repo.
Given that almost all the library has been written in Python, you will not have to 
compile but the Cdens utility, which should be as easy as :

    gcc density.c -lm -Ofast density.c

We recommend creating a path variable to call the different scripts. You can
append the next line to your .bashrc

    A2MD=/path/to/scripts/

If you experience any problem, contact us!
NOTE: We don't have a Mac, so we unknow which sort of problems Mac users will experience.

### Command Line Tools

We'll begin by learning how to use the scripts contained in the *scripts* folder. There
you can observe three scripts:

- **a2mdrun.py** : operations with a2md models.
- **a2mdpredict.py** : operations with a2md model prediction using neural networks.
- **a2mdutils.py**: utilities to work with the different file formats.
- **a2mdcompare.py**: programs to evaluate electron density comparissons.

To understand the different scripts, we'll run some case studies:

#### Visualizing an electron density from a QM calulation

We'll generate a DX file from a QM calculation contained at the 
./a2mdtest/a2mdtest/gdb_000214/ folder

    Cdens -i gdb_000214.wfn -d gdb_000214.dx -r 0.2 -s 3.0
    
The result is a .dx called gdb_000214.dx, which can be read, among others, in
Chimera. You should be visualizing the density of benzene right now.

#### Fit an A2MD model to benzene

A2MD module consists on linear models of electron density which are fitted
by using a "sophisticated" linear regression to quantum reference electron density.
The first step is to obtain a sample of the electron density to fit our model. To do
so, we'll use again (and for the last time) our C utility Cdens.

    Cdens -i gdb_000214.wfn -u gdb_000214.surface1.csv -r 0.25 -s 3.0
    
This will provide a csv file containing three columns of coordinates and a fourth
column with the electron density in electron/(Bohr^3) units (**BTW: we always use atomic
units; at least internally**). 

Our fitting method consists on sampling different surfaces at the core-valence region. 
We like the isosurfaces 0.25, 0.1, and 0.01.

    Cdens -i gdb_000214.wfn -u gdb_000214.surface1.csv -r 0.25 -s 10
    Cdens -i gdb_000214.wfn -u gdb_000214.surface2.csv -r 0.10 -s 10
    Cdens -i gdb_000214.wfn -u gdb_000214.surface3.csv -r 0.01 -s 10
    ls *.surface?.csv | sed "/^x,/d"
    ls *.surface?.csv | sed "s/,/  /g"
    cat *.surface?.csv > gdb_000214.csv
    rm *surface?.csv
    
 Next step is getting a mol2 file in the same coordinates that the resulting electron
 density. This might sound obvious, but it is not: some qm programs like G09
 change the file to what they call "Standard Orientation". So the coordinates
 of the WFN might not coincide with the coordinates of input of your G09 file.
 We recommend to check coordinates with the WFN. And remember : **WFN file 
 coordinate units are Bohrs, and Mol2 files and G09 inputs are Angstroms**.
 
 In our case, we already have a well oriented Mol2 file called "gdb_000214.mol2".
 In the future, we'll see some options we provide to ease this step.
 
 Now, you can call A2MDrun script to fit an A2MD model declared with the information
 of the Mol2 file and the electron density sample.
 
    python3 ${A2MD}/a2mdrun.py fit gdb_000214.mol2 gdb_000214.csv   

If everything run smoothly, you should have a *gdb_000214.ppp* file in your workfolder.

#### Understanding A2MD/fit

If you type 
    
    python3 ${A2MD}/a2mdrun.py fit --help
    
you will get the following message:

    Usage: a2mdrun.py fit [OPTIONS] NAME SAMPLE
    
    Options:
      --opt_mode TEXT                 either restricted, unrestricted or
                                      semirestricted
    
      --scheme TEXT                   either default, harmonic, extended, spheric
      --regularization_constant FLOAT
                                      defines penalty on coefficient norm        
      --output TEXT                   file where to store the output parameters  
      --cluster TEXT                  use rbf to clusterize by distance signature
      --verbose INTEGER               0 for no output, 1 for error, 2 for info   
      --help                          Show this message and exit.
    
There are many options, all related to some of the critical steps of the a2md fitting
process. Let's explore them:

- First of all, we have to understand our input: the Mol2 file: 

        Mol2 -> (elements, coordinates, atom types, charges, topology)
    
    In A2MD, we use three of these pieces of information: elements, coordinates, 
    charges and topology. From these data, A2MD builds an internal representation
    of the molecule that will be parametrized according to different criteria.
    
- For each molecule, a specific model of electron density is required according 
to its coordinates, atomic species, and bonding patterns. We provide different 
parametrization schemes with different features. Just to enumerate them:
    1. Default: a model containing exponential and exponential-gaussian terms, as exposed in 
    Cuevas and Pacios, 2020.
    2. Harmonic: a new model that we are currently building, based on spherical harmonics.
    3. Extended: a model containing parameters for sulphur species. Under testing.
    4. Spherical: classic spherical Pacios 92 model.
    
    More technical details about each of the sets of parameters will be introduced later. And you will check how
    easy is to create your own.    
    
- Once each molecule has been given a "naïve model of electron density", a possible middle step
is to set up symmetry restraints, so equal or equivalent atoms are given the same electron density
parameters during optimization. This option is set up by the "--cluster" option. By now, only a symmetry clustering 
technnique has been applied, based on gaussian environments of atoms (RBF). Some more techniques may 
arrive in the near future.
- The optimization technique that is being used right now is a simple-but-sophisticated least squares with
Lagrange multipliers and regularization pennalties. The regularization pennalties are used to keep the model
as simple as possible for a given density dataset (otherwise, things can go south in very unexpected ways. 
Simple is always better than intrincate). But the value of the pennalties is arbitrary, so we provide the option
--regularization_constant so the user can play with the constant. In general lines, we don't advise to move it from
its reference value (1e-4).
- The lagrange multipliers are used to constraint electron density to match the number of electrons of the
molecule. But there are multiple ways to do it. By now, we offer two solid options: a "restricted" optimization
in which the atomic charges of the Mol2 file that you provided are used to restraint the electron density functions
corresponding to each of the atoms; and an "unrestricted" optimization, in which only a restraint for global atomic 
charge is set. Which is better? It depends. If you performed some detailed charge analysis (for instance, NPA) in your 
reference calculation, then it is better to use it. But if you skipped it, then it will be better to use a global
charge restraint.

#### Visualizing the result

So what can we do with the electron density model? We are developping new functionalities. But let's start
with something easy: visualizing the resulting electron density model.

    python3 ${A2MD}/a2mdrun.py write-dx --output=gdb_000214.a2md.dx --expand=3.0 --res=0.25 --kind=density gdb_000214.mol2 gdb_000214.ppp   

This should return an electron density model that is generally similar to the high-level reference
electron density, even if it contains way fewer functions.

#### Setting up a QM calculation

Now let's make a step back: how do we get a reference calculation? By now, you might have guessed that we need an
AIM Wavefunction File, which is a representation of the wavefunction following a fixed format file. Although you
may already have experience with QM software, we have created a script to convert Mol2 files to Gaussian/ORCA 
inputs that would allow you to get the wavefunction file.

For instance, in the case of benzene, you can use:

    python3 ${A2MD}/a2mdutils.py prepare-qm --charge=0 --multiplicity=1 --wfn=True --population=NPA --basis="6-311++G(d,p)" --method=MP2 --nprocs=1 --program=g09 gdb_000214.mol2 
    
The resulting gdb_000214.g09.input file can be run in Gaussian09.

    g09 < gdb_000214.g09.input > gdb_000214.g09.output 

Since G09 is a commercial software, we also provide the possibility of using ORCA, which is free. But some of the stuff
is way more difficult to do in ORCA, since you will need to perform a separated population analysis, and to make some
modifications to the WFN files. If you want to use ORCA to run your models, contact us. 

What is the best setup (theory level, population 
technique, etc) for the models? We have mostly played with MP2/6-311++G(d,p) level and Natural Population Analysis. 
Results can be checked in our article. But building this open library would not make sense if we didn't encourage 
you to play with other theory levels and set ups in general.

We are open to include more software options! 

#### Be careful with rotations!!

Gaussian09 usually moves the molecule coordinates to allign its inertia axis with the x,y,z axis. So, unless told 
the opposite (symmetry=None) the resulting density will be in a new coordinates axis, and the density will not match.

Since we like the idea of using an objetive system of reference, we provide a utility to update the original mol2 file
to the WFN coordinates. The script can also read QM charges from a Gaussian09 output if they were obtained either using
NPA or MK population analysis.

    python3 ${A2MD}/a2mdutils.py update-mol2 gdb_000214.mol2 gdb_000214.wfn gdb_000214.g09.output gdb_000214.r.mol2
    
#### Skipping QM calculations: using machine learning

To avoid making a QM calculation, we have developped a new technique based on neural networks to generate optimized
electron density models for organic molecules. If you wan to know the details, you can read our latest article (Cuevas
and Pacios 2020).

Developing a new model of electron density should be as easy as:

    python3 ${A2MD}/a2mdpredict.py --predictor=a2mdc --output=gdb_000214.a2mdnet.ppp gdb_000214.mol2
    
As you can notice, it should have taken less than 1 second, and most models we have developped until now seem 
reasonable. But be careful! This new tool is still experimental, and inputs lying very far away from the training sets
could led to unreasonable results.

#### Studying molecular features: the electrostatic potential

We are still expanding the features that can be studied using our representation of electron density. But
one of the implemented ones is the electrostatic potential calculated from a continuum representation of charge.

You can study the EP of your fitted model using:

    python3 ${A2MD}/a2mdrun.py write-dx --output=gdb_000214.a2md.ep.dx --expand=3.0 --res=0.25 --kind=ep gdb_000214.mol2 gdb_000214.ppp

#### Exploring the scripts

You can always use the --help flag to find out what each of the scripts functions does. As A2MD is still growing,
we are constantly adding new functions. 


### Library

### Cdens

