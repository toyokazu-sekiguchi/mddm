# Description
Code for cosmological evolution with decaying dark matter and MCMC parameter estimation based on recent cosmographic data.

# Prerequisites (vesion number usued in development)
Python3 (3.6), SWIG (3.0.12)

## Python modules
numpy (1.18.1),scipy (1.4.1), h5py (2.10.0), emcee (3.0.2), tqdm (4.43.0), getdist (1.1.0)

# Installation
Compilation is required in order to build a python interface for HyRec (written in C).
1. Git clone source file via `https://github.com/toyokazu-sekiguchi/mddm.git`.
2. Go to `mddm/HyRec`.
3. Edit `INC_PY` and `LIB_PY` in `Makefile` appropriately, so that correct `Python.h` and `libpython3.*` are incorporated. Then `Make pyrec`.
4. Go back to the parent directory. 

# Usage
There are two stages in the analysis: MCMC and 

## stage 1: MCMC (including calculation for fiducial model)
Basic usage:
`python3 driver.py test.ini`

`test.ini` specifies a variety of parameters and consists of three sections:
1. [OUTPUT]
2. [FIDUCIAL COSMOLOGY]
3. [DDM SETUP]
4. [LIKELIHOODS]
5. [MCMC]

### roles of python files:
* const.py: definition of units and constants
* mdd.py: calculation of cosmological background evolution. 

## stage 2: Postprocessing
Basic usage:
`python3 post.py dist.ini`

`dist.ini` specifies 

# Notes
## Cosmological assumptions
* Flatness is assumed.
* Neutrinos are assumed to consist of three mass eigenstates.
* 4He abundance Yp(\omega_b, N_\nu) is fitted with BBN.dat taken from CLASS which are originally obtained using the PArthENoPE code (http://parthenope.na.infn.it).
* Recombination history is computed based on HyReC (https://pages.jh.edu/~yalihai1/hyrec/hyrec.html).
* BBN calculation is 
* Restarting functionarity is for the time being not supported and will be implemented in future.
