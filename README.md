# Description
Code for cosmological evolution with decaying dark matter and MCMC parameter estimation based on recent cosmographic data.

# Prerequisites (vesion number usued in development)
Python3 (3.6), SWIG (3.0.12)

## Python modules
numpy (1.18.1),scipy (1.4.1), h5py (2.10.0), emcee (3.0.2), tqdm (4.43.0), getdist (1.1.0)

# Installation
1. Git clone source file via `https://github.com/toyokazu-sekiguchi/mddm.git`.
2. Go to `mddm/HyRec`.
3. Edit `INC_PY` and `LIB_PY` in `Makefile` appropriately, so that correct `Python.h` and `libpython3.*` are incorporated. Then `Make pyrec`.
4. Go back to the parent directory. 

# Usage
There are two stages in the analysis: MCMC and 

## stage 1: MCMC (including calculation for fiducial model)
Basic usage:
`python3 driver.py test.ini`

`test.ini` includes a variety of parameters and consists of three sections: 


# Notes
* Recombination history is computed based on HyReC (Al\"imoud & Hirata) https://pages.jh.edu/~yalihai1/hyrec/hyrec.html
* Restarting functionarity is for the time being not supported and will be implemented in future.
