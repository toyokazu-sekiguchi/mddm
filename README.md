# Description
Code for cosmological evolution with decaying dark matter and MCMC parameter estimation based on recent cosmographic data.

# Prerequisites (vesion number usued in development)
* Python3 (3.6, 3.7)
* SWIG (3.0.12)
* C compiler (gcc Apple LLVM version 10.0.1/Intel icc 18.0.1)

## Python modules
* numpy (1.18.1)
* scipy (1.4.1)
* h5py (2.10.0)
* emcee (3.0.2)
* tqdm (4.43.0)
* getdist (1.1.0)

# Installation
Compilation is required in order to build a python interface for HyRec (written in C).
1. Git clone source file via `https://github.com/toyokazu-sekiguchi/mddm.git`.
2. Go to `mddm/HyRec`.
3. Edit `INC_PY` and `LIB_PY` in `Makefile` appropriately, so that correct `Python.h` and `libpython3.*` are incorporated. Then `Make pyrec`.
4. Go back to the parent directory. 

# Usage
There are two stages in the analysis: MCMC and postprocessing.

## stage 1: MCMC (including calculation for fiducial model)
Basic usage:
`python3 driver.py test.ini`

`test.ini` specifies a variety of parameters and consists of five sections:
* [OUTPUT]
 - `root`: chacters specifying output prefix.
* [FIDUCIAL COSMOLOGY]
 - `ob`, `odm`, `ol`: density parameters $\Omega_i h^2$ of baryon, dark matter (assuming no decay) and cosmological constant. Note that actual value of $\Omega_{dm} h^2$ and henceforth $h = \sqrt{\sum_i \Omega_i h^2}$ differs from the input value when decay rate is finite.
 - `decay_rate`: decay rate of decaying dark matter in units of 1/Gyr.
 - `mratio`: mass ratio of massive decay-product to decaying particle $m_1/m_0$.
 - `nnu`: effective number of neutrinos. The total number of neutrinos are enhanced by this factor (temperature is fixed to the standard value i.e. $T_\nu = (4/11)^{1/3} T_\gamma$.
 - `mnu`: sum of neutrino mass in units of eV.
 - `neutrino_hierarchy`: flag for neutrino mass hierarchy. 1 for normal, 0 for degenerate and -1 for inverted ones.
* [DDM SETUP]
 -
* [LIKELIHOODS]
 - a
 - `ob`, `odm`, `ol`: density parameters $\Omega_ih^2$ of baryon, dark matter and cosmological constant.a
 - `ob`, `odm`, `ol`: density parameters $\Omega_ih^2$ of baryon, dark matter and cosmological constant.a
 - `ob`, `odm`, `ol`: density parameters $\Omega_ih^2$ of baryon, dark matter and cosmological constant.
* [MCMC]



### roles of python files:
* const.py: definition of units and constants
* mdd.py: calculation of cosmological background evolution. 
* likelihoods.py: calculation of likelihood function incorporating recent BAO (arXiv:1607.03155, arXiv:1801.03062, arXiv:1702.00176), direct Hubble measurement (arXiv:2001.03624) and CMB \Theta_* (arXiv:1807.06209).
* mcmc.py: MCMC analysis based on Affine Invariant MCMC sampler (emcee). Parallelization is supported based on the multiprocessing python module.
* driver.py: main function

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
