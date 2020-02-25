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
  - `root`: Chacters specifying output prefix.
* [FIDUCIAL COSMOLOGY]
  - `ob`, `odm`, `ol`: Density parameters $\Omega_i h^2$ of baryon, dark matter (assuming no decay) and cosmological constant. Note that actual value of $\Omega_{dm} h^2$ and henceforth $h = \sqrt{\sum_i \Omega_i h^2}$ differs from the input value when decay rate is finite.
  - `decay_rate`: Decay rate of decaying dark matter in units of 1/Gyr.
  - `mratio`: Mass ratio of massive decay-product to decaying particle $m_1/m_0$.
  - `nnu`: Effective number of neutrinos. The total number of neutrinos are enhanced by this factor (temperature is fixed to the standard value i.e. $T_\nu = (4/11)^{1/3} T_\gamma$.
  - `mnu`: Sum of neutrino mass in units of eV.
  - `neutrino_hierarchy`: Flag for neutrino mass hierarchy. 1 for normal, 0 for degenerate and -1 for inverted ones.
* [DDM SETUP]
Background evolution is computed based on iteration method using the following parameters.
  - `num_a`: Number of scale factor bin. Empirically `30` is recommended.  
  - `max_it`: Maximum iteration number. This is not relevant because usually convergence is achieved within one iteration.
  - `tol`: Tolerance parameter for the error from the true evolution measured by the current energy density of massless decay-product. Setting to `1e-10` works.
* [LIKELIHOODS]
  - `use_BAO`,`use_H0`, `use_CMB`: Flags for whether data is incorporated in likelihood calculation. They should be either `true` of `false`
* [MCMC]
  - `ob`, `odm`, `ol`, `decay_rate`, `mratio`, `nnu`, `mnu`: When each parameter is varied in the parameter estimation, three numbers should be given in order: lower limit, upper limit, initial fluctuations. When left as blank, corresponding parameter is fixed to the fiducial value.
  - `nwalkers`: Number of walkers in affine invariant MCMC sampler. This should be at least twice the number of varied parameters.
  - `nsteps`: Number of steps for MCMC analysis
  - `parallel`: If `true`, parallelization is implemented in the MCMC calculation.

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
