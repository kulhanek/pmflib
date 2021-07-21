# PMFLib - A Toolkit for Free Energy Calculations
PMFLib implements free energy calculation methods employing biased molecular dynamics (MD) simulations. The following methods are available:
* Adaptive Biasing Force (ABF) and combination of Umbrella Sampling and Adaptive Biasing Force (US-ABF)
* Constrained Dynamics (CST), also known as Blue Moon
* Direct and Well-tempered Metadynamics (MTD)
* Adaptive Biasing Potential (ABP)

Some methods provide mean forces, which need to be integrated to get free energies. PMFLib offers the following integration techniques:
* Reverse Finite Differences (RDF)
* Radial Basis Functions (RBF)
* Gaussian Process Regression (GPR) with hyperparameter optimization

In addition to the free energy, some methods can provide enthalpy and entropy from the same biased MD simulation.

Supplementary  methods include:
* Monitoring of collective variables
* Restrained Dynamics

Sampling can be improved by Multiple-walker Approach (MWA), which uses server/client architecture with fully asynchronous communications.
MWA is available for ABF, ABP, and MTD.

Further details: [https://pmflib.ncbr.muni.cz](https://pmflib.ncbr.muni.cz)

## Building Package
This package can be build:
* [independently](https://github.com/kulhanek/pmflib-build)
* as an extension of [sander](https://github.com/kulhanek/sander-pmf-build)
