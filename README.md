[![DOI](https://zenodo.org/badge/878002821.svg)](https://doi.org/10.5281/zenodo.15781961)

# Integrated estimation of the spatial population density surface using semi-continuous sampling data
Charles Edwards, Nokuthaba Sibanda and Marie-Julie Roux

Methods in Ecology and Evolution (2025)

### This repository contains the data and code necessary to reproduce the results of the above publication.

Raw data are found in the `data` directory as comma delimited files.

Code and formatted input data, including initial values, are found in the `code` directory. Within the `code` directory are subdirectories for each species and model run. With reference to Table S1 in the Supplementary Information, the following models are provided:

Exponential models:

E1: `exp`

E2: `exp_reg`

E3: `exp_hyb`

E4: `exp_hyb_reg`

Weibull models:

W1: `wei`

W2: `wei_reg`

W3: `wei_hyb`

W4: `wei_hyb_reg`

Each model is given the suffix `sdd_*` for the survey only model, and `idd_*` for the integrated model. Makefiles are provided for each model run. 

Output tables and figures can be generated using the Makefile in the `code/scripts` directory.
