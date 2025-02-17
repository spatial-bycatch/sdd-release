# Integrated spatial density estimation -- model release

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