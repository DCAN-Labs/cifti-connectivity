# `src` folder

This folder contains all of the scripts used by the `cifti_conn_wrapper.py` wrapper. There should be at least 16 files in this folder. The first is this `README`. The other 15 are listed below:

## Files to run compiled scripts

All of these files were created by the MATLAB compiler and must be in this folder for `cifti_conn_wrapper.py` to run.

### Python Functions Source Script

1. `input_validation.py`

### Compiled MATLAB Scripts

1. `cifti_conn_matrix_for_wrapper`
1. `cifti_conn_pairwise_corr_for_wrapper`
1. `cifti_conn_template_for_wrapper`

### BASH Scripts to Run Compiled MATLAB Scripts

1. `run_cifti_conn_matrix_for_wrapper.sh`
1. `run_cifti_conn_pairwise_corr_for_wrapper.sh`
1. `run_cifti_conn_template_for_wrapper.sh`

## Other files

Although these files are not necessary for `cifti_conn_wrapper.py` to run, they are included in case changes need to be made to the MATLAB scripts, which would then need to be re-compiled.

### Uncompiled MATLAB Scripts

1. `cifti_conn_matrix_for_wrapper.m`
1. `cifti_conn_pairwise_corr_for_wrapper.m`
1. `cifti_conn_template_for_wrapper.m`

### Other MATLAB Compiler Inputs

1. `compile.sh`
1. `isthisanoutlier.m`

### Other MATLAB Compiler Outputs

1. `auto_generated_readme.txt`
1. `mccExcludedFiles.log`
1. `requiredMCRProducts.txt`
