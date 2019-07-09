# Connectivity Matrices 

Scripts to generate correlation matrices from BOLD dense or parcellated time series data.

## Installation

Clone this repository and save it somewhere on the Linux system that you want to use it from.

## Dependencies
- [Python v3.5.2](https://www.python.org/downloads/release/python-352) or greater
- [MathWorks MATLAB Runtime Environment (MRE) version 9.1 (2016b)](https://www.mathworks.com/products/compiler/matlab-runtime.html) or greater
  
## Purpose
- Run `cifti_conn_matrix` on your dtseries or ptseries to generate a correlation matrix (in Fischer Z by default) or all the greyordinates/parcellations to all other greyordinate/parcellations.
- Run `cifti_conn_template` to build a template connectivity matrix of a list of subjects. (It adds all the connectivity matrices one by one, then divides by the number of subjects.)
- Run `cifti_conn_pairwise_corr` to generate a correlation of correlation martices.  This compares the connectivity matrix of each individual to the template, and provides a vector where each element is the correlation of connectivity to that greyorindate/parcellation.

## Usage

Any of the three matrix, template, or pairwise_corr functions can be completed by running `cifti_conn_wrapper.py` from within the directory cloned from this repo. `cifti_conn_wrapper.py` requires two positional arguments and can take many optional arguments. 

### Required Positional Arguments
 
1. `series_file` is one argument, a path to the dense timeseries (`dtseries`) or parcellated timeseries (`ptseries`) `.conc` file. A `.conc` file is a text file with paths to each file being examined.
2. `tr` takes the repetition time (time interval between frames of a scan) for your data as a floating-point number.
3. `scripts_to_run` is one or more argument(s), the name(s) of each script that you want to run.

Those arguments must be given in that order: `series_file` first, `tr` second, and `scripts_to_run` last. For example, here is a valid basic call of this wrapper:
```
python3 cifti_conn_wrapper.py ./raw/group_ptseries.conc 2.5 cifti_conn_matrix
```

This wrapper can run any combination of the three scripts in any order. To run multiple scripts, list all of their names in the order that you want the wrapper to run them. Here is an example call which will run all three scripts in order:
```
python3 cifti_conn_wrapper.py ./raw/group_ptseries.conc 2.5 cifti_conn_matrix cifti_conn_template cifti_conn_pairwise_corr
```

For more usage information, call this script with the `--help` command: `python3 cifti_conn_wrapper.py --help`

### Server-Dependent Arguments
 
If these arguments are excluded, then by default the wrapper will use hardcoded paths which are only valid on the RUSHMORE or Exacloud servers. If this script is run on a different server or locally, then these arguments are required.

`--mre_dir` takes one argument, a valid path to the MRE compiler runtime directory. Example: `--mre_dir /usr/home/code/Matlab2016bRuntime/v91`

`--wb_command` takes one argument, a valid path to the workbench command file. Example: `--wb_command /usr/local/home/wb_command`

### Optional Arguments: Run Modes

These arguments can be included without a value.

`--keep_conn_matrices` will make the wrapper keep the `dconn/pconn` files after creating them. Otherwise, it will delete the `d/pconns` after adding them to the average `d/pconn`. The `d/pconn` files are needed to run `cifti_conn_pairwise_corr`, which compares the `d/pconn` files to the average file created by `cifti_conn_template`. So if this flag is excluded and `cifti_conn_pairwise_corr` is run, then the `d/pconn` files will be kept until the `pairwise_corr` script finishes.

`--beta8` will run a beta version to reduce file size. Include this argument to reduce floating point precision and discard lower triangle of matrix. Exclude it to leave the same.  If included, this will produce 8Gb `.dconns`. Otherwise, this will make 32Gb `.dconns`. This option does nothing for `ptseries`.

`--make_conn_conc` will make a list of connectivity matrices created by this wrapper. By default, the wrapper will not make a list.

`--remove_outliers` will remove outliers from the BOLD signal if this flag is included. Otherwise (by default), frames will only be censored by the FD threshold.

`--suppress_warnings` will prevent the wrapper from asking user for confirmation if the `.dconn` files created by the wrapper will exceed a certain threshold. By default, the wrapper will warn the user if it will create files totaling over 100 GB. This argument does nothing for `ptseries`.

### Optional Arguments: Numerical Values

`--smoothing_kernel` takes the smoothing kernel as one floating-point number. Only include this argument if smoothing will be used. Smoothing on ptseries is not supported. The default value is 'none'.

`--minutes` takes the minutes limit as one floating-point number. The minutes limit is the number of minutes to be used to generate the correlation matrix. The default minutes limit of 'none' will make an `allframesbelowFDX` `.dconn` file. Subjects will have differing numbers of time points that go into each `.dconn`.

`--fd_threshold` takes the motion threshold (maximum amount of acceptable motion between frames of a scan) for your data as a floating-point number. The default value is `0.2`.

### Optional Arguments: Directories and Files

Each of the arguments below accepts one value, a valid file or directory path. Each argument has a default value, but the user can optionally specify a different path.
 
`--input` takes a path to the directory containing the input data. The default value of this argument is `./raw/`.
 
`--output` takes a path to the directory which will be filled with output data. The default value of this argument is `./data/`.
 
`--motion` takes the name of a `.conc` file pointing to FNL motion mat files for each dt or ptseries. The default value is `./raw/ADHD_DVARS_group_motion.conc`, a list of paths to files called `Power_2014FDonly.mat`. Type `none` to use no motion file.
 
`--left` takes the `.conc` file name of subjects' left midthicknessfile. This argument is only needed for smoothing. If this flag is included but no filename is given, then `ADHD_DVARS_group_left_midthickness_surfaces.conc` will be used as the default name. If only a filename is given, then the wrapper will look for the file in the `--input` directory. However, this argument also accepts absolute paths. 

`--right` takes the `.conc` file name of subjects' right midthicknessfile. This argument is only needed for smoothing. If this flag is included but no filename is given, then `ADHD_DVARS_group_right_midthickness_surfaces.conc` will be used as the default name. If only a filename is given, then the wrapper will look for the file in the `--input` directory. However, this argument also accepts absolute paths.

`--template` takes the full path of a template file created by `cifti_conn_template`. If `cifti_conn_pairwise_corr` is run before running `cifti_conn_template`, this file should already exist at the specified path. By default, the wrapper will build the name of this file by combining the values of the `series_file`, `--motion`, `--fd_threshold`, and `--minutes` arguments. The wrapper will then either create or look for a file by that name in the `--output` directory if needed. This argument is unnecessary for `cifti_conn_matrix` to run.

`--additional_mask` takes the path to an additional mask on top of the FD threshold. The mask should be a `.mat` file with `0`s and `1`s where `0`s are frames to be discarded and `1`s are frames to be used to make your matrix. This mask can be used to extract rests between blocks of task. This vector of ones and zeros should be the same length as your dtseries. The default value is "none".

### Advanced Usage Example
```
python3 cifti_conn_wrapper.py ./raw/dtseries_file.conc 3.0 cifti_conn_template cifti_conn_pairwise_corr --beta8 --keep_conn_matrices --motion none --template <Name of template file to be made> --output <Directory where output data will be placed> --make_conn_conc
```

## Expected Naming Convention for Input Imaging Data 
For dtseries, subjects are expected to have the naming convention `XXXXX_Atlas.dtseries.nii`. Smoothing will create an additional dtseries called `XXXX_Atlas_SMOOTHED.dtseries.nii`.

## Outputs
Within the output folder, here is what the outputs will look like:
- The output of `cifti_conn_matrix` will look like `./data/XXXXXX-X_FNL_preproc_Gordon_subcortical.ptseries.nii_5_minutes_of_data_at_FD_0.2.pconn.nii`
- The output of `cifti_conn_template` will look like `dtseries_AVG.dconn.nii`
- The output of `cifti_conn_pairwise_corr` will look like `./data/XXXXXX-X_FNL_preproc_Gordon_subcortical.ptseries.nii_5_minutes_of_data_at_FD_0.2.pconn.nii`

## Explanation of Process

### cifti_conn_matrix

This code build a connectivity matrix from BOLD times series data. The code has the option of using motion censoring (highly encouraged), to ensure that the connectivity matrix is accurate. It takes these arguments from the wrapper: `series_file`, `time_series`, `motion_file`, `fd_threshold`, `tr`, `minutes_limit`, `smoothing_kernel`, `left`, `right`, and `beta8`.

The `time_series` argument passed to `cifti_conn_matrix` is either 'dtseries' or 'ptseries'. The wrapper infers it from the `series_file`.

Options for `cifti_conn_matrix` which have not yet been added to the MATLAB scripts or the wrapper include the following:
- `removeoutliers` should be either 0 or 1.  If 1, outliers are removed from the BOLD signal.  If 0 frames will only be censored by the FD threshold.
- `additionalmask` should take one file name. Addional mask on top of the FD threshold.  This mask can be used to extract rests between blocks of task. This vector of 1s and zeros should be the same length as your dtseries.
- `make_dconn_conc` should be either 0 or 1. - make a list of connectivity matrices that were created by this script.

### cifti_conn_template
This code builds a template d/pconn from a list of d/ptseries.  If the d/pconn exists, if will load it instead of making it anew (calls `cifti_conn_matrix`). It takes the same arguments from the wrapper as `cifti_conn_matrix`, as well as `keep_conn_matrices`.

### cifti_conn_pairwise_corr
This compares the connectivity matrix of each individual to the template (see cifti_conn_template), and provides a vector where each element is the correlation of connectivity to that greyorindate/parcellation. It takes these arguments from the wrapper: `template`, `p_or_dconn`, `conn_matrices`, and `keep_conn_matrices`.

The `p_or_dconn` argument passed to `cifti_conn_pairwise_corr` is simply the `time_series` argument reformatted: `dtseries` becomes `dconn`, and `ptseries` becomes `pconn`. The `conn_matrices` argument is a path to a `.conc` list of connectivity matrices to use in the analysis. The wrapper builds the name of the `conn_matrices` file by combining the `series_file`, `time_series`, `minutes`, and `fd_threshold` arguments.

## Known Issues
  - The code added a pause function each time it computes outliers.  The pause function was added because sometimes the outlier file was read before it was done being written.
  - The outlier file persists after creation (currently saved as `dtseries_temp.txt`.  This file need to be deleted as part of regular clean-up, but has not been implemeneted.

## Recent Updates
  - v1: added to make a "all-frames" dconn (using "none" as minutes limit). This will make dconns from subjects with differing amounts of time series data, but for subjects with a considerable amounts of data, this becomes less of an issue.
  	-  added an option to keep dconns in cifti_conn_template.m rather then automatically deleting them.

### 4/26/2019
- v1: add an option to remove outliers from the bold signal.
- Added an option to provide an additional mask in addition  to the FD mask provided.  This allows one to use periods of rest between tasks.
- Added an option to not make a list of dconns after running.  When running this code in parallel often it would created many small files.

### 7/3/2019
- Wrote Python CLI wrapper to run all three scripts, wrote documentation for the wrapper, and reorganized directory structure

## Feature Requests
 - https://trello.com/b/hQkyYits/robo-science
 
