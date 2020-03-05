# Connectivity Matrices

Run `cifti_conn_wrapper.py` to generate correlation matrices from BOLD dense or parcellated time series data.

## Prerequisites

### Installation

Clone this repository and save it somewhere on the Linux system that you want to use it from.

### Dependencies

- [Python v3.5.2](https://www.python.org/downloads/release/python-352) or greater
- MathWorks [MATLAB Runtime Environment (MRE) version 9.1 (2016b)](https://www.mathworks.com/products/compiler/matlab-runtime.html)
- Washington University [Workbench Command (wb_command)](https://github.com/Washington-University/workbench)

### Requirements For `.conc` Files

This wrapper uses `.conc` text files for several of its inputs. All of them should only contain valid file paths. Although relative paths can be used, they will only work when running the wrapper from the right directory, so using absolute paths is recommended. Each line of each `.conc` file should have one path and nothing else.  

All `.conc` files used as input to this wrapper should have the same number of file paths, and therefore the same number of lines. That is required because the wrapper assumes that each path represents the same subject session as the path at the same line in every other input `.conc` file.

None of the time series files listed in the `.conc` files should have exactly the same base file name. If any do, then the wrapper will overwrite one of their outputs with the other one's outputs.

If you want to run this wrapper on only one subject session, you do not need to use a `.conc` file. Instead, you can use the file paths that would have gone into the `.conc` files directly. For example, you can use your `*.ptseries.nii` file path for the `series-file` argument, the `*power_2014_FD_only.mat` file for the `--motion` argument, etc.

### Expected Naming Convention for Input Imaging Data

This wrapper processes dense (`dtseries`) or parcellated (`ptseries`) time series data. For dtseries, each subject session is expected to have one file following the naming convention `XXXXX_Atlas.dtseries.nii`. Smoothing will create an additional dtseries called `XXXX_Atlas_SMOOTHED.dtseries.nii`.

## Usage

### Modes

- Run `matrix` mode on your `dtseries` or `ptseries` to generate a correlation matrix (in Fisher-Z by default) of all the greyordinates/parcellations to all other greyordinate/parcellations.

- Run `template` mode to build an average template connectivity matrix of a list of subjects. It adds all the connectivity matrices one by one, then divides by the number of subjects.

- Run `pairwise_corr` mode to make a correlation of correlation martices. This compares the connectivity matrix of each individual to the template, and provides a vector where each element is the correlation of connectivity to that greyorindate/parcellation.

`cifti_conn_wrapper.py` can run any of these modes. It requires 4 positional arguments, and accepts many optional arguments. Each mode corresponds to one compiled MATLAB file in the `./src/` directory.

### Required Positional Arguments

These arguments must be given in order:

1. `series-file` takes one path to the `.conc` file with a list of paths to each dense (`dtseries`) or parcellated (`ptseries`) timeseries file.

2. `tr` takes the repetition time (time interval between frames of a scan) for your data as a floating-point number.

3. `output` takes a path to the directory which will be filled with output data.

4. `scripts-to-run` takes one or more argument(s), the name(s) of each mode that you want to run.

For more usage information, call this script with the `--help` command: `python3 cifti_conn_wrapper.py --help`

### Examples

#### Simplest Usage

Since only the four positional argument are technically required, this is a valid command:

```sh
series_conc=./raw/group_ptseries.conc
tr=2.5
output=./data/
mode=matrix
python3.5 cifti_conn_wrapper.py ${series_conc} ${tr} ${output} ${mode}
```

However, running with no optional arguments will not do any motion correction or smoothing. It will also include subject sessions with any amount of good data.

#### Common Use Case

Here is a common use case of the `cifti_conn_wrapper` in `matrix` mode:

```sh
wb_path=/home/code/workbench/bin_linux64/wb_command
fd_thresh=0.3
mre=/home/code/MATLAB_MCR/v91
motion=raw/round1_motion.conc
python3.5 cifti_conn_wrapper.py ${series_conc} ${tr} ${output} matrix --motion ${motion} --mre-dir ${mre} --wb-command ${wb_path} --minutes 4 --fd-threshold ${fd_thresh}
```

This case includes the 4 required arguments, as well as the `fd-threshold`, `minutes`, `motion`, `mre-dir`, and `wb-command` optional arguments.

#### Multiple Modes

This wrapper can run any combination of the three modes in any order. To run multiple modes, list all of their names in the order that you want the wrapper to run them. Here is an example command which will run all three scripts in order:

```sh
python3.5 cifti_conn_wrapper.py ./raw/group_ptseries.conc 4 ./data/ matrix template pairwise_corr
```

## Options

There are three kinds of options: Server-dependent flags used by all 3 run modes, flags used by `matrix` and `template` mode, and flags used by `template` and `pairwise_corr` mode.

### Which Modes Use Which Flags

|                        | matrix | template | pairwise_corr |
|------------------------|:------:|:--------:|:-------------:|
| `--mre-dir`            | Y      | Y        | Y             |
| `--wb-command`         | Y      | Y        | Y             |
| `--additional-mask`    | Y      | Y        |               |
| `--beta8`              | Y      | Y        |               |
| `--dtseries`           | Y      | Y        |               |
| `--fd-threshold`       | Y      | Y        |               |
| `--left` and `--right` | Y      | Y        |               |
| `--make-conn-conc`     | Y      | Y        |               |
| `--minutes`            | Y      | Y        |               |
| `--motion`             | Y      | Y        |               |
| `--remove-outliers`    | Y      | Y        |               |
| `--smoothing-kernel`   | Y      | Y        |               |
| `--suppress-warnings`  | Y      | Y        |               |
| `--vector-censor`      | Y      | Y        |               |
| `--keep-conn-matrices` |        | Y        | Y             |
| `--template`           |        | Y        | Y             |

### Server-Dependent Arguments for All 3 Modes

These arguments apply to all 3 run modes. If they are excluded, then by default the wrapper will use hardcoded paths which are only valid on the RUSHMORE or Exacloud servers. So these arguments are required if this script is run on a different server or locally. However, if there is already a `wb_command` file in the user's BASH PATH variable, then the script will use that.

- `--mre-dir` takes one argument, a valid path to the MATLAB Runtime Environment directory. Example: `--mre-dir /usr/home/code/Matlab2016bRuntime/v91`

- `--wb-command` takes one argument, a valid path to the workbench command file. Example: `--wb-command /usr/local/home/wb_command`

### Optional Arguments for Modes Creating Matrix or Template

These optional arguments apply only to the first 2 run modes, `matrix` and `template.`

#### Numerical Value Flags

- `--smoothing-kernel` takes the smoothing kernel as one floating-point number. Include this argument to use smoothing, but not otherwise. If the `series_file` has `ptseries` files, then smoothing will use the `--dtseries` argument. By default, the wrapper will not do smoothing.

- `--minutes` takes the minutes limit as one floating-point number. The minutes limit is the minimum number of minutes of good data necessary for a subject session to be included in the correlation matrix. By default, the wrapper will have no minutes limit and make an `allframesbelowFDX` matrix. Subjects will have different numbers of time points that in each connectivity matrix.

- `--fd-threshold` takes the motion threshold distinguishing good from bad data. This floating-point number is the maximum amount of acceptable motion between frames of a scan. Raising the FD threshold excludes more data by setting a more stringent quality requirement. The default value is `0.2`.

#### Runtime Option Flags

These arguments can be included without a value.

- `--remove-outliers` will remove outliers from the BOLD signal. By default, frames will only be censored by the FD threshold.

- `--make-conn-conc` will make a list of connectivity matrices created by this wrapper. By default, the wrapper will not make a list. Running `pairwise_corr` mode will only work if there is already a list of connectivity matrices. Create one by running `matrix` or `template` mode with the `--make-conn-conc` flag. That will automatically build the `.conc` file name and save that file in the `output` folder.

- `--beta8` will run a beta version to reduce file size. Include this argument to reduce floating point precision and discard lower triangle of matrix. Exclude it to leave the same.  If included, this will produce 8Gb `.dconns`. Otherwise, this will make 33Gb `.dconns`. This option does nothing for `ptseries`.

- `--suppress-warnings` will prevent the wrapper from asking user for confirmation if the `.dconn` files created by the wrapper will exceed a certain threshold. By default, the wrapper will warn the user if it will create files totaling over 100 GB. This argument does nothing for `ptseries`.

- `--vector-censor` allows the user to choose the vector censoring/removal option for motion data. The default option , `combined` removal, will combine `Power_2014_FD_only` vector censoring (`fd`) with outlier vector censoring (`outlier`). The options are `outlier`, `fd`, `combined`, and `frame` removal. This argument only affects motion correction, so it does nothing unless the `--motion` flag is included.

#### Directory and File Path Flags

Each of the arguments below accepts one value, a valid file or directory path. Each argument has a default value, but the user can optionally specify a different path.

- `--motion` takes the name of a `.conc` text file listing paths pointing to FNL motion mat files for each dt or ptseries. This flag is necessary for motion correction, since by default the wrapper will not do any motion correction.

- `--left` and `--right` take the `.conc` file name of subjects' left and right midthickness files, respectively. These arguments are only needed for smoothing. If these flags are included without filenames is given, then `ADHD_DVARS_group_left`/`right_midthickness_surfaces.conc` will be used as the default name. If only a filename is given, then the wrapper will look for the file in the `--input` directory. However, this argument also accepts absolute paths.

- `--dtseries` takes the path to 1 .conc file with a list of `.dtseries.nii` file paths. If the `series-file` has a list of paths to .ptseries.nii files, then a dtseries .conc file is still needed for outlier detection and removal. If `--dtseries` is excluded, then this script will try to find a dtseries `.conc` file in the same location as the `series-file` argument.

- `--additional-mask` takes the path to an additional mask on top of the FD threshold. The mask should be a `.txt` file with `0`s and `1`s where `0`s are frames to be discarded and `1`s are frames to be used to make your matrix. This mask can be used to extract rests between blocks of task. This vector of `0`s and `1`s should be the same length as your dtseries. By default, no additional mask is used.

### Optional Arguments for Template-Creation and Pairwise-Correlation Modes

These arguments only apply to the second and last run modes, `template` and `pairwise_corr`.

- `--keep-conn-matrices` will make the wrapper keep the `dconn`/`pconn` files after creating them. Otherwise, it will delete the `d`/`pconns` after adding them to the average `d`/`pconn`. The `d`/`pconn` files are needed to run the `pairwise_corr` mode, which compares the `d`/`pconn` files to the average file created by the `template` mode. So if this flag is excluded and `pairwise_corr` is run, then the `d`/`pconn` files will be kept until `pairwise_corr` finishes.

- `--template` takes the full path to a template file created by running `template` mode. If `pairwise_corr` mode is run before `template` mode, then the `--template` file should already exist at the specified path. By default, the wrapper will build the name of this file by combining the values of the `series-file`, `--motion`, `--fd-threshold`, and `--minutes` arguments. The wrapper will then create, or look for, a file by that name in the `output` directory if needed.

## Outputs

Within the output folder, here is what the outputs will look like:

- The output of `cifti_conn_matrix` will look like `./data/XXXXXX-X_FNL_preproc_Gordon_subcortical.ptseries.nii_5_minutes_of_data_at_FD_0.2.pconn.nii`
- The output of `cifti_conn_template` will look like `./data/dtseries_AVG.dconn.nii`
- The output of `cifti_conn_pairwise_corr` will look like `./data/XXXXXX-X_FNL_preproc_Gordon_subcortical.ptseries.nii_5_minutes_of_data_at_FD_0.2.pconn.nii`

## Explanation of Process

### cifti_conn_matrix_for_wrapper

This code build a connectivity matrix from BOLD times series data. The code has the option of using motion censoring (highly encouraged), to ensure that the connectivity matrix is accurate. It takes these arguments from the wrapper: `series-file`, `time-series`, `motion-file`, `fd-threshold`, `tr`, `minutes-limit`, `smoothing-kernel`, `left`, `right`, and `beta8`.

The `time-series` argument passed to `cifti_conn_matrix` is either 'dtseries' or 'ptseries'. The wrapper infers it from the `series-file`.

To avoid conflating multiple files with the same name listed in the input `.conc`, for each connectivity matrix, this wrapper generates a random hash composed of the first character of each folder name in its `.d`/`ptseries.nii` file's path. The wrapper also appends the first 5 characters of the `.d`/`ptseries.nii` file name. This gives each output connectivity matrix a unique but consistent filename.

### cifti_conn_template_for_wrapper

This code builds a template d/pconn from a list of d/ptseries. If the d/pconn exists, if will load it instead of making it anew (calls `cifti_conn_matrix_for_wrapper`). It takes the same arguments from the wrapper as `cifti_conn_matrix_for_wrapper`, as well as `keep-conn-matrices`.

### cifti_conn_pairwise_corr_exaversion

This compares the connectivity matrix of each individual to the template (see `cifti_conn_template_for_wrapper`), and provides a vector where each element is the correlation of connectivity to that greyorindate/parcellation. It takes these arguments from the wrapper: `template`, `p`-or-`dconn`, `make-conn-conc`, and `keep-conn-matrices`.

The `p`-or-`dconn` argument passed to `cifti_conn_pairwise_corr_exaversion` is simply the `time-series` argument reformatted: `dtseries` becomes `dconn`, and `ptseries` becomes `pconn`. By default, the wrapper builds the name of the `matrices-conc` file by combining the `series-file`, `time-series`, `minutes`, and `fd-threshold` arguments.

## Known Issues

- The code added a pause function each time it computes outliers.  The pause function was added because sometimes the outlier file was read before it was done being written.
- The outlier file persists after creation (currently saved as `dtseries_temp.txt`.  This file need to be deleted as part of regular clean-up, but has not been implemeneted.

## Updates

- v1: added to make a "all-frames" dconn (using "none" as minutes limit). This will make dconns from subjects with differing amounts of time series data, but for subjects with a considerable amounts of data, this becomes less of an issue.
  - added an option to keep dconns in cifti_conn_template.m rather then automatically deleting them.

### 4/26/2019

- v1: add an option to remove outliers from the bold signal.
- Added an option to provide an additional mask in addition  to the FD mask provided.  This allows one to use periods of rest between tasks.
- Added an option to not make a list of dconns after running.  When running this code in parallel often it would created many small files.

### 7/3/2019

- Wrote Python CLI wrapper to run all three scripts, wrote documentation for the wrapper, and reorganized directory structure

### 8/23/2019

- Fixed a bug where the wrapper would not run unless connectivity matrix list `.conc` file already exited, even though the wrapper includes the functionality to make that list.

### 11/1/2019

- Fixed a bug where `cifti-conn-matrix` would treat two different files with the same name in different folders as if they were the same file, skipping over one. Also, updated README to reflect that `--remove-outliers`, `--additional-mask`, and `--make-conn-conc` were added.

### 11/13/2019

- Added the `--dtseries` parameter.
- Updated description of how the hash appended to each connectivity matrix filename is generated.

### 3/5/2020

- Added the `--vector-censor` parameter.

## Feature Requests

- [Robo-Science Trello Board](https://trello.com/b/hQkyYits/robo-science)
