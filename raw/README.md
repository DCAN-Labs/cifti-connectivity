# `raw` folder

By default, this folder should contain all of the raw data used by the `cifti_conn_wrapper.py` wrapper. There should be 5 files in this folder, including this `README`.

## Files belonging in this folder
1. Group motion `.conc` file
1. `dtseries` or `ptseries` time series `.conc` file
1. Left midthickness file
1. Right midthickness file

If these files are kept in a different folder than this one, that will need to be specified when running `cifti_conn_wrapper.py` by using non-default values for the `series_file`, `--motion`, `--left`, and `--right` arguments.