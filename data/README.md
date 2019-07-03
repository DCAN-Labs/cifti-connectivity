# `data` folder

By default, this folder will contain all of the files created by the `cifti_conn_wrapper.py` wrapper. The number of files in this folder will depend on which scripts are run (`matrix`, `template`, or `pairwise_corr`) from the wrapper, and whether the `d`/`pconn` files are kept using the `--keep_conn_matrices` argument.

To create the output files in a different folder than this one, use a non-default value for the `--output` argument.