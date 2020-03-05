#! /usr/bin/env python3

"""
CIFTI connectivity wrapper
Greg Conan: conan@ohsu.edu
Created 2019-06-18
Updated 2020-03-05
"""

##################################
#
# Command-line wrapper for CIFTI connectivity matrix scripts.
# This wrapper can run any one or multiple of the three following scripts:
# 1. Run `cifti_conn_matrix` on your dtseries or ptseries to generate a
#    correlation matrix or all the greyordinates/parcellations to all other
#    greyordinate/parcellations.
# 2. Run `cifti_conn_template` to build a template connectivity matrix of a
#    list of subjects.
# 3. Run `cifti_conn_pairwise_corr` to generate a correlation of correlation
#    matrices.
# The functions which call these scripts are at the end of this file.
#
##################################

### Imports and functionality for defining constants

# Standard imports
import argparse
from datetime import datetime
import glob
import os
import socket
import subprocess
import sys

# Get path to directory containing this file
try:
    PWD = os.path.dirname(os.path.abspath(sys.argv[0]))
except OSError as e:
    sys.exit("{} {} cannot find the path to the directory it is in."
                .format(e, sys.argv[0]))

# Function to get default input (raw) and wrapped scripts (src) directories
def local_path_to(parent, path):
    """
    :param parent: String naming the directory containing the file at path
    :param path: String naming the relative path to a file
    :return: String which is an absolute path to {pwd}/parent/path
    """
    return os.path.join(PWD, parent, path)

# Verify that file with input validation functions can be accessed
input_validation_filepath = local_path_to("src", "input_validation.py")
if not os.access(input_validation_filepath, os.R_OK):
    sys.exit("Cannot access necessary input validation functions from {}"
             .format(input_validation_filepath))

# Local custom imports
from src.input_validation import (CHOICES_TO_RUN,
                                  none_or_valid_float_value_as_string,
                                  valid_output_dir,
                                  valid_readable_dir,
                                  valid_readable_file,
                                  valid_readable_file_or_none,
                                  validate_cli_args,
                                  WARNING_IF_DCONN_SIZE_EXCEEDS,
                                  will_delete_conn_matrices_later)

### Constants

# Default names of data files used by wrapped scripts
midthicknessfile = "ADHD_DVARS_group_{}_midthickness_surfaces.conc"
GROUP_LEFT_CONC = local_path_to("raw", midthicknessfile.format("left"))
GROUP_MOTION_CONC = local_path_to("raw", "ADHD_DVARS_group_motion.conc")
GROUP_RIGHT_CONC = local_path_to("raw", midthicknessfile.format("right"))

# Hardcoded parts of default directory paths for MRE dir and wb_command
EXA_DIR = "/home/exacloud/lustre1/fnl_lab/"
UTIL_DIR = "code/external/utilities/"
MRE_DIR = "Matlab2016bRuntime/v91"

### Functions


def main():
    cli_args = get_cli_args()

    def timestamp_when(step, completion):
        """
        Get the time and date when a step of this script started/finished
        :param step: String naming the step which either started or finished
        :param completion: String which is either 'started' or 'finished'
        :return: String with the date and time when step started or finished
        """
        return("CIFTI connectivity matrix wrapper {} running {} script at {}"
               .format(completion, step, datetime.now().strftime(
                           "%H:%M:%S on %b %d, %Y"
                       )))

    # Run all scripts the user said to run, in the order the user gave them,
    # and print when each script started
    timestamps = []
    for script in cli_args.scripts:
        timestamps.append(timestamp_when(script, "started"))
        print(timestamps[-1])

        # Call script
        globals()["cifti_conn_" + script](cli_args)

    # Tell user how long each script took to run by printing all timestamps
    timestamps.append(timestamp_when(script, "finished"))
    print("\n".join(timestamps))


def get_cli_args():
    """
    Get and validate all args from command line using argparse.
    :return: Namespace containing all validated inputted command line arguments
    """
    # Options for --vector-censor parameter
    censoring_options = ("combined", "fd", "outlier", "frame")

    # Create arg parser
    parser = argparse.ArgumentParser(
        description="Wrapper to generate connectivity matrices."
    )

    # Required arguments: Time series filename, repetition time of scan,
    # script(s) to run, and output directory
    parser.add_argument(
        "series_file",
        type=valid_readable_file,
        help=("Name of dense or parcellated timeseries .conc file (i.e. text "
              "file with paths to each file being examined).")
    )
    parser.add_argument(
        "tr",
        type=none_or_valid_float_value_as_string,
        help=("Specify the repetition time (time interval between frames of a "
              "scan) for your data.")
    )
    parser.add_argument(
        "output",
        type=valid_output_dir,
        help=("Location to save all output files to. This must be a valid "
              "path to a directory, but if the directory does not exist yet, "
              "then this wrapper will create it at that path.")
    )
    parser.add_argument(
        "scripts",
        nargs="+",
        choices=CHOICES_TO_RUN,
        help=("Select which script(s) this wrapper should run. It will run each "
              "script in the order that you enter them. Script options: "
              + ", ".join(CHOICES_TO_RUN))
    )

    # Optional: Get name of additional mask to put on top of FD threshold
    parser.add_argument(
        "-mask",
        "--additional-mask",
        type=valid_readable_file_or_none,
        default="none",
        dest="mask",
        help=("Additional mask on top of the FD threshold. The mask should be "
              "a .mat file with 0s and 1s where 0 are frames to be discarded "
              "and 1s are frames to be used to make your matrix. This mask "
              "can be used to extract rests between blocks of task. This "
              "vector of ones and zeros should be the same length as your "
              "dtseries.")
    )

    # Option: Specify whether to use beta file size reduction for dconns
    parser.add_argument(
        "-b",
        "--beta8",
        action="store_const",
        const="1",
        default="0",
        help=("Beta version to reduce file size. Include this flag to reduce "
              "floating point precision and discard lower triangle of matrix. "
              "Exclude it to leave the same.  If included, this will produce "
              "8Gb .dconns. Otherwise, this will make 33Gb .dconns. This "
              "option does nothing for ptseries.")
    )

    # Optional: Select vector censoring/removal mode
    parser.add_argument(
        "-censor",
        "--vector-censor",
        choices=censoring_options,
        default=censoring_options[0],
        help=("Choose vector censoring/removal option for motion data. The "
              "default option, {}, will combine Power_2014_FD_only vector "
              "censoring ('fd') with outlier vector censoring ('outlier'). The"
              " options are {}. This argument only affects motion correction."
              .format(censoring_options[0], censoring_options))
    )

    # Optional: Specify whether to make list of conn matrices made by wrapper
    parser.add_argument(
        "-conc",
        "--make-conn-conc",
        action="store_const",
        const="1",
        default="0",
        help=("Make a list of connectivity matrices created by this wrapper. "
              "By default, the wrapper will not make a list.")
    )

    # Optional: .dtseries file for outlier detection
    parser.add_argument(
        "-dt",
        "--dtseries",
        type=valid_readable_file,
        help=("Path to a .conc file with a list of .dtseries.nii file paths. "
              "If the series_file has a list of paths to .ptseries.nii files, "
              "then a dtseries .conc file is still needed for outlier "
              "detection and removal. If this argument is excluded and "
              "outlier detection/removal will be done, then this script will "
              "try to find a .dtseries.nii file in the same location as the "
              "series_file argument.")
    )

    # Optional: Get frame displacement motion threshold
    default_fd_threshold = "0.2"
    parser.add_argument(
        "-fd",
        "--fd-threshold",
        type=none_or_valid_float_value_as_string,
        default=default_fd_threshold,
        dest="fd",
        help=("Specify motion threshold (maximum amount of acceptable motion "
              "between frames of a scan) for your data. Default value is "
              + default_fd_threshold)
    )

    # Optional: Specify whether to keep or to delete dconn/pconn files after
    # creating them.
    parser.add_argument(
        "-keep",
        "--keep-conn-matrices",
        action="store_const",
        const="1",
        default="0",
        help=("If this flag is included, then the wrapper will keep the "
              "dconn/pconn after creating them. Otherwise, it will delete the "
              "d/pconns after adding them to the average d/pconn.")
    )

    # Optional: Get name of left surface .conc file
    parser.add_argument(
        "-l",
        "--left",
        type=valid_readable_file_or_none,
        nargs="?",
        const=GROUP_LEFT_CONC,
        default="none",
        help=(".conc file name of subjects' left midthicknessfile. Include "
              "this flag if smoothing will be run. If this flag is included"
              "but no filename is given, then {} will be used as the default "
              "name.".format(GROUP_LEFT_CONC))
    )

    # Optional: Get minutes limit
    parser.add_argument(
        "-min",
        "--minutes",
        default="none",
        type=none_or_valid_float_value_as_string,
        help=("Specify the number of minutes to be used to generate the "
              "correlation matrix. The default minutes limit of 'none' "
              "will make an 'allframesbelowFDX' .dconn. Subjects will have "
              "differing numbers of time points that go into each .dconn")
    )

    # Optional: Get name of .conc file listing FNL motion mat files
    parser.add_argument(
        "-m",
        "--motion",
        type=valid_readable_file_or_none,
        nargs="?",
        const=GROUP_MOTION_CONC,
        default="none",
        help=("Path to .conc file that points to FNL motion mat files for each"
              "dt or ptseries. By default motion correction will not be done. "
              "If this flag is included without a filename, then motion "
              "correction will be done using the default motion file(s) at "
              "{}. To use your own motion file, type it after this flag as an "
              "argument.".format(GROUP_MOTION_CONC))
    )

    # Optional: Give path to MATLAB Runtime Environment (MRE) directory
    parser.add_argument(
        "-mre",
        "--mre-dir",
        type=valid_readable_dir,
        default=get_default_mre_dir(),
        help=("Path to directory containing MATLAB Runtime Environment (MRE)"
              "version 9.1. This is used to run compiled MATLAB scripts. This "
              "argument must be a valid path to an existing folder.")
    )

    # Optional: Specify whether to remove outliers from BOLD signal
    parser.add_argument(
        "-outliers",
        "--remove-outliers",
        action="store_const",
        const="1",
        default="0",
        help=("If this flag is included, then outliers will be removed from "
              "the BOLD signal. Otherwise, frames will only be censored by "
              "the FD threshold.")
    )

    # Optional: Get name of right surface .conc file
    parser.add_argument(
        "-r",
        "--right",
        type=valid_readable_file_or_none,
        nargs="?",
        const=GROUP_RIGHT_CONC,
        default="none",
        help=("Path to subjects' right midthicknessfile .conc file. Include "
              "this flag if smoothing will be run. If this flag is included "
              "but no filename is given, then {} will be used as the default "
              "name.".format(GROUP_RIGHT_CONC))
    )

    # Optional: Get smoothing kernel for dconns
    parser.add_argument(
        "-smooth",
        "--smoothing-kernel",
        default="none",
        type=none_or_valid_float_value_as_string,
        help=("Specify smoothing kernel. If this flag is excluded, then no "
              "smoothing kernel will be used. Smoothing on ptseries is not "
              "supported.")
    )

    # Optional: Get name of template file to create
    parser.add_argument(
        "-t",
        "--template",
        type=valid_readable_file,
        help="Path to the template file to be created."
    )

    # Optional: Ignore warnings about output .dconn file sizes
    parser.add_argument(
        "-warn-off",
        "--suppress-warnings",
        action="store_true",
        help=("By default, the wrapper will ask user for confirmation if "
              "the .dconn files created by the wrapper will exceed {} GB. "
              "Include this argument to ignore this warning. Irrelevant for "
              "ptseries.".format(WARNING_IF_DCONN_SIZE_EXCEEDS))
    )

    # Optional: Specify path to wb_command
    parser.add_argument(
        "-wb",
        "--wb-command",
        type=valid_readable_file,
        default=get_default_wb_command(),
        help="Path to workbench command file called 'wb_command'."
    )

    return validate_cli_args(parser.parse_args(), parser)


def get_default_mre_dir():
    """
    :return: String which is the absolute path to the MRE directory if on the 
             Exacloud or Rushmore servers
    """
    return {"rushmore": "/mnt/max/shared/" + UTIL_DIR + MRE_DIR,
            "exacloud": EXA_DIR + UTIL_DIR + MRE_DIR,
            "other": None}[get_server_host()]


def get_server_host():
    """
    :return: String naming the server: "rushmore", "exacloud", or "other".
    """
    host = socket.gethostname()
    if "exa" in host:
        host = "exacloud"
    elif host != "rushmore":
        host = "other"
    return host


def get_default_wb_command():
    """
    :return: String, absolute path to wb_command if on Exacloud or Rushmore
    """
    return {"rushmore": "/mnt/max/software/workbench/bin_linux64/wb_command",
            "exacloud": (EXA_DIR + UTIL_DIR +
                         "workbench-1.3.2/bin_rh_linux64/wb_command"),
            "other": None}[get_server_host()]


def get_conc_file_paths(cli_args):
    """
    Build an incomplete path to .conc file(s) from input parameters, then use
    it to return a list of path(s) to either .conc file(s) or a matrix
    :param cli_args: argparse namespace with all command-line arguments
    :return: List of strings which are paths to .conc files or matrix files
    """

    def get_paths(cli_args, fd_or_none):
        """
        Build incomplete path depending on whether .conc file was named using
        the FD threshold or the absent motion file (hence the string 'none')
        """
        return os.path.join(cli_args.output, "{}_{}conn_of_{}*{}.conc".format(
            os.path.basename(cli_args.series_file), cli_args.time_series[0], 
            cli_args.time_series, fd_or_none
        ))

    if os.path.splitext(cli_args.series_file)[1] != ".conc":
        paths_list = [cli_args.series_file]

    else:  # Validate the path to .conc file(s) to be returned
        paths = get_paths(cli_args, cli_args.fd)
        if not next(glob.iglob(paths), None):
            paths = get_paths(cli_args, "none")
            if not next(glob.iglob(paths), None):
                print("No .conc file(s) at {}".format(paths))
        paths_list = glob.glob(paths)
    return paths_list


def get_matrix_or_template_parameters(cli_args):
    """
    cifti_conn_matrix and cifti_conn_template both have the same required
    parameters, with only a few exceptions. This function returns a list of
    all parameters required by both scripts.
    :param cli_args: Full argparse namespace containing all CLI args, including
                     all parameters used in both matrix and template mode.
    :return: A list of all parameters required by matrix and template scripts.
    """
    return [cli_args.mre_dir,
            cli_args.wb_command,
            cli_args.series_file,
            cli_args.time_series,
            cli_args.motion,
            cli_args.fd,
            cli_args.tr,
            cli_args.minutes,
            cli_args.smoothing_kernel,
            cli_args.left,
            cli_args.right,
            cli_args.beta8,
            cli_args.remove_outliers,
            cli_args.mask,
            cli_args.make_conn_conc,
            cli_args.output,
            cli_args.dtseries,
            cli_args.vector_censor]


def get_script_filename(index):
    """
    :param index: Integer showing which script named in CHOICES_TO_RUN to use:
                  0 for matrix, 1 for template, and 2 for pairwise_corr
    :return: String naming the file to run that script
    """
    return "run_cifti_conn_{}_for_wrapper.sh".format(CHOICES_TO_RUN[index])


def cifti_conn_matrix(cli_args):
    """
    Run cifti_conn_matrix script to generate a correlation matrix of all
    greyordinates/parcellations to all greyordinates/parcellations
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    subprocess.check_call((local_path_to("src", get_script_filename(0)),
                           *get_matrix_or_template_parameters(cli_args)))


def cifti_conn_template(cli_args):
    """
    Run cifti_conn_template script to build a template connectivity matrix of
    a list of subjects.
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    # If user is running pairwise_corr, keep conn matrices until that finishes
    keep_conn_matrices = cli_args.keep_conn_matrices
    if will_delete_conn_matrices_later(cli_args):
        print("Warning: CIFTI connectivity matrix files will not be deleted "
              "until {} finishes running.".format(CHOICES_TO_RUN[2]))
        keep_conn_matrices = "1"

    # Call cifti_conn_template script
    subprocess.check_call((local_path_to("src", get_script_filename(1)),
                           *get_matrix_or_template_parameters(cli_args),
                           keep_conn_matrices, cli_args.template))


def has_only_conn_paths(file_path):
    """
    :param file_path: String which is a valid path to a readable file
    :return: True if file_path points to a file with only .conn.nii file paths,
             otherwise False
    """
    with open(file_path) as infile:
        all_conn_paths = True
        next_line = infile.readline()
        while all_conn_paths and next_line:
            if next_line.strip()[-8:] != "conn.nii":
                all_conn_paths = False
            next_line = infile.readline()
    return all_conn_paths


def cifti_conn_pairwise_corr(cli_args):
    """
    Run cifti_conn_pairwise_corr script to generate a correlation of
    correlation matrices.
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    # Run pairwise_corr on each .conc file listing conn matrices made by
    # by a previous step of the script
    for conc_file in ([cli_args.series_file] if
                      has_only_conn_paths(cli_args.series_file) else
                      get_conc_file_paths(cli_args)):

        # Call cifti_conn_pairwise_corr script
        subprocess.check_call((local_path_to("src", get_script_filename(2)),
                               cli_args.mre_dir,
                               cli_args.wb_command,
                               cli_args.template,
                               cli_args.time_series[0] + "conn",
                               conc_file,
                               cli_args.keep_conn_matrices,
                               cli_args.output))

        # If user said not to keep .conc file but pairwise_corr used it,
        # then try to delete it after pairwise_corr finishes
        if will_delete_conn_matrices_later(cli_args):
            try:
                os.remove(conc_file)
            except FileNotFoundError as e:
                print("Could not find and remove connectivity matrix file at "
                      + e.filename)


if __name__ == '__main__':
    main()
