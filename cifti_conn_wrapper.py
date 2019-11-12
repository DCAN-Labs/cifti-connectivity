#! /usr/bin/env python3

"""
CIFTI connectivity wrapper
Greg Conan: conan@ohsu.edu
Created 2019-06-18
Last Updated 2019-11-12
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

import argparse
from datetime import datetime
from glob import iglob
import os
from socket import gethostname
import subprocess
import sys

### Constants

# Names of scripts to run from this one
CHOICES_TO_RUN = ["cifti_conn_matrix", "cifti_conn_template",
                  "cifti_conn_pairwise_corr"]

# Default directories for input data, output data, and wrapped scripts (src)
try:
    PWD = os.path.dirname(os.path.abspath(__file__))
    assert os.access(os.path.join(PWD, "cifti_conn_wrapper.py"), os.R_OK)
except (OSError, AssertionError):
    PWD = os.getcwd()
DEFAULT_INPUT = os.path.join(PWD, "raw")
DEFAULT_OUTPUT = os.path.join(PWD, "data")
DEFAULT_SOURCE = os.path.join(PWD, "src")

# Default names of data files used by wrapped scripts
GROUP_LEFT_CONC = os.path.join(DEFAULT_INPUT,
                      "ADHD_DVARS_group_left_midthickness_surfaces.conc")
GROUP_MOTION_CONC = os.path.join(DEFAULT_INPUT, "ADHD_DVARS_group_motion.conc")
GROUP_RIGHT_CONC = os.path.join(DEFAULT_INPUT,
                      "ADHD_DVARS_group_right_midthickness_surfaces.conc")

# MATLAB Runtime Environment (MRE) directories which depend on the host server
MRE_EXACLOUD = ("/home/exacloud/lustre1/fnl_lab/code/external/utilities/"
                "Matlab2016bRuntime/v91")
MRE_RUSHMORE = ("/mnt/max/shared/code/external/utilities/Matlab2016bRuntime/"
                "v91")

# Default names of wrapped scripts to run
SCRIPT_MATRIX = "run_cifti_conn_matrix_for_wrapper.sh"
SCRIPT_PAIRWISE_CORR = "run_cifti_conn_pairwise_corr_exaversion.sh"
SCRIPT_TEMPLATE = "run_cifti_conn_template_for_wrapper.sh"

# If .dconn files' total size exceeds this many gigabytes, then warn the user
WARNING_IF_DCONN_SIZE_EXCEEDS = 100

# Workbench commands which depend on the host server
WB_EXACLOUD = ("/home/exacloud/lustre1/fnl_lab/code/external/utilities/"
               "workbench-1.3.2/bin_rh_linux64/wb_command")
WB_RUSHMORE = "/mnt/max/software/workbench/bin_linux64/wb_command"


### Functions


def main():

    cli_args = get_cli_args()

    def print_timestamp_when(completion, step):
        """
        Print message with the date and time when step started or finished
        :param completion: String which is either 'started' or 'finished'
        :param step: String naming the step which either started or finished
        :return: N/A
        """
        print("CIFTI connectivity matrix wrapper {} running {} at {}"
              .format(completion, step, datetime.now().strftime(
                  "%H:%M:%S on %b %d, %Y")))

    # Run all scripts the user said to run, in the order the user gave them
    for script in cli_args.scripts_to_run:
        print_timestamp_when("started", script)

        # Call script
        globals()[script](cli_args)

        print_timestamp_when("finished", script)


def get_cli_args():
    """
    Get and validate all args from command line using argparse.
    :return: Namespace containing all validated inputted command line arguments
    """
    # Create arg parser
    parser = argparse.ArgumentParser(
        description="Wrapper to generate connectivity matrices."
    )

    # Required arguments: Time series filename, repetition time of scan, and
    # script(s) to run
    parser.add_argument(
        "series_file",
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
        "scripts_to_run",
        nargs="+",
        choices=CHOICES_TO_RUN,
        help=("Select which script(s) this wrapper should run. It will run each "
              "script in the order that you enter them. Script options: "
              + ", ".join(CHOICES_TO_RUN))
    )

    # Optional: Get name of additional mask to put on top of FD threshold
    parser.add_argument(
        "-a",
        "--additional_mask",
        default="none",
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

    # Optional: Specify whether to make list of conn matrices made by wrapper
    parser.add_argument(
        "-c",
        "--make_conn_conc",
        action="store_const",
        const="1",
        default="0",
        help=("Make a list of connectivity matrices created by this wrapper. "
              "By default, the wrapper will not make a list.")
    )

    # Optional: .dtseries file for outlier detection
    parser.add_argument(
        "-d",
        "--dtseries",
        help=("Path to a .conc file with a list of .dtseries.nii file paths. "
              "If the series_file has a list of paths to .ptseries.nii files, "
              "then a dtseries .conc file is still needed for outlier "
              "detection and removal. If this argument is excluded, then this "
              "script will try to find a .dtseries.nii file in the same "
              "location as the series_file argument.")
    )

    # Optional: Give path to MATLAB Runtime Environment (MRE) directory
    parser.add_argument(
        "-e",
        "--mre_dir",
        help=("Path to directory containing MATLAB Runtime Environment (MRE)"
              "version 9.1 or newer. This is used to run compiled MATLAB "
              "scripts. This argument must be a valid path to an existing "
              "folder.")
    )

    # Optional: Get frame displacement motion threshold
    default_fd_threshold = "0.2"
    parser.add_argument(
        "-f",
        "--fd_threshold",
        default=default_fd_threshold,
        type=none_or_valid_float_value_as_string,
        help=("Specify motion threshold (maximum amount of acceptable motion "
              "between frames of a scan) for your data. Default value is "
              + default_fd_threshold)
    )

    # Optional: Get path to directory containing HCP CIFTI conn scripts
    parser.add_argument(
        "-i",
        "--input",
        help=("Directory containing all HCP conn data files. By default, this "
              "will be the directory containing the time series file.")
    )

    # Optional: Specify whether to keep or to delete dconn/pconn files after
    # creating them.
    parser.add_argument(
        "-k",
        "--keep_conn_matrices",
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
        nargs="?",
        const=GROUP_LEFT_CONC,
        default="none",
        help=(".conc file name of subjects' left midthicknessfile. Include "
              "this flag if smoothing will be run. If this flag is included"
              "but no filename is given, then {} will be used as the default "
              "name.".format(GROUP_LEFT_CONC))
    )

    # Optional: Get name of .conc file listing FNL motion mat files
    parser.add_argument(
        "-m",
        "--motion",
        nargs="?",
        const=GROUP_MOTION_CONC,
        default="none",
        help=("Name of .conc file that points to FNL motion mat files for each"
              "dt or ptseries. By default motion correction will not be done. "
              "If this flag is included without a filename, then motion "
              "correction will be done using the default motion file(s) at {}."
              ". To use your own motion file, type it after this flag as an "
              "argument.".format(GROUP_MOTION_CONC))
    )

    # Optional: Specify output folder
    parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_OUTPUT,
        help=("Location to save all output files to. By default, this will be "
              + DEFAULT_OUTPUT)
    )

    # Optional: Get minutes limit
    parser.add_argument(
        "-p",
        "--minutes",
        default="none",
        type=none_or_valid_float_value_as_string,
        help=("Specify the number of minutes to be used to generate the "
              "correlation matrix. The default minutes limit of 'none' "
              "will make an 'allframesbelowFDX' .dconn. Subjects will have "
              "differing numbers of time points that go into each .dconn")
    )

    # Optional: Get name of right surface .conc file
    parser.add_argument(
        "-r",
        "--right",
        nargs="?",
        const=GROUP_RIGHT_CONC,
        default="none",
        help=(".conc file name of subjects' right midthicknessfile. Include "
              "this flag if smoothing will be run. If this flag is included "
              "but no filename is given, then {} will be used as the default "
              "name.".format(GROUP_RIGHT_CONC))
    )

    # Optional: Get smoothing kernel for dconns
    parser.add_argument(
        "-s",
        "--smoothing_kernel",
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
        help="File path and name of template file to be created."
    )

    # Optional: Ignore warnings about output .dconn file sizes
    parser.add_argument(
        "-u",
        "--suppress_warnings",
        action="store_true",
        help=("By default, the wrapper will ask user for confirmation if "
              "the .dconn files created by the wrapper will exceed {} GB. "
              "Include this argument to ignore this warning. Irrelevant for "
              "ptseries.".format(WARNING_IF_DCONN_SIZE_EXCEEDS))
    )

    # Optional: Specify whether to remove outliers from BOLD signal
    parser.add_argument(
        "-v",
        "--remove_outliers",
        action="store_const",
        const="1",
        default="0",
        help=("If this flag is included, then outliers will be removed from "
              "the BOLD signal. Otherwise, frames will only be censored by "
              "the FD threshold.")
    )

    # Optional: Specify path to wb_command
    parser.add_argument(
        "-w",
        "--wb_command",
        help="Path to workbench command file called 'wb_command'."
    )

    return validate_cli_args(parser.parse_args(), parser)


def none_or_valid_float_value_as_string(str_to_check):
    """
    Unless a string is "none", tries to convert it to a float and back to check
    that it represents a valid float value. Throws ValueError if type
    conversion fails. This function is only needed because the MATLAB scripts
    take some arguments either as a float value in string form or as "none".
    :param str_to_check: string to validate
    :return: string which is either "none" or represents a valid float value
    """
    if str_to_check != "none":
        str_to_check = str(float(str_to_check))
    return str_to_check


def validate_cli_args(cli_args, parser):
    """
    Check that all command line arguments will allow this script to work.
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: Validated command-line arguments argparse namespace
    """
    # Validate input and output directories
    cli_args.input = validate_and_get_input_dir(cli_args, parser)
    cli_args.output = validate_and_get_output_dir(cli_args.output, parser)

    # Validate time series, motion, left surface, and right surface files
    # For all of them except time_series, "none" is a valid value
    if cli_args.series_file == "none" or cli_args.dtseries == "none":
        parser.error("series_file and --dtseries must be valid file paths, "
                     "not 'none'.")
    for file_arg in ["series_file", "motion", "left", "right",
                     "additional_mask"]:
        cli_args = validate_readable_file(file_arg, cli_args, parser)

    # Add whether time series is dense or parcellated to CLI args namespace
    setattr(cli_args, "time_series", get_valid_series_type(cli_args, parser))

    # Validate extra .dtseries .conc file path argument
    if not cli_args.dtseries:
        if cli_args.time_series == "dtseries":
            cli_args.dtseries = cli_args.series_file
        else:
            dtseries_file = cli_args.series_file.replace("dts", "pts")
            if "pts" in cli_args.series_file:
                cli_args.dtseries = dtseries_file
            else:
                parser.error("Cannot find dtseries file. Please enter the "
                             "path to an existing .conc file as the "
                             "--dtseries argument.")
    cli_args = validate_readable_file("dtseries", cli_args, parser)

    # Return cli_args with workbench command, MRE dir, and template file added
    return add_wb_command_mre_dir_and_template_to(cli_args, parser)


def get_valid_series_type(cli_args, parser):
    """
    Infer whether series_file is dense or parcellated by reading series_file
    provided by user. Raise a parser error if the series_file is neither or is
    otherwise invalid.
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: "dtseries" or "ptseries"
    """
    try:
        if cli_args.series_file[-11:] == "tseries.nii":
            time_series = cli_args.series_file[-12:-4]
        else:
            with open(cli_args.series_file, "r") as series_file:
                time_series = series_file.readline().split(".")[-2]

        # If size of .dconn files to make exceeds threshold, then warn user
        if time_series == "dtseries":
            if ((not cli_args.suppress_warnings) and
                    (CHOICES_TO_RUN[0] in cli_args.scripts_to_run or
                     CHOICES_TO_RUN[1] in cli_args.scripts_to_run)):
                warn_user_about_dconn_size(cli_args)

        # If invalid time series is given, then tell user and crash
        elif time_series != "ptseries":
            parser.error("Time series file must contain only a list of "
                         ".ptseries.nii or .dtseries.nii file paths.")
    except OSError:
        parser.error("Could not read time series file paths from "
                     + cli_args.series_file)
    except IndexError:
        parser.error("Series file must be a path to a .dtseries.nii file or "
                     "to a .ptseries.nii file, or a path to a .conc file "
                     "which only has a list of paths to those .nii files.")

    return time_series


def will_make_conn_conc(cli_args):
    """
    Check whether this wrapper will actually make a conn_matries .conc file
    :param cli_args: argparse namespace with all command-line arguments
    :return: True if this wrapper will create a .conc file with a list of
    connectivity matrices made by the wrapper, and False otherwise.
    """
    return ((CHOICES_TO_RUN[0] in cli_args.scripts_to_run or CHOICES_TO_RUN[1]
             in cli_args.scripts_to_run) and cli_args.make_conn_conc == "1")


def validate_and_get_input_dir(cli_args, parser):
    """
    Validate the directory containing all connectivity scripts and files
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if invalid dir path
    :return: Path to valid input directory as a string
    """
    input_dir = cli_args.input

    # If an input directory path is invalid, infer one from time series file
    if input_dir is None or not os.access(input_dir, os.R_OK):
        try:
            input_path = os.path.dirname(cli_args.series_file)
            assert os.path.isdir(input_path)
            input_dir = os.path.abspath(input_path)

        # If the input directory cannot be inferred from the series file, then
        # raise a parser error
        except (OSError, AssertionError):
            err_msg = ("--input argument {} is not a readable directory."
                       .format(input_dir))
            if not os.access(cli_args.series_file, os.R_OK):
                err_msg = ("{} This is because the time series file {} cannot "
                           "be read.".format(err_msg, cli_args.series_file))
            parser.error(err_msg)

    # Ensure that input folder path is formatted correctly with trailing slash
    if input_dir[-1] != "/":
        input_dir += "/"

    return input_dir


def validate_and_get_output_dir(output, parser):
    """
    Validate that output directory exists and is writeable, or create one, and
    then return valid path to that output directory
    :param output: File path of output directory to validate
    :param parser: argparse ArgumentParser to raise error if path is invalid
    :return: Valid path to output directory
    """
    try:
        output = os.path.abspath(output)
        if not os.path.isdir(output):
            os.makedirs(output)
        assert os.access(output, os.W_OK)
    except (OSError, TypeError):
        parser.error("Cannot make output folder at " + output)
    except AssertionError:
        parser.error("Cannot write to output folder at " + output)

    # Ensure that output directory is formatted properly with trailing slash
    if output[-1] != "/":
        output += "/"

    return output


def validate_readable_file(file_arg_name, cli_args, parser):
    """
    Ensure that an attribute of the argparse namespace either is "none" or
    points to a readable file. If neither, then try to find the file in the
    --input directory, and update the namespace accordingly if the file is
    found. Otherwise, raise an argparse error.
    :param file_arg_name: Filename attribute of cli_args to validate
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if path is invalid
    :return: cli_args, but with validated file_arg_name attribute
    """
    # Get the value of the file_arg_name attribute of cli_args object
    file_arg = getattr(cli_args, file_arg_name)

    # Unless file_arg is "none" or a readable file, try to find file_arg
    # in the cli_args.input directory, and update cli_args object accordingly
    if file_arg != "none" and not os.access(file_arg, os.R_OK):
        new_file_arg = os.path.join(cli_args.input, os.path.basename(file_arg))
        if os.access(new_file_arg, os.R_OK):
            setattr(cli_args, file_arg_name, new_file_arg)

        # If file_arg is still not found, raise a parser error
        else:
            msg = "Cannot find readable {} file at {}".format(file_arg_name,
                                                               file_arg)
            if file_arg != new_file_arg:
                msg = "{} or at {}".format(msg, new_file_arg)
            parser.error(msg)

    return cli_args


def warn_user_about_dconn_size(cli_args):
    """
    If the size of all of the .dconn files exceeds a certain threshold, then
    tell the user and make the user approve of .dconn creation to continue
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    # Open series file again to count its lines
    with open(cli_args.series_file, "r") as series_file:

        # Figure out how many .dconn files will be made (1 per line in series file,
        # but only 1 if .dconn files are deleted after being added to rolling avg)
        number_of_dconns_to_make = 1
        if (cli_args.keep_conn_matrices == "1" or not
                will_delete_conn_matrices_later(cli_args)):
            while series_file.readline():
                number_of_dconns_to_make += 1

    # Calculate how many gigabytes each .dconn file will be
    if cli_args.beta8 == "1":
        gb_per_dconn = 8
    else:
        gb_per_dconn = 33
    gb_of_files_created = number_of_dconns_to_make * gb_per_dconn

    # Warn user and ask for confirmation if file size exceeds threshold
    if gb_of_files_created > WARNING_IF_DCONN_SIZE_EXCEEDS:
        check = None
        while check.lower() != "y":
            check = input(
                "Warning: This will create about {} GB of .dconn files. Enter "
                "'y' to continue or 'n' to quit: ".format(gb_of_files_created)
            )
            if check.lower() == "n":
                sys.exit(1)


def get_matrices_concname(cli_args):
    """
    Create list of connectivity matrices from command-line argument inputs,
    validate it, add it to argparse namespace, and return the namespace
    :param cli_args: argparse namespace with all command-line arguments
    :return: Path to .conc file with list of connectivity matrices
    """
    # Get series_file name without its path
    series = os.path.basename(cli_args.series_file)

    # Re-format p/d series to create filename of matrices .conc list
    if cli_args.time_series == "ptseries":
        p_or_d = "_pconn_of_ptseries"
    else:
        p_or_d = "_dconn_of_dtseries"

    print(cli_args.fd_threshold)

    # Get parts of conn matrices list filename about minutes and FD threshold
    if cli_args.minutes != "none":
        minutes_and_fd = ("{}_minutes_of_data_at_FD_{}".format(
                          cli_args.minutes, cli_args.fd_threshold))
    else:
        minutes_and_fd = "_all_frames_at_FD_{}".format(cli_args.fd_threshold)

    return "".join((cli_args.output, series, p_or_d, minutes_and_fd, ".conc"))


def add_wb_command_mre_dir_and_template_to(cli_args, parser):
    """
    Validate and get workbench command, MRE directory, and template file. The
    first two of those depend on the server that this script is called from.
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if path is invalid
    :return: cli_args, but with validated wb_command, mre_dir, and template
    attributes
    """
    # If the user gave a wb_command, validate it; otherwise try to get a
    # default one, and ask user for one if on an unknown server
    if cli_args.wb_command:
        cli_args = validate_readable_file("wb_command", cli_args, parser)
    else:
        cli_args.wb_command = get_wb_command()

    # If no MRE dir was provided, use a default depending on the host server
    if not cli_args.mre_dir:
        host = gethostname()
        if host == "rushmore":
            cli_args.mre_dir = MRE_RUSHMORE
        elif "exa" in host:
            cli_args.mre_dir = MRE_EXACLOUD
        else:
            parser.error("Please enter a MATLAB Runtime Environment directory "
                         "using the --mre_dir flag.")
        if not os.access(cli_args.mre_dir, os.R_OK):
            parser.error("Cannot read MATLAB Runtime Environment directory at "
                         + cli_args.mre_dir)

    # If no template file name was given, then make one
    if not cli_args.template:
        cli_args.template = get_template_file_path(cli_args)

    # If pairwise_corr is being run before template, validate that template
    # file exists
    for script in cli_args.scripts_to_run:
        if script == CHOICES_TO_RUN[1]:
            break
        elif script == CHOICES_TO_RUN[2]:
            cli_args = validate_readable_file("template", cli_args, parser)
            cli_args.template = os.path.abspath(cli_args.template)
            break

    return cli_args


def get_wb_command():
    """
    Returns a valid path to a workbench command file
    :return: valid workbench command filepath
    """
    # If wb_command is already in BASH PATH, then get it and format it properly
    try:
        wb_command = str(subprocess.check_output((
            "which", "wb_command"
        )).decode("utf-8").strip())

    # Otherwise, use hardcoded wb_command depending on host server, or get
    # wb_command from user if host server is not known
    except subprocess.CalledProcessError:

        # If host server is RUSHMORE or EXACLOUD, then get its wb_command
        host = gethostname()
        if host == "rushmore":
            wb_command = WB_RUSHMORE
        elif "exa" in host:
            wb_command = WB_EXACLOUD

        # If on an unknown host, either quit if user says to, or get wb_command
        # from user and validate it
        else:
            print("No workbench command found on {}.".format(host))
            command_found = False
            while not command_found:
                wb_command = input(
                    "Please enter a valid path to a workbench command, or "
                    "enter 'q' to quit: "
                )
                if wb_command == "q":
                    sys.exit(1)
                elif os.access(wb_command, os.R_OK):
                    command_found = True
                else:
                    print("Error: Workbench command not found at {}.".format(
                              wb_command))

    return wb_command


def get_template_file_path(cli_args):
    """
    Build and return name of file created by cifti_conn_template
    :param cli_args: argparse namespace with all command-line arguments
    :return: String name of file created by cifti_conn_template
    """
    # Get name of series file without its directory path
    series_file_name = os.path.basename(cli_args.series_file)

    # Build filename depending on minutes_limit and motion_file
    if cli_args.minutes == "none":
        fd = '_all_frames_at_FD_'
    else:
        fd = '_FD_'
    if cli_args.motion == "none":
        motion_and_smooth = cli_args.motion + '_and_smoothing_'
    elif cli_args.minutes == "none":
        motion_and_smooth = cli_args.fd_threshold + '_and_smoothing_'
    else:
        motion_and_smooth = ("{}_at_{}_minutes_remaining_smoothing_".format(
                             cli_args.fd_threshold, cli_args.minutes))

    # Return the full file path and name of template file, assuming that
    # template file is in cli_args.output directory
    return "".join((cli_args.output, series_file_name, fd, motion_and_smooth,
                    cli_args.smoothing_kernel, '_AVG.',
                    cli_args.time_series[0], "conn.nii"))


def get_conc_file_paths(cli_args):
    """
    Build and return an incomplete path to one .conc file from input parameters
    """
    paths = os.path.join(cli_args.output, "{}_{}conn_of_{}*{}.conc".format(
        os.path.basename(cli_args.series_file), cli_args.time_series[0], 
        cli_args.time_series, cli_args.fd_threshold
    ))
    if not next(iglob(paths), None):
        paths = os.path.join(cli_args.output, "{}_{}conn_of_{}*{}.conc".format(
            os.path.basename(cli_args.series_file), cli_args.time_series[0], 
            cli_args.time_series, "none"
        ))
        if not next(iglob(paths), None):
            print("No .conc file(s) at {}".format(paths))
    return paths


def get_matrix_or_template_parameters(cli_args):
    """
    cifti_conn_matrix and cifti_conn_template both have the same required
    parameters, with only a few exceptions. This function returns a list of
    all parameters required by both scripts.
    :param cli_args: Full argparse namespace containing all CLI args, including
    all necessary parameters for cifti_conn_matrix and cifti_conn_template.
    :return: A list of all parameters required by matrix and template scripts.
    """
    return([
        cli_args.wb_command,
        cli_args.series_file,
        cli_args.time_series,
        cli_args.motion,
        cli_args.fd_threshold,
        cli_args.tr,
        cli_args.minutes,
        cli_args.smoothing_kernel,
        cli_args.left,
        cli_args.right,
        cli_args.beta8,
        cli_args.remove_outliers,
        cli_args.additional_mask,
        cli_args.make_conn_conc,
        cli_args.output,
        cli_args.dtseries
    ])


def will_delete_conn_matrices_later(cli_args):
    """
    Checks whether this wrapper should keep connectivity matrices until done
    running cifti_conn_pairwise corr and then delete them afterwards.
    :param cli_args: argparse namespace with all command-line arguments
    :return: True if user said not to keep connectivity matrices even though
    they are needed to run cifti_conn_pairwise_corr, and False otherwise.
    """
    return (cli_args.keep_conn_matrices == "0"
            and CHOICES_TO_RUN[2] in cli_args.scripts_to_run)


def cifti_conn_matrix(cli_args):
    """
    Run cifti_conn_matrix script to generate a correlation matrix of all
    greyordinates/parcellations to all greyordinates/parcellations
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    subprocess.check_call([os.path.join(DEFAULT_SOURCE, SCRIPT_MATRIX),
                           cli_args.mre_dir]
                          + get_matrix_or_template_parameters(cli_args))


def cifti_conn_template(cli_args):
    """
    Run cifti_conn_template script to build a template connectivity matrix of
    a list of subjects.
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    # If user is running pairwise_corr, keep conn matrices until that finishes
    if will_delete_conn_matrices_later(cli_args):
        print("Warning: CIFTI conn matrix files will not be deleted until {} "
              " finishes running.".format(CHOICES_TO_RUN[2]))
        cli_args.keep_conn_matrices = "1"

    # Get required parameters for cifti_conn_template script
    parameters = get_matrix_or_template_parameters(cli_args)
    parameters.append(cli_args.keep_conn_matrices)
    parameters.append(cli_args.template)

    # Call cifti_conn_template script
    subprocess.check_call([os.path.join(DEFAULT_SOURCE, SCRIPT_TEMPLATE),
                           cli_args.mre_dir] + parameters)


def cifti_conn_pairwise_corr(cli_args):
    """
    Run cifti_conn_pairwise_corr script to generate a correlation of
    correlation matrices.
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    # Reformat time series variable to pass to cifti_conn_pairwise_corr script
    p_or_d = cli_args.time_series[0] + "conn"

    # Call cifti_conn_pairwise_corr script
    for conc_file in iglob(get_conc_file_paths(cli_args)):
        print(conc_file)
        subprocess.check_call((
            os.path.join(DEFAULT_SOURCE, SCRIPT_PAIRWISE_CORR),
            cli_args.mre_dir,
            cli_args.wb_command,
            cli_args.template,
            p_or_d,
            conc_file,
            cli_args.keep_conn_matrices
        ))

        # If user said not to keep conn matrices but pairwise_corr used them,
        # then try to delete the conn matrices after pairwise_corr finishes
        if will_delete_conn_matrices_later(cli_args):
            try:
                with open(conc_file) as infile:
                    for line in infile:
                        print("Deleting " + line)
                        os.unlink(line)
            except FileNotFoundError as e:
                print("Could not find file at " + e.filename)


if __name__ == '__main__':
    main()
