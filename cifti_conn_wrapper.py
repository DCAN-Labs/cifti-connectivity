#! /usr/bin/env python3

"""
CIFTI connectivity wrapper
Greg Conan: conan@ohsu.edu
Created 2019-06-18
Last Updated 2020-02-07
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
import glob
import os
import socket
import subprocess
import sys

### Constants

# Names of scripts to run from this one
CHOICES_TO_RUN = ["matrix", "template", "pairwise_corr"]

# Get path to directory containing this file
try:
    PWD = os.path.dirname(os.path.abspath(sys.argv[0]))
except OSError as e:
    sys.exit("{} {} cannot find the path to the directory it is in."
                .format(e, sys.argv[0]))

# Default directories for input data and wrapped scripts (src)
DEFAULT_INPUT = os.path.join(PWD, "raw")
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

    def timestamp_when(step, completion):
        """
        Get the time and date when a step of this script started/finished
        :param step: String naming the step which either started or finished
        :param completion: String which is either 'started' or 'finished'
        :return: message with the date and time when step started or finished
        """
        return("CIFTI connectivity matrix wrapper {} running {} script at {}"
               .format(completion, step, datetime.now().strftime(
                  "%H:%M:%S on %b %d, %Y")))

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
    # Create arg parser
    parser = argparse.ArgumentParser(
        description="Wrapper to generate connectivity matrices."
    )

    # Required arguments: Time series filename, repetition time of scan,
    # script(s) to run, and output directory
    parser.add_argument(
        "series_file",
        type=readable_file,
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
        type=readable_file_or_none,
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
        type=readable_file,
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

    # Optional: Get path to directory containing HCP CIFTI conn scripts
    parser.add_argument(
        "-in",
        "--input",
        type=readable_dir,
        help=("Directory containing all HCP conn data files. By default, this "
              "will be the directory containing the time series file.")
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
        type=readable_file_or_none,
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
        type=readable_file_or_none,
        nargs="?",
        const=GROUP_MOTION_CONC,
        default="none",
        help=("Name of .conc file that points to FNL motion mat files for each"
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
        type=readable_dir,
        help=("Path to directory containing MATLAB Runtime Environment (MRE)"
              "version 9.1. This is used to run compiled MATLAB "
              "scripts. This argument must be a valid path to an existing "
              "folder.")
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
        type=readable_file_or_none,
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
        type=readable_file,
        help="File path and name of template file to be created."
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
        type=readable_file,
        help="Path to workbench command file called 'wb_command'."
    )

    return validate_cli_args(parser.parse_args(), parser)


def readable_file(path):
    """
    Throw argparse exception unless parameter is a valid readable filename 
    string. This is used instead of argparse.FileType("r") because the latter 
    leaves an open file handle, which has caused problems.
    :param path: Parameter to check if it represents a valid filename
    :return: String representing a valid filename
    """
    try:
        assert os.access(path, os.R_OK)
        return os.path.abspath(path)
    except (AssertionError, OSError, TypeError):
        raise argparse.ArgumentTypeError("Cannot read file at {}".format(path))


def readable_file_or_none(path):
    """
    Throw argparse exception unless path either is "none" or points to a 
    readable file.
    :param path: Parameter to check if it represents a valid filename
    :return: String representing a valid filename, or the string "none"
    """
    return path if path == "none" else readable_file(path)


def readable_dir(path):
    """
    Throw argparse exception unless path is a string representing a valid path 
    to a readable directory.
    :param path: Parameter to check if it represents a valid filename
    :return: String representing a valid path to an existing readable directory
    """
    if os.path.isdir(path):
        return readable_file(path)
    else:
        raise argparse.ArgumentTypeError(("{} is not a readable "
                                          "directory").format(path))


def none_or_valid_float_value_as_string(str_to_check):
    """
    Unless a string is "none", tries to convert it to a float and back to check
    that it represents a valid float value. Throws ValueError if type
    conversion fails. This function is only needed because the MATLAB scripts
    take some arguments either as a float value in string form or as "none".
    :param str_to_check: string to validate
    :return: string which is either "none" or represents a valid float value
    """
    return str_to_check if str_to_check == "none" else str(float(str_to_check))


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

    # Add whether time series is dense or parcellated to CLI args namespace
    setattr(cli_args, "time_series", get_valid_series_type(cli_args, parser))

    # If the user did not give a wb_command, then try to get a default one, or
    # ask user for one if on an unknown server
    if not cli_args.wb_command:
        cli_args.wb_command = get_valid_wb_command()

    # Validate other arguments
    cli_args.dtseries = validate_and_get_dtseries(cli_args, parser)
    cli_args.mre_dir = validate_and_get_mre_dir(cli_args.mre_dir, parser)
    cli_args.template = validate_and_get_template(cli_args)

    # Validate that all .conc file inputs have an equal number of lines
    validate_concs_same_length([k for k, v in vars(cli_args).items()
                                if isinstance(v, str) and v[-5:] == ".conc"],
                               cli_args, parser)
    return cli_args


def validate_and_get_input_dir(cli_args, parser):
    """
    Validate the directory containing all connectivity scripts and files
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if invalid dir path
    :return: Path to valid input directory as a string
    """
    # If an input directory path is invalid, infer one from time series file
    if not cli_args.input:
        try:
            return format_with_trailing_slash(readable_dir(os.path.dirname(
                cli_args.series_file)))

        # If the input directory cannot be inferred from the series file, then
        # raise a parser error
        except (OSError, AssertionError):
            err_msg = ("--input argument {} is not a readable directory."
                       .format(input_dir))
            if not os.access(cli_args.series_file, os.R_OK):
                err_msg = ("{} This is because the time series file {} cannot "
                           "be read.".format(err_msg, cli_args.series_file))
            parser.error(err_msg)


def format_with_trailing_slash(path):
    """
    :param path: String representing a file path
    :return: path, but with one "/" at the end
    """
    return path if path[-1] == "/" else path + "/"


def validate_and_get_output_dir(output, parser):
    """
    Validate that output directory exists and is writeable, or create one, and
    then return valid path to that output directory
    :param output: String representing file path of output directory to check
    :param parser: argparse.ArgumentParser to raise error if path is invalid
    :return: Valid path to output directory
    """
    try:  
        # Output must be a directory, not a file
        if os.path.isfile(output):  
            raise OSError()
        else:
            os.makedirs(os.path.abspath(output), exist_ok=True)

        # Output must be writeable and formatted with a trailing slash
        assert os.access(output, os.W_OK)
        return format_with_trailing_slash(output)

    except (OSError, TypeError):
        parser.error("Cannot make output folder at " + output)
    except AssertionError:
        parser.error("Cannot write to output folder at " + output)


def warn_user_about_dconn_size(cli_args):
    """
    If the size of all of the .dconn files exceeds a certain threshold, then
    tell the user and make the user approve of .dconn creation to continue
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    # Open series file again to count its lines
    with open(cli_args.series_file, "r") as series_file:

        # Figure out how many .dconn files to make (1 per line in series 
        # file, but only 1 if .dconn files are deleted after being added  
        # to the rolling average)
        number_of_dconns_to_make = 1
        if (cli_args.keep_conn_matrices == "1" or not
                will_delete_conn_matrices_later(cli_args)):
            while series_file.readline():
                number_of_dconns_to_make += 1

    # Calculate how many gigabytes each .dconn file will be
    gb_per_dconn = 8 if cli_args.beta8 == "1" else 33
    gb_of_files_created = number_of_dconns_to_make * gb_per_dconn

    # Warn user and ask for confirmation if file size exceeds threshold
    if gb_of_files_created > WARNING_IF_DCONN_SIZE_EXCEEDS:
        check = "y"
        while check.lower() != "y":
            check = input(
                "Warning: This will create about {} GB of .dconn files. Enter "
                "'y' to continue or 'n' to quit: ".format(gb_of_files_created)
            )
            if check.lower() == "n":
                sys.exit(1)


def will_delete_conn_matrices_later(cli_args):
    """
    Checks whether this wrapper should keep connectivity matrices until done
    running cifti_conn_pairwise corr and then delete them afterwards.
    :param cli_args: argparse namespace with all command-line arguments
    :return: True if user said not to keep connectivity matrices even though
    they are needed to run cifti_conn_pairwise_corr, and False otherwise.
    """
    return (cli_args.keep_conn_matrices == "0"
            and CHOICES_TO_RUN[2] in cli_args.scripts)


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
        elif cli_args.series_file[-8:] == "conn.nii":
            time_series = cli_args.series_file[-9] + "tseries"
        else:
            with open(cli_args.series_file, "r") as series_file:
                time_series = series_file.readline().split(".")[-2]

        # If size of .dconn files to make exceeds threshold, then warn user
        if time_series == "dtseries":
            if ((not cli_args.suppress_warnings) 
                and os.path.splitext(cli_args.series_file)[1] != ".nii" 
                and (CHOICES_TO_RUN[0] in cli_args.scripts
                     or CHOICES_TO_RUN[1] in cli_args.scripts)):
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
    

def get_valid_wb_command():
    """
    Try to find a valid workbench command file by checking the BASH path, by
    checking default locations on known servers, or finally by asking the user
    :return: String representing a valid path to the wb_command file
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
        host = socket.gethostname()
        if host == "rushmore":
            wb_command = WB_RUSHMORE
        elif "exa" in host:
            wb_command = WB_EXACLOUD

        # If on an unknown host, either quit if user says to, or get wb_command
        # from user and validate it
        else:
            print("No workbench command found on {}.".format(host))
            while not os.access(wb_command, os.R_OK):
                wb_command = input("Please enter a valid path to a workbench "
                                   "command, or 'q' to end the program: ")
                if wb_command == "q":
                    sys.exit(1)

    return readable_file(wb_command)


def validate_and_get_dtseries(cli_args, parser):
    """
    Validate the path to the .dtseries file used for outlier removal
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if invalid dir path
    :return: Path to valid .dtseries.nii file as a string, or the "none" string
    """
    if cli_args.remove_outliers != "1":
        dtseries = "none"
    elif cli_args.dtseries:
        dtseries = cli_args.dtseries
    elif cli_args.time_series == "dtseries":
        dtseries = cli_args.series_file
    else:
        dtseries = cli_args.series_file.replace("ptseries", "dtseries")
    return readable_file_or_none(dtseries)


def validate_and_get_mre_dir(path, parser):
    """
    Validate and get path to Matlab Runtime Environment directory. The default
    path depends on the server that this file is called from.
    :param path: None, or string representing a valid path to MRE directory 
    :param parser: argparse.ArgumentParser to raise error if anything's invalid
    :return: String representing valid path to MRE directory
    """
    # If no MRE dir was provided, use a default depending on the host server
    if not path:
        host = socket.gethostname()
        if host == "rushmore":
            path = MRE_RUSHMORE
        elif "exa" in host:
            path = MRE_EXACLOUD
        else:
            parser.error("Please enter a MATLAB Runtime Environment directory "
                         "using the --mre_dir flag.")
        if not os.access(path, os.R_OK):
            parser.error("Cannot read MATLAB Runtime Environment directory at "
                         + path)
    return readable_file(path)


def validate_and_get_template(cli_args):
    """
    Validate and get path to template file, either to be imported (if 
    pairwise_corr will be run before template) or made (if template will be run 
    before pairwise_corr).
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if path is invalid
    :return: cli_args, but with validated wb_command, mre_dir, and template
    attributes
    """
    # If no template file name was given, then make one
    if not cli_args.template:
        template = get_template_file_path(cli_args)

    # If pairwise_corr will be run before template, validate that template
    # file exists; otherwise it's unneeded / this script will make the template
    for script in cli_args.scripts:
        if script == CHOICES_TO_RUN[1]:
            break
        elif script == CHOICES_TO_RUN[2]:
            template = readable_file(cli_args.template)
            break
    return template


def get_template_file_path(cli_args):
    """
    Build and return the name of the file created by cifti_conn_template
    :param cli_args: argparse namespace with all command-line arguments
    :return: String name of file created by cifti_conn_template
    """
    # Get name of series file without its directory path
    series_file_name = os.path.basename(cli_args.series_file)

    # Build filename parts about minutes_limit and motion_file
    fd_str = '_all_frames_at_FD_' if cli_args.minutes == "none" else '_FD_'
    if cli_args.motion == "none":
        motion_and_smooth = cli_args.motion + '_and_smoothing_'
    elif cli_args.minutes == "none":
        motion_and_smooth = cli_args.fd + '_and_smoothing_'
    else:
        motion_and_smooth = ("{}_at_{}_minutes_remaining_smoothing_".format(
                             cli_args.fd, cli_args.minutes))

    # Return the full file path and name of template file, assuming that
    # template file is in cli_args.output directory
    return "".join((cli_args.output, series_file_name, fd_str, 
                    motion_and_smooth, cli_args.smoothing_kernel, '_AVG.',
                    cli_args.time_series[0], "conn.nii"))


def validate_concs_same_length(conc_file_args, cli_args, parser):
    """
    Throw argparse error unless all .conc files have the same number of lines
    :param conc_file_args: List of names of arguments in cli_args which are
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    valid paths to .conc files which each contain a list of file paths
    :return: N/A
    """
    # Get the line number of every .conc file
    prev_file_lines = -1
    for conc_arg in conc_file_args:
        conc_path = getattr(cli_args, conc_arg)
        with open(conc_path) as conc_file:
            lines = 0
            for line in conc_file:
                lines += 1

        # If each file has the same number of lines as the previous one, then
        # keep checking; otherwise raise a parser error
        if prev_file_lines == -1 or lines == prev_file_lines:
            prev_file_lines = lines
            prev_conc = conc_path
        else:
            parser.error("These .conc files have different numbers of lines:\n"
                         "{}\n{}\n\nPlease use .conc files which all have the "
                         "same number of lines.".format(prev_conc, conc_path))


def get_conc_file_paths(cli_args):
    """
    Build an incomplete path to .conc file(s) from input parameters, then use
    it to return a list of path(s) to either .conc file(s) or a matrix
    :param cli_args: argparse namespace with all command-line arguments
    :return: List of strings which are paths to .conc files or matrix files
    """

    def get_paths(cli_args, fd_or_none):
        """
        Builds incomplete path depending on whether .conc file was named using
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
    all necessary parameters for cifti_conn_matrix and cifti_conn_template.
    :return: A list of all parameters required by matrix and template scripts.
    """
    return([
        cli_args.mre_dir,
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
        cli_args.dtseries
    ])


def cifti_conn_matrix(cli_args):
    """
    Run cifti_conn_matrix script to generate a correlation matrix of all
    greyordinates/parcellations to all greyordinates/parcellations
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    subprocess.check_call((
        os.path.join(DEFAULT_SOURCE, SCRIPT_MATRIX),
        *get_matrix_or_template_parameters(cli_args)
    ))


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
        print("Warning: CIFTI conn matrix files will not be deleted until {} "
              "finishes running.".format(CHOICES_TO_RUN[2]))
        keep_conn_matrices = "1"

    # Call cifti_conn_template script
    subprocess.check_call((
        os.path.join(DEFAULT_SOURCE, SCRIPT_TEMPLATE),
        *get_matrix_or_template_parameters(cli_args),
        keep_conn_matrices,
        cli_args.template
    ))


def cifti_conn_pairwise_corr(cli_args):
    """
    Run cifti_conn_pairwise_corr script to generate a correlation of
    correlation matrices.
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    # Run pairwise_corr on each .conc file listing all of conn matrices made by
    # by a previous step of the script
    for conc_file in get_conc_file_paths(cli_args):

        # Call cifti_conn_pairwise_corr script
        subprocess.check_call((
            os.path.join(DEFAULT_SOURCE, SCRIPT_PAIRWISE_CORR),
            cli_args.mre_dir,
            cli_args.wb_command,
            cli_args.template,
            cli_args.time_series[0] + "conn",
            conc_file,
            cli_args.keep_conn_matrices,
            cli_args.output
        ))

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
