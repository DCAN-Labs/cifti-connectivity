#! /usr/bin/env python3

"""
CIFTI conn wrapper
Greg Conan: conan@ohsu.edu
Created 2019-06-18
Last Updated 2019-07-05
"""

##################################
#
# Wrapper for CIFTI conn scripts that can be run from the command line.
# This wrapper can run any one or multiple of the three following scripts:
# 1. Run `cifti_conn_matrix` on your dtseries or ptseries to generate a
#    correlation matrix or all the greyordinates/parcellations to all other
#    greyordinate/parcellations.
# 2. Run `cifti_conn_template` to build a template connectivity matrix of a
#    list of subjects.
# 3. Run `cifti_conn_pairwise_corr` to generate a correlation of correlation
#    matrices.
#
##################################

import argparse
import os
from pathlib import Path
from socket import gethostname
import subprocess
import sys

# Constants

# Names of scripts to run from this one
CHOICES_TO_RUN = ["cifti_conn_matrix", "cifti_conn_template",
                  "cifti_conn_pairwise_corr"]

# Default directories for input data, output data, and wrapped scripts (src)
DEFAULT_INPUT = "./raw/"
DEFAULT_OUTPUT = "./data/"
DEFAULT_SOURCE = "./src/"

# Default names of data files used by wrapped scripts
GROUP_LEFT_CONC = (DEFAULT_INPUT
                   + "ADHD_DVARS_group_left_midthickness_surfaces.conc")
GROUP_MOTION_CONC = DEFAULT_INPUT + "ADHD_DVARS_group_motion.conc"
GROUP_RIGHT_CONC = (DEFAULT_INPUT
                    + "ADHD_DVARS_group_right_midthickness_surfaces.conc")

# MATLAB Runtime Environment (MRE) directories which depend on the host server
MRE_EXACLOUD = ("/home/exacloud/lustre1/fnl_lab/code/external/utilities/"
                "Matlab2016bRuntime/v91")
MRE_RUSHMORE = ("/mnt/max/shared/code/external/utilities/Matlab2016bRuntime/"
                "v91")

# Default names of wrapped scripts to run
SCRIPT_MATRIX = "run_cifti_conn_matrix_for_wrapper.sh"
SCRIPT_PAIRWISE_CORR = "run_cifti_conn_pairwise_corr_exaversion.sh"
SCRIPT_TEMPLATE = "run_cifti_conn_template_for_wrapper.sh"

# If .dconn files' size exceeds this many gigabytes, then warn the user
WARNING_IF_DCONN_SIZE_EXCEEDS = 300

# Workbench commands which depend on the host server
WB_EXACLOUD = ("/home/exacloud/lustre1/fnl_lab/code/external/utilities/"
               "workbench-1.3.2/bin_rh_linux64/wb_command")
WB_RUSHMORE = "/mnt/max/software/workbench/bin_linux64/wb_command"


def main():

    cli_args = get_cli_args()

    # Run all scripts the user said to run in the order the user gave them
    for script in cli_args.scripts_to_run:

        # Tell user when script started
        print("\nCIFTI conn wrapper started running " + script + " at:")
        subprocess.call("date")

        # Call script
        globals()[script](cli_args)

        # Tell user when script finished
        print("\nCIFTI conn wrapper finished running " + script + " at:")
        subprocess.call("date")


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

    # Optional: Get path to directory containing HCP CIFTI conn scripts
    parser.add_argument(
        "-i",
        "--input",
        default=DEFAULT_INPUT,
        help=("Directory containing all HCP conn data files. By default, this "
              "will be " + DEFAULT_INPUT + " But if that does not exist or "
              "is empty, the default will instead be the directory containing "
              "the time series file.")
    )

    # Optional: Get name of .conc file listing FNL motion mat files
    parser.add_argument(
        "-m",
        "--motion",
        default=GROUP_MOTION_CONC,
        help=("Name of .conc file that points to FNL motion mat files for each"
              "dt or ptseries. Default name is " + GROUP_MOTION_CONC)
    )

    # Optional: Get name of template file to create
    parser.add_argument(
        "-t",
        "--template",
        help="File path and name of template file to be created."
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

    # Optional: Get name of left surface .conc file
    parser.add_argument(
        "-l",
        "--left",
        nargs="?",
        const=GROUP_LEFT_CONC,
        default="none",
        help=(".conc file name of subjects' left midthicknessfile. Include "
              "this flag if smoothing will be run. If this flag is included"
              "but no filename is given, then " + GROUP_LEFT_CONC
              + " will be used as the default name.")
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
              "but no filename is given, then " + GROUP_RIGHT_CONC
              + " will be used as the default name.")
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
              "8Gb .dconns. Otherwise, this will make 32Gb .dconns. This "
              "option does nothing for ptseries.")
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

    # Optional: Specify output folder
    parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_OUTPUT,
        help=("Location to save all output files to. By default, this will be "
              "the directory specified in --input.")
    )

    # Optional: Give path to MATLAB Runtime Environment (MRE) directory
    parser.add_argument(
        "--mre_dir",
        help=("Path to directory containing MATLAB Runtime Environment (MRE)"
              "version 9.1 or newer. This is used to run compiled MATLAB "
              "scripts. This argument must be a valid path to an existing "
              "folder.")
    )

    # Optional: Specify path to wb_command
    parser.add_argument(
        "-w",
        "--wb_command",
        help="Path to workbench command file."
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

    # Optional: Specify whether to remove outliers from BOLD signal
    parser.add_argument(
        "--remove_outliers",
        action="store_const",
        const="1",
        default="0",
        help=("If this flag is included, then outliers will be removed from "
              "the BOLD signal. Otherwise, frames will only be censored by "
              "the FD threshold.")
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

    # Optional: Ignore warnings about output .dconn file sizes
    parser.add_argument(
        "--suppress_warnings",
        action="store_true",
        help=("By default, the wrapper will ask user for confirmation if "
              "the .dconn files created by the wrapper will exceed "
              + str(WARNING_IF_DCONN_SIZE_EXCEEDS) + " gigabytes. Include "
              "this argument to ignore this warning. Irrelevant for ptseries.")
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
    if cli_args.series_file == "none":
        parser.error("series_file must be a valid file path and not 'none'.")
    else:
        for file_arg in ["series_file", "motion", "left", "right"]:
            cli_args = validate_readable_file(file_arg, cli_args, parser)

    # Infer whether series_file is dense or parcellated by reading series_file
    # provided by user
    try:
        with open(cli_args.series_file, "r") as series_file:
            time_series = series_file.readline().split(".")[-2]

            # If size of .dconn files to make exceeds threshold, then warn user
            if time_series == "dtseries":
                if (not cli_args.suppress_warnings and
                        (CHOICES_TO_RUN[0] in cli_args.scripts_to_run or
                         CHOICES_TO_RUN[1] in cli_args.scripts_to_run)):
                    warn_user_about_dconn_size(series_file, cli_args)

            # If invalid time series is given, then tell user and crash
            elif time_series != "ptseries":
                parser.error("Time series file must contain only a list of "
                             ".ptseries.nii or .dtseries.nii file paths.")
    except OSError:
        print("Error: Could not read time series file paths from "
              + cli_args.series_file)
        sys.exit(1)

    # Add whether time series is dense or parcellated to CLI args namespace
    cli_args.__setattr__("time_series", time_series)

    # Add list of connectivity matrices to CLI args namespace if pairwise_corr
    # is being run
    if CHOICES_TO_RUN[2] in cli_args.scripts_to_run:
        cli_args.__setattr__("conn_matrices", get_conn_matrices_list(cli_args))
        cli_args = validate_readable_file("conn_matrices", cli_args, parser)

    # Return cli_args with workbench command, MRE dir, and template file added
    return add_wb_command_mre_dir_and_template_to(cli_args, parser)


def validate_and_get_input_dir(cli_args, parser):
    """
    Validate the directory containing all connectivity scripts and files
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if invalid dir path
    :return: Path to valid input directory as a string
    """
    input_dir = cli_args.input

    # If a input directory path is invalid, infer one from time series file
    if not os.access(input_dir, os.R_OK):
        try:
            input_path = Path(cli_args.series_file).resolve().input
            assert input_path.is_dir()
            input_dir = str(input_path)
        except (OSError, AssertionError):
            parser.error(input_dir + " is not an existing directory.")

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
        Path(output).mkdir(exist_ok=True, parents=True)
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
        new_file_arg = cli_args.input + Path(file_arg).name
        if os.access(new_file_arg, os.R_OK):
            setattr(cli_args, file_arg_name, new_file_arg)

        # If file_arg is still not found, raise a parser error
        else:
            msg = ("Cannot find readable " + file_arg_name + " file at "
                   + file_arg)
            if file_arg != new_file_arg:
                msg += (" or at " + new_file_arg)
            parser.error(msg)

    return cli_args


def warn_user_about_dconn_size(series_file, cli_args):
    """
    If the size of all of the .dconn files exceeds a certain threshold, then
    tell the user and make the user approve of .dconn creation to continue
    :param series_file: Currently-open file containing the names of all .dconn
    files to create.
    :param cli_args: argparse namespace with all command-line arguments
    :param beta8: Argument specifying whether to use file size compression
    :return: N/A
    """
    # Figure out how many .dconn files will be made (1 per line in series file,
    # but only 1 if .dconn files are deleted after being added to rolling avg)
    number_of_dconns_to_make = 1
    if (cli_args.keep_conn_matrices == "1" or
            will_delete_conn_matrices_later(cli_args)):
        while series_file.readline():
            number_of_dconns_to_make += 1

    print(cli_args.beta8)

    # Calculate how many gigabytes each .dconn file will be
    if cli_args.beta8 == "1":
        gb_per_dconn = 8
    else:
        gb_per_dconn = 33
    gb_of_files_created = number_of_dconns_to_make * gb_per_dconn

    # Warn user and ask for confirmation if file size exceeds threshold
    if gb_of_files_created > WARNING_IF_DCONN_SIZE_EXCEEDS:
        check = None
        while check != "y":
            check = input("Warning: This will create about "
                          + str(gb_of_files_created) + " GB of .dconn files. "
                          "Press y to continue or n to quit: ")
            if check == "n":
                sys.exit(1)


def get_conn_matrices_list(cli_args):
    """
    Create list of connectivity matrices from command-line argument inputs,
    validate it, add it to argparse namespace, and return the namespace
    :param cli_args: argparse namespace with all command-line arguments
    :return: Path to .conc file with list of connectivity matrices
    """
    # Get series_file name without its path
    series = Path(cli_args.series_file).name

    # Re-format p/d series to create filename of conn_matrices list
    if cli_args.time_series == "ptseries":
        p_or_d = "_pconn_of_ptseries"
    else:
        p_or_d = "_dconn_of_dtseries"

    # Get parts of conn_matrices_list filename about minutes and FD threshold
    if cli_args.minutes != "none":
        minutes_and_fd = (cli_args.minutes + "_minutes_of_data_at_FD_"
                          + cli_args.fd_threshold + ".conc")
    else:
        minutes_and_fd = "_all_frames_at_FD_" + cli_args.fd_threshold

    return cli_args.output + series + p_or_d + minutes_and_fd + ".conc"


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
        cli_args = validate_readable_file(cli_args, "wb_command", parser)
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
            cli_args.template = str(Path(cli_args.template).resolve())
            break

    return cli_args


def get_wb_command():
    """
    Returns a valid path to a workbench command file
    :return: valid workbench command filepath
    """
    # If wb_command is already in BASH PATH, then get it and format it properly
    try:
        wb_command = subprocess.check_output([
            "which", "wb_command"
        ]).decode("utf-8").strip()

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
            print("No workbench command found on " + host + ".")
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
                    print("Error: Workbench command not found at "
                          + wb_command + ".", end=" ")

    return wb_command


def get_template_file_path(cli_args):
    """
    Build and return name of file created by cifti_conn_template
    :param cli_args: argparse namespace with all command-line arguments
    :return: String name of file created by cifti_conn_template
    """
    # Get name of series file without its directory path
    series_file_name = str(Path(cli_args.series_file).name)

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
        motion_and_smooth = (cli_args.fd_threshold + '_at_' + cli_args.minutes
                             + '_minutes_remaining_smoothing_')

    # Return the full file path and name of template file, assuming that
    # template file is in cli_args.output directory
    return (cli_args.output + series_file_name + fd + motion_and_smooth
            + cli_args.smoothing_kernel + '_AVG.' + cli_args.time_series[0]
            + "conn.nii")


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
        cli_args.output
    ])


def will_delete_conn_matrices_later(cli_args):
    """
    Checks whether this wrapper should keep connectivity matrices until done
    running cifti_conn_pairwise corr and then delete them afterwards.
    :param cli_args: argparse namespace with all command-line arguments
    :return: True if user said not to keep connectivity matrices even though
    they are needed to run cifti_conn_pairwise_corr, and False otherwise.
    """
    return (
        cli_args.keep_conn_matrices == "0"
        and CHOICES_TO_RUN[2] in cli_args.scripts_to_run
    )


def cifti_conn_matrix(cli_args):
    """
    Run cifti_conn_matrix script to generate a correlation matrix of all
    greyordinates/parcellations to all greyordinates/parcellations
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    # Get all parameters to pass to matrix script, including its options
    parameters = get_matrix_or_template_parameters(cli_args)

    # Call matrix script
    subprocess.check_call([DEFAULT_SOURCE + SCRIPT_MATRIX, cli_args.mre_dir]
                          + parameters)


def cifti_conn_template(cli_args):
    """
    Run cifti_conn_template script to build a template connectivity matrix of
    a list of subjects.
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    # If user is running pairwise_corr, keep conn matrices until that finishes
    if will_delete_conn_matrices_later(cli_args):
        print("Warning: CIFTI conn matrix files will not be deleted until "
              + CHOICES_TO_RUN[2] + " finishes running.")
        cli_args.keep_conn_matrices = "1"

    # Get required parameters for cifti_conn_template script
    parameters = get_matrix_or_template_parameters(cli_args)
    parameters.append(cli_args.keep_conn_matrices)
    parameters.append(cli_args.template)

    # Call cifti_conn_template script
    subprocess.check_call([DEFAULT_SOURCE + SCRIPT_TEMPLATE, cli_args.mre_dir]
                          + parameters)


def cifti_conn_pairwise_corr(cli_args):
    """
    Run cifti_conn_pairwise_corr script to generate a correlation of
    correlation matrices.
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    # Reformat time series variable to pass to cifti_conn_pairwise_corr script
    p_or_d = cli_args.time_series[0] + "conn"

    # Get parameters to call cifti_conn_pairwise_corr script
    parameters = [
        cli_args.wb_command,
        cli_args.template,
        p_or_d,
        cli_args.conn_matrices,
        cli_args.keep_conn_matrices,
    ]

    # Call cifti_conn_pairwise_corr script
    subprocess.check_call([DEFAULT_SOURCE + SCRIPT_PAIRWISE_CORR,
                           cli_args.mre_dir] + parameters)

    # If user said not to keep conn matrices but pairwise_corr used them,
    # then try to delete the conn matrices after pairwise_corr finishes
    if will_delete_conn_matrices_later(cli_args):
        try:
            with open(cli_args.conn_matrices) as infile:
                for line in infile:
                    print("Deleting " + line)
                    Path(line.strip()).unlink()
        except FileNotFoundError as e:
            print("Could not find file at " + e.filename)


if __name__ == '__main__':
    main()