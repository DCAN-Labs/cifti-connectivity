#! /usr/bin/env python3

"""
CIFTI connectivity wrapper
Greg Conan: conan@ohsu.edu
Created 2019-06-18
Updated 2020-03-03
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

# Imports
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

# To get default directories for input data (raw) and wrapped scripts (src)
def local_path_to(directory, path):
    return os.path.join(PWD, directory, path)

# Default names of data files used by wrapped scripts
midthicknessfile = "ADHD_DVARS_group_{}_midthickness_surfaces.conc"
GROUP_LEFT_CONC = local_path_to("raw", midthicknessfile.format("left"))
GROUP_MOTION_CONC = local_path_to("raw", "ADHD_DVARS_group_motion.conc")
GROUP_RIGHT_CONC = local_path_to("raw", midthicknessfile.format("right"))

# MATLAB Runtime Environment (MRE) directory depending on the host server
MRE_EXACLOUD = ("/home/exacloud/lustre1/fnl_lab/code/external/utilities/"
                "Matlab2016bRuntime/v91")
MRE_RUSHMORE = ("/mnt/max/shared/code/external/utilities/"
                "Matlab2016bRuntime/v91")

# If .dconn files' total size exceeds this many gigabytes, then warn user
WARNING_IF_DCONN_SIZE_EXCEEDS = 100

# Workbench command file depending on the host server
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
               .format(completion, step, datetime.now().strftime("%H:%M:%S on"
                       " %b %d, %Y")))

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
        type=valid_readable_dir,
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
        type=valid_readable_file_or_none,
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
        type=valid_readable_file,
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
        type=valid_readable_file,
        help="Path to workbench command file called 'wb_command'."
    )

    return validate_cli_args(parser.parse_args(), parser)


def valid_readable_file(path):
    """
    Throw argparse exception unless parameter is a valid readable filename 
    string. This is used instead of argparse.FileType("r") because the latter 
    leaves an open file handle, which has caused problems.
    :param path: Parameter to check if it represents a valid filename
    :return: String representing a valid filename
    """
    return validate(path, lambda x: os.access(x, os.R_OK),
                    os.path.abspath, "Cannot read file at {}")


def validate(path, is_real, make_valid, err_msg, prepare=None):
    """
    Parent/base function used by different type validation functions. Raises an
    argparse.ArgumentTypeError if the input path is somehow invalid.
    :param path: String to check if it represents a valid path 
    :param is_real: Function which returns true if and only if path is real
    :param make_valid: Function which returns a fully validated path
    :param err_msg: String to show to user to tell them what is invalid
    :param prepare: Function to create something at path before validation
    :return: path, but fully validated as pointing to the right file or dir
    """
    try:
        if prepare:
            prepare(path)
        assert is_real(path)
        return make_valid(path)
    except (OSError, TypeError, AssertionError, ValueError, 
            argparse.ArgumentTypeError):
        raise argparse.ArgumentTypeError(err_msg.format(path))


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


def valid_output_dir(path):
    """
    Try to create a directory to write files into at the given path, and throw 
    argparse exception if that fails
    :param path: String which is a valid (not necessarily real) folder path
    :return: String which is a validated absolute path to real writeable folder
    """
    return validate(path, lambda x: os.path.isdir(x) and os.access(x, os.W_OK),
                    lambda y: format_with_trailing_slash(valid_readable_file(y)),
                    "Cannot make output directory at {}", 
                    lambda z: os.makedirs(z, exist_ok=True))


def format_with_trailing_slash(path):
    """
    :param path: String representing a file path
    :return: path, but with one "/" at the end
    """
    return path if path[-1] == "/" else path + "/"


def valid_readable_file_or_none(path):
    """
    Throw argparse exception unless path either is "none" or points to a 
    readable file.
    :param path: Parameter to check if it represents a valid filename
    :return: String representing a valid filename, or the string "none"
    """
    return path if path == "none" else valid_readable_file(path)


def valid_readable_dir(path):
    """
    Throw argparse exception unless path is a string representing a valid path 
    to a readable directory.
    :param path: Parameter to check if it represents a valid filename
    :return: String representing a valid path to an existing readable directory
    """
    return validate(path, os.path.isdir, valid_readable_file,
                    "Cannot read directory at {}")


def validate_cli_args(cli_args, parser):
    """
    Check that all command line arguments will allow this script to work.
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: Validated command-line arguments argparse namespace
    """
    # Validate file path and directory arguments
    for arg in ("time_series", "dtseries", "template", "mre_dir"):
        cli_args = try_to_validate_file_arg(
            arg, globals()["validate_and_get_" + arg], cli_args, parser
        )

    # If the user did not give a wb_command, then try to get a default one, or
    # ask user for one if on an unknown server
    if not cli_args.wb_command:
        try:
            cli_args.wb_command = str(subprocess.check_output(
                ("which", "wb_command")
            ).decode("utf-8").strip())

        # Otherwise, use hardcoded wb_command depending on host server, or get
        # wb_command from user if host server is not known
        except subprocess.CalledProcessError:
            cli_args.wb_command = get_valid_server_arg(
                WB_RUSHMORE, WB_EXACLOUD, "workbench command"
            )

    # Validate that all .conc file inputs have an equal number of lines
    validate_concs_same_length([k for k, v in vars(cli_args).items()
                                if isinstance(v, str) and v[-5:] == ".conc"],
                               cli_args, parser)
    return cli_args


def try_to_validate_file_arg(arg_name, validate_fn, cli_args, parser):
    """
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    """
    try:
        if not getattr(cli_args, arg_name, None):
            setattr(cli_args, arg_name, validate_fn(cli_args))
        return cli_args
    except argparse.ArgumentTypeError as e:
        parser.error("Invalid {} path argument: {}".format(arg_name, e))


def validate_and_get_mre_dir(cli_args):
    """
    :param cli_args: argparse namespace with all command-line arguments
    :return: String which is a valid path to a (MRE) directory
    """
    return cli_args.mre_dir if cli_args.mre_dir else get_valid_server_arg(
        MRE_RUSHMORE, MRE_EXACLOUD, "MATLAB Runtime Environment directory"
    )


def validate_and_get_time_series(cli_args):
    """
    Infer whether series_file is dense or parcellated by reading series_file
    provided by user. Raise a parser error if the series_file is neither or is
    otherwise invalid.
    :param cli_args: argparse namespace with all command-line arguments
    :return: "dtseries" or "ptseries"
    """
    try:
        if cli_args.series_file[-11:] == "tseries.nii":
            time_series = cli_args.series_file[-12:-4]
        elif cli_args.series_file[-8:] == "conn.nii":
            time_series = cli_args.series_file[-9:-4] # + "tseries"
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

        # If user gave d/pconns, ensure that run mode is pairwise_corr
        elif time_series in ("dconn", "pconn"):
            if cli_args.scripts != ["pairwise_corr"]:
                err_msg = ("{} input arguments are only valid for running in "
                           "pairwise_corr mode.".format(time_series))
                raise argparse.ArgumentTypeError(err_msg)

        # If invalid time series is given, then tell user and crash
        elif time_series != "ptseries":
            err_msg = ("Time series file must either be a file with the "
                       "extensions .ptseries.nii, .dtseries.nii, .pconn.nii, "
                       "and .dconn.nii, or a .conc file with a list of paths "
                       "to files with those extensions.")
            raise argparse.ArgumentTypeError(err_msg)

    except OSError:
        raise argparse.ArgumentTypeError("Could not read time series file "
                                         "paths from " + cli_args.series_file)
    except IndexError:
        err_msg = ("Series file must be a path to a .nii file, or to a "
                   ".conc file which only has a list of .nii file paths.")
        raise argparse.ArgumentTypeError(err_msg)
    return time_series


def validate_and_get_dtseries(cli_args):
    """
    Validate the path to the .dtseries file used for outlier removal
    :param cli_args: argparse namespace with all command-line arguments
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
    return valid_readable_file_or_none(dtseries)


def get_valid_server_arg(if_rushmore, if_exacloud, param_name):
    """
    Validate and get a given path which varies depending on the server that
    this file is called from.
    :param if_rushmore: String which is a valid path on the Rushmore server
    :param if_exacloud: String which is a valid path on the Exacloud server
    :param param_name: String naming what the path to get actually means
    :return: String representing a valid path to server-dependent file or dir
    """
    host = socket.gethostname()
    if host == "rushmore":
        result = if_rushmore
    elif "exa" in host:
        result = if_exacloud

    # If on an unknown host, quit if user says to, or get valid path from user
    else:
        result = None
        print("No {} found on {} server.".format(param_name, host))
        result = get_required_user_input(
            ("Please enter a valid {} path, or enter 'n' to end the program: "
             .format(param_name)),
            ("Terminating {} because {} could not be found."
             .format(sys.argv[0], param_name)), result,
            lambda x: not os.access(x, os.R_OK)
        )
    return result


def get_required_user_input(input_msg, quit_msg, check_var, check_cond):
    """
    Keep asking user for input until they either give valid input or quit
    :param input_msg: String telling user what to input
    :param quit_msg: String to show user upon terminating the script
    :param check_var: String which will be the user's input
    :param check_cond: Function which accepts user input to quit or continue
    :return: check_var from user input, but validated
    """
    while check_cond(check_var):
        check_var = input(input_msg)
        if check_var.lower() == "n":
            sys.exit(quit_msg)
    return check_var


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
        warning = ("Warning: This will create about {} GB of .dconn files. Enter "
                    "'y' to continue or 'n' to quit: ".format(gb_of_files_created))
        get_required_user_input(warning, "Terminating script.", None, 
                                lambda x: x.lower() != "y" )


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
            template = valid_readable_file(cli_args.template)
            break
    return template


def get_template_file_path(cli_args):
    """
    Build and return the name of the file created by cifti_conn_template
    :param cli_args: argparse namespace with all command-line arguments
    :return: String name of file created by cifti_conn_template
    """
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
    return os.path.join(cli_args.output, "".join((
        os.path.basename(cli_args.series_file), fd_str, motion_and_smooth, 
        cli_args.smoothing_kernel, '_AVG.', cli_args.time_series[0], "conn.nii"
    )))


def validate_concs_same_length(conc_file_args, cli_args, parser):
    """
    Throw argparse error unless all .conc files have the same number of lines
    :param conc_file_args: List of names of .conc file arguments in cli_args
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    valid paths to .conc files which each contain a list of file paths
    :return: N/A
    """
    # Get the number of lines in every .conc file
    prev_file_lines = -1
    for conc_arg in conc_file_args:
        conc_path = getattr(cli_args, conc_arg)
        with open(conc_path) as conc_file:
            lines = 0
            while conc_file.readline():
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
    all necessary parameters for cifti_conn_matrix and cifti_conn_template.
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
            cli_args.dtseries]


def get_script_filename(index):
    """
    :param index: Integer showing which script named in CHOICES_TO_RUN to use
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
        print("Warning: CIFTI conn matrix files will not be deleted until {} "
              "finishes running.".format(CHOICES_TO_RUN[2]))
        keep_conn_matrices = "1"

    # Call cifti_conn_template script
    subprocess.check_call((local_path_to("src", get_script_filename(1)),
                           *get_matrix_or_template_parameters(cli_args),
                           keep_conn_matrices,
                           cli_args.template))


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
    