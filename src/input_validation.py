#! /usr/bin/env python3

"""
Input validation functions source 
Greg Conan: conan@ohsu.edu
Created 2020-03-04
Updated 2021-03-29
"""

##################################
#
# Input validation functions for CIFTI connectivity wrapper:
# cifti_conn_wrapper.py
#
##################################

### Imports
import argparse
import os
import subprocess
import sys

### Constants

# Names of scripts to run from main script
CHOICES_TO_RUN = ["matrix", "template", "pairwise_corr"]

# If .dconn files' total size exceeds this many gigabytes, then warn user
WARNING_IF_DCONN_SIZE_EXCEEDS = 100

### Functions (listed in the order that cifti_conn_wrapper calls them)


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
    for arg in ("time_series", "dtseries", "template",
                "mre_dir", "wb_command"):
        cli_args = try_to_validate_file_arg(
            arg, globals()["validate_and_get_" + arg], cli_args, parser
        )

    # If user is smoothing, then ensure that midthickness surfaces were given
    if getattr(cli_args, "smoothing_kernel", "none") != "none":
        for arg in ("left", "right"):
            file_arg = getattr(cli_args, arg, None)
            if file_arg == "none" or not os.path.exists(file_arg):
                parser.error("Invalid --{0} argument. --{0} is "
                             "needed for smoothing.".format(arg))

    # Validate that all .conc file inputs have an equal number of lines
    validate_concs_same_length([k for k, v in vars(cli_args).items()
                                if isinstance(v, str) and v[-5:] == ".conc"],
                               cli_args, parser)
    return cli_args


def try_to_validate_file_arg(arg_name, validate_fn, cli_args, parser):
    """
    :param arg_name: String naming the cli_arg to try to validate
    :param validate_fn: Function to validate the argument in cli_args
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: cli_args, but with a valid arg_name attribute/argument
    """
    try:
        if not getattr(cli_args, arg_name, None):
            setattr(cli_args, arg_name, validate_fn(cli_args))
        return cli_args
    except argparse.ArgumentTypeError as e:
        parser.error("Invalid {} path argument: {}".format(arg_name, e))


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
            time_series = cli_args.series_file[-9:-4]
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
        get_required_user_input(warning, "Terminating script.",
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


def validate_and_get_mre_dir(cli_args):
    """
    :param cli_args: argparse namespace with all command-line arguments
    :return: String which is a valid path to a (MRE) directory
    """
    return (cli_args.mre_dir if cli_args.mre_dir else
            get_valid_server_arg("MATLAB Runtime Environment directory"))


def validate_and_get_wb_command(cli_args):
    """
    :param cli_args: argparse namespace with all command-line arguments
    :return: String which is a valid path to a workbench command executable
    """
    if not cli_args.wb_command:
        try:
            cli_args.wb_command = str(subprocess.check_output(
                ("which", "wb_command")
            ).decode("utf-8").strip())
        except subprocess.CalledProcessError:
            cli_args.wb_command = get_valid_server_arg("workbench command")
    return cli_args.wb_command


def get_valid_server_arg(param_name):
    """
    Validate and get a given path which varies depending on the server that
    this file is called from. Asks user to quit or to provide a valid path
    :param param_name: String naming what the path to get actually means
    :return: String representing a valid path to server-dependent file or dir
    """
    return get_required_user_input(
        ("No {0} found on your server. Please enter a valid {0} "
            "path, or enter 'n' to end the program: ".format(param_name)),
        ("Terminating {} because {} could not be found."
         .format(sys.argv[0], param_name)),
        lambda x: not os.access(x, os.R_OK)
    )


def get_required_user_input(input_msg, quit_msg, check_cond):
    """
    Keep asking user for input until they either give valid input or quit
    :param input_msg: String telling user what to input
    :param quit_msg: String to show user upon terminating the script
    :param check_cond: Function which accepts user input to quit or continue
    :return: User input, but validated, unless user chooses to quit
    """
    check_var = input(input_msg)
    while check_cond(check_var):
        if check_var.lower() == "n":
            sys.exit(quit_msg)
        check_var = input(input_msg)
    return check_var


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