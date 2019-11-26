function cifti_conn_template_for_wrapper(wb_command, dt_or_ptseries_conc_file, series, motion_file, FD_threshold, TR, minutes_limit, smoothing_kernel, left_surface_file, right_surface_file, bit8, remove_outliers, additional_mask, make_dconn_conc, output_directory, dtseries_conc, keep_conn_matrices, template_file)
% wb_command = workbench command path
% dt_or_ptseries_conc_file = dense timeseries or parcellated timeseries conc
% file (i.e. text file with paths to each file being included
% series = 'dtseries' or 'ptseries'; specify if files are dense or parcellated
% motion_file = conc file that points to FNL motion mat files for each dt or
% ptseries; Can also be a text file with a vectors of 1 and 0s;set to 'none' if you
% would like to include all frames
% FD threshold = specify motion threshold (e.g. 0.2); set to 0 if you'd
% like to include all frames, but make sure motion file is also set to
% 'none'
% TR = The TR of your study
% minutes_limit= the amount of remaining data you are requring to include a
% subject. Note: all subjects should conform to your minutes limit, if not
% the code will stop and tell you.
% smoothing_kernal = do you want to smooth your data? if so put in the
% kernal here, if not put 'none'
% left_surface_file = if you are smoothing you'll need concs of each subjects midthickness file
% right_surface_file = if you are smoothing you'll need concs of each subjects midthickness file
% bit8 = if you want to reduce the size of the file choose bit8 = 1, if not
% zero (of note as of 11/20/2017 having this on for this code has not been
% tested - df)
% keep_conn_matrices = set to 1 to keep dconn.nii, set to 0 to delete them.
% removeoutliers. 0 or 1.  If 1, outliers are removed from the BOLD signal.  If 0 frames will only be censored by the FD threshold.
% Addional mask on top of the FD threshold.  This mask can be used to extract rests between blocks of task. This vector of 1s and zeros should be the same length as your dtseries.
% make_dconn_conc - make a list of connectivity matrices that were created by this script. If this is 'none', no list will be created; otherwise it should be a valid path.
%  NOTE: this will delete dconns that have been previously made if the name
%  matches the expected output.

%% 3/12/18 -RH 
%Added an option to set minutes limit to none.  This will use all frames
%below a specific threshold.
%Added an arguement to keep .dconns/pconns after being made.

%% Change strings to numbers if not already numbers  (i.e.if using the shell command)
if ~exist('output_directory', 'var')
    output_directory = [];
end

%Check to make sure that TR, bit8, keep_conn_matrices, and FD are numbers
FD_threshold = make_numeric(FD_threshold, 'FD_threshold');
TR = make_numeric(TR, 'TR');
bit8 = make_numeric(bit8, 'bit8');
keep_conn_matrices = make_numeric(keep_conn_matrices, 'keep_conn_matrices');

%Check to make sure that minutes limit and smoothing kernels are number
%(unless you've set them to 'none')
smoothing_kernel = make_numeric_or_none(smoothing_kernel, 'smoothing_kernel', 'No smoothing kernel supplied');
minutes_limit = make_numeric_or_none(minutes_limit, 'minutes_limit', 'Warning minutes limit is set to none.  Subjects might have different number of data point in individual matrices');

%% Load concatenated paths (i.e. paths to ciftis)
%check to see if there one subject or a a list of subjects in conc file.
conc = strsplit(dt_or_ptseries_conc_file, '.');
conc = char(conc(end));
if strcmpi('conc', conc)
    A = importdata(dt_or_ptseries_conc_file);
else
    A = {dt_or_ptseries_conc_file};
end

for i = 1:length(A)
    if ~exist(A{i})
        NOTE = ['Series ' num2str(i) ' does not exist']
        return
    end
end
disp('All series files exist continuing ...')

%% Combine dtseries for a given wave using d(p)tseries_cifti.conc
if series == 'ptseries'
    suffix = 'ptseries.nii';
    suffix2 = 'pconn.nii';
elseif series == 'dtseries'
    suffix = 'dtseries.nii';
    suffix2 = 'dconn.nii';
else
    'series needs to be "ptseries" or "dtseries"';
    return
end

%% Load Motion vector paths
if is_none(motion_file)
    disp('No motion files, will use all frames to generate matrices')
    B = {motion_file};
else
    if strcmpi('conc', conc)
        B = importdata(motion_file);
    else
        B = {motion_file};
    end
    
    % Check that all motion files in conc file exist
    for i = 1:length(B)
        if ~exist(B{i})
            NOTE = ['motion file ' num2str(i) ' does not exist']
            return
        end
    end
    disp('All motion files exist continuing ...')
end

%% Check to see if surface files exist
%check to make sure there is a list of subjects in conc file.
conc = strsplit(dt_or_ptseries_conc_file, '.');
conc = char(conc(end));
if is_none(smoothing_kernel)
    Note = ['No smoothing, continuing...'];
    if strcmpi('conc', conc)
        left = left_surface_file;
        right = right_surface_file;
    else
        disp(['Need to input a conc file with several subjects or ' ...
              'dtseries to average'])
        return
    end
else
    C = get_surface_files('left', left_surface_file, conc);
    D = get_surface_files('right', right_surface_file, conc);
    left = C{1};
    right = D{1};
end

%% Generate Average Matrix
% This section and the next one were rewritten by Greg Conan on 2019-11-06 
% to trim redundant code and modularize it into separate functions.
% Make first dconn.nii then rename

% Map a matrix filepath to whether that matrix has been added to the avg
already_added = containers.Map('KeyType', 'char', 'ValueType', 'logical');

% Get dconn paths and file names of matrix, dconns, motion, output conc
matrix_filename = get_matrix_filename(A, 1, smoothing_kernel, suffix);
dconn_paths = get_dconn_paths(FD_threshold, minutes_limit, ...
    matrix_filename, motion_file, output_directory, suffix2);
motion_B = get_motion_file(motion_file, B, 1);
output_conc = get_output_conc(dt_or_ptseries_conc_file, output_directory);

% Make average matrix unless one already exists
paths = dir(dconn_paths);
if isempty(paths)
    dconn_paths = char(cifti_conn_matrix_for_wrapper( ...
        wb_command, A{1}, series, motion_B, num2str(FD_threshold), ...
        TR, minutes_limit, smoothing_kernel, left, right, bit8, ...
        remove_outliers, additional_mask, make_dconn_conc, ...
        output_directory, dtseries_conc));
else
    dconn_paths = char(get_already_existing_matrix(paths, already_added));
end

% Rename average dconn, do cifti math, clean variables, then do next
rename_avg_dconn(dconn_paths, output_conc, keep_conn_matrices, suffix2);
already_added(dconn_paths) = 1;

%% Build rest of template
for i = 2:length(A)
    
    % Get matrix filename, dconn paths, and motion file
    matrix_filename = get_matrix_filename(A, i, smoothing_kernel, suffix);
    dconn_paths = char(get_dconn_paths(FD_threshold, minutes_limit, ...
        matrix_filename, motion_file, output_directory, suffix2));
    motion_B = get_motion_file(motion_file, B, i);
    
    % Make average matrix unless one already exists
    paths = dir(dconn_paths);
    if isempty(paths)
        dconn_paths = char(cifti_conn_matrix_for_wrapper( ...
            wb_command, A{i}, series, motion_B, num2str(FD_threshold), ...
            TR, minutes_limit, smoothing_kernel, left, right, bit8, ...
            remove_outliers, additional_mask, make_dconn_conc, ...
            output_directory, dtseries_conc));
    else
        dconn_paths = char(get_already_existing_matrix(paths, ...
                                                       already_added));
    end
    
    % Rename average dconn, do cifti math, clean variables, then do next
    rename_avg_dconn(dconn_paths, output_conc, keep_conn_matrices, ...
                     suffix2);
    cifti_math_and_cleanup(dconn_paths, output_conc, ...
        keep_conn_matrices, suffix2, wb_command);
    already_added(dconn_paths) = 1;
end

%% Generate Mean by dividing by N
N = num2str(length(A)); % number of observations

%mean = s1/N
cmd = [wb_command ' -cifti-math " s1/' N '" ' output_conc '_AVG.' ...
       suffix2 ' -var s1 ' output_conc '_AVG.' suffix2];
execute_and_clear(cmd);

%% Rename file to display chosen options
cmd = ['mv ' output_conc '_AVG.' suffix2 ' ' template_file];
system(cmd);
disp('Done making average template');
end


function num_var = make_numeric(num_var, var_name)
    % Convert num_var to a number if it is not already a number or 'none'
    if is_none(num_var) && ~isnumeric(num_var)
        disp([var_name ' passed in as a string. Converting to numeric.'])
        num_var = str2num(num_var);
    end
end


function num_var = make_numeric_or_none(num_var, var_name, varargin)
    % Convert num_var to a number if it is not already a number or 'none'
    if is_none(num_var)
        if nargin > 0
            disp(varargin)
        end
    else
        num_var = make_numeric(num_var, var_name);
    end
end


function surface_files = get_surface_files(LR, surface_conc, conc)
    % Get all paths to surface files
    if strcmpi('conc', conc)
        surface_files = importdata(surface_conc); 
    else
        surface_files = {surface_conc};
    end
    
    % Check to make sure surface files exist
    for i = 1:length(surface_files)
        if ~exist(surface_files{i}, 'file')
            disp(['Subject ' LR ' surface ' num2str(i) ' does not exist'])
            return
        end
    end
    disp(['All ' LR ' surface files for smoothing exist. Continuing ...'])
end


function matrix_filename = get_matrix_filename(A, ix, smoothing, suffix)
    % Get name of connectivity matrix file
    if is_none(smoothing)
        A_smoothed = A{ix};
    else
        A_smoothed = [A{ix}(1:length(A{ix})-13) '_SMOOTHED_' ...
                      num2str(smoothing) '.' suffix];
    end
    [~, matrix_filename, matrix_ext] = fileparts(char(A_smoothed));
    matrix_filename = [matrix_filename matrix_ext];
end


function dconn_paths = get_dconn_paths(fd, minutes, matrix, ...
        motion_file, output_dir, ext)
    % Get path to dconns
    if is_none(motion_file)
        dconn_paths = char(strcat(output_dir, matrix, ...
            '_all_frames_at_FD_', {motion_file}, '*.', {ext}));
    else
        if is_none(minutes)
            dconn_paths = char(strcat(output_dir, matrix, ...
                '_all_frames_at_FD_', {num2str(fd)}, '*.', {ext}));
        else
            dconn_paths = char(strcat(output_dir, matrix, '_', ...
                {num2str(minutes)}, '_minutes_of_data_at_FD_', ...
                {num2str(fd)}, '*.', {ext}));
        end
    end
end


function motion_file = get_motion_file(motion_file, B, ix)
    % Get the path to a motion file, or 'none' if there is no motion file
    if ~is_none(motion_file)
        motion_file = B{ix};
    end
end


function output_conc = get_output_conc(input_conc, output_directory)
    % Get filename of output .conc file listing paths to conn matrices
    [~, series_conc_name, series_conc_ext] = fileparts(char(input_conc));
    output_conc = [output_directory series_conc_name series_conc_ext]
end


function matrix_file = get_already_existing_matrix(paths, matrices_added)
    % Get the name of the matrix file which would have been created, but
    % has already been created. Handle 2 matrices which have the same name
    % but were in different folders separately.
    for m = 1:size(paths)
        matrix_file = [paths(m).folder filesep paths(m).name];
        if isKey(matrices_added, matrix_file)
            clear matrix_file
        else
            break
        end
    end
    if exist(matrix_file)
        disp(['conn.nii already exists for this subject. Renaming to '...
              'AVG for first subject.'])
    else
        disp(['Error: Mismatch in  number of matrix files: ' matrix_file])
        return
    end
end


function rename_avg_dconn(paths, output_conc, keep_conn_matrices, suffix2)
    % Add the word "AVG" to average matrix to not overwrite the originals
    if ~keep_conn_matrices
        cmd = ['mv ' paths ' ' output_conc '_AVG.' suffix2];
    else
        cmd = ['cp ' paths ' ' output_conc '_AVG.' suffix2];
    end
    system(cmd);
end


function cifti_math_and_cleanup(dconn_paths, output_conc, ...
                                keep_conn_matrices, suffix2, wb_command)
    % Run cifti math on average matrices, and then delete dconn files 
    % unless keep_conn_matrices is true
    cmd = [wb_command ' -cifti-math  " M1 + M2 " ' output_conc '_AVG.' ...
           suffix2 ' -var M1 ' output_conc '_AVG.' suffix2 ' -var M2 ' ...
           dconn_paths];
    execute_and_clear(cmd);
    if keep_conn_matrices == 0
        delete(dconn_paths);
    else
    end
end


function execute_and_clear(cmd)
    % Given a string that can be executed on the command line, execute it
    % and then clear the variable holding that string.
    tic;
    system(cmd);
    toc;
    clear cmd
end


function result = is_none(param)
    % Return 1 if param is a string saying none; otherwise return 0
    result = strcmpi(param,'none');
end