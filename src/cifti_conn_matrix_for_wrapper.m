function outfile = cifti_conn_matrix_for_wrapper(wb_command, ...
    dt_or_ptseries_conc_file, series, motion_file, FD_threshold, TR, ...
    minutes_limit, smoothing_kernel, left_surface_file, ...
    right_surface_file, bit8, remove_outliers, additional_mask, ...
    make_dconn_conc, output_directory, dtseries_conc, vector_censor)
% This script is a modularized version of cifti_conn_matrix_exaversion, 
% except that it also accepts and uses these parameters: remove_outliers, 
% additional_mask, and make_dconn_conc.

%% Parameter definitions
% dt_or_ptseries_conc_file = dense timeseries or parcellated timeseries
%     conc file (i.e. text file with paths to each file being examined
% series = 'dtseries' if files are parcellated files or 'ptseries' if not
% motion_file = conc file that points to FNL motion mat files for each dt 
%     or ptseries (note: if no motion file to be used, type 'none')
% FD threshold = specify motion threshold (e.g. 0.2)
% TR = Repetition time of your data
% minutes_limit = specify the number of minutes to be used to generate the
%     correlation matrix
% smoothing_kernal = specify smoothing kernal (note: if no smoothing file 
%     to be used, type 'none')
% bit8 = set to 1 if you want to make the outputs smaller in 8bit

%% Input validation
orig_conc_folders = split(dt_or_ptseries_conc_file, filesep);
orig_conc_file = char(orig_conc_folders(end));
no_output = 0;
if ~exist('output_directory', 'var')
    output_directory = [];
    no_output = 1;
    disp('All series files exist. Continuing ...')
else
    output_directory = char(output_directory);
end
if isempty(output_directory)
    no_output = 1;
end

% Convert number variables to cnumbers (only within compiled code)
% Check to ensure that input variables which should be numbers are numbers
bit8 = make_numeric(bit8, 'Beta 8');
FD_threshold = make_numeric(FD_threshold, 'FD threshold');
TR = make_numeric(TR, 'TR');
minutes_limit = make_numeric(minutes_limit, 'Minutes limit');
remove_outliers = make_numeric(remove_outliers, 'Outlier removal');
smoothing_kernel = make_numeric(smoothing_kernel, 'Smoothing kernel');

if (~strcmpi(smoothing_kernel, 'none')) && (strcmpi(series, 'ptseries'))
    disp('Check your settings. Smoothing not allowed on ptseries.');
    return
end

%% Load concatenated paths (i.e. paths to ciftis)
[A, conc] = get_paths_from_conc(dt_or_ptseries_conc_file);
if ~strcmpi(dtseries_conc, 'none')
    [dtseries_E, ~] = get_paths_from_conc(dtseries_conc);
else
    dtseries_E = 'none';
end
    
%% Set file extensions
if strcmpi(series, 'ptseries')
    suffix = 'ptseries.nii';
    suffix2 = 'pconn.nii';
    
    % Set default precision
    % Connectivity file sizes are NOT reduced for ptseries
    output_precision = ' ';
    
elseif strcmpi(series, 'dtseries')
    suffix = 'dtseries.nii';
    suffix2 = 'dconn.nii';
    if exist('bit8', 'var')
        bit8 = make_numeric(bit8, 'Beta 8');
        
        % Set default precision. 
        if bit8  % Connectivity file sizes are reduced for dtseries
            output_precision = [' -cifti-output-datatype INT8 ' ...
                                '-cifti-output-range -1.0 1.0'];
        else
            output_precision = ' ';
        end
    else
        output_precision = ' ';
    end
else
    'Series needs to be "ptseries" or "dtseries"';
    return
end

%% preallocate memory if you want to do smoothing
if ~strcmpi(smoothing_kernel, 'none')
    disp(['Smoothing kernel: ' smoothing_kernel])
    A_smoothed = num2cell(zeros(size(A)));
    C = get_surface_files('left', left_surface_file, conc);
    D = get_surface_files('right', right_surface_file, conc);
end

%% Generate Motion Vectors and correlation matrix
if strcmpi(motion_file, 'none') % run this if no motion censoring
    disp('No motion files, will use all frames to generate matrices')
    for i = 1:length(A)
        [orig_cifti_filename, input_directory, ...
            ~] = get_folder_params(A{i}, filesep);
        if no_output
            output_directory = char(input_directory);
        end
       
        % Run minutes limit calculation
        if strcmpi(smoothing_kernel, 'none')
            left = 'none';
            right = 'none';            
        else
            left = C{i};
            right = D{i};
        end
        min_file_end = '_cifti_censor_FD_vector_All_Good_Frames.txt';
        [A_smoothed{i}, outfile] = minutes_limit_calculation(...
            FD_threshold, input_directory, left, '_all_frames_at_', ...
            min_file_end, motion_file, orig_cifti_filename, 'none', ...
            output_directory, output_precision, right, ...
            smoothing_kernel, suffix, suffix2, wb_command);
    end
else % use motion censoring
    % rename temp file so it isn't overwritten when running in parallel
    stdev_temp_filename=[output_directory orig_conc_file '_temp.txt'];
    B = import_conc_or_file(motion_file, conc);
    
    % Verify that motion files exist and their length matches timeseries'
    if length(A) ~= length(B)
        disp(['Length of motion .conc file does not match length of ' ...
              'series .conc file.'])
    end
    verify_paths_exist(B, 'Motion');
    
    %% Generate motion vector for subject and generate correlation matrix
    for i = 1:length(B)
        motion_exten = strsplit(B{i}, '.');
        motion_exten = char(motion_exten(end));
        
        % Get input folder parameters
        [orig_cifti_filename, input_directory, ...
            ~] = get_folder_params(A{i}, filesep);
        if no_output
            output_directory = char(input_directory);
        end
        other_motion_mask = ~strcmp('mat', motion_exten);
        if strcmpi(dtseries_E, 'none')
            outlier_rmv_file = 'none';
        else
            outlier_rmv_file = dtseries_E{i};
        end
        
        % use an external mask (.txt) instead of calculating the mask here
        [orig_motion_filename, inputB_directory, ...
            ~] = get_folder_params(B{i}, filesep);
        if other_motion_mask 
            FDvec = B{i};
            FDvec = remove_outliers_if(remove_outliers, FDvec, ...
                outlier_rmv_file, stdev_temp_filename, wb_command, 30); 
            
        else  % Use power 2014 motion
            FDvec = get_FDvec(B{i}, FD_threshold, vector_censor);
            if no_output
                output_directory = char(inputB_directory);
            end
            FDvec = add_mask_to_FDvec(FDvec, additional_mask);
            FDvec = remove_outliers_if(remove_outliers, FDvec, ...
                outlier_rmv_file, stdev_temp_filename, wb_command, 30);
        end  % [input_directory orig_cifti_filename] was instead of rmv
        
        % Assign variables for minutes limit calculation
        if strcmpi(minutes_limit, 'none')            
            min_cmd_part = '_all_frames_at_';
            min_file_end = '_cifti_censor_FD_vector_All_Good_Frames.txt';
            save_out_file(output_directory, orig_motion_filename, ...
                          FD_threshold, min_file_end, FDvec);
        else
            % Get the number of good minutes in your data
            good_frames_idx = find(FDvec == 1);
            good_minutes = (length(good_frames_idx)*TR)/60;
            
            % If there is less than 30 seconds of good data, then do not
            % generate the correlation matrix
            if good_minutes < 0.5
                disp(['Subject ' num2str(i) ' has less than 30 ' ...
                      'seconds of good data.'])
            
            % if there is not enough data for your subject, just generate 
            % the matrix with all available frames
            elseif minutes_limit > good_minutes 
                min_cmd_part = '_all_frames_at_';
                min_file_end = ['_cifti_censor_FD_vector_All_Good_' ...
                                'Frames.txt'];
                save_out_file(output_directory, orig_motion_filename, ...
                              FD_threshold, min_file_end, FDvec);
            
            % if there is enough data, match the amount of data used for
            % your subject matrix to the minutes_limit
            elseif minutes_limit <= good_minutes
                
                % Number of good frames to randomly pull
                rand_good_frames = sort(randperm(length( ...
                    good_frames_idx), round(minutes_limit*60/TR)));
                FDvec_cut = zeros(length(FDvec), 1);
                ones_idx = good_frames_idx(rand_good_frames);
                
                % The new vector that should match good frames with the
                % minutes limit
                FDvec_cut(ones_idx) = 1; 
                min_cmd_part = ['_' num2str(minutes_limit) ...
                                '_minutes_of_data_at_'];
                min_file_end = ['_cifti_censor_FD_vector' min_cmd_part ...
                                num2str(FD_threshold) '_threshold.txt'];
                save_out_file(output_directory, orig_motion_filename, ...
                              FD_threshold, min_file_end, FDvec_cut);
            else
                disp(['Something is wrong about the number of good ' ...
                      'minutes calculation.'])
            end
        end
        
        % Run minutes limit calculation
        if strcmpi(minutes_limit, 'none') || good_minutes >= 0.5
            if strcmpi(smoothing_kernel, 'none')
                left = 'none';
                right = 'none';
            else
                left = C{i};
                right = D{i};
            end
            [A_smoothed{i}, outfile] = minutes_limit_calculation(...
                FD_threshold, input_directory, left, min_cmd_part, ...
                min_file_end, motion_file, orig_cifti_filename, ...
                orig_motion_filename, output_directory, ...
                output_precision, right, smoothing_kernel, suffix, ...
                suffix2, wb_command);
            if strcmpi(smoothing_kernel, 'none')
                clear A_smoothed
            end
        end        
    end
end

%% Generate final conc file for other code (e.g. cifti_conn_pairwise_corr)
% Greg Conan initially modularized everything from here down on 2019-11-11
if make_dconn_conc
    disp(['Done making (p or d)conn.nii files for subjects. Making ' ...
          'output .conc files now.'])
            
    %% Save all paths into different conc files
    if strcmpi(smoothing_kernel, 'none')
        smoothing = 'none';
    else  % use smoothed data
        smoothing = A_smoothed;
    end
    if strcmpi(motion_file, 'none')
        motion_file_obj = 'none';
        stdev = 'none';
    else
        motion_file_obj = B;
        stdev = stdev_temp_filename;
    end    
    
    % Sort paths into separate lists based on how their total amount of
    % good minutes compares to the minutes_limit threshold
    [dconn_paths_all_frames, dconn_paths_all_frames_at_thresh, ...
     dconn_paths_all_frames_at_thresh_min_lim, ...
     subjectswithoutenoughdata] = sort_paths(A, motion_file_obj, ...
        dtseries_E, FD_threshold, minutes_limit, motion_file, ...
        no_output, output_directory, remove_outliers, smoothing, ...
        stdev, suffix2, TR, wb_command, additional_mask, vector_censor);
    
    % Save all of those path lists into .conc files
    finish_and_save_concs(A, FD_threshold, minutes_limit, ...
        motion_file, no_output, output_directory, ...
        dt_or_ptseries_conc_file, series(1), smoothing_kernel, ...
        dconn_paths_all_frames, dconn_paths_all_frames_at_thresh, ...
        dconn_paths_all_frames_at_thresh_min_lim, ...
        subjectswithoutenoughdata)
end
disp('Done making output concs.')
end


%% Other functions used by the main function


function num_var = make_numeric(num_var, var_name)
    % Convert num_var to a number if it is not already a number or 'none'
    if ~strcmpi(num_var, 'none') && ~isnumeric(num_var)
        disp([var_name ' passed in as a string. Converting to numeric.'])
        num_var = str2num(num_var);
    end
end


function [paths, conc] = get_paths_from_conc(conc_file)
    % Check if conc_file is 1 subject or a list of them; verify all subjects
    conc = strsplit(conc_file, '.');
    conc = char(conc(end));
    paths = import_conc_or_file(conc_file, conc);
    verify_paths_exist(paths, 'Subject series');
end


function verify_paths_exist(paths, to_what)
    % Validate that all paths point to real files
    for i = 1:length(paths)
        if ~exist(paths{i}, 'file')
            disp([to_what ' file ' num2str(i) ' does not exist.'])
            return
        end
    end
    disp([to_what ' files all exist. Continuing...'])
end


function from_files = import_conc_or_file(path, conc)
    % Import data from a file, or from a .conc list of files
    if strcmpi('conc', conc)
        from_files = importdata(path);
    else
        from_files = {path};
    end
end


function surface_files = get_surface_files(LR, surface_conc, conc)
    % Get and validate all paths to surface files
    surface_files = import_conc_or_file(surface_conc, conc);
    verify_paths_exist(surface_files, ['Subject ' LR ' surface']);
end


function [cifti, in_dir, folders] = get_folder_params(A_i, sep)
    % Return the path to the original cifti file, the path to the input
    % directory, and an object holding all parts of A_i's file path
    folders = split(A_i, sep);
    cifti = char(folders(end));
    folder_input = join(folders(1:end-1), sep);
    in_dir = [char(folder_input) sep];
end


function execute_display_and_clear(cmd, msg)
    % Given a string that can be executed on the command line, execute it
    % and then clear the variable holding that string. Also display msg.
    tic;
    disp(msg)
    system(cmd);
    toc;
    clear cmd
end


function [A_smoothed_i, outfile] = minutes_limit_calculation(...
        FD_threshold, input_dir, left, min_cmd_part, min_file_end, ...
        motion_file, orig_cifti_file, orig_motion_file, output_dir, ...
        precision, right, smoothing, suffix, suffix2, wb_command)
    % Calculate the number of good minutes and get the output file name
    if strcmpi(motion_file, 'none')
       min_cmd_motion = [min_cmd_part 'FD_' motion_file]; 
       min_file_end = '.txt';
       weights = '';
    else 
       min_cmd_motion = [min_cmd_part 'FD_' num2str(FD_threshold)];
    end
    if strcmpi(smoothing, 'none')
        outfile = [char(output_dir) char(orig_cifti_file) ...
                   min_cmd_motion '.' suffix2];
        A_smoothed_i = 'none';
        
    else % Smooth
        A_smoothed_i = [char(output_dir) char(orig_cifti_file( ...
            1:length(char(orig_cifti_file))-13)) '_SMOOTHED_' ...
            num2str(smoothing) '.' suffix];
        outfile = [A_smoothed_i min_cmd_motion '.' suffix2];
        if ~exist(A_smoothed_i, 'file')
            cmd = [wb_command ' -cifti-smoothing ' char(input_dir) ...
                   char(orig_cifti_file) ' ' num2str(smoothing) ' ' ...
                   num2str(smoothing) ' COLUMN ' char(output_dir) ...
                   char(orig_cifti_file(1:length(char(orig_cifti_file) ...
                   )-13)) '_SMOOTHED_' num2str(smoothing) '.' suffix ...
                   ' -left-surface ' left '  -right-surface ' right];
            execute_display_and_clear(cmd, cmd);
        else %smoothed series already exists
            disp('Smoothed series already created for this subject')
        end
    end
    
    % Check to see if the file already exists
    if ~exist(outfile, 'file')
        if strcmpi(smoothing, 'none')
            if ~exist('weights', 'var')
                weights = [' -weights ' char(output_dir) ...
                           char(orig_motion_file) '_' ...
                           num2str(FD_threshold) min_file_end];
            end
            cmd = [wb_command ' ' precision ...
                   ' -cifti-correlation ' char(input_dir) ...
                   char(orig_cifti_file) ' ' char(output_dir) ...
                   char(orig_cifti_file) min_cmd_motion '.' suffix2 ...
                   weights ' -fisher-z'];
        else
            if ~exist('weights', 'var')
                weights = [' -weights ' char(output_dir) ...
                           char(orig_motion_file) '_' ...
                           num2str(FD_threshold) min_file_end];
            end
            cmd = [wb_command ' ' precision ' -cifti-correlation ' ...
                   A_smoothed_i ' ' A_smoothed_i min_cmd_motion '.' ...
                   suffix2 weights ' -fisher-z'];
        end
        execute_display_and_clear(cmd, ['Running wb_command cifti-' ...
            'correlation. This may take a few minutes.']);
    else
        disp([outfile ' already exists'])
    end
end


function [all_frames, at_thresh, at_thr_min_lim, without] = sort_paths(...
        A, B, dtseries_E, FD_threshold, minutes_limit, motion_file, ...
        no_output, output_directory, remove_outliers, smoothed, ...
        stdev_temp_filename, suffix2, TR, wb_command, mask, vector_censor)
    % Build and return 4 lists of paths to connectivity matrix files, such
    % that each list contains matrices with a certain amount of good data 
    % compared to the minutes_limit and FD_threshold:
    without = [];        % All subject sessions with < 30 sec of good data
    at_thr_min_lim = []; % Other frames at FD_threshold and minutes_limit
    at_thresh = [];      % Other frames at FD_threshold
    all_frames = [];     % Other frames regardless of minutes_limit & FD
    
    if no_output
        output_directory = char(input_directory);
    end
    
    % Sort all connectivity matrix file paths
    for i = 1:length(A)
        
        % Input validation
        [orig_cifti_path, input_directory, ...
            ~] = get_folder_params(A{i}, filesep);
        if strcmpi(smoothed, 'none')
            A_dir = strcat(char(output_directory), char(orig_cifti_path));
        else
            A_dir = smoothed{i};
        end
        
        % Sort subject sessions 
        if strcmpi(motion_file, 'none')
            all_frames = [all_frames; strcat({A_dir}, ...
                '_all_frames_at_FD_', motion_file, '.', {suffix2})];
        else
            FDvec = get_FDvec(B{i}, FD_threshold, vector_censor);
            if strcmpi(dtseries_E, 'none')
                dt = dtseries_E;
            else
                dt = dtseries_E{i};
            end
            FDvec = additional_frame_removal(mask, dt, FDvec, ...
                        input_directory, orig_cifti_path, ...
                        remove_outliers, stdev_temp_filename, wb_command);
            [at_thresh, at_thr_min_lim, without] = collect_conns( ...
                A_dir, suffix2, FD_threshold, FDvec, minutes_limit, ...
                TR, at_thresh, at_thr_min_lim, without);
        end
    end
end


function FDvec = additional_frame_removal(additional_mask, dt, FDvec, ...
    input_directory, orig_cifti, remove_outliers, stdev_temp_filename, ...
    wb_command)
    % additional frame removal depending on additional_mask
    if strcmpi(additional_mask, 'none')
        FDvec = remove_frames_from_FDvec([input_directory orig_cifti], ...
                    FDvec, stdev_temp_filename, wb_command, 30);
    else
        FDvec = add_mask_to_FDvec(FDvec, additional_mask);
        FDvec = remove_outliers_if(remove_outliers, FDvec, dt, ...
                                   stdev_temp_filename, wb_command, 30);
    end 
end


function FDvec = get_FDvec(file_with_FDvec, fd, vector_censor)
    % Get FDvec from the file at file_with_FDvec
    load(file_with_FDvec)
    allFD = zeros(1, length(motion_data)); % motion_data is from file
    for j = 1:length(motion_data)
        allFD(j) = motion_data{j}.FD_threshold;
    end
    FDidx = find(allFD == fd);
    
    % Added by Greg Conan 2020-03-04 for different vector_censor options
    vector_field = [vector_censor '_removal'];
    if isfield(motion_data{FDidx}, vector_field)
        FDvec = motion_data{FDidx}.(vector_field);
    else
        disp(['Error: ' file_with_FDvec ' does not include ' ...
              vector_field ' censoring option.']) 
        return
    end
    FDvec = abs(FDvec-1);
end


function FDvec = remove_outliers_if(remove_outliers, FDvec, filename, ...
    stdev_temp_filename, wb_command, wait_time)
    % Remove outliers from external mask if user said to; otherwise pass
    if remove_outliers && ~strcmpi(filename, 'none')
        disp(['Removing outliers using .dtseries file ' filename])
        FDvec = remove_frames_from_FDvec(filename, FDvec, ...
            stdev_temp_filename, wb_command, wait_time);
        
    else  % exist('remove_outliers','var') == 1 && remove_outliers == 0;
        disp(['Motion censoring performed on FD alone. '...
              'Frames with outliers in BOLD std dev not removed']);
    end
end


function FDvec = add_mask_to_FDvec(FDvec, additional_mask)
    % Add to FDvec to account for additional_mask
    if strcmpi(additional_mask, 'none')
        disp(['No additional mask supplied. Using full time series. ' ...
              'Frames will be excluded only by FD (unless removal of ' ...
              'outliers is indicated).'])
    else
        % Load .txt file with 0s and 1s. 0s are frames to be discarded,
        % and 1s are frames to make your matrix.
        additionalvec = load(additional_mask); 
        FDvec_temp = FDvec & (additionalvec);
        FDvec = FDvec_temp;
    end
end


function FDvec = remove_frames_from_FDvec(cifti_file, FDvec, ...
        stdev, wb_command, wait_time)
    % additional frame removal based on Outliers command: isoutlier with
    % "median" method.
    cmd = [wb_command ' -cifti-stats ' char(cifti_file) ...
           ' -reduce STDEV > ' stdev];
    system(cmd);
    if wait_time > 0
        disp(['Waiting ' num2str(wait_time) ' seconds for writing of temp ' ...
              'file before reading. (not an error)'])
        pause(wait_time);
    end
    clear cmd
    STDEV_file=load(stdev); % load stdev of .nii file.
    FDvec_keep_idx = find(FDvec==1); % find frames kept from the FD mask
    
    % find outlier
    Outlier_file=isthisanoutlier(STDEV_file(FDvec_keep_idx), 'median'); 
    Outlier_idx=find(Outlier_file==1); % find outlier indices
    FDvec(FDvec_keep_idx(Outlier_idx)) = 0; % set outliers to zero in FDvec
    clear STDEV_file FDvec_keep_idx Outlier_file Outlier_idx
end


function finish_and_save_concs(A, FD_threshold, minutes_limit, ...
        motion_file, no_output, output_directory, series_path, p_or_d, ...
        smoothing_kernel, all_frames, at_thresh, at_thr_min_lim, without)
    % Given lists of paths to connectivity matrices, save them into their
    % respective .conc files
    if no_output
        [~, input_directory, ~] = get_folder_params(A{1}, filesep);
        output_directory =char(input_directory);
    end
    conc_names = {'all_frames', 'at_thresh', 'at_thr_min_lim', 'without'};
    concs = {all_frames, at_thresh, at_thr_min_lim, without};
    for i=1:4
        conc_path = get_conc_path(series_path, p_or_d, ...
            FD_threshold, conc_names(i), minutes_limit, ...
            smoothing_kernel, motion_file, output_directory);
        make_conn_conc(conc_path, concs(i));
    end
end


function [at_thresh, at_thresh_min_lim, without] = collect_conns( ...
    A_i, ext, FD_threshold, FDvec, min_lim, TR, ...
    at_thresh, at_thresh_min_lim, without)
    % Get the names of connectivity matrices which meet FD threshold or 
    % lack enough subject data
    fd = num2str(FD_threshold);
    if strcmpi(min_lim, 'none')
        at_thresh = [at_thresh; strcat({A_i}, '_all_frames_at_FD_', ...
                     {fd}, '.', {ext})];
    else
        % Sort pconn based on how its good minutes compare to the limit
        good_frames_idx = find(FDvec == 1);
        good_minutes = (length(good_frames_idx)*TR)/60;
        if good_minutes < 0.5
            without = [without; strcat({A_i}, '_lessthan30sec_', ...
                                       {fd}, '.', {ext})];
        elseif good_minutes < min_lim
            at_thresh = [at_thresh; strcat({A_i}, ...
                '_all_frames_at_FD_', {fd}, '.', {ext})];
        elseif good_minutes >= min_lim
            at_thresh_min_lim = [at_thresh_min_lim; strcat({A_i}, '_', ...
                {num2str(min_lim)}, '_minutes_of_data_at_FD_', {fd}, ...
                 '.', {ext})];
        else
            disp('Something is wrong generating .conc files.')
        end
    end
end


function save_out_file(out_dir, orig_motion, fd, min_file_end, FDvec)
    % Build output file name and save it
    fileID = fopen([char(out_dir) char(orig_motion) '_' num2str(fd) ...
                    min_file_end], 'w');
    fprintf(fileID, '%1.0f\n', FDvec);
    fclose(fileID);
end


function conc_path = get_conc_path(series_path, p_or_d, FD_threshold, ...
        minutes_condition, minutes_limit, smoothing, motion, output_dir)
    % Build and return the path to one .conc file using input parameters
    [~, series, ext] = fileparts(series_path);
    conc_prefix = [series ext '_' p_or_d 'conn_of_' p_or_d 'tseries'];
    fd = num2str(FD_threshold);
    minutes_condition = char(minutes_condition);
    
    if strcmpi(minutes_condition, 'all_frames')
        min = '_all_frames_at_FD_';
    elseif strcmpi(minutes_condition, 'at_thresh_min_lim')
        min = ['_' num2str(minutes_limit) '_minutes_of_data_at_FD_'];
    elseif strcmpi(minutes_condition, 'without')
        min = '_lessthan30sec_';
    else
        min = '_all_frames_at_FD_';
    end
    
    if strcmpi(motion, 'none')
        conc_extra = [min motion];
    else
        if strcmpi(smoothing, 'none')
            conc_extra = [min fd];
        else
            conc_extra = ['_SMOOTHED_' smoothing min fd];
        end
    end
    conc_path = [output_dir conc_prefix conc_extra '.conc'];
end


function make_conn_conc(concname, paths)
    % Save .conc file listing paths to all connectivity matrices
    paths = paths{1};
    if ~isempty(paths)
        fileID = fopen(char(concname), 'w');
        [nrows, ~] = size(paths);
        for row = 1:nrows
            fprintf(fileID, '%s\n', char(paths{row}));
        end
        fclose(fileID);
    end
end