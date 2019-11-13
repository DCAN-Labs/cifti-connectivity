function outfile = cifti_conn_matrix_for_wrapper(wb_command, dt_or_ptseries_conc_file,series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernel,left_surface_file, right_surface_file, bit8,remove_outliers, additional_mask, make_dconn_conc, output_directory, dtseries_conc)
% This script is a modularized version of cifti_conn_matrix_exaversion, except that it
% accepts and uses these parameters: remove_outliers, additional_mask, make_dconn_conc

%dt_or_ptseries_conc_file = dense timeseries or parcellated timeseries conc file (i.e. text file with paths to each file being examined
%series = 'dtseries' or 'ptseries'; specify if files are dense or parcellated
%motion_file = conc file that points to FNL motion mat files for each dt or ptseries (note: if no motion file to be used, type 'none')
%FD threshold = specify motion threshold (e.g. 0.2)
%TR = TR of your data
%minutes_limit = specify the number of minutes to be used to generate the correlation matrix
%smoothing_kernal = specify smoothing kernal (note: if no smoothing file to be used, type 'none')


% NEED TO ADD Write "Outlier frames" to file.

%Uncomment for debugging
% function  cifti_conn_matrix_exaversion
% dt_or_ptseries_conc_file='/mnt/max/shared/data/study/ABCD/ABCD_SIEMENS_SEFMNoT2/cub-sub-NDARINVLWRKNUN1/ses-baselineYear1Arm1/HCP_release_20161027_v1.1/cub-sub-NDARINVLWRKNUN1/MNINonLinear/Results/cub-sub-NDARINVLWRKNUN1_FNL_preproc_v2_Atlas.dtseries.nii';
% series='dtseries';
% motion_file= '/mnt/max/shared/data/study/ABCD/ABCD_SIEMENS_SEFMNoT2/cub-sub-NDARINVLWRKNUN1/ses-baselineYear1Arm1/HCP_release_20161027_v1.1/cub-sub-NDARINVLWRKNUN1/FNL_preproc_v2/analyses_v2/motion/power_2014_FD_only.mat'; %'none';
% FD_threshold=0.35;
% TR=0.8;
% minutes_limit=3;
% smoothing_kernal= '1.7'; %'none'
% left_surface_file='/mnt/max/shared/data/study/ABCD/ABCD_SIEMENS_SEFMNoT2/cub-sub-NDARINVLWRKNUN1/ses-baselineYear1Arm1/HCP_release_20161027_v1.1/cub-sub-NDARINVLWRKNUN1/MNINonLinear/fsaverage_LR32k/cub-sub-NDARINVLWRKNUN1.L.midthickness.32k_fs_LR.surf.gii';
% right_surface_file='/mnt/max/shared/data/study/ABCD/ABCD_SIEMENS_SEFMNoT2/cub-sub-NDARINVLWRKNUN1/ses-baselineYear1Arm1/HCP_release_20161027_v1.1/cub-sub-NDARINVLWRKNUN1/MNINonLinear/fsaverage_LR32k/cub-sub-NDARINVLWRKNUN1.R.midthickness.32k_fs_LR.surf.gii';
% bit8 = set to 1 if you want to make the outputs smaller in 8bit

orig_conc_folders = split(dt_or_ptseries_conc_file, filesep);
orig_conc_file = char(orig_conc_folders(end));

no_output = 0;
if exist('output_directory','var') == 0
    output_directory =[];
    no_output = 1;
    disp('All series files exist continuing ...')
else
    output_directory = char(output_directory);
end

if isempty(output_directory)
    no_output = 1;
end
%convert number variables to cnumbers (only within compiled code);
%Check to make sure that FD_threshold is a number
if isnumeric(FD_threshold)
else
    FD_threshold = str2num(FD_threshold);
end

%Check to make sure that TR is a number
if isnumeric(TR)
else
    TR = str2num(TR);
end

%Check to make sure that  minutes limit is a number (unless you've set it
%to 'none')
if strcmpi(minutes_limit,'none')
    minutes_limit= 'none';
elseif isnumeric(minutes_limit)
    disp('minutes limit is passed in as numeric.')
else
    disp('minutes limit is passed in as string. converting to numeric')
    minutes_limit = str2num(minutes_limit);
end

%Check to make sure that bit8 is a number
if isnumeric(bit8)
else
    bit8 = str2num(bit8);
end

if exist('remove_outliers','var')
    if isnumeric(remove_outliers)
    else
        remove_outliers = str2num(remove_outliers);
    end
end

if strcmpi(smoothing_kernel,'none')
    smoothing_kernel= 'none';
elseif isnumeric(smoothing_kernel)
else
    smoothing_kernel = str2num(smoothing_kernel);
end

if (strcmpi(smoothing_kernel,'none')~=1) && (strcmpi(series,'ptseries'))
    disp('Check your settings. Smoothing not allowed on ptseries.');
    return
else
end

%% Load concatenated paths (i.e. paths to ciftis)
[A, conc] = get_paths_from_conc(dt_or_ptseries_conc_file);
[dtseries_E, ~] = get_paths_from_conc(dtseries_conc);

%% set file extensions
if strcmpi(series, 'ptseries')
    suffix = 'ptseries.nii';
    suffix2 = 'pconn.nii';
    output_precision = ' '; %set default precision. Connectivity file sizes are NOT reduced for ptseries
elseif strcmpi(series, 'dtseries')
    suffix = 'dtseries.nii';
    suffix2 = 'dconn.nii';
    if exist('bit8','var') == 1
        if isnumeric(bit8)==1
        else
            bit8 = str2num(bit8);
        end
        if bit8 == 1
            output_precision = ' -cifti-output-datatype INT8 -cifti-output-range -1.0 1.0'; %set default precision. Connectivity file sizes are reduced for dtseries
        else
            output_precision = ' ';
        end
    else
        output_precision = ' ';
    end
else
    'series needs to be "ptseries" or "dtseries"';
    return
end

disp(smoothing_kernel)
%% prealocate memory if you want to do smoothing
if strcmpi(smoothing_kernel, 'none')
else
    A_smoothed = num2cell(zeros(size(A)));
    % moved to line 178 (use motion censoring) Anders Perrone - 20171107
    %stdev_temp_filename=[char(dt_or_ptseries_conc_file) '_temp.txt']; %rename temp file so that it can't be overwritten when the script is run in parallel.
    if strcmpi('conc',conc) == 1
        C = importdata(left_surface_file); %check to make sure surface files exist
        D = importdata(right_surface_file); %check to make sure surface files exist
    else
        C = {left_surface_file};
        D = {right_surface_file};
    end
    for i = 1:length(C)
        if exist(C{i}) == 0
            NOTE = ['Subject left surface ' num2str(i) ' does not exist']
            return
        else
        end
    end
    disp('All left surface files for smoothing exist continuing ...')
    for i = 1:length(D)
        if exist(D{i}) == 0
            NOTE = ['Subject right surface ' num2str(i) ' does not exist']
            return
        else
        end
    end
    disp('All right surface files for smoothing exist continuing ...')
end

%% Generate Motion Vectors and correlation matrix
hashes = cell(1, length(A)); % Preallocate list mapping input ix to hash
if strcmpi(motion_file,'none') % run this if no motion censoring
    disp('No motion files, will use all frames to generate matrices')
    for i = 1:length(A)
        [orig_cifti_filename, input_directory, ...
            folders] = get_folder_params(A{i}, filesep);
        if no_output
            output_directory =char(input_directory);
        end
       
        % Run minutes limit calculation
        hash = get_hash(folders);
        hashes{i} = hash;
        if strcmpi(smoothing_kernel, 'none')
            left = 'none';
            right = 'none';            
        else
            left = C{i};
            right = D{i};
        end
        min_file_end = ['_cifti_censor_FD_vector_All_Good_Frames_' ...
            hash '.txt'];
        [A_smoothed{i}, outfile] = minutes_limit_calculation(...
            FD_threshold, input_directory, left, '_all_frames_at_', ...
            min_file_end, motion_file, orig_cifti_filename, 'none', ...
            hash, output_directory, output_precision, right, ...
            smoothing_kernel, suffix, suffix2, wb_command);
    end
else %use motion censoring
    stdev_temp_filename=[output_directory orig_conc_file '_temp.txt']; %rename temp file so that it can't be overwritten when the script is run in parallel.
    if strcmp('conc',conc) == 1
        B = importdata(motion_file);
    else
        B = {motion_file};
    end
    
    %Make sure length of motion files matches length of timeseries
    if length(A)==length(B)
    else
        disp('length of motion conc file does not match length of series con file')
    end
    
    % Check that all motion files in conc file exist
    for i = 1:length(B)
        if exist(B{i}) == 0
            NOTE = ['motion file ' num2str(i) ' does not exist']
            return
        else
        end
    end
    disp('All motion files exist continuing ...')
    
    %% Generate motion vector for subject and generate correlation matrix
    for i = 1:length(B)
        motion_exten = strsplit(B{i}, '.');
        motion_exten = char(motion_exten(end));
        
        % Get input folder parameters and make hash
        [orig_cifti_filename, input_directory, ...
            folders] = get_folder_params(A{i}, filesep);
        hash = get_hash(folders);
        hashes{i} = hash;
        if no_output
            output_directory = char(input_directory);
        end
        
        if strcmp('mat',motion_exten) == 1
            other_motion_mask = 0; 
        else
            other_motion_mask = 1;  
        end
        if other_motion_mask %use an external mask rather than calculating the mask here (e.g. .txt file.
            FDvec=B{i};
            if exist('remove_outliers','var') == 0 %remove outliers from external mask
                disp('Removal outliers not specified.  It will be performed by default.')
                FDvec = remove_frames_from_FDvec(dtseries_E{i}, ...
                    FDvec, stdev_temp_filename, wb_command, 30);  
                
            elseif remove_outliers == 1 %remove outliers from external mask
                disp('Removing outliers using .dtseries file...')
                FDvec = remove_frames_from_FDvec(dtseries_E{i}, ...
                    FDvec, stdev_temp_filename, wb_command, 30);    
            else % exist('remove_outliers','var') == 1 && remove_outliers == 0;
                disp('Motion censoring performed on FD alone. Frames with outliers in BOLD std dev not removed');
            end
        else  %Use power 2014 motion
            FDvec = get_FDvec(B{i}, FD_threshold);
            [orig_motion_filename, inputB_directory, ...
                ~] = get_folder_params(B{i}, filesep);
            if no_output
                output_directory = char(inputB_directory);
            end
            if strcmpi(additional_mask,'none')==1
                disp('No additional mask supplied. Using full time series. Frames will be excluded only by FD (unless removal of outliers is indicated).')
            else
                additionalvec = load(additional_mask); %Load .mat file with 0s and 1s.  where 0 are frames to be discarded and 1s are frames to be used to make your matrix.
                FDvec_temp = FDvec & (additionalvec);
                FDvec = FDvec_temp;
            end
            if exist('remove_outliers','var') == 0 || remove_outliers == 1
                disp('Removal outliers not specified.  It will be performed by default.')
                FDvec = remove_frames_from_FDvec([input_directory ...
                    orig_cifti_filename], FDvec, stdev_temp_filename, ...
                    wb_command, 30);
            else % exist('remove_outliers','var') == 1 && remove_outliers == 0;
                disp('Motion censoring performed on FD alone. Frames with outliers in BOLD std dev not removed');
            end
        end
        
        % Assign variables for minutes limit calculation
        if strcmpi(minutes_limit,'none')            
            min_cmd_part = '_all_frames_at_';
            min_file_end = ['_cifti_censor_FD_vector_All_Good_Frames_' hash '.txt'];
            
            fileID = fopen([char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) min_file_end],'w');
            fprintf(fileID,'%1.0f\n',FDvec);
            fclose(fileID);
        else 
            good_frames_idx = find(FDvec == 1);
            good_minutes = (length(good_frames_idx)*TR)/60; % number of good minutes in your data
            
            if good_minutes < 0.5 % if there is less than 30 seconds, don't generate the correlation matrix
                subject_has_too_few_frames = ['Subject ' num2str(i) ' has less than 30 seconds of good data']
                
            elseif minutes_limit > good_minutes % if there is not enough data for your subject, just generate the matrix with all available frames
                min_cmd_part = '_all_frames_at_';
                min_file_end = ['_cifti_censor_FD_vector_All_Good_Frames_' hash '.txt'];
                
                fileID = fopen([char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) min_file_end],'w');
                fprintf(fileID,'%1.0f\n',FDvec);
                fclose(fileID);
                               
            elseif minutes_limit <= good_minutes % if there is enough data, match the amount of data used for your subject matrix to the minutes_limit
                good_frames_needed = round(minutes_limit*60/TR); %number of good frames to randomly pull
                rand_good_frames = sort(randperm(length(good_frames_idx),good_frames_needed));
                FDvec_cut = zeros(length(FDvec),1);
                ones_idx = good_frames_idx(rand_good_frames);
                FDvec_cut(ones_idx) = 1; % the new vector that should match good frames with the minutes limit
 
                min_cmd_part = ['_' num2str(minutes_limit) '_minutes_of_data_at_'];
                min_file_end = ['_cifti_censor_FD_vector' min_cmd_part num2str(FD_threshold) '_threshold' '_' hash '.txt'];
                fileID = fopen([char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) min_file_end],'w');
                fprintf(fileID,'%1.0f\n', FDvec_cut);
                fclose(fileID);
            else
                disp('something is wrong with regard to number of good minutes calc')
            end
        end
        
        % Run minutes limit calculation
        if strcmpi(minutes_limit,'none') || good_minutes >= 0.5
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
                orig_motion_filename, hash, output_directory, ...
                output_precision, right, smoothing_kernel, suffix, ...
                suffix2, wb_command);
            
            if strcmpi(smoothing_kernel, 'none')
                clear A_smoothed
            end
        end        
    end
end

%% GENERATE FINAL CONC text file here to use in other code (e.g. cifti_conn_pairwise_corr.m)
% Greg Conan modularized everything from here down on 2019-11-11
if make_dconn_conc
    disp('Done making (p or d)conn.nii for subjects.  Making output .concs')
            
    %% Get all paths to put into different conc files
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
    [dconn_paths_all_frames, dconn_paths_all_frames_at_thresh, ...
     dconn_paths_all_frames_at_thresh_min_lim, ...
     subjectswithoutenoughdata] = sort_paths(A, motion_file_obj, ...
        FD_threshold, hashes, minutes_limit, motion_file, no_output, ...
        output_directory, smoothing, stdev, suffix2, TR, wb_command, ...
        additional_mask);
    
    %% Save all of those paths into .conc files
    finish_and_save_concs(A, FD_threshold, minutes_limit, ...
        motion_file, no_output, output_directory, ...
        dt_or_ptseries_conc_file, series(1), smoothing_kernel, ...
        dconn_paths_all_frames, dconn_paths_all_frames_at_thresh, ...
        dconn_paths_all_frames_at_thresh_min_lim, ...
        subjectswithoutenoughdata)
end
disp('Done making output concs.')
end


function [paths, conc] = get_paths_from_conc(conc_file)
    %check to see if there 1 subject or a list of subjects in conc file.
    conc = strsplit(conc_file, '.');
    conc = char(conc(end));
    if strcmp('conc',conc) == 1
        paths = importdata(conc_file);
    else
        paths = {conc_file};
    end

    % Validate that all paths in .conc file point to real files
    for i = 1:length(paths)
        if exist(paths{i}) == 0
            NOTE = ['Subject Series ' num2str(i) ' does not exist']
            return
        else
        end
    end
    disp('All series files exist continuing ...')
end


function [cifti, in_dir, folders] = get_folder_params(A_i, sep)
    % Return the path to the original cifti file, the path to the input
    % directory, and an object holding all parts of A_i's file path
    folders = split(A_i, sep);
    cifti = char(folders(end));
    folder_input = join(folders(1:end-1), sep);
    in_dir = [char(folder_input) sep];
end


function hash = get_hash(folders)
    % Return a string composed of the first character in every element of
    % the 'folders' array except the last element, which has its first 5 
    % characters at the end. This consistently but uniquely identifies
    % each file based on its original location and ID.
    folders_len = length(folders)-2;
    hash = blanks(folders_len);
    for d=1:folders_len
        hash(1, d) = folders{d+1}(1);
    end
    hash = [hash folders{end}(1:5)];
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


function not_exist = matrix_not_already_exist(to_make)
    % Return 1 if to_make does not exist, and 0 if it does
    not_exist = (exist(to_make) == 0);
end


function [A_smoothed_i, outfile] = minutes_limit_calculation(...
        FD_threshold, input_directory, left, min_cmd_part, ...
        min_file_end, motion_file, orig_cifti_filename, ...
        orig_motion_filename, hash, output_dir, output_precision, ...
        right, smoothing_kernel, suffix, suffix2, wb_command)
    % Calculate minutes_limit
    if strcmpi(motion_file, 'none')
       min_cmd_motion = [min_cmd_part 'FD_' motion_file]; 
       min_file_end = [hash '.txt'];
       weights = '';
    else
       min_cmd_motion = [min_cmd_part 'FD_' num2str(FD_threshold)];
    end
    if strcmpi(smoothing_kernel,'none')
        outfile = [char(output_dir) char(orig_cifti_filename) min_cmd_motion '_' hash '.' suffix2];
        A_smoothed_i = 'none';
    else % Smooth
        A_smoothed_i = [char(output_dir) char(orig_cifti_filename(1:length(char(orig_cifti_filename))-13)) '_SMOOTHED_' num2str(smoothing_kernel) '_' hash '.' suffix];
        outfile = [A_smoothed_i min_cmd_motion '_' hash '.' suffix2];
        if exist(A_smoothed_i) == 0
            Column = 'COLUMN';
            cmd = [wb_command ' -cifti-smoothing ' char(input_directory) char(orig_cifti_filename) ' ' num2str(smoothing_kernel) ' ' num2str(smoothing_kernel) ' ' Column ' ' char(output_dir) char(orig_cifti_filename(1:length(char(orig_cifti_filename))-13)) '_SMOOTHED_' num2str(smoothing_kernel)  '_' hash '.' suffix ' -left-surface ' left '  -right-surface ' right];
            execute_display_and_clear(cmd, cmd);
        else %smoothed series already exists
            disp('Smoothed series already created for this subject')
        end
    end
    if matrix_not_already_exist(outfile) % check to see if the file already exists
        if strcmpi(smoothing_kernel,'none')
            if ~exist('weights', 'var')
                weights = [' -weights ' char(output_dir) char(orig_motion_filename) '_' num2str(FD_threshold) min_file_end];
            end
            cmd = [wb_command ' ' output_precision ' -cifti-correlation ' char(input_directory) char(orig_cifti_filename) ' ' char(output_dir) char(orig_cifti_filename) min_cmd_motion '_' hash '.' suffix2 weights ' -fisher-z'];
        else
            if ~exist('weights', 'var')
                weights = [' -weights ' char(output_dir) char(orig_motion_filename) '_' num2str(FD_threshold) min_file_end];
            end
            cmd = [wb_command ' ' output_precision ' -cifti-correlation ' A_smoothed_i ' ' A_smoothed_i min_cmd_motion '_' hash '.' suffix2 weights ' -fisher-z'];
        end
        execute_display_and_clear(cmd, 'Running wb_command cifti-correlation.  This may take a few minutes.');
    else
        disp([outfile ' already exists'])
    end
end


function [all_frames, at_thresh, at_thr_min_lim, without] = sort_paths(...
        A, B, FD_threshold, hashes, minutes_limit, motion_file, ...
        no_output, output_directory, smoothed, stdev_temp_filename, ...
        suffix2, TR, wb_command, additional_mask)
    % Build and return 4 lists of paths to connectivity matrix files, such
    % that each list contains matrices with a certain amount of data 
    % compared to the minutes_limit and FD_threshold:
    all_frames = [];     % All frames regardless of minutes_limit & FD
    at_thresh = [];      % All frames at FD_threshold
    at_thr_min_lim = []; % All frames at FD_threshold and minutes_limit
    without = [];        % Subjects with < 30 seconds of data
    for i = 1:length(A)
        hash = hashes{i};
        [orig_cifti_filename, input_directory, ~] = get_folder_params(A{i}, filesep);
        if no_output
            output_directory = char(input_directory);
        end
        if strcmpi(smoothed, 'none')==1
            A_dir = strcat(char(output_directory),char(orig_cifti_filename));
        else
            A_dir = smoothed{i};
        end
        if strcmpi(motion_file,'none')==1
            all_frames = [all_frames; strcat({A_dir}, '_all_frames_at_FD_', {num2str(FD_threshold)}, '_', hash, '.', {suffix2})];
        else
            FDvec = get_FDvec(B{i}, FD_threshold);
            FDvec = additional_frame_removal(additional_mask, FDvec, ...
                        input_directory, orig_cifti_filename, ...
                        stdev_temp_filename, wb_command);
            
            [at_thresh, at_thr_min_lim, without] = collect_conns( ...
                A_dir, suffix2, FD_threshold, FDvec, hash, ...
                minutes_limit, TR, at_thresh, at_thr_min_lim, without);
        end
    end
end


function FDvec = additional_frame_removal(additional_mask, FDvec, ...
    input_directory, orig_cifti_filename, stdev_temp_filename, wb_command)
    % additional frame removal depending on additional_mask
    if strcmpi(additional_mask, 'none')==1
        FDvec = remove_frames_from_FDvec([input_directory ...
                orig_cifti_filename], FDvec, ...
            stdev_temp_filename, wb_command, 30);
    else
        FDvec = add_mask_to_FDvec(FDvec, additional_mask);
        if exist('remove_outliers','var') == 0 || remove_outliers == 1
            disp('Removal outliers not specified.  It will be performed by default.')
            FDvec = remove_frames_from_FDvec([input_directory ...
                orig_cifti_filename], FDvec, ...
                stdev_temp_filename, wb_command, 30);
            
        elseif exist('remove_outliers','var') == 1 && remove_outliers == 0
            disp('Motion censoring performed on FD alone. Frames with outliers in BOLD std dev not removed');
        end
    end
end


function FDvec = get_FDvec(B_i, fd)
    % Get FDvec from the file at B_i
    load(B_i)
    allFD = zeros(1,length(motion_data)); % motion_data comes from load(B_i)
    for j = 1:length(motion_data)
        allFD(j) = motion_data{j}.FD_threshold;
    end
    FDidx = find(allFD == fd);
    FDvec = motion_data{FDidx}.frame_removal;
    FDvec = abs(FDvec-1);
end


function FDvec = remove_frames_from_FDvec(cifti_file, FDvec, ...
        stdev, wb_command, wait)
    % additional frame removal based on Outliers command: isoutlier with "median" method.
    cmd = [wb_command ' -cifti-stats ' char(cifti_file) ' -reduce STDEV > ' stdev];
    system(cmd);
    if wait > 0
        disp(['waiting ' num2str(wait) ' seconds for writing of temp file before reading. (not an error)'])
        pause(wait);
    end
    clear cmd
    STDEV_file=load(stdev); % load stdev of .nii file.
    FDvec_keep_idx = find(FDvec==1); %find the kept frames from the FD mask
    Outlier_file=isthisanoutlier(STDEV_file(FDvec_keep_idx),'median'); %find outlier
    Outlier_idx=find(Outlier_file==1); %find outlier indices
    FDvec(FDvec_keep_idx(Outlier_idx))=0; %set outliers to zero within FDvec
    clear STDEV_file FDvec_keep_idx Outlier_file Outlier_idx
end


function [at_thresh, at_thresh_min_lim, without] = collect_conns( ...
    A_i, ext, FD_threshold, FDvec, hash, min_lim, TR, ...
    at_thresh, at_thresh_min_lim, without)
    % Get the names of connectivity matrices which meet FD threshold or 
    % lack enough subject data
    fd = num2str(FD_threshold);
    if strcmpi(min_lim,'none')==1
        at_thresh = [at_thresh; strcat({A_i}, '_all_frames_at_FD_', {fd}, '_', hash, '.', {ext})];
    else
        good_frames_idx = find(FDvec == 1);
        good_minutes = (length(good_frames_idx)*TR)/60;
        if good_minutes < 0.5
            without = [without; strcat({A_i}, '_lessthan30sec_', {fd}, '_', hash, '.', {ext})];
        elseif min_lim > good_minutes
            at_thresh = [at_thresh; strcat({A_i}, '_all_frames_at_FD_', {fd}, '_', hash, '.', {ext})];
        elseif min_lim <= good_minutes
            at_thresh_min_lim = [at_thresh_min_lim; strcat({A_i}, '_', {num2str(min_lim)}, '_minutes_of_data_at_FD_', {fd}, '_', hash, '.', {ext})];
        else
            disp('something is wrong generating conc files')
        end
    end
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
    conc_names = {'all_frames', 'at_thresh', 'at_thr_min_lim', ...
                  'without'};
    concs = {all_frames, at_thresh, at_thr_min_lim, without};
    for i=1:4
        conc_path = get_conc_path(series_path, p_or_d, ...
            FD_threshold, conc_names(i), minutes_limit, ...
            smoothing_kernel, motion_file, output_directory);
        make_conn_conc(conc_path, concs(i));
    end
end


function FDvec = add_mask_to_FDvec(FDvec, additional_mask)
    % Add to FDvec to account for additional_mask
    if strcmpi(additional_mask,'none')==1
        disp('No additional mask supplied. Using full time series. Frames will be excluded only by FD (unless removal of outliers is indicated).')
    else
        additionalvec = load(additional_mask); %Load .mat file with 0s and 1s.  where 0 are frames to be discarded and 1s are frames to be used to make your matrix.
        FDvec_temp = FDvec & (additionalvec);
        FDvec = FDvec_temp;
    end
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
    if strcmpi(motion, 'none')==1
        conc_extra = [min motion];
    else
        if strcmpi(smoothing, 'none')==1
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
        concname = char(concname);
        fileID = fopen(concname,'w');
        [nrows, ~] = size(paths);
        for row = 1:nrows
            fprintf(fileID, '%s\n', char(paths{row}));
        end
        fclose(fileID);
    end
end