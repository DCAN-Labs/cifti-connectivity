function cifti_conn_template_for_wrapper(wb_command, dt_or_ptseries_conc_file, series, motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal, left_surface_file, right_surface_file, bit8, remove_outliers, additional_mask, make_conn_conc, output_directory, keep_conn_matrices, template_file)
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
% make_dconn_conc - make a list of connectivity matrices that were created by this script.
% template_file = File path and name of template file to be created.
%  NOTE: this will delete dconns that have been previously made if the name
%  matches the expected output.

%% 3/12/18 -RH 
%Added an option to set minutes limit to none.  This will use all frames
%below a specific threshold.
%Added an arguement to keep .dconns/pconns after being made.

% Old hardcoded wb_command paths
%wb_command = '/Applications/workbench/bin_macosx64/wb_command'; 
%wb_command = 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/bin/wb_command';
%wb_command = '/home/exacloud/lustre1/fnl_lab/code/external/utilities/workbench-9253ac2/bin_rh_linux64/wb_command';

%% Change strings to numbers if not already numbers  (i.e.if using the shell command)
%Check to make sure that smoothing kernal is a number
if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
    smoothing_kernal= 'none';
elseif isnumeric(smoothing_kernal)==1
    disp('kernal is passed in as numeric')
else
    disp('kernal is  passed in as string. converting to numeric')
    smoothing_kernal = str2num(smoothing_kernal);
end
disp(smoothing_kernal);

if exist('output_directory','var') == 0
    output_directory = [];
end

%Check to make sure that FD_threshold is a number
if isnumeric(FD_threshold)==1
else
    FD_threshold = str2num(FD_threshold);
end

%Check to make sure that TR is a number
if isnumeric(TR)==1
else
    TR = str2num(TR);
end

%Check to make sure that  minutes limit is a number (unless you've set it
%to 'none')
if strcmp(minutes_limit,'none')==1 || strcmp(minutes_limit,'None')==1 || strcmp(minutes_limit,'NONE')==1
    disp('Warning minutes limit is set to none.  Subjects might have different number of data point in individual matrices');
    minutes_limit= 'none';
elseif isnumeric(minutes_limit)==1
    disp('minutes limit is passed in as numeric.')
else
    disp('minutes limit is passed in as string. converting to numeric')
    minutes_limit = str2num(minutes_limit);
end

%Check to make sure that bit8 is a number
if isnumeric(bit8)==1
else
    bit8 = str2num(bit8);
end

%Check to make sure that keep_conn_matrices is a number
if isnumeric(keep_conn_matrices)==1
else
    keep_conn_matrices = str2num(keep_conn_matrices);
end
%% Load concatenated paths (i.e. paths to ciftis)
%check to see if there one subject or a a list of subjects in conc file.

conc = strsplit(dt_or_ptseries_conc_file, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    A = importdata(dt_or_ptseries_conc_file);
else
    A = {dt_or_ptseries_conc_file};
end

for i = 1:length(A)
    if exist(A{i}) == 0
        NOTE = ['Series ' num2str(i) ' does not exist']
        return
    else
    end
end
disp('All series files exist continuing ...')

%% Combine dtseries for a given wave using d(p)tseries_cifti.conc
if series == 'ptseries'
    suffix = 'ptseries.nii';
    suffix2 = 'pconn.nii';
    suffix3 = 'pconn';
elseif series == 'dtseries'
    suffix = 'dtseries.nii';
    suffix2 = 'dconn.nii';
    suffix3 = 'dconn';
else
    'series needs to be "ptseries" or "dtseries"';
    return
end

%% Load Motion vector paths
if strcmp(motion_file,'none')==1 || strcmp(motion_file,'None')==1
    'No motion files, will use all frames to generate matrices'
else
    if strcmp('conc',conc) == 1
        B = importdata(motion_file);
    else
        B = {motion_file};
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
end

%% Check to see if surface files exist
%check to make sure there is a list of subjects in conc file.
conc = strsplit(dt_or_ptseries_conc_file, '.');
conc = char(conc(end));
if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
    Note = ['No smoothing, continuing...'];
    if strcmp('conc',conc) == 1
    else
        'Need to input a conc file with several subjects or dtseries to average'
        return
    end
else
    if strcmp('conc',conc) == 1
        C = importdata(left_surface_file); %check to make sure surface files exist
        D = importdata(right_surface_file); %check to make sure surface files exist
    else
        'Need to input a conc file with several subjects or dtseries to average'
        return
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
%% Generate Average Matrix
%make first dconn.nii then rename
if strcmp(minutes_limit,'none')==1 || strcmp(minutes_limit,'None')==1 || strcmp(minutes_limit,'NONE')==1
    if strcmp(motion_file,'none')==1 || strcmp(motion_file,'None')==1 || strcmp(motion_file,'NONE')==1
        if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
            [~, matrix_filename, matrix_ext] = fileparts(char({A{1}}));
            dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_all_frames_at_FD_', {motion_file}, '.', {suffix2}))];
            if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                cifti_conn_matrix_for_wrapper(wb_command, A{1},series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file,right_surface_file,bit8, remove_outliers, additional_mask, , output_directory)
                
                % rename out put to average
                %dconn_paths_all_frames = [char(strcat({A{1}}, '_all_frames_at_FD_', {motion_file}, '.', {suffix2}))];
                if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                    NOTE = ['Exiting...file ' dconn_paths_none_minutes_at_thresh_min_lim ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                    return
                else
                end
            else
                disp('conn.nii exist for this subject yet already.  Renaming to AVG for first subject.')
            end
            
            [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
            output_conc = [output_directory series_conc_name series_conc_ext];
            if keep_conn_matrices == 0
                cmd = ['mv ' dconn_paths_none_minutes_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            else
                cmd = ['cp ' dconn_paths_none_minutes_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            end
            system(cmd);
        else
            A_smoothed = [A{1}(1:length(A{1})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
            [~, matrix_filename, matrix_ext] = fileparts(char(A_smoothed));
            dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_all_frames_at_FD_', {motion_file}, '.', {suffix2}))];
            if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                cifti_conn_matrix_for_wrapper(wb_command, A{1},series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,C{1}, D{1},bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                %rename output to average
                if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                    NOTE = ['Exiting...file ' dconn_paths_none_minutes_at_thresh_min_lim ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                    return
                else
                end
            else
                disp('conn.nii exist for this subject yet already.  Renaming to AVG for first subject.')
            end
            [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
            output_conc = [output_directory series_conc_name series_conc_ext];
            if keep_conn_matrices == 0
                cmd = ['mv ' dconn_paths_none_minutes_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            else
                cmd = ['cp ' dconn_paths_none_minutes_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            end
            system(cmd);
        end
        
    else % motion correction
        if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
            [~, matrix_filename, matrix_ext] = fileparts(char({A{1}}));
            dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_all_frames_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
            if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                cifti_conn_matrix_for_wrapper(wb_command, A{1},series,B{1}, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file, right_surface_file,bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                % rename out put to average
                %dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat({A{1}}, '_all_frames_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                    NOTE = ['Exiting...file ' dconn_paths_none_minutes_at_thresh_min_lim ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                    return
                else
                end
            else
                disp('conn.nii exist for this subject yet already.  Renaming to AVG for first subject.')
            end
            [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
            output_conc = [output_directory series_conc_name series_conc_ext];
            if keep_conn_matrices == 0
            cmd = ['mv ' dconn_paths_none_minutes_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            else
            cmd = ['cp ' dconn_paths_none_minutes_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            end
            system(cmd); 
            
        else
            A_smoothed = [output_directory A{1}(1:length(A{1})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
            [~, matrix_filename, matrix_ext] = fileparts(char(A_smoothed));
            dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_all_frames_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
            if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                cifti_conn_matrix_for_wrapper(wb_command, A{1},series,B{1}, FD_threshold, TR, minutes_limit, smoothing_kernal,C{1}, D{1},bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                % rename out put to average
                %A_smoothed = [A{1}(1:length(A{1})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                %dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat(A_smoothed, '_all_frames_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                    NOTE = ['Exiting...file ' dconn_paths_none_minutes_at_thresh_min_lim ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                    return
                else
                end
                
            else
                disp('conn.nii exist for this subject yet already.  Renaming to AVG for first subject.')
            end
            [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
            output_conc = [output_directory series_conc_name series_conc_ext];
            if keep_conn_matrices == 0
                cmd = ['mv ' dconn_paths_none_minutes_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            else
                cmd = ['cp ' dconn_paths_none_minutes_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            end
            system(cmd);
            
        end
    end
else %minutes limit is a number
    
    if strcmp(motion_file,'none')==1 || strcmp(motion_file,'None')==1 || strcmp(motion_file,'NONE')==1
        if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
            [~, matrix_filename, matrix_ext] = fileparts(char({A{1}}));
            dconn_paths_all_frames = [char(strcat(output_directory, matrix_filename, matrix_ext, '_all_frames_at_FD_', {motion_file}, '.', {suffix2}))];
            if exist(dconn_paths_all_frames) == 0
                disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                cifti_conn_matrix_for_wrapper(wb_command, A{1},series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file,right_surface_file,bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                
                % rename out put to average
                %dconn_paths_all_frames = [char(strcat({A{1}}, '_all_frames_at_FD_', {motion_file}, '.', {suffix2}))];
                if exist(dconn_paths_all_frames) == 0
                    NOTE = ['Exiting...file ' dconn_paths_all_frames ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                    return
                else
                end
            else
                disp('conn.nii exist for this subject yet already.  Renaming to AVG for first subject.')
            end
            [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
            output_conc = [output_directory series_conc_name series_conc_ext];
            if keep_conn_matrices == 0
                cmd = ['mv ' dconn_paths_all_frames ' ' output_conc '_AVG.' suffix2];
            else
                cmd = ['cp ' dconn_paths_all_frames ' ' output_conc '_AVG.' suffix2];
            end
            system(cmd);
            
        else %smoothing
            
            A_smoothed = [A{1}(1:length(A{1})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
            [~, matrix_filename, matrix_ext] = fileparts(char(A_smoothed));
            dconn_paths_all_frames = [char(strcat(output_directory, matrix_filename, matrix_ext, '_all_frames_at_FD_', {motion_file}, '.', {suffix2}))];
            if exist(dconn_paths_all_frames) == 0
                disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                cifti_conn_matrix_for_wrapper(wb_command, A{1},series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,C{1}, D{1},bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                %rename output to average
                %A_smoothed = [A{1}(1:length(A{1})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                %dconn_paths_all_frames = [char(strcat(A_smoothed, '_all_frames_at_FD_', {motion_file}, '.', {suffix2}))];
                if exist(dconn_paths_all_frames) == 0
                    NOTE = ['Exiting...file ' dconn_paths_all_frames ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                    return
                else
                end
            else
                disp('conn.nii exist for this subject yet already.  Renaming to AVG for first subject.')
                
            end
            [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
            output_conc = [output_directory series_conc_name series_conc_ext];
            if keep_conn_matrices == 0
                cmd = ['mv ' dconn_paths_all_frames ' ' output_conc '_AVG.' suffix2];
            else
                cmd = ['cp ' dconn_paths_all_frames ' ' output_conc '_AVG.' suffix2];
            end
            system(cmd);
        end
    else %use motion correction
        if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
            [~, matrix_filename, matrix_ext] = fileparts(char({A{1}}));
            dconn_paths_all_frames_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_', {num2str(minutes_limit)}, '_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
            if exist(dconn_paths_all_frames_at_thresh_min_lim) == 0
                disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                
                cifti_conn_matrix_for_wrapper(wb_command, A{1},series,B{1}, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file, right_surface_file,bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                
                % rename out put to average
                %dconn_paths_all_frames_at_thresh_min_lim = [char(strcat({A{1}}, '_', {num2str(minutes_limit)}, '_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                if exist(dconn_paths_all_frames_at_thresh_min_lim) == 0
                    NOTE = ['Exiting...file ' dconn_paths_all_frames_at_thresh_min_lim ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                    return
                else
                end
            else
                disp('conn.nii exist for this subject yet already.  Renaming to AVG for first subject.')
                
            end
            [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
            output_conc = [output_directory series_conc_name series_conc_ext];
            if keep_conn_matrices == 0
                cmd = ['mv ' dconn_paths_all_frames_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            else
                cmd = ['cp ' dconn_paths_all_frames_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            end
            system(cmd);
            
        else
            A_smoothed = [A{1}(1:length(A{1})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
            [~, matrix_filename, matrix_ext] = fileparts(char(A_smoothed));
            dconn_paths_all_frames_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_', {num2str(minutes_limit)}, '_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
            if exist(dconn_paths_all_frames_at_thresh_min_lim) == 0
                disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                cifti_conn_matrix_for_wrapper(wb_command, A{1},series,B{1}, FD_threshold, TR, minutes_limit, smoothing_kernal,C{1}, D{1},bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                % rename out put to average
                %A_smoothed = [A{1}(1:length(A{1})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                %dconn_paths_all_frames_at_thresh_min_lim = [char(strcat(A_smoothed, '_', {num2str(minutes_limit)}, '_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                if exist(dconn_paths_all_frames_at_thresh_min_lim) == 0
                    NOTE = ['Exiting...file ' dconn_paths_all_frames_at_thresh_min_lim ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                    return
                else
                end
                disp('conn.nii exist for this subject yet already.  Renaming to AVG for first subject.')
            else
            end
            [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
            output_conc = [output_directory series_conc_name series_conc_ext];
            if keep_conn_matrices == 0
                cmd = ['mv ' dconn_paths_all_frames_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            else
                cmd = ['cp ' dconn_paths_all_frames_at_thresh_min_lim ' ' output_conc '_AVG.' suffix2];
            end
            system(cmd);
        end
    end
end %minutes limit

for i = 2:length(A) %build rest of template
    if strcmp(minutes_limit,'none')==1 || strcmp(minutes_limit,'None')==1 || strcmp(minutes_limit,'NONE')==1
        if strcmp(motion_file,'none')==1 || strcmp(motion_file,'None')==1 || strcmp(motion_file,'NONE')==1
            if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                [~, matrix_filename, matrix_ext] = fileparts(char({A{i}}));
                dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_none_minutes_of_data_at_FD_', {motion_file}, '.', {suffix2}))];
                if exist(dconn_paths_all_frames_at_thresh_min_lim) == 0
                    disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                    cifti_conn_matrix_for_wrapper(wb_command, A{i},series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file, right_surface_file,bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                    
                    % add image to average
                    %dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat({A{i}}, '_none_minutes_of_data_at_FD_', {motion_file}, '.', {suffix2}))];
                    if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                        NOTE = ['Exiting...file ' dconn_paths_all_frames ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                        return
                    else
                    end
                else
                    disp('conn.nii exist for this subject yet already.  Adding to the average dconn.')
                end
                [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
                output_conc = [output_directory series_conc_name series_conc_ext];
                cmd = [wb_command ' -cifti-math  " M1 + M2 " ' output_conc '_AVG.' suffix2 ' -var M1 ' output_conc '_AVG.' suffix2 ' -var M2 '  dconn_paths_none_minutes_at_thresh_min_lim];
                tic;
                system(cmd);
                toc;
                clear cmd
                if keep_conn_matrices == 0
                    delete(dconn_paths_none_minutes_at_thresh_min_lim)
                else
                end
                
            else %smooth data
                A_smoothed = [A{i}(1:length(A{i})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                [~, matrix_filename, matrix_ext] = fileparts(char(A_smoothed));
                dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_none_minutes_of_data_at_FD_', {motion_file}, '.', {suffix2}))];
                if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                    disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                    cifti_conn_matrix_for_wrapper(wb_command, A{i},series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,C{i}, D{i},bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                    
                    % add image to average
                    %A_smoothed = [A{i}(1:length(A{i})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                    %dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat(A_smoothed, '_none_minutes_of_data_at_FD_', {motion_file}, '.', {suffix2}))];
                    if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                        NOTE = ['Exiting...file ' dconn_paths_none_minutes_at_thresh_min_lim ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                        return
                    else
                    end
                else
                    disp('conn.nii exist for this subject yet already.  Adding to the average dconn.')
                end
                
                [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
                output_conc = [output_directory series_conc_name series_conc_ext];
                cmd = [wb_command ' -cifti-math  " M1 + M2 " ' output_conc '_AVG.' suffix2 ' -var M1 ' output_conc '_AVG.' suffix2 ' -var M2 '  dconn_paths_none_minutes_at_thresh_min_lim];
                tic;
                system(cmd);
                toc;
                clear cmd
                if keep_conn_matrices == 0
                    delete(dconn_paths_none_minutes_at_thresh_min_lim)
                else
                end
                
            end
            
        else %use motion correction
            if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                [~, matrix_filename, matrix_ext] = fileparts(char({A{i}}));
                dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_all_frames_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                    
                    disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                    cifti_conn_matrix_for_wrapper(wb_command, A{i},series,B{i}, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file, right_surface_file,bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                    % add image to average
                    %dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat({A{i}}, '_none_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                    
                    if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                        NOTE = ['Exiting...file ' dconn_paths_none_minutes_at_thresh_min_lim ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                        return
                    else
                    end
                else
                    disp('conn.nii exist for this subject yet already.  Adding to the average dconn.')
                end
                [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
                output_conc = [output_directory series_conc_name series_conc_ext];
                cmd = [wb_command ' -cifti-math  " M1 + M2 " ' output_conc '_AVG.' suffix2 ' -var M1 ' output_conc '_AVG.' suffix2 ' -var M2 '  dconn_paths_none_minutes_at_thresh_min_lim];
                tic;
                system(cmd);
                toc;
                clear cmd
                if keep_conn_matrices == 0
                    delete(dconn_paths_none_minutes_at_thresh_min_lim)
                else
                end
                
            else %smooth data
                A_smoothed = [A{i}(1:length(A{i})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                [~, matrix_filename, matrix_ext] = fileparts(char(A_smoothed));
                dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_all_frames_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                    disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                    cifti_conn_matrix_for_wrapper(wb_command, A{i},series,B{i}, FD_threshold, TR, minutes_limit, smoothing_kernal,C{i}, D{i},bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                    % add image to average
                    %A_smoothed = [A{i}(1:length(A{i})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                    %dconn_paths_none_minutes_at_thresh_min_lim = [char(strcat(A_smoothed, '_none_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                    if exist(dconn_paths_none_minutes_at_thresh_min_lim) == 0
                        NOTE = ['Exiting...file ' dconn_paths_none_minutes_at_thresh_min_lim ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                        return
                    else
                    end
                else
                    disp('conn.nii exist for this subject yet already.  Adding to the average dconn.')
                end
                [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
                output_conc = [output_directory series_conc_name series_conc_ext];
                cmd = [wb_command ' -cifti-math  " M1 + M2 " ' dt_or_ptseries_conc_file '_AVG.' suffix2 ' -var M1 ' dt_or_ptseries_conc_file '_AVG.' suffix2 ' -var M2 '  dconn_paths_none_minutes_at_thresh_min_lim];
                tic;
                system(cmd);
                toc;
                clear cmd
                if keep_conn_matrices == 0
                    delete(dconn_paths_none_minutes_at_thresh_min_lim)
                else
                end
                
            end
        end
    else %minutes limit is a number and not 'none'.
        if strcmp(motion_file,'none')==1 || strcmp(motion_file,'None')==1
            if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                [~, matrix_filename, matrix_ext] = fileparts(char({A{i}}));
                dconn_paths_all_frames = [char(strcat(output_directory, matrix_filename, matrix_ext, '_all_frames_at_FD_', {motion_file}, '.', {suffix2}))];
                if exist(dconn_paths_all_frames) == 0
                    disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                    
                    cifti_conn_matrix_for_wrapper(wb_command, A{i},series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file, right_surface_file,bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                    
                    % add image to average
                    %dconn_paths_all_frames = [char(strcat({A{i}}, '_all_frames_at_FD_', {motion_file}, '.', {suffix2}))];
                    if exist(dconn_paths_all_frames) == 0
                        NOTE = ['Exiting...file ' dconn_paths_all_frames ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                        return
                    else
                    end
                else
                    disp('conn.nii exist for this subject yet already.  Adding to the average dconn.')
                end
                
                [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
                output_conc = [output_directory series_conc_name series_conc_ext];
                cmd = [wb_command ' -cifti-math  " M1 + M2 " ' output_conc '_AVG.' suffix2 ' -var M1 ' output_conc '_AVG.' suffix2 ' -var M2 '  dconn_paths_all_frames];
                tic;
                system(cmd);
                toc;
                clear cmd
                if keep_conn_matrices == 0
                    delete(dconn_paths_all_frames)
                else
                end
                
                
            else
                A_smoothed = [A{i}(1:length(A{i})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                [~, matrix_filename, matrix_ext] = fileparts(char(A_smoothed));
                dconn_paths_all_frames = [char(strcat(output_directory, matrix_filename, matrix_ext, '_all_frames_at_FD_', {motion_file}, '.', {suffix2}))];
                if exist(dconn_paths_all_frames) == 0
                    disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                    
                    cifti_conn_matrix_for_wrapper(wb_command, A{i},series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,C{i}, D{i},bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                    
                    % add image to average
                    %A_smoothed = [A{i}(1:length(A{i})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                    %dconn_paths_all_frames = [char(strcat(A_smoothed, '_all_frames_at_FD_', {motion_file}, '.', {suffix2}))];
                    if exist(dconn_paths_all_frames) == 0
                        NOTE = ['Exiting...file ' dconn_paths_all_frames ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                        return
                    else
                    end
                else
                    disp('conn.nii exist for this subject yet already.  Adding to the average dconn.')   
                end
                [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
                output_conc = [output_directory series_conc_name series_conc_ext];
                cmd = [wb_command ' -cifti-math  " M1 + M2 " ' output_conc '_AVG.' suffix2 ' -var M1 ' output_conc '_AVG.' suffix2 ' -var M2 '  dconn_paths_all_frames];
                tic;
                system(cmd);
                toc;
                clear cmd
                if keep_conn_matrices == 0
                    delete(dconn_paths_all_frames)
                else
                end
                
            end
            
        else %use motion correction
            if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                [~, matrix_filename, matrix_ext] = fileparts(char({A{i}}));
                dconn_paths_all_frames_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_', {num2str(minutes_limit)}, '_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                if exist(dconn_paths_all_frames_at_thresh_min_lim) == 0
                    disp('conn.nii does not exist for this subject yet.  Running code to create matrix.')
                    
                    cifti_conn_matrix_for_wrapper(wb_command, A{i},series,B{i}, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file, right_surface_file,bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                    
                    % add image to average
                    %dconn_paths_all_frames_at_thresh_min_lim = [char(strcat({A{i}}, '_', {num2str(minutes_limit)}, '_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                    if exist(dconn_paths_all_frames_at_thresh_min_lim) == 0
                        NOTE = ['Exiting...file ' dconn_paths_all_frames_at_thresh_min_lim ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                        return
                    else
                    end
                else
                    disp('conn.nii exist for this subject yet already.  Adding to the average dconn.')
                end
                
                [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
                output_conc = [output_directory series_conc_name series_conc_ext];
                cmd = [wb_command ' -cifti-math  " M1 + M2 " ' output_conc '_AVG.' suffix2 ' -var M1 ' output_conc '_AVG.' suffix2 ' -var M2 '  dconn_paths_all_frames_at_thresh_min_lim];
                tic;
                system(cmd);
                toc;
                clear cmd
                if keep_conn_matrices == 0
                    delete(dconn_paths_all_frames_at_thresh_min_lim)
                else
                end
                
            else
                if exist(dconn_paths_all_frames_at_thresh_min_lim) == 0
                    A_smoothed = [A{i}(1:length(A{i})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                    [~, matrix_filename, matrix_ext] = fileparts(char(A_smoothed));
                    dconn_paths_all_frames_at_thresh_min_lim = [char(strcat(output_directory, matrix_filename, matrix_ext, '_', {num2str(minutes_limit)}, '_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                    cifti_conn_matrix_for_wrapper(wb_command, A{i},series,B{i}, FD_threshold, TR, minutes_limit, smoothing_kernal,C{i}, D{i},bit8, remove_outliers, additional_mask, make_conn_conc, output_directory)
                    % add image to average
                    %A_smoothed = [A{i}(1:length(A{i})-13) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                    %dconn_paths_all_frames_at_thresh_min_lim = [char(strcat(A_smoothed, '_', {num2str(minutes_limit)}, '_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '.', {suffix2}))];
                    if exist(dconn_paths_all_frames_at_thresh_min_lim) == 0
                        NOTE = ['Exiting...file ' dconn_paths_all_frames_at_thresh_min_lim ' does not exist, likely because subject does not meet criteron or there was an error making their connectivity matrix']
                        return
                    else
                    end
                else
                    disp('conn.nii exist for this subject yet already.  Adding to the average dconn.')
                end
                [~, series_conc_name, series_conc_ext] = fileparts(char(dt_or_ptseries_conc_file));
                output_conc = [output_directory series_conc_name series_conc_ext];
                cmd = [wb_command ' -cifti-math  " M1 + M2 " ' output_conc '_AVG.' suffix2 ' -var M1 ' output_conc '_AVG.' suffix2 ' -var M2 '  dconn_paths_all_frames_at_thresh_min_lim];
                tic;
                system(cmd);
                toc;
                clear cmd
                if keep_conn_matrices == 0
                    delete(dconn_paths_all_frames_at_thresh_min_lim)
                else
                end
                
            end
        end
    end
end

%% Generate Mean by dividing by N

N = num2str(length(A)); % number of observations

%mean = s1/N
cmd = [wb_command ' -cifti-math " s1/' N '" ' output_conc '_AVG.' suffix2 ' -var s1 ' output_conc '_AVG.' suffix2];
tic;
system(cmd);
toc;
clear cmd

%% Rename file to display chosen options
cmd = ['mv ' output_conc '_AVG.' suffix2 ' ' template_file];
system(cmd);
disp('Done making average template');
end
