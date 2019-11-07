function outfile = cifti_conn_matrix_for_wrapper(wb_command, dt_or_ptseries_conc_file,series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file, right_surface_file, bit8,remove_outliers, additional_mask, make_dconn_conc, output_directory, dtseries_conc)

% This should be the same as cifti_conn_matrix_exaversion, except that it
% accepts and uses these parameters: remove_outliers, additional_mask, make_dconn_conc

%dt_or_ptseries_conc_file = dense timeseries or parcellated timeseries conc file (i.e. text file with paths to each file being examined
%series = 'dtseries' or 'ptseries'; specify if files are dense or parcellated
%motion_file = conc file that points to FNL motion mat files for each dt or ptseries (note: if no motion file to be used, type 'none')
%FD threshold = specify motion threshold (e.g. 0.2)
%TR = TR of your data
%minutes_limit = specify the number of minutes to be used to generate the correlation matrix
%smoothing_kernal = specify smoothing kernal (note: if no smoothing file to be used, type 'none')
% example:
% cifti_conn_matrix('/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/ADHD_DVARS_group_dtseries.conc','dtseries','/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/ADHD_DVARS_group_motion.conc',0.2,2.5,5,2.55,'/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/ADHD_DVARS_group_right_midthickness_surfaces.conc','/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/ADHD_DVARS_group_right_midthickness_surfaces.conc');
% example without smoothing:
% cifti_conn_matrix('/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/ADHD_DVARS_group_dtseries.conc','dtseries','/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/ADHD_DVARS_group_motion.conc',0.2,2.5,5,'none','none','none');


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

orig_conc_folders = split(dt_or_ptseries_conc_file,'/');
orig_conc_file = char(orig_conc_folders(end));
orig_conc_folder_input = join(orig_conc_folders(1:end-1),'/');
orig_conc_folder_input = [char(orig_conc_folder_input) '/'];

no_output = 0;
if exist('output_directory','var') == 0
    output_directory =[];
    no_output = 1;    disp('All series files exist continuing ...')

else
    output_directory = char(output_directory);
end
if isempty(output_directory)
    no_output = 1;
end
%convert number variables to cnumbers (only within compiled code);
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

if exist('remove_outliers','var') == 1
    if isnumeric(remove_outliers)==1
    else
        remove_outliers = str2num(remove_outliers);
    end
end


% if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
% smoothing_kernal= 'none';
% else
% end

if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
    smoothing_kernal= 'none';
elseif isnumeric(smoothing_kernal)==1
else
    smoothing_kernal = str2num(smoothing_kernal);
end


if (strcmp(smoothing_kernal,'none')~=1) && (strcmp(series,'ptseries')==1)
    disp('Check your settings. Smoothing not allowed on ptseries.');
    return
else
end

if exist('make_dconn_conc','var') == 1
    if isnumeric(make_dconn_conc)==1
    else
        make_dconn_conc = str2num(make_dconn_conc);
    end
end

%%
%HARDCODE WARNING
%wb_command =
%'/home/exacloud/lustre1/fnl_lab/code/external/utilities/workbench/bin_rh_linux64/wb_command'; %old version
%wb_command = '/home/exacloud/tempwork/fnl_lab/code/external/utilities/workbench-1.3.2/bin_rh_linux64/wb_command';%exacloud path
%wb_command = '/usr/local/bin/wb_command';
%wb_command = 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/bin/wb_command'; % workbench command path
%wb_command = '/Applications/workbench/bin_macosx64/wb_command'; % workbench command path

%% Load concatenated paths (i.e. paths to ciftis)
[A, conc] = get_paths_from_conc(dt_or_ptseries_conc_file);
[dtseries_E, ~] = get_paths_from_conc(dtseries_conc);

%% set file extensions
disp(series)
if strcmp(series, 'ptseries')
    suffix = 'ptseries.nii';
    suffix2 = 'pconn.nii';
    suffix3 = 'pconn';
    output_precision = ' '; %set default precision. Connectivity file sizes are NOT reduced for ptseries
elseif strcmp(series, 'dtseries')
    suffix = 'dtseries.nii';
    suffix2 = 'dconn.nii';
    suffix3 = 'dconn';
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

disp(smoothing_kernal)
%% prealocate memory if you want to do smoothing
if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
else
    A_smoothed = num2cell(zeros(size(A)));
    % moved to line 178 (use motion censoring) Anders Perrone - 20171107
    %stdev_temp_filename=[char(dt_or_ptseries_conc_file) '_temp.txt']; %rename temp file so that it can't be overwritten when the script is run in parallel.
    if strcmp('conc',conc) == 1
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

if strcmp(motion_file,'none')==1 || strcmp(motion_file,'None')==1 % run this if no motion censoring
    'No motion files, will use all frames to generate matrices'
    for i = 1:length(A)
        folders = split(A{i},'/');
        orig_cifti_filename = char(folders(end));
        folder_input = join(folders(1:end-1),'/');
        input_directory = [char(folder_input) '/'];
        if no_output
            output_directory =char(input_directory);
        end
        if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
            outfile = [char(output_directory) char(orig_cifti_filename) '_all_frames_at_FD_' motion_file '_' char(folders(end-1)) '.' suffix2];
        else
            A_smoothed{i} = [char(output_directory) char(orig_cifti_filename(1:length(char(orig_cifti_filename))-13)) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
            outfile = [A_smoothed{i} '_all_frames_at_FD_' motion_file '_' char(folders(end-1)) '.' suffix2];
            if exist(A_smoothed{i}) == 0 %check to see if smoothing file does not exist
                clear smoothedfile_name
                Column = 'COLUMN';
                cmd = [wb_command ' -cifti-smoothing ' char(input_directory) char(orig_cifti_filename) ' ' num2str(smoothing_kernal) ' ' num2str(smoothing_kernal) ' ' Column ' ' char(output_directory) char(orig_cifti_filename(1:length(char(orig_cifti_filename))-13)) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix  ' -left-surface ' C{i} '  -right-surface ' D{i}];
                disp('207')
                disp(cmd)
                system(cmd);
                clear cmd
            else %smoothed series already exists
                disp('Smoothed series already created for this subject')
            end
        end
        if matrix_not_already_exist(outfile) % make sure the matrix doesn't already exist
            %%%%
            if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                
                cmd = [wb_command output_precision ' -cifti-correlation ' char(input_directory) char(orig_cifti_filename) ' ' char(output_directory) char(orig_cifti_filename) '_all_frames_at_FD_' motion_file  '_' char(folders(end-1)) '.' suffix2  ' -fisher-z'];
                %cmd = [wb_command ' -cifti-correlation ' A{i} ' ' A{i} '_all_frames_at_FD_' motion_file '.' suffix2  ' -fisher-z'];
                tic;
                disp('Running wb_command cifti-correlation.  This may take a few minutes.')
                system(cmd);
                toc;
                clear cmd
            else
                cmd = [wb_command output_precision ' -cifti-correlation ' A_smoothed{i} ' ' A_smoothed{i} '_all_frames_at_FD_' motion_file  '_' char(folders(end-1)) '.' suffix2  ' -fisher-z'];
                %cmd = [wb_command ' -cifti-correlation ' A{i} ' ' A{i} '_all_frames_at_FD_' motion_file '.' suffix2  ' -fisher-z'];
                tic;
                disp('Running wb_command cifti-correlation.  This may take a few minutes.')
                system(cmd);
                toc;
                clear cmd
            end
        else
            if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                subject_conn_exists = [char(input_directory) char(orig_cifti_filename) '_all_frames_at_FD_' motion_file  '_' char(folders(end-1)) '.' suffix2 ' already exists']
            else
                subject_conn_exists = [A_smoothed{i} '_all_frames_at_FD_' motion_file  '_' char(folders(end-1)) '.' suffix2 ' already exists']
            end
        end
    end
else %use motion cenosoring
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
            
            if strcmp('mat',motion_exten) == 1
            other_motion_mask =0; 
            else
             other_motion_mask =1;  
            end
        if other_motion_mask ==1 %use an external mask rather than calculating the mask here (e.g. .txt file.
            FDvec=B{i};
            if exist('remove_outliers','var') == 0 || remove_outliers == 1 %remove outliers from external mask
                disp('Removal outliers not specified.  It will be performed by default.')
                %% additional frame removal based on Outliers command: isoutlier with "median" method.
                cmd = [wb_command ' -cifti-stats ' dtseries_E{i} ' -reduce STDEV > ' stdev_temp_filename];      % Changed A{i} to new dtseries .conc file parameter
                system(cmd);
                clear cmd
                disp('waiting 10 seconds for writing of temp file before reading. Line:258 (not an error)')
                pause(10);  % to give time to write temp file, failed to load with pause(1)
                STDEV_file=load(stdev_temp_filename); % load stdev of .nii file.
                FDvec_keep_idx = find(FDvec==1); %find the kept frames from the FD mask
                Outlier_file=isthisanoutlier(STDEV_file(FDvec_keep_idx),'median'); %find outlier
                Outlier_idx=find(Outlier_file==1); %find outlier indices
                FDvec(FDvec_keep_idx(Outlier_idx))=0; %set outliers to zero within FDvec
                clear STDEV_file FDvec_keep_idx Outlier_file Outlier_idx
                
            else exist('remove_outliers','var') == 1 && remove_outliers == 0;
                disp('Motion censoring performed on FD alone. Frames with outliers in BOLD std dev not removed');
            end
        else
            %Use power 2014 motion
            load(B{i})
            folders = split(A{i},'/');
            orig_cifti_filename = char(folders(end));
            folder_input = join(folders(1:end-1),'/');
            input_directory = [char(folder_input) '/'];
            if no_output
                output_directory = char(input_directory);
            end
            foldersB = split(B{i},'/');
            orig_motion_filename = char(foldersB(end));
            folderB_input = join(foldersB(1:end-1),'/');
            inputB_directory = [char(folderB_input) '/'];
            if no_output
                output_directory = char(inputB_directory);
            end
            allFD = zeros(1,length(motion_data));
            for j = 1:length(motion_data)
                allFD(j) = motion_data{j}.FD_threshold;
            end
            FDidx = find(round(allFD,3) == round(FD_threshold,3));
            FDvec = motion_data{FDidx}.frame_removal;
            FDvec = abs(FDvec-1);
            
            if strcmp(additional_mask,'none')==1 || strcmp(additional_mask,'None')==1 || strcmp(additional_mask,'NONE')==1
                disp('No additional mask supplied. Using full time series. Frames will be excluded only by FD (unless removal of outliers is indicated).')
            else
                additionalvec = load(additional_mask); %Load .mat file with 0s and 1s.  where 0 are frames to be discarded and 1s are frames to be used to make your matrix.
                FDvec_temp = FDvec & (additionalvec);
                FDvec = FDvec_temp;
            end
            
            if exist('remove_outliers','var') == 0 || remove_outliers == 1
                disp('Removal outliers not specified.  It will be performed by default.')
                
                %% additional frame removal based on Outliers command: isoutlier with "median" method.
                cmd = [wb_command ' -cifti-stats ' char(input_directory) char(orig_cifti_filename) ' -reduce STDEV > ' stdev_temp_filename]
                system(cmd);
                clear cmd
                disp('waiting 60 seconds for writing of temp file before reading. Line:258 (not an error)')
                pause(30);  % to give time to write temp file, failed to load with pause(1)
                STDEV_file=load(stdev_temp_filename); % load stdev of .nii file.
                FDvec_keep_idx = find(FDvec==1); %find the kept frames from the FD mask
                Outlier_file=isthisanoutlier(STDEV_file(FDvec_keep_idx),'median'); %find outlier
                Outlier_idx=find(Outlier_file==1); %find outlier indices
                FDvec(FDvec_keep_idx(Outlier_idx))=0; %set outliers to zero within FDvec
                clear STDEV_file FDvec_keep_idx Outlier_file Outlier_idx
                
            else exist('remove_outliers','var') == 1 && remove_outliers == 0;
                disp('Motion censoring performed on FD alone. Frames with outliers in BOLD std dev not removed');
            end
        end
        
        %%
        if strcmp(minutes_limit,'none')==1 || strcmp(minutes_limit,'None')==1 || strcmp(minutes_limit,'NONE')==1
            
            fileID = fopen([char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) '_cifti_censor_FD_vector_All_Good_Frames.txt'],'w');
            fprintf(fileID,'%1.0f\n',FDvec);
            fclose(fileID);
            
            if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                outfile = [char(output_directory) char(orig_cifti_filename) '_all_frames_at_FD_' num2str(FD_threshold) '_' char(folders(end-1)) '.' suffix2];
            else
                A_smoothed{i} = [char(output_directory) char(orig_cifti_filename(1:length(char(orig_cifti_filename))-13)) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                outfile = [A_smoothed{i} '_all_frames_at_FD_' num2str(FD_threshold) '_' char(folders(end-1)) '.' suffix2];
                if exist(A_smoothed{i}) ==0 %check to see if smoothed file exists yet.
                    Column = 'COLUMN';
                    cmd = [wb_command ' -cifti-smoothing ' char(input_directory) char(orig_cifti_filename) ' ' num2str(smoothing_kernal) ' ' num2str(smoothing_kernal) ' ' Column ' ' char(output_directory) char(orig_cifti_filename(1:length(char(orig_cifti_filename))-13)) '_SMOOTHED_' num2str(smoothing_kernal)  '_' char(folders(end-1)) '.' suffix  ' -left-surface ' C{i} '  -right-surface ' D{i}]
                    disp('320')
                    disp(cmd)
                    system(cmd);
                    clear cmd
                else %smoothed series already exists
                    disp('Smoothed series already created for this subject')
                end
            end
            
            if matrix_not_already_exist(outfile) % if the matrix already exists skip      add folder_input as param?
                if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                    cmd = [wb_command output_precision ' -cifti-correlation ' char(input_directory) char(orig_cifti_filename) ' ' char(output_directory) char(orig_cifti_filename) '_all_frames_at_FD_' num2str(FD_threshold) '_' char(folders(end-1)) '.' suffix2 ' -weights ' char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) '_cifti_censor_FD_vector_All_Good_Frames.txt -fisher-z'];
                    %cmd = [wb_command ' -cifti-correlation ' A{i} ' ' A{i} '_all_frames_at_FD_' num2str(FD_threshold) '.' suffix2 ' -weights ' B{i} '_' num2str(FD_threshold) '_cifti_censor_FD_vector_All_Good_Frames.txt -fisher-z'];
                    tic;
                    disp('Running wb_command cifti-correlation.  This may take a few minutes.')
                    system(cmd);
                    toc;
                    clear cmd
                else
                    cmd = [wb_command output_precision ' -cifti-correlation ' A_smoothed{i} ' ' A_smoothed{i} '_all_frames_at_FD_' num2str(FD_threshold)  '_' char(folders(end-1)) '.' suffix2 ' -weights ' char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) '_cifti_censor_FD_vector_All_Good_Frames.txt -fisher-z'];
                    tic;
                    disp('Running wb_command cifti-correlation.  This may take a few minutes.')
                    system(cmd);
                    toc;
                    clear cmd
                end
            else
                if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                    subject_conn_exists = [char(output_directory) char(orig_cifti_filename) '_all_frames_at_FD_' num2str(FD_threshold)  '_' char(folders(end-1)) '.' suffix2 ' already exists']
                else
                    subject_conn_exists = [A_smoothed{i} '_all_frames_at_FD_' num2str(FD_threshold)  '_' char(folders(end-1)) '.' suffix2 ' already exists']
                end
            end
            
        else %Start minutes limit calculation
            good_frames_idx = find(FDvec == 1);
            good_minutes = (length(good_frames_idx)*TR)/60; % number of good minutes in your data
            
            if good_minutes < 0.5 % if there is less than 30 seconds, don't generate the correlation matrix
                subject_has_too_few_frames = ['Subject ' num2str(i) ' has less than 30 seconds of good data']
                
            elseif minutes_limit > good_minutes % if there is not enough data for your subject, just generate the matrix with all available frames
                fileID = fopen([char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) '_cifti_censor_FD_vector_All_Good_Frames.txt'],'w');
                fprintf(fileID,'%1.0f\n',FDvec);
                fclose(fileID);
                
                if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                    outfile = [char(output_directory) char(orig_cifti_filename) '_all_frames_at_FD_' num2str(FD_threshold) '_' char(folders(end-1)) '.' suffix2];
                else
                    A_smoothed{i} = [char(output_directory) char(orig_cifti_filename(1:length(char(orig_cifti_filename))-13)) '_SMOOTHED_' num2str(smoothing_kernal) '.' suffix];
                    outfile = [A_smoothed{i} '_all_frames_at_FD_' num2str(FD_threshold) '_' char(folders(end-1)) '.' suffix2];
                    if exist(A_smoothed{i}) ==0 %check to see if smoothed file exists yet.
                        Column = 'COLUMN';
                        cmd = [wb_command ' -cifti-smoothing ' char(input_directory) char(orig_cifti_filename) ' ' num2str(smoothing_kernal) ' ' num2str(smoothing_kernal) ' ' Column ' ' char(output_directory) char(orig_cifti_filename(1:length(char(orig_cifti_filename))-13)) '_SMOOTHED_' num2str(smoothing_kernal)  '_' char(folders(end-1)) '.' suffix  ' -left-surface ' C{i} '  -right-surface ' D{i}];
                        disp('373')
                        disp(cmd)
                        system(cmd);
                        clear cmd
                    else %smoothed series already exists
                        disp('Smoothed series already created for this subject')
                    end
                end
                
                if matrix_not_already_exist(outfile) % if the matrix already exists skip
                    if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                        cmd = [wb_command output_precision ' -cifti-correlation ' char(input_directory) char(orig_cifti_filename) ' ' char(output_directory) char(orig_cifti_filename) '_all_frames_at_FD_' num2str(FD_threshold)  '_' char(folders(end-1)) '.' suffix2 ' -weights ' char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) '_cifti_censor_FD_vector_All_Good_Frames.txt -fisher-z'];
                        %cmd = [wb_command ' -cifti-correlation ' A{i} ' ' A{i} '_all_frames_at_FD_' num2str(FD_threshold) '.' suffix2 ' -weights ' B{i} '_' num2str(FD_threshold) '_cifti_censor_FD_vector_All_Good_Frames.txt -fisher-z'];
                        tic;
                        disp('Running wb_command cifti-correlation.  This may take a few minutes.')
                        system(cmd);
                        toc;
                        clear cmd
                    else
                        cmd = [wb_command output_precision ' -cifti-correlation ' A_smoothed{i} ' ' A_smoothed{i} '_all_frames_at_FD_' num2str(FD_threshold) '.' suffix2 ' -weights ' char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) '_cifti_censor_FD_vector_All_Good_Frames'  '_' char(folders(end-1)) '.txt -fisher-z'];
                        tic;
                        disp('Running wb_command cifti-correlation.  This may take a few minutes.')
                        system(cmd);
                        toc;
                        clear cmd
                    end
                    
                else
                    if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                        subject_conn_exists = [char(output_directory) char(orig_cifti_filename) '_all_frames_at_FD_' num2str(FD_threshold)  '_' char(folders(end-1)) '.' suffix2 ' already exists']
                    else
                        subject_conn_exists = [A_smoothed{i} '_all_frames_at_FD_' num2str(FD_threshold)  '_' char(folders(end-1)) '.' suffix2 ' already exists']
                    end
                end
                
                
            elseif minutes_limit <= good_minutes % if there is enough data, match the amount of data used for your subject matrix to the minutes_limit
                good_frames_needed = round(minutes_limit*60/TR); %number of good frames to randomly pull
                rand_good_frames = sort(randperm(length(good_frames_idx),good_frames_needed));
                FDvec_cut = zeros(length(FDvec),1);
                ones_idx = good_frames_idx(rand_good_frames);
                FDvec_cut(ones_idx) = 1; % the new vector that should match good frames with the minutes limit
                
                fileID = fopen([char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) '_cifti_censor_FD_vector_' num2str(minutes_limit) '_minutes_of_data_at_' num2str(FD_threshold) '_threshold.txt'],'w');
                fprintf(fileID,'%1.0f\n',FDvec_cut);
                fclose(fileID);
                if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                    outfile = [char(output_directory) char(orig_cifti_filename) '_' num2str(minutes_limit) '_minutes_of_data_at_FD_' num2str(FD_threshold) '_' char(folders(end-1)) '.' suffix2]
                else % Smooth
                    A_smoothed{i} = [char(output_directory) char(orig_cifti_filename(1:length(char(orig_cifti_filename))-13)) '_SMOOTHED_' num2str(smoothing_kernal) '_' char(folders(end-1)) '.' suffix];
                    outfile = [A_smoothed{i} '_' num2str(minutes_limit) '_minutes_of_data_at_FD_' num2str(FD_threshold) '_' char(folders(end-1)) '.' suffix2]
                    if exist(A_smoothed{i}) == 0
                        Column = 'COLUMN';
                        cmd = [wb_command ' -cifti-smoothing ' char(input_directory) char(orig_cifti_filename) ' ' num2str(smoothing_kernal) ' ' num2str(smoothing_kernal) ' ' Column ' ' char(output_directory) char(orig_cifti_filename(1:length(char(orig_cifti_filename))-13)) '_SMOOTHED_' num2str(smoothing_kernal)  '_' char(folders(end-1)) '.' suffix  ' -left-surface ' C{i} '  -right-surface ' D{i}];
                        disp('426')
                        disp(cmd)
                        system(cmd);
                        clear cmd
                    else %smoothed series already exists
                        disp('Smoothed series already created for this subject')
                    end
                end
                if matrix_not_already_exist(outfile) % check to see if the file already exists
                    if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                        cmd = [wb_command output_precision ' -cifti-correlation ' char(input_directory) char(orig_cifti_filename) ' ' char(output_directory) char(orig_cifti_filename) '_' num2str(minutes_limit) '_minutes_of_data_at_FD_' num2str(FD_threshold) '_' char(folders(end-1)) '.' suffix2 ' -weights ' char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) '_cifti_censor_FD_vector_' num2str(minutes_limit) '_minutes_of_data_at_' num2str(FD_threshold) '_threshold.txt -fisher-z'];
                        %cmd = [wb_command ' -cifti-correlation ' A{i} ' ' A{i} '_' num2str(minutes_limit) '_minutes_of_data_at_FD_' num2str(FD_threshold) '.' suffix2 ' -weights ' B{i} '_' num2str(FD_threshold) '_cifti_censor_FD_vector_' num2str(minutes_limit) '_minutes_of_data_at_' num2str(FD_threshold) '_threshold.txt -fisher-z'];
                        tic;
                        disp('Running wb_command cifti-correlation.  This may take a few minutes.')
                        system(cmd);
                        toc;
                        clear cmd
                    else
                        cmd = [wb_command output_precision ' -cifti-correlation ' A_smoothed{i} ' ' A_smoothed{i} '_' num2str(minutes_limit) '_minutes_of_data_at_FD_' num2str(FD_threshold) '_' char(folders(end-1)) '.' suffix2 ' -weights ' char(output_directory) char(orig_motion_filename) '_' num2str(FD_threshold) '_cifti_censor_FD_vector_' num2str(minutes_limit) '_minutes_of_data_at_' num2str(FD_threshold) '_threshold'  '_' char(folders(end-1)) '.txt -fisher-z'];
                        %cmd = [wb_command ' -cifti-correlation ' A{i} ' ' A{i} '_' num2str(minutes_limit) '_minutes_of_data_at_FD_' num2str(FD_threshold) '.' suffix2 ' -weights ' B{i} '_' num2str(FD_threshold) '_cifti_censor_FD_vector_' num2str(minutes_limit) '_minutes_of_data_at_' num2str(FD_threshold) '_threshold.txt -fisher-z'];
                        tic;
                        disp('Running wb_command cifti-correlation.  This may take a few minutes.')
                        system(cmd);
                        toc;
                        clear cmd
                    end
                else
                    if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
                        subject_conn_exists = [char(output_directory) char(orig_cifti_filename) '_' num2str(minutes_limit) '_minutes_of_data_at_FD_' num2str(FD_threshold)  '_' char(folders(end-1)) '.' suffix2 ' already exists']
                    else
                        subject_conn_exists = [A_smoothed{i} '_' num2str(minutes_limit) '_minutes_of_data_at_FD_' num2str(FD_threshold)  '_' char(folders(end-1)) '.' suffix2 ' already exists']
                    end
                    
                end
                
            else
                'something is wrong with regard to number of good minutes calc'
            end
        end
    end
end

if make_dconn_conc == 1
    disp('Done making (p or d)conn.nii for subjects.  Making output .concs')
    
    %% GENERATE FINAL CONC text file here to use in other code (e.g. cifti_conn_pairwise_corr.m)
    
    if strcmp('conc',conc) == 1
        concname = [dt_or_ptseries_conc_file '_' suffix3 '_of_' series]; %prefix name
    else
        concname = strsplit(dt_or_ptseries_conc_file, '/');
        concname = char(concname(end));
        concname = [concname '_' suffix3 '_of_' series]; %prefix name
    end
    
    %different conc files that potentially could be produced
    subjectswithoutenoughdata = [];
    dconn_paths_all_frames = [];
    dconn_paths_all_frames_at_thresh=[];
    dconn_paths_all_frames_at_thresh_min_lim =[];
    
    %% Generate conc files
    if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
        for i = 1:length(A)
            folders = split(A{i},'/');
            orig_cifti_filename = char(folders(end));
            folder_input = join(folders(1:end-1),'/');
            input_directory = [char(folder_input) '/'];
            if no_output
                output_directory =char(input_directory);
            end
            if strcmp(motion_file,'none')==1 || strcmp(motion_file,'None')==1
                dconn_paths_all_frames = [dconn_paths_all_frames; strcat(char(output_directory),char(orig_cifti_filename), '_all_frames_at_FD_', {FD_threshold}, '.', {suffix2})];
            else
                load(B{i})
                allFD = zeros(1,length(motion_data));
                for j = 1:length(motion_data)
                    allFD(j) = motion_data{j}.FD_threshold;
                end
                FDidx = find(allFD == FD_threshold);
                FDvec = motion_data{FDidx}.frame_removal;
                FDvec = abs(FDvec-1);
                
                %% additional frame removal based on Outliers command: isoutlier with "median" method.
                cmd = [wb_command ' -cifti-stats ' char(input_directory) char(orig_cifti_filename) ' -reduce STDEV > ' stdev_temp_filename];
                system(cmd);
                disp('waiting 60 seconds for writing of temp file before reading. Line:460 (not an error)')
                pause(30);
                STDEV_file=load(stdev_temp_filename); % load stdev of .nii file.
                FDvec_keep_idx = find(FDvec==1); %find the kept frames from the FD mask
                Outlier_file=isthisanoutlier(STDEV_file(FDvec_keep_idx),'median'); %find outlier
                Outlier_idx=find(Outlier_file==1); %find outlier indices
                FDvec(FDvec_keep_idx(Outlier_idx))=0; %set outliers to zero within FDvec
                clear STDEV_file FDvec_keep_idx Outlier_file Outlier_idx
                
                if strcmp(minutes_limit,'none')==1 || strcmp(minutes_limit,'None')==1 || strcmp(minutes_limit,'NONE')==1
                    dconn_paths_all_frames_at_thresh = [dconn_paths_all_frames_at_thresh; strcat(char(output_directory),char(orig_cifti_filename), '_all_frames_at_FD_', {num2str(FD_threshold)},  '_', char(folders(end-1)), '.', {suffix2})];
                else
                    good_frames_idx = find(FDvec == 1);
                    good_minutes = (length(good_frames_idx)*TR)/60;
                    
                    if good_minutes < 0.5
                        subjectswithoutenoughdata = [subjectswithoutenoughdata; strcat(char(output_directory),char(orig_cifti_filename), '_lessthan30sec_', {num2str(FD_threshold)}, '_', char(folders(end-1)), '.', {suffix2})];
                    elseif minutes_limit > good_minutes
                        dconn_paths_all_frames_at_thresh = [dconn_paths_all_frames_at_thresh; strcat(char(output_directory),char(orig_cifti_filename), '_all_frames_at_FD_', {num2str(FD_threshold)}, '_', char(folders(end-1)), '.', {suffix2})];
                    elseif minutes_limit <= good_minutes
                        dconn_paths_all_frames_at_thresh_min_lim = [dconn_paths_all_frames_at_thresh_min_lim; strcat(char(output_directory),char(orig_cifti_filename), '_', {num2str(minutes_limit)}, '_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '_', char(folders(end-1)), '.', {suffix2})];
                    else
                        'something is wrong generating conc files'
                    end
                end
                
            end
        end
    else %use smoothed data
        for i = 1:length(A)
            folders = split(A{i},'/');
            orig_cifti_filename = char(folders(end));
            folder_input = join(folders(1:end-1),'/');
            input_directory = [char(folder_input) '/'];
            if no_output
                output_directory =char(input_directory);
            end
            if strcmp(motion_file,'none')==1 || strcmp(motion_file,'None')==1
                dconn_paths_all_frames = [dconn_paths_all_frames; strcat({A_smoothed{i}}, '_all_frames_at_FD_', {FD_threshold}, '_', char(folders(end-1)), '.', {suffix2})];
            else
                load(B{i})
                allFD = zeros(1,length(motion_data));
                for j = 1:length(motion_data)
                    allFD(j) = motion_data{j}.FD_threshold;
                end
                FDidx = find(allFD == FD_threshold);
                FDvec = motion_data{FDidx}.frame_removal;
                FDvec = abs(FDvec-1);
                
                %% additional frame removal based on Outliers command: isoutlier with "median" method.
                cmd = [wb_command ' -cifti-stats ' char(input_directory) char(orig_cifti_filename) ' -reduce STDEV > ' stdev_temp_filename];
                system(cmd);
                disp('waiting 60 seconds for writing of temp file before reading. Line:505 (not an error)')
                pause(30);
                clear cmd
                STDEV_file=load(stdev_temp_filename); % load stdev of .nii file.
                FDvec_keep_idx = find(FDvec==1); %find the kept frames from the FD mask
                Outlier_file=isthisanoutlier(STDEV_file(FDvec_keep_idx),'median'); %find outlier
                Outlier_idx=find(Outlier_file==1); %find outlier indices
                FDvec(FDvec_keep_idx(Outlier_idx))=0; %set outliers to zero within FDvec
                clear STDEV_file FDvec_keep_idx Outlier_file Outlier_idx
                if strcmp(minutes_limit,'none')==1 || strcmp(minutes_limit,'None')==1 || strcmp(minutes_limit,'NONE')==1
                    dconn_paths_all_frames_at_thresh = [dconn_paths_all_frames_at_thresh; strcat({A_smoothed{i}}, '_all_frames_at_FD_', {num2str(FD_threshold)}, '_', char(folders(end-1)), '.', {suffix2})];
                else
                    good_frames_idx = find(FDvec == 1);
                    good_minutes = (length(good_frames_idx)*TR)/60;
                    
                    if good_minutes < 0.5
                        subjectswithoutenoughdata = [subjectswithoutenoughdata; strcat({A_smoothed{i}}, '_lessthan30sec_', {num2str(FD_threshold)}, '_', char(folders(end-1)), '.', {suffix2})];
                    elseif minutes_limit > good_minutes
                        dconn_paths_all_frames_at_thresh = [dconn_paths_all_frames_at_thresh; strcat({A_smoothed{i}}, '_all_frames_at_FD_', {num2str(FD_threshold)}, '_', char(folders(end-1)), '.', {suffix2})];
                    elseif minutes_limit <= good_minutes
                        dconn_paths_all_frames_at_thresh_min_lim = [dconn_paths_all_frames_at_thresh_min_lim; strcat({A_smoothed{i}}, '_', {num2str(minutes_limit)}, '_minutes_of_data_at_FD_', {num2str(FD_threshold)}, '_', char(folders(end-1)), '.', {suffix2})];
                    else
                        'something is wrong generating conc files'
                    end
                end
            end
            
        end
        
    end
    
    if strcmp(smoothing_kernal,'none')==1 || strcmp(smoothing_kernal,'None')==1 || strcmp(smoothing_kernal,'NONE')==1
        for i = 1:length(A)
            folders = split(A{i},'/');
            orig_cifti_filename = char(folders(end));
            folder_input = join(folders(1:end-1),'/');
            input_directory = [char(folder_input) '/'];
            if no_output
                output_directory =char(input_directory);
            end
            if strcmp(motion_file,'none')==1 || strcmp(motion_file,'None')==1
                fileID = fopen([concname '_all_frames_at_FD_' motion_file '.conc'],'w');
                [nrows,ncols] = size(dconn_paths_all_frames);
                for row = 1:nrows
                    fprintf(fileID,'%s\n' ,dconn_paths_all_frames{row,:});
                end
                fclose(fileID);
            else
                load(B{i})
                allFD = zeros(1,length(motion_data));
                for j = 1:length(motion_data)
                    allFD(j) = motion_data{j}.FD_threshold;
                end
                FDidx = find(allFD == FD_threshold);
                FDvec = motion_data{FDidx}.frame_removal;
                FDvec = abs(FDvec-1);
                
                if strcmp(additional_mask,'none')==1 || strcmp(additional_mask,'None')==1 || strcmp(additional_mask,'NONE')==1
                    disp('No additional mask supplied. Using full time series. Frames will be excluded only by FD (unless removal of outliers is indicated).')
                else
                    additionalvec = load(additional_mask); %Load .mat file with 0s and 1s.  where 0 are frames to be discarded and 1s are frames to be used to make your matrix.
                    FDvec_temp = FDvec & (additionalvec);
                    FDvec = FDvec_temp;
                end
                if exist('remove_outliers','var') == 0 || remove_outliers == 1
                    disp('Removal outliers not specified.  It will be performed by default.')
                    %% additional frame removal based on Outliers command: isoutlier with "median" method.
                    cmd = [wb_command ' -cifti-stats ' char(input_directory) char(orig_cifti_filename) ' -reduce STDEV > ' stdev_temp_filename];
                    system(cmd);
                    disp('waiting 60 seconds for writing of temp file before reading. Line:558 (not an error)')
                    pause(30);
                    clear cmd
                    STDEV_file=load(stdev_temp_filename); % load stdev of .nii file.
                    FDvec_keep_idx = find(FDvec==1); %find the kept frames from the FD mask
                    Outlier_file=isthisanoutlier(STDEV_file(FDvec_keep_idx),'median'); %find outlier
                    Outlier_idx=find(Outlier_file==1); %find outlier indices
                    FDvec(FDvec_keep_idx(Outlier_idx))=0; %set outliers to zero within FDvec
                    clear STDEV_file FDvec_keep_idx Outlier_file Outlier_idx
                    
                else exist('remove_outliers','var') == 1 && remove_outliers == 0;
                    disp('Motion censoring performed on FD alone. Frames with outliers in BOLD std dev not removed');
                end
                
                good_frames_idx = find(FDvec == 1);
                good_minutes = (length(good_frames_idx)*TR)/60;
                
                if strcmp(minutes_limit,'none')==1 || strcmp(minutes_limit,'None')==1 || strcmp(minutes_limit,'NONE')==1
                    fileID = fopen([concname '_all_frames_at_FD_' num2str(FD_threshold) '.conc'],'w');
                    [nrows,ncols] = size(dconn_paths_all_frames_at_thresh);
                    for row = 1:nrows
                        fprintf(fileID,'%s\n' ,dconn_paths_all_frames_at_thresh{row,:});
                    end
                    
                else
                    
                    if good_minutes < 0.5
                        fileID = fopen([concname '_lessthan30sec_' num2str(FD_threshold) '.conc'],'w');
                        [nrows,ncols] = size(subjectswithoutenoughdata);
                        for row = 1:nrows
                            fprintf(fileID,'%s\n' ,subjectswithoutenoughdata{row,:});
                        end
                        fclose(fileID);
                    elseif minutes_limit > good_minutes
                        fileID = fopen([concname '_all_frames_at_FD_' num2str(FD_threshold) '.conc'],'w');
                        [nrows,ncols] = size(dconn_paths_all_frames_at_thresh);
                        for row = 1:nrows
                            fprintf(fileID,'%s\n' ,dconn_paths_all_frames_at_thresh{row,:});
                        end
                        fclose(fileID);
                    elseif minutes_limit <= good_minutes
                        fileID = fopen([concname '_' num2str(minutes_limit) '_minutes_of_data_at_FD_' num2str(FD_threshold) '.conc'],'w');
                        [nrows,ncols] = size(dconn_paths_all_frames_at_thresh_min_lim);
                        for row = 1:nrows
                            fprintf(fileID,'%s\n' ,dconn_paths_all_frames_at_thresh_min_lim{row,:});
                        end
                        fclose(fileID);
                    end
                end
            end
        end
    else %do use smoothed data
        for i = 1:length(A)
            folders = split(A{i},'/');
            orig_cifti_filename = char(folders(end));
            folder_input = join(folders(1:end-1),'/');
            input_directory = [char(folder_input) '/'];
            if no_output
                output_directory =char(input_directory);
            end
            if strcmp(motion_file,'none')==1 || strcmp(motion_file,'None')==1
                fileID = fopen([concname '_all_frames_at_FD_' motion_file '.conc'],'w');
                [nrows,ncols] = size(dconn_paths_all_frames);
                for row = 1:nrows
                    fprintf(fileID,'%s\n' ,dconn_paths_all_frames{row,:});
                end
                fclose(fileID);
            else
                load(B{i})
                allFD = zeros(1,length(motion_data));
                for j = 1:length(motion_data)
                    allFD(j) = motion_data{j}.FD_threshold;
                end
                FDidx = find(allFD == FD_threshold);
                FDvec = motion_data{FDidx}.frame_removal;
                FDvec = abs(FDvec-1);
                
                %% additional frame removal based on Outliers command: isoutlier with "median" method.
                cmd = [wb_command ' -cifti-stats ' char(input_directory) char(orig_cifti_filename) ' -reduce STDEV > ' stdev_temp_filename];
                system(cmd);
                disp('waiting 60 seconds for writing of temp file before reading. Line:627 (not an error)')
                pause(30);
                clear cmd
                STDEV_file=load(stdev_temp_filename); % load stdev of .nii file.
                FDvec_keep_idx = find(FDvec==1); %find the kept frames from the FD mask
                Outlier_file=isthisanoutlier(STDEV_file(FDvec_keep_idx),'median'); %find outlier
                Outlier_idx=find(Outlier_file==1); %find outlier indices
                FDvec(FDvec_keep_idx(Outlier_idx))=0; %set outliers to zero within FDvec
                clear STDEV_file FDvec_keep_idx Outlier_file Outlier_idx
                if strcmp(minutes_limit,'none')==1 || strcmp(minutes_limit,'None')==1 || strcmp(minutes_limit,'NONE')==1
                    fileID = fopen([concname '_SMOOTHED_' num2str(smoothing_kernal) '_all_frames_at_FD_' num2str(FD_threshold) '.conc'],'w');
                    [nrows,ncols] = size(dconn_paths_all_frames_at_thresh);
                    for row = 1:nrows
                        fprintf(fileID,'%s\n' ,dconn_paths_all_frames_at_thresh{row,:});
                    end
                else
                    good_frames_idx = find(FDvec == 1);
                    good_minutes = (length(good_frames_idx)*TR)/60;
                    
                    if good_minutes < 0.5
                        fileID = fopen([concname '_SMOOTHED_' num2str(smoothing_kernal) '_lessthan30sec_' num2str(FD_threshold) '.conc'],'w');
                        [nrows,ncols] = size(subjectswithoutenoughdata);
                        for row = 1:nrows
                            fprintf(fileID,'%s\n' ,subjectswithoutenoughdata{row,:});
                        end
                        fclose(fileID);
                    elseif minutes_limit > good_minutes
                        fileID = fopen([concname '_SMOOTHED_' num2str(smoothing_kernal) '_all_frames_at_FD_' num2str(FD_threshold) '.conc'],'w');
                        [nrows,ncols] = size(dconn_paths_all_frames_at_thresh);
                        for row = 1:nrows
                            fprintf(fileID,'%s\n' ,dconn_paths_all_frames_at_thresh{row,:});
                        end
                        fclose(fileID);
                    elseif minutes_limit <= good_minutes
                        fileID = fopen([concname '_SMOOTHED_' num2str(smoothing_kernal) '_' num2str(minutes_limit) '_minutes_of_data_at_FD_' num2str(FD_threshold) '.conc'],'w');
                        [nrows,ncols] = size(dconn_paths_all_frames_at_thresh_min_lim);
                        for row = 1:nrows
                            fprintf(fileID,'%s\n' ,dconn_paths_all_frames_at_thresh_min_lim{row,:});
                        end
                        fclose(fileID);
                    end
                end
            end
        end
    end
else
end
disp('Done making output concs.')
end


function not_exist = matrix_not_already_exist(to_make)
    
    not_exist = (exist(to_make) == 0);
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

