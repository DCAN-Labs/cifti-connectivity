function  cifti_conn_pairwise_corr_exaversion(wb_command, pconn_dconn_template,pconn_or_dconn, subject_conn_conc, keep_conn_matrices, output_dir)

%pconn_dconn_template = template correlation matrix to compare data to
%pconn_or_dconn = specify if your data are parcellated or dense 
%subject_conn_conc = conc file that has full paths to your cifti connectivity matrices

%pconn_dconn_template='/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii';
%pconn_or_dconn='dconn';
%subject_conn_conc='/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/cub-sub-NDARINVLWRKNUN1_FNL_preproc_v2_Atlas.dtseries.nii_dconn_of_dtseries_SMOOTHED_1.7_3_minutes_of_data_at_FD_0.35.conc';

%%
%wb_command = '/home/exacloud/lustre1/fnl_lab/code/external/utilities/workbench-9253ac2/bin_rh_linux64/wb_command';%exacloud path
%wb_command = '/home/exacloud/tempwork/fnl_lab/code/external/utilities/workbench-1.3.2/bin_rh_linux64/wb_command';
%wb_command = 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/bin/wb_command'; % workbench command path
%wb_command = '/Applications/workbench/bin_macosx64/wb_command'; % workbench command path
%% make sure template exists (i.e. paths to template)

if ~exist(pconn_dconn_template, 'file')
    disp(pconn_dconn_template);
    disp('Template does not exist.')
    return
else
    [~,pconn_dconn_template_basename,extension] = fileparts(pconn_dconn_template);
    [~, pconn_dconn_template_basename, ~] = fileparts(pconn_dconn_template_basename)
    % pconn_dconn_template_basename = [filename extension];
end

%% make sure files connectivity matrices in conc exist
%conc = strsplit(subject_conn_conc, '/');
disp(subject_conn_conc)
%conc = char(conc(end));
%disp(conc)
A = import_conc_or_file(subject_conn_conc, 'pconn or dconn .conc or file');

%% Generate workbench command to do the pairwise correlation for all matrices in conc file

if pconn_or_dconn == 'pconn'
    suffix = 'pscalar.nii';
    disp(pconn_or_dconn);
elseif pconn_or_dconn == 'dconn'
    suffix = 'dscalar.nii';
    disp(pconn_or_dconn);
else
    'matrices needs to be "pconn" or "dconn"';
    return
end

[template_path, template_name, template_ext] = fileparts(pconn_dconn_template);
[~, template_name, ~] = fileparts(template_name);
scalar_paths = [];
for i = 1:length(A)
    [~, A_i_basename, A_i_ext] = fileparts(A{i});  % set scalar name
    out_conn_path_pt1 = [output_dir filesep A_i_basename A_i_ext '_to_']
    dconn_vs_atlas = [out_conn_path_pt1 pconn_dconn_template_basename ...
                      '.' suffix]; 
    if ~exist(dconn_vs_atlas) % make sure the matrix doesn't already exist
        cmd = [wb_command ' -cifti-pairwise-correlation ' pconn_dconn_template ' ' A{i} ' ' dconn_vs_atlas] % A{i} '_to_' template_name '.' suffix]
        execute_display_and_clear(cmd, cmd);
    else
        dconn_vs_atlas = [A{i} '_to_' pconn_dconn_template_basename '.dscalar.nii already exists']
    end
    
    %add option to delete matrices after making scalars
    msg = strcat('matrix: ', (A{i}));
    if str2num(keep_conn_matrices) == 0
        execute_display_and_clear(['rm -f ' A{i}], ['removing ' msg])
    else
       disp(['keeping ' msg])
    end
    
    %scalar_paths = [scalar_paths; strcat({pwd}, '/', {subject_conn_conc}, '_', {num2str(i)}, '_to_', {pconn_dconn_template}, '.', {suffix})]; %for scalar conc file
    scalar_paths = [scalar_paths; strcat(out_conn_path_pt1, {template_name}, '.', {suffix})]; %for scalar conc file 
end
Summary_of_pairwise = [num2str(i),'  total scalars made from subject list.'];
disp(Summary_of_pairwise);
%% GENERATE FINAL CONC text file here to use in other code

%concname = [pwd '/' subject_conn_conc '_vs_' pconn_dconn_template];
[~, new_scalar_base, new_scalar_ext] = fileparts(subject_conn_conc);
new_scalar_path = [output_dir filesep new_scalar_base new_scalar_ext ...
                   '_vs_' template_name '_scalar.conc'];
disp(['Saving to ' new_scalar_path])
fileID = fopen(new_scalar_path, 'w');
nrows = length(A);
for row = 1:nrows
    fprintf(fileID, '%s\n', scalar_paths{row,:});
end
fclose(fileID);
%clear concname scalar_paths
end


function paths = import_conc_or_file(path, file_type)
    % Import data from a file, or from a .conc list of files
    if is_conc(path)
        paths = importdata(path);
    else
        paths = {path};
    end

    % Verify that all paths point to real files
    for i = 1:length(paths)
        if ~exist(paths{i}, 'file')
            disp(['The ' file_type ' file ' num2str(i) ' does not exist.'])
            return
        end
    end
    disp(['All ' file_type ' files exist. Continuing ...'])
end


function conc_bool = is_conc(file_path)
    % Return 1 if file_path is a .conc file and 0 otherwise
    conc = strsplit(file_path, '.');
    conc = char(conc(end));
    conc_bool = strcmpi(conc, 'conc');
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