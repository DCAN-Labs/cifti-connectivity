function  cifti_conn_pairwise_corr_exaversion(wb_command, pconn_dconn_template,pconn_or_dconn, subject_conn_conc, keep_conn_matrices)

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

if exist(pconn_dconn_template)==0
    disp(pconn_dconn_template);
    NOTE = ['Template does not exist']
    return
else
    [~,filename,extension] = fileparts(pconn_dconn_template);
    pconn_dconn_template_basename = [filename extension];
end

%% make sure files connectivity matrices in conc exist
%conc = strsplit(subject_conn_conc, '/');
conc=subject_conn_conc;
disp(conc)
%conc = char(conc(end));
%disp(conc)
A = importdata(conc);

for i = 1:length(A)
    disp(A{i});
    if exist(A{i}) == 0
        NOTE = ['matrix ' i ' does not exist']
        return
    else
    end
end
'All matrix files exist continuing ...'

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

scalar_paths = [];
for i = 1:length(A)
    dconn_vs_altas = [A{i} '_to_' pconn_dconn_template_basename '.dscalar.nii']; % set scalar name
    if exist(dconn_vs_altas) == 0 % make sure the matrix doesn't already exist
        [template_path,template_name,template_ext] = fileparts(pconn_dconn_template);
        cmd = [wb_command ' -cifti-pairwise-correlation ' pconn_dconn_template ' ' A{i} ' ' A{i} '_to_' template_name '.' suffix];
        disp(cmd)
        tic;
        system(cmd);
        toc;
        clear cmd
    else
        dconn_vs_altas = [A{i} '_to_' pconn_dconn_template_basename '.dscalar.nii already exists']
    end
    
    %add option to delete matrices after making scalars
    if str2num(keep_conn_matrices) ==0
      cmd = ['rm -f ' A{i}];
      disp(strcat('removing matrix: ',(A{i})))
      system(cmd);
      clear cmd
    else
       disp(strcat('keeping matrix: ',(A{i})))
    end
    
    %scalar_paths = [scalar_paths; strcat({pwd}, '/', {subject_conn_conc}, '_', {num2str(i)}, '_to_', {pconn_dconn_template}, '.', {suffix})]; %for scalar conc file
    scalar_paths = [scalar_paths; strcat(A{i}, '_to_', {template_name}, '.', {suffix})]; %for scalar conc file 
end
Summary_of_pairwise = [num2str(i),'  total scalars made from subject list.'];
disp(Summary_of_pairwise);
%% GENERATE FINAL CONC text file here to use in other code ()

%concname = [pwd '/' subject_conn_conc '_vs_' pconn_dconn_template];
[template_path,template_name,template_ext] = fileparts(pconn_dconn_template)
disp(subject_conn_conc)
concname = [subject_conn_conc '_vs_' template_name];
disp([concname '_scalar.conc'])
fileID = fopen([concname '_scalar.conc'],'w');
nrows = length(A);
for row = 1:nrows
    fprintf(fileID,'%s\n' ,scalar_paths{row,:});
end
fclose(fileID);
%clear concname scalar_paths
