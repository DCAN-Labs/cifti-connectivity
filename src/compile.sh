#! /usr/bin/env bash


#rm -rf cifti_conn_template_exaversion
#rm -rf cifti_conn_matrix_exaversion

#rm -rf ${FNL_preproc_dir}/filtered_movement_regressors
#rm -rf ${FNL_preproc_dir}/analyses_v2
#rm -rf ${FNL_preproc_dir}/specific_filtered_movement_regressors


#/mnt/max/software/MATLAB/R2016b/bin/mcc -v -m -R -singleCompThread -R -nodisplay -o cifti_conn_matrix_for_wrapper cifti_conn_matrix_for_wrapper.m -a /mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/isthisanoutlier.m
#/mnt/max/software/MATLAB/R2016b/bin/mcc -v -m -R -singleCompThread -R -nodisplay -o cifti_conn_template_for_wrapper cifti_conn_template_for_wrapper.m -a /mnt/max/shared/code/internal/utilities/cifti_conn_wrapper/src/cifti_conn_matrix_for_wrapper.m -a /mnt/max/shared/code/internal/utilities/cifti_conn_wrapper/src/isthisanoutlier.m
#/mnt/max/software/MATLAB/R2016b/bin/mcc -v -m -R -singleCompThread -R -nodisplay -o cifti_conn_template_exaversion cifti_conn_template_exaversion.m -a /mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/cifti_conn_matrix_exaversion.m -a /mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/isthisanoutlier.m
#/mnt/max/software/MATLAB/R2016b/bin/mcc -v -m -R -singleCompThread -R -nodisplay -o cifti_conn_matrix_exaversion cifti_conn_matrix_exaversion.m -a /mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/isthisanoutlier.m
/mnt/max/software/MATLAB/R2016b/bin/mcc -v -m -R -singleCompThread -R -nodisplay -o cifti_conn_pairwise_corr_for_wrapper cifti_conn_pairwise_corr_for_wrapper.m -a /mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/isthisanoutlier.m 
#/mnt/max/software/MATLAB/R2016b/bin/mcc -v -m -R -singleCompThread -R -nodisplay -o cifti_conn_template_parallel cifti_conn_template_parallel.m -a /mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/cifti_conn_matrix_exaversion.m -a /mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/isthisanoutlier.m 

#${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o filtered_movement_regressors           ${FNL_preproc_dir}/filtered_movement_regressors.m           -a ${HCP_Mat_Path} -a ${FNL_preproc_dir}/scripts
#${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o FNL_preproc_Matlab                     ${FNL_preproc_dir}/FNL_preproc_Matlab.m                     -a ${HCP_Mat_Path} -a ${FNL_preproc_dir}/scripts -a ${framewise_disp_path}
#${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o analyses_v2                            ${FNL_preproc_dir}/analyses_v2.m                            -a ${HCP_Mat_Path} -a ${FNL_preproc_dir}/scripts -a ${framewise_disp_path}
#${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o specific_filtered_movement_regressors  ${FNL_preproc_dir}/specific_filtered_movement_regressors.m  -a ${HCP_Mat_Path} -a ${FNL_preproc_dir}/scripts

#add MCR_CACHE_ROOT to all run scripts for Exahead processing
sed -i '/exe_dir=`dirname "$0"`/a if [ ! -d $TMPDIR/$USER ]; then\n    mkdir $TMPDIR/$USER\nfi\nexport MCR_CACHE_ROOT=$TMPDIR/$USER' run_cifti_conn_matrix_for_wrapper.sh

