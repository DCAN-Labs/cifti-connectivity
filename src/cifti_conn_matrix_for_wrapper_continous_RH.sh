#! /bin/sh

# enter code directory
#cd "$( dirname "${BASH_SOURCE[0]}" )"

## Matlab command and usage

#cifti_conn_matrix_for_wrapper_continous(wb_command(1), dt_or_ptseries_conc_file(2), series(3), motion_file(4), FD_threshold(5), TR(6), minutes_limit(7), smoothing_kernel(8), left_surface_file(9), right_surface_file(10), bit8(11), remove_outliers(12), additional_mask_conc(13), make_dconn_conc(14), output_directory(15), dtseries_conc(16), use_continous_minutes(17), memory_limit_value(18)

X="addpath('/home/faird/shared/code/internal/utilities/cifti_connectivity/src/'); cifti_conn_matrix_for_wrapper_continous('${1}', '${2}', '${3}','${4}', '${5}', '${6}', '${7}', '${8}', '${9}', '${10}', '${11}', '${12}', '${13}', '${14}', '${15}', '${16}', '${17}', '${18}')"

#Hermosillo R. 11/1/2022

###########################################################################################
#matlab_exec=/home/exacloud/lustre1/fnl_lab/code/external/GUIs/MATLAB/R2018a/matlab
#matlab_exec=/home/exacloud/lustre1/fnl_lab/code/external/GUIs/MATLAB/R2016b/matlab
matlab_exec=/panfs/roc/msisoft/matlab/R2019a/bin/matlab

RandomHash=`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 16`
Tempmatlabcommand="matlab_command""$RandomHash"".m"

if [ -f "matlab_command""$RandomHash"".m" ]
then
	#echo "matlab_command.m found removing â¦"
	rm -fR "matlab_command""$RandomHash"".m"
fi

#echo ${X} 
echo ${X} > "matlab_command""$RandomHash"".m"
cat "matlab_command""$RandomHash"".m"
${matlab_exec} -nodisplay -nosplash < "matlab_command""$RandomHash"".m"
rm -f "matlab_command""$RandomHash"".m"
