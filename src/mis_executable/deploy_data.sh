

# global
Rexe=/usr/local/bin/Rscript
transform_matrix_exe=/RAID_32T/fbh/at_ze/jobs/mis_executable/transform_matrix.R
crdata_subpath=/outs/filtered_feature_bc_matrix
ssodata_subpath='/Solo.out/Gene/filtered'
utdata_subpath=/count
pipekit=('CR' 'SSo' 'UtIIST')


# $1: mis_path
# $2: suppress
# $3: suppress

echo [deploy_data.sh] is running...

mis_path=$1
data_fullpath=${mis_path}/data
mis=$(basename "$mis_path")

if [ -d ${data_fullpath} ];then
	echo [deploy_data.sh] ${data_fullpath} is exist, erasing and redeploy.
	rm ${data_fullpath} -r
fi

mkdir ${data_fullpath}

for pipe in ${pipekit[@]}
do
	if [[ ${pipe} == ${2} ]] || [[ ${pipe} == ${3} ]];then
		continue
	fi
	pipe_path=${mis_path}/${pipe}
	if [ -d ${pipe_path} ];then
		mkdir ${data_fullpath}/${pipe}
		for ref_g in `ls ${pipe_path}`
		do
			mkdir ${data_fullpath}/${pipe}/${ref_g}
			ref_g_path=${pipe_path}/${ref_g}
			mat_wth_path=${data_fullpath}/${pipe}/${ref_g}/mat_wth
			mat_lg_path=${data_fullpath}/${pipe}/${ref_g}/mat_lg
			mkdir ${mat_wth_path}
			mkdir ${mat_lg_path}
			
			if [ ${pipe} = 'CR' ];then
				data_path=${ref_g_path}${crdata_subpath}
				crfmt_path=${data_fullpath}/${pipe}/${ref_g}/crfmt
				mkdir ${crfmt_path}
				cp ${data_path}/* ${crfmt_path}
				${Rexe} ${transform_matrix_exe} crfmt2mat_wth ${crfmt_path} ${mat_wth_path}/${mis}_wth.tsv
				gzip -r ${mat_wth_path}
			elif [ ${pipe} = 'SSo' ];then
				data_path=${ref_g_path}${ssodata_subpath}
				crfmt_path=${data_fullpath}/${pipe}/${ref_g}/crfmt
				mkdir ${crfmt_path}
				cp ${data_path}/* ${crfmt_path}
				${Rexe} ${transform_matrix_exe} crfmt2mat_wth ${crfmt_path} ${mat_wth_path}/${mis}_wth.tsv
				gzip -r ${mat_wth_path}
			elif [ ${pipe} = 'UtIIST' ];then
				data_path=${ref_g_path}${utdata_subpath}
				cp ${data_path}/* ${mat_wth_path}
			fi
			
			#gunzip -r ${mat_wth_path}
			#${Rexe} ${transform_matrix_exe} mat_wth2mat_lg ${mat_wth_path}/${mis}_wth.tsv ${mat_lg_path}/${mis}_lg.tsv
			#gzip -r ${mat_wth_path}
			#gzip -r ${mat_lg_path}
		done
	fi
done