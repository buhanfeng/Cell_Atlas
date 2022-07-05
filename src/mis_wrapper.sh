

# global
Rexe=/usr/local/bin/Rscript
pyexe=/RAID_32T/fbh/tools/miniconda3/bin/python
shexe=/usr/bin/sh


# $0
# $1: project_root
# $2: mission_list
# $3: executable file
# $4: parameter, XX,XX,XX
project_root=$1
config_file=$2
executable_file=$3
parameter=`echo $4 | sed 's/,/ /g'`

echo [${0}] is running ... pwd is forced
echo [${0}] project_root set to: $project_root
echo [${0}] config_file set to: $config_file
echo [${0}] executable_file set to: $executable_file
echo [${0}] parameter set to: $parameter

bash_env=`pwd`
echo [${0}] env path is: $bash_env

for line in `cat ${config_file}`
do
	cd ${bash_env}
	line=$line |tr -d '\r\n'
	line=$line |tr -d '\n'
    echo [$0] precessing mis: $line
	mis_path=${project_root}/${line}
	ty="${executable_file##*.}"
	if [ ${ty}x == 'Rx' ];then
		echo "execute cmd: ${Rexe} ${executable_file} ${parameter}"
		${Rexe} ${executable_file} ${mis_path} ${parameter}
	elif [ ${ty}x == 'shx' ];then
		echo "execute cmd: ${shexe} ${executable_file} ${parameter}"
		${shexe} ${executable_file} ${mis_path} ${parameter}
	elif [ ${ty}x == 'pyx' ];then
		echo "execute cmd: ${pyexe} ${executable_file} ${parameter}"
		${pyexe} ${executable_file} ${mis_path} ${parameter}
	fi
done