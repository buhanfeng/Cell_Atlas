

# global



# $1: mis_path
# $2: pipe
# $3: ref_g
# $4: data_type

echo [pool_data.sh] is running...

mis_path=$1
pipe=$2
ref_g=$3
data_type=$4

mis=$(basename "$mis_path")
root=$(dirname "$mis_path")

des=${root}/${pipe}_${ref_g}_${data_type}`date -d yesterday +%Y%m%d`

if [ -d ${des} ];then
    echo 'warning:' ${des} is exist, please check.
else
    mkdir ${des}
fi

cp ${mis_path}/data/${pipe}/${ref_g}/${data_type}/${mis}* ${des}


