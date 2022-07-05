

# global
ssodata_subpath='/Solo.out/Gene/filtered'
pipekit=('SSo')


# $1: mis_path
# $2: suppress
# $3: suppress

echo [gzip_SSo.sh] is running...

mis_path=$1
mis=$(basename "$mis_path")

pipe_path=${mis_path}/SSo

echo ${mis_path}

if [ -d ${pipe_path} ];then
    for ref_g in `ls ${pipe_path}`
    do
        ref_g_path=${pipe_path}/${ref_g}
        data_path=${ref_g_path}${ssodata_subpath}
        echo ${data_path}
        gzip -r ${data_path}
    done
fi
