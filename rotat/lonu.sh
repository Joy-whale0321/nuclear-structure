#!/bin/bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
# Additional commands for my local environment
export SPHENIX=/sphenix/u/jzhang1
export MYINSTALL=$SPHENIX/install

# Setup MYINSTALL to local directory and run sPHENIX setup local script
# to adjust PATH, LD LIBRARY PATH, ROOT INCLUDE PATH, etc
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

echo "sPHENIX environment setup finished"

# 定义Python脚本路径
PYTHON_SCRIPT="/sphenix/user/jzhang1/nuclear-structure/rotat/nucleonposition_withrot.py"

# 定义第一个参数
Nuclear_system="NeNe"

condornum=$2

for EQMD_struct1 in {1..5}; do
    # 内层循环从 EQMD_struct1 的值开始，避免重复组合
    for EQMD_struct2 in $(seq $EQMD_struct1 5); do
        python $PYTHON_SCRIPT --sys $Nuclear_system --runnum $condornum --struct1 EQMD_struct1 --struct2 EQMD_struct2
    done
done




