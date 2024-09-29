#!/bin/bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
# Additional commands for my local environment
export SPHENIX=/sphenix/u/jzhang1
export MYINSTALL=$SPHENIX/install

# Setup MYINSTALL to local directory and run sPHENIX setup local script
# to adjust PATH, LD LIBRARY PATH, ROOT INCLUDE PATH, etc
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

echo "sPHENIX environment setup finished"

# 定义第一个参数
for EQMD_struct1 in {1..5}; do
    # 内层循环从 EQMD_struct1 的值开始，避免重复组合
    for EQMD_struct2 in $(seq $EQMD_struct1 5); do
        hadd "/sphenix/user/jzhang1/nuclear-structure/output/method1/NeEQMD/NeNe_EQMD${EQMD_struct1}${EQMD_struct2}_10k.root" "/sphenix/user/jzhang1/nuclear-structure/output/method1/NeEQMD/NeNe_EQMD${EQMD_struct1}${EQMD_struct2}_run*root"

        rm -rf "/sphenix/user/jzhang1/nuclear-structure/output/method1/NeEQMD/NeNe_EQMD${EQMD_struct1}${EQMD_struct2}_run*root"
    done
done




