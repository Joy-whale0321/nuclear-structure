# main function

import random
import math
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import pandas as pd
import time
import ROOT
import argparse
import os

from scipy.spatial.transform import Rotation as R

import genOOstruct.py
from readdocx import extract_xyz_from_docx
from rotatmatrix import get_rotmatrix
from calepsilon import calculate_Qn

parser = argparse.ArgumentParser(description='Process nuclear system.')
parser.add_argument('--sys', type=str, default="NeNe", help='Input the nuclear system')
parser.add_argument('--runnum', type=int, default=0, help='Input the condor number')
parser.add_argument('--struct1', type=int, default=1, help='NeNe_EQMD_struct')
parser.add_argument('--struct2', type=int, default=1, help='NeNe_EQMD_struct')
args = parser.parse_args()

Nuclear_system = args.sys
runnum = args.runnum
print(f'The nuclear system is: {Nuclear_system}')

EQMD_struct1 = args.struct1
EQMD_struct2 = args.struct2

# output_dir = "/sphenix/user/jzhang1/nuclear-structure/output/method2/NeEQMD/"
# os.system(f"mkdir -p {output_dir}")
output_dir = "/mnt/e/git-repo/nuclear-structure/rotat/"
root_filename = f"{output_dir}{Nuclear_system}_EQMD{EQMD_struct1}{EQMD_struct2}_run{runnum}.root"
print(f"Creating ROOT file: {root_filename}")

root_file = ROOT.TFile(root_filename, "RECREATE")
tree_pos = ROOT.TTree("nucleon_position_origin", "nucleon position origin data")
tree_rot = ROOT.TTree("nucleon_position_rotate", "nucleon position rotate data")

x1, y1, z1 = np.zeros(1, dtype=float), np.zeros(1, dtype=float), np.zeros(1, dtype=float)
x2, y2, z2 = np.zeros(1, dtype=float), np.zeros(1, dtype=float), np.zeros(1, dtype=float)

tree_pos.Branch("x1", x1, "x1/D")
tree_pos.Branch("y1", y1, "y1/D")
tree_pos.Branch("z1", z1, "z1/D")
tree_pos.Branch("x2", x2, "x2/D")
tree_pos.Branch("y2", y2, "y2/D")
tree_pos.Branch("z2", z2, "z2/D")

tree_rot.Branch("x1", x1, "x1/D")
tree_rot.Branch("y1", y1, "y1/D")
tree_rot.Branch("z1", z1, "z1/D")
tree_rot.Branch("x2", x2, "x2/D")
tree_rot.Branch("y2", y2, "y2/D")
tree_rot.Branch("z2", z2, "z2/D")

c2_2_array = []
epsilon2_Q2_array = []

nevents = 1
epsilon2_array = np.zeros(nevents)

for events_i in range(nevents):    
    if Nuclear_system == "OO":
        cluster_type = "square"  # 可选 "tetrahedron，square"
        cluster_origins = get_cluster_origins(cluster_type)

        nucleons_group1 = generate_nucleon_positions(cluster_origins)  # 前16个核子 
        nucleons_group2 = generate_nucleon_positions(cluster_origins)  # 后16个核子
        
    elif Nuclear_system == "NeNe":
        nucleons_group1 = extract_xyz_from_docx(EQMD_struct1)
        nucleons_group2 = extract_xyz_from_docx(EQMD_struct2)

    else:
        raise ValueError(f"未知的 cluster_type: {cluster_type}")

    # 得到旋转矩阵，将 Sympy 矩阵转换为 NumPy 矩阵，方便数值计算; 对每组核子应用旋转
    rotation_matrix1 = get_rotmatrix(2)
    rotation_matrix2 = get_rotmatrix(2)
 
    R1_np = np.array(rotation_matrix1).astype(np.float64)
    R2_np = np.array(rotation_matrix2).astype(np.float64)

    # 旋转！
    nucleons_group1_rotated = np.dot(nucleons_group1, R1_np.T)
    nucleons_group2_rotated = np.dot(nucleons_group2, R2_np.T)
    
    # 获取 group1 和 group2 中的核子数量
    num_group1 = nucleons_group1_rotated.shape[0]
    num_group2 = nucleons_group2_rotated.shape[0]
    nucleon_totalnumber = nucleons_group2_rotated.shape[0]

    distances_squared = np.zeros((num_group1, num_group2))
    status_vector = np.zeros(2*nucleon_totalnumber)
    R_nucleon = 0.85
    # 计算每对核子的距离平方，只使用x和y坐标（忽略z坐标）
    for group_i in range(num_group1):
        for group_j in range(num_group2):
            delta_x = nucleons_group1_rotated[group_i, 0] - nucleons_group2_rotated[group_j, 0]
            delta_y = nucleons_group1_rotated[group_i, 1] - nucleons_group2_rotated[group_j, 1]
            distances_squared[group_i, group_j] = delta_x**2 + delta_y**2

            if distances_squared[group_i, group_j] < (2*R_nucleon)**2:
                status_vector[group_i - 1]  = 1  # 标记 group1 中核子 i 的状态为 1
                status_vector[group_j - 1 + nucleon_totalnumber] = 1  # 标记 group2 中核子 j 的状态为 1

    # 提取参与碰撞核子信息; 计算二阶偏心度
    participant_x = []
    participant_y = []
    for cal_i in range(len(status_vector)):
        if status_vector[cal_i] == 1:
            if cal_i < nucleon_totalnumber:
                npart_x, npart_y = nucleons_group1_rotated[cal_i, 0], nucleons_group1_rotated[cal_i, 1]                
            else:
                npart_x, npart_y = nucleons_group2_rotated[cal_i - nucleon_totalnumber, 0], nucleons_group2_rotated[cal_i - nucleon_totalnumber, 1]        
            participant_x.append(npart_x)
            participant_y.append(npart_y)

    # 将参与者的坐标转换为 numpy 数组，方便后续计算  
    participant_x = np.array(participant_x)
    participant_y = np.array(participant_y)
    participant_phi = np.arctan2(participant_y, participant_x)

    # 计算eccentricity
    if len(participant_x) > 0:    
        x2_avg = np.mean(participant_x ** 2)
        y2_avg = np.mean(participant_y ** 2)
        xy_avg = np.mean(participant_x * participant_y)  
        epsilon_2 = np.sqrt((x2_avg - y2_avg) ** 2 + 4 * xy_avg ** 2) / (x2_avg + y2_avg)

        epsilon2_array[events_i] = epsilon_2
    else:
        epsilon2_array[events_i] = 1.1

    # 调用calculate_Qn函数
    Q2_event = calculate_Qn(participant_phi, 2)
    M_event = len(participant_phi)

    c2_2_event = (abs(Q2_event) ** 2 - M_event) / (M_event * (M_event-1))
    epsilon2_Q2_event = np.sqrt(c2_2_event)

    c2_2_array.append(c2_2_event)
    epsilon2_Q2_array.append(epsilon2_Q2_event)

    # 填充数据
    for i in range(nucleons_group1.shape[0]):  # 遍历每一行
        x1[0] = nucleons_group1[i, 0]  
        y1[0] = nucleons_group1[i, 1]  
        z1[0] = nucleons_group1[i, 2]  
        x2[0] = nucleons_group2[i, 0]  
        y2[0] = nucleons_group2[i, 1]  
        z2[0] = nucleons_group2[i, 2]  
        tree_pos.Fill()

    for i in range(nucleons_group1_rotated.shape[0]):  # 遍历每一行
        x1[0] = nucleons_group1_rotated[i, 0]  
        y1[0] = nucleons_group1_rotated[i, 1]  
        z1[0] = nucleons_group1_rotated[i, 2]  
        x2[0] = nucleons_group2_rotated[i, 0]  
        y2[0] = nucleons_group2_rotated[i, 1]  
        z2[0] = nucleons_group2_rotated[i, 2]  
        tree_rot.Fill()

# Q-cummulant to epsilon
if len(c2_2_array) > 0:
    c2_2_ensemble_average = np.mean(c2_2_array)
    epsilon2_Q2_ave = np.sqrt(c2_2_ensemble_average)
else:
    print("No valid events found for eccentricity calculation.")

hist_epsilon2_xyz = ROOT.TH1D("epsilon2_xyz", "epsilon2_xyz", 300, -1, 2)
for epsilon2_xyz_i in epsilon2_array:
    hist_epsilon2_xyz.Fill(epsilon2_xyz_i)

hist_epsilon2_Q2 = ROOT.TH1D("epsilon2_Q2", "epsilon2_Q2", 300, -1, 2)
for epsilon2_Q2_i in epsilon2_Q2_array:
    hist_epsilon2_Q2.Fill(epsilon2_Q2_i)

hist_epsilon2_xyz.Write()
hist_epsilon2_Q2.Write()
tree_pos.Write()
tree_rot.Write()

root_file.Close()
