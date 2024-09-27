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

from scipy.spatial.transform import Rotation as R

from readdocx import extract_xyz_from_docx
from rotatmatrix import get_rotmatrix


def get_cluster_origins(cluster_type):
    if cluster_type == "square":
        return np.array([
            [1.5, 1.5, 0],    # cluster A
            [1.5, -1.5, 0],   # cluster B
            [-1.5, -1.5, 0],  # cluster C
            [-1.5, 1.5, 0]    # cluster D
        ])
    elif cluster_type == "tetrahedron":
        return np.array([
            [3.5/(2*math.sqrt(2)), 3.5/(2*math.sqrt(2)), 3.5/(2*math.sqrt(2))],    # cluster A
            [-3.5/(2*math.sqrt(2)), 3.5/(2*math.sqrt(2)), 3.5/(2*math.sqrt(2))],   # cluster B
            [3.5/(2*math.sqrt(2)),-3.5/(2*math.sqrt(2)), 3.5/(2*math.sqrt(2))],  # cluster C
            [3.5/(2*math.sqrt(2)),3.5/(2*math.sqrt(2)), -3.5/(2*math.sqrt(2))]    # cluster D
        ])
    else:
        raise ValueError("Unknown cluster type. Please choose 'square' or 'tetrahedron'.")

def generate_nucleon_positions(cluster_origins):
    num_nucleons = 16
    nucleon_radius = 0.85

    positions = []
    d_min = 2 * nucleon_radius  # 两个核子之间的最小距离，取为两倍核子半径

    # 定义R高斯分布的参数
    mu = 0  # 均值
    sigma = 1.23  # 标准差
    normalize_guassR = 1 / (sigma * np.sqrt(2 * np.pi))

    for i in range(num_nucleons):
        while True:
            r = np.random.normal(loc=mu, scale=sigma)
            theta = np.arccos(np.random.uniform(-1, 1))
            phi = np.random.uniform(0, 2 * np.pi)

            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)

            new_position = np.array([x, y, z])

            cluster_number = i // 4
            new_position += cluster_origins[cluster_number]

            if not positions:
                positions.append(new_position)
                break

            # 检查新生成的核子与已有核子之间的距离
            distances = [np.linalg.norm(new_position - pos) for pos in positions]
            
            # if len(distances) > 0:
            #     print(f"核子 {i} 与已有核子的距离: {distances}")
            
            if all(distance >= d_min for distance in distances):
                positions.append(new_position)
                break

    return np.array(positions)

start_time = time.time()

parser = argparse.ArgumentParser(description='Process nuclear system.')
parser.add_argument('--sys', type=str, default="OO", help='Input the nuclear system')
parser.add_argument('--runnum', type=int, default=0, help='Input the condor number')
args = parser.parse_args()

Nuclear_system = args.sys
runnum = args.runnum
print(f'The nuclear system is: {Nuclear_system}')

output_dir = "/sphenix/user/jzhang1/nuclear-structure/output/method1/"
root_filename = f"{output_dir}{Nuclear_system}_run{runnum}.root"
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

nevents = 1000
epsilon_2_array = np.zeros(nevents)

for events_i in range(nevents):    
    if Nuclear_system == "OO":
        cluster_type = "square"  # 可选 "tetrahedron，square"
        cluster_origins = get_cluster_origins(cluster_type)

        nucleons_group1 = generate_nucleon_positions(cluster_origins)  # 前16个核子 
        nucleons_group2 = generate_nucleon_positions(cluster_origins)  # 后16个核子
        
    elif Nuclear_system == "NeNe":
        nucleons_group1 = extract_xyz_from_docx()
        nucleons_group2 = extract_xyz_from_docx()

    else:
        raise ValueError(f"未知的 cluster_type: {cluster_type}")

    # print("Nucleons Group 1:")
    # print(nucleons_group1)

    # 得到旋转矩阵，将 Sympy 矩阵转换为 NumPy 矩阵，方便数值计算; 对每组核子应用旋转
    rotation_matrix1 = get_rotmatrix(1)
    rotation_matrix2 = get_rotmatrix(1)
 
    R1_np = np.array(rotation_matrix1).astype(np.float64)
    R2_np = np.array(rotation_matrix2).astype(np.float64)

    # 旋转！
    nucleons_group1_rotated = np.dot(nucleons_group1, R1_np.T)
    nucleons_group2_rotated = np.dot(nucleons_group2, R2_np.T)
    # nucleons_group1_rotated = nucleons_group1
    # nucleons_group2_rotated = nucleons_group2
    
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
    if len(participant_x) > 0:    
        x2_avg = np.mean(participant_x ** 2)
        y2_avg = np.mean(participant_y ** 2)
        xy_avg = np.mean(participant_x * participant_y)  
        epsilon_2 = np.sqrt((x2_avg - y2_avg) ** 2 + 4 * xy_avg ** 2) / (x2_avg + y2_avg)

        epsilon_2_array[events_i] = epsilon_2

    else:
        epsilon_2_array[events_i] = 1.1

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

end_time3 = time.time()
print(f"代码运行时间: {end_time3 - start_time} 秒")

histogram = ROOT.TH1D("epsilon_2_hist", "Epsilon_2 Distribution", 300, -1, 2)
for epsilon_2 in epsilon_2_array:
    histogram.Fill(epsilon_2)

histogram.Write()
tree_pos.Write()
tree_rot.Write()

root_file.Close()
