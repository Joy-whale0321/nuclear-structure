import random
import math
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import pandas as pd
import time
import ROOT


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

# 从高斯分布中生成样本
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


root_file = ROOT.TFile("/mnt/e/git-repo/nuclear-structure/output/tetrahedron_wo_overlap_100k.root", "RECREATE")

# 旋转 and 旋转矩阵
alpha, beta, gamma = sp.symbols('alpha beta gamma')

R_x = sp.Matrix([
    [1, 0, 0],
    [0, sp.cos(alpha), -sp.sin(alpha)],
    [0, sp.sin(alpha), sp.cos(alpha)]
])

R_y = sp.Matrix([
    [sp.cos(beta), 0, sp.sin(beta)],
    [0, 1, 0],
    [-sp.sin(beta), 0, sp.cos(beta)]
])

R_z = sp.Matrix([
    [sp.cos(gamma), -sp.sin(gamma), 0],
    [sp.sin(gamma), sp.cos(gamma), 0],
    [0, 0, 1]
])
R_combined = R_z * R_y * R_x

start_time = time.time()

nevents = 100000
epsilon_2_array = np.zeros(nevents)
# 初始化一个空的 DataFrame 来存储所有run的数据
all_runs_data = pd.DataFrame()

for events_i in range(nevents):
    cluster_type = "tetrahedron"  # 可选 "tetrahedron，square"
    cluster_origins = get_cluster_origins(cluster_type)

    # 2个O原子核 两组，每组16个核子
    nucleons_group1 = generate_nucleon_positions(cluster_origins) # 前16个核子 
    nucleons_group2 = generate_nucleon_positions(cluster_origins) # 后16个核子

    # 随机生成 0 到 2π 之间的 6 个数
    random_angles = np.random.uniform(0, 2 * np.pi, 6)

    alpha1_val, beta1_val, gamma1_val = random_angles[:3]
    alpha2_val, beta2_val, gamma2_val = random_angles[3:]

    R1_evaluated = R_combined.subs({alpha: alpha1_val, beta: beta1_val, gamma: gamma1_val})
    R2_evaluated = R_combined.subs({alpha: alpha2_val, beta: beta2_val, gamma: gamma2_val})

    # 将 Sympy 矩阵转换为 NumPy 矩阵，方便数值计算; 对每组核子应用旋转
    R1_np = np.array(R1_evaluated).astype(np.float64)
    R2_np = np.array(R2_evaluated).astype(np.float64)

    nucleons_group1_rotated = np.dot(nucleons_group1, R1_np.T)
    nucleons_group2_rotated = np.dot(nucleons_group2, R2_np.T)
    # nucleons_group1_rotated = nucleons_group1
    # nucleons_group2_rotated = nucleons_group2
    
    # 获取 group1 和 group2 中的核子数量
    num_group1 = nucleons_group1_rotated.shape[0]
    num_group2 = nucleons_group2_rotated.shape[0]

    distances_squared = np.zeros((num_group1, num_group2))
    status_vector = np.zeros(32)
    R_nucleon = 0.85
    # 计算每对核子的距离平方，只使用x和y坐标（忽略z坐标）
    for group_i in range(num_group1):
        for group_j in range(num_group2):
            delta_x = nucleons_group1_rotated[group_i, 0] - nucleons_group2_rotated[group_j, 0]
            delta_y = nucleons_group1_rotated[group_i, 1] - nucleons_group2_rotated[group_j, 1]
            distances_squared[group_i, group_j] = delta_x**2 + delta_y**2

            if distances_squared[group_i, group_j] < (2*R_nucleon)**2:
                status_vector[group_i - 1]  = 1  # 标记 group1 中核子 i 的状态为 1
                status_vector[group_j + 15] = 1  # 标记 group2 中核子 j 的状态为 1

    # 提取参与碰撞核子信息; 计算二阶偏心度
    participant_x = []
    participant_y = []
    for cal_i in range(len(status_vector)):
        if status_vector[cal_i] == 1:
            if cal_i < 16:
                npart_x, npart_y = nucleons_group1[cal_i, 0], nucleons_group1[cal_i, 1]                
            else:
                npart_x, npart_y = nucleons_group2[cal_i - 16, 0], nucleons_group2[cal_i - 16, 1]        
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
                
    # # 将 group1 和 group2 转换为 DataFrame，并命名列 (x, y, z)
    # df_group1 = pd.DataFrame(nucleons_group1_rotated, columns=['x', 'y', 'z'])
    # df_group2 = pd.DataFrame(nucleons_group2_rotated, columns=['x', 'y', 'z'])
    # df_group1['group'] = 'group1'
    # df_group2['group'] = 'group2'
    # df_group1['run'] = f'run_{events_i+1}'  # 运行编号
    # df_group2['run'] = f'run_{events_i+1}'  # 运行编号

    # # 合并 group1 和 group2 的数据
    # df_combined = pd.concat([df_group1, df_group2], ignore_index=True)

    # # 在每个 run 数据末尾添加一行 NaN 作为空行
    # empty_row = pd.DataFrame([[np.nan] * df_combined.shape[1]], columns=df_combined.columns)
    # df_combined = pd.concat([df_combined, empty_row], ignore_index=True)

    # # 追加到所有 runs 的数据
    # all_runs_data = pd.concat([all_runs_data, df_combined], ignore_index=True)

# all_runs_data.to_excel('/mnt/c/Users/12896/Desktop/togpt/nucleons_tetrahedron_runs.xlsx', index=False)

end_time = time.time()
print(f"代码运行时间: {end_time - start_time} 秒")

histogram = ROOT.TH1D("epsilon_2_hist", "Epsilon_2 Distribution", 300, -1, 2)

for epsilon_2 in epsilon_2_array:
    histogram.Fill(epsilon_2)

histogram.Write()

root_file.Close()
