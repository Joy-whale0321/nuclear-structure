import random
import math
import sympy as sp
import numpy as np
from shapely.geometry import Polygon
from scipy.spatial.transform import Rotation as R

import time
import ROOT

def quaternion_to_u1_u2_u3(q):
    q_w, q_x, q_y, q_z = q

    # 还原 u1
    u1 = q_y**2 + q_z**2

    # 还原 u2
    u2 = np.arctan2(q_w, q_x) / (2 * np.pi) + 0.5

    # 还原 u3
    u3 = np.arctan2(q_y, q_z) / (2 * np.pi) + 0.5

    return u1, u2, u3

root_file = ROOT.TFile("/mnt/e/git-repo/nuclear-structure/output/rot_test.root", "RECREATE")

# 定义旋转角度符号变量
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

histogram1 = ROOT.TH1D("u1", "u1", 1500, -1, 2)
histogram2 = ROOT.TH1D("u2", "u2", 1500, -1, 2)
histogram3 = ROOT.TH1D("u3", "u3", 1500, -1, 2)

start_time = time.time()

nevents = 100000
for events_i in range(nevents): 
    angles = [random.uniform(0, 2 * math.pi) for _ in range(3)]

    alpha1_val, beta1_val, gamma1_val = angles[:3]

    R1_evaluated = R_combined.subs({alpha: alpha1_val, beta: beta1_val, gamma: gamma1_val})

    quaternion = R.from_matrix(R1_evaluated).as_quat()

    u1, u2, u3 = quaternion_to_u1_u2_u3(quaternion)
    
    histogram1.Fill(u1)
    histogram2.Fill(u2)
    histogram3.Fill(u3)

histogram1.Write()
histogram2.Write()
histogram3.Write()

root_file.Close()

end_time3 = time.time()
print(f"代码运行时间: {end_time3 - start_time} 秒")