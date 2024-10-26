# check the rotation method is ok?

import random
import math
import sympy as sp
import numpy as np
from scipy.spatial.transform import Rotation as R
import time
import ROOT

def quaternion_to_u1_u2_u3(q):
    q_w, q_x, q_y, q_z = q

    u1 = q_y**2 + q_z**2

    u2 = np.arctan2(q_w, q_x) / (2 * np.pi) + 0.5

    u3 = np.arctan2(q_y, q_z) / (2 * np.pi) + 0.5

    return u1, u2, u3

root_file = ROOT.TFile("/mnt/e/git-repo/nuclear-structure/output/rot_test2.root", "RECREATE")

histogram1 = ROOT.TH1D("u1", "u1", 1500, -1, 2)
histogram2 = ROOT.TH1D("u2", "u2", 1500, -1, 2)
histogram3 = ROOT.TH1D("u3", "u3", 1500, -1, 2)

hist_alpha = ROOT.TH1D("alpha", "alpha", 500, -np.pi, np.pi)
hist_beta  = ROOT.TH1D("beta",  "beta",  500, -np.pi, np.pi)
hist_gamma = ROOT.TH1D("gamma", "gamma", 500, -np.pi, np.pi)

start_time = time.time()

nevents = 100000
for events_i in range(nevents): 
    
    random_rotation = R.random()
    quaternion = random_rotation.as_quat()

    u1, u2, u3 = quaternion_to_u1_u2_u3(quaternion)
    
    rotation = R.from_quat([quaternion[1], quaternion[2], quaternion[3], quaternion[0]])

    histogram1.Fill(u1)
    histogram2.Fill(u2)
    histogram3.Fill(u3)

    # 将旋转对象转换为欧拉角 (空间固定 ZYX 顺序, 按照 'xyz' 对应 fixed ZYX)
    euler_angles = rotation.as_euler('xyz', degrees=False)
    alpha, beta, gamma = euler_angles
    
    hist_alpha.Fill(alpha)
    hist_beta.Fill(beta)
    hist_gamma.Fill(gamma)

histogram1.Write()
histogram2.Write()
histogram3.Write()

hist_alpha.Write()
hist_beta.Write()
hist_gamma.Write()

root_file.Close()

end_time3 = time.time()
print(f"代码运行时间: {end_time3 - start_time} 秒")
