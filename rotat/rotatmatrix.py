import sympy as sp
import numpy as np
from scipy.spatial.transform import Rotation as R

def get_rotmatrix(rot_method=3):

    if rot_method == 1:
        # euler anguler
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

        random_angles = np.random.uniform(0, 2 * np.pi, 3)
        alpha_val, beta_val, gamma_val = random_angles[:3]

        rotation_matrix = R_combined.subs({alpha: alpha_val, beta: beta_val, gamma: gamma_val})

    elif rot_method == 2:
        # 随机生成四元数，转换成旋转矩阵
        u1, u2, u3 = np.random.uniform(0, 1, 3)

        q_w = np.sqrt(1 - u1) * np.sin(2 * np.pi * u2)
        q_x = np.sqrt(1 - u1) * np.cos(2 * np.pi * u2)
        q_y = np.sqrt(u1) * np.sin(2 * np.pi * u3)
        q_z = np.sqrt(u1) * np.cos(2 * np.pi * u3)

        rotation_matrix = np.array([
            [1 - 2*(q_y**2 + q_z**2), 2*(q_x*q_y - q_w*q_z), 2*(q_x*q_z + q_w*q_y)],
            [2*(q_x*q_y + q_w*q_z), 1 - 2*(q_x**2 + q_z**2), 2*(q_y*q_z - q_w*q_x)],
            [2*(q_x*q_z - q_w*q_y), 2*(q_y*q_z + q_w*q_x), 1 - 2*(q_x**2 + q_y**2)]
        ])

    elif rot_method == 3:
        # 随机生成一个四元数表示的旋转
        r = R.random()
        rotation_matrix = r.as_matrix()

    else:
        raise ValueError("rot_method must be 1, 2, or 3")
    
    return rotation_matrix