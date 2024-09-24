import random
import math
import sympy as sp
import numpy as np
from shapely.geometry import Polygon

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


# 定义顶点
point_A = sp.Matrix([1, 1, 0])
point_B = sp.Matrix([1, -1, 0])
point_C = sp.Matrix([-1, -1,0])
point_D = sp.Matrix([-1, 1, 0])

# 旋转顶点
rotated_A = R_combined * point_A
rotated_B = R_combined * point_B
rotated_C = R_combined * point_C
rotated_D = R_combined * point_D

# 投影顶点
projected_A = rotated_A.extract([0, 1], [0])  # 取第0列的前两行
projected_B = rotated_B.extract([0, 1], [0])  # 取第0列的前两行
projected_C = rotated_C.extract([0, 1], [0])  # 取第0列的前两行
projected_D = rotated_D.extract([0, 1], [0])  # 取第0列的前两行

# 生成6个随机角度(0到2π之间的弧度)，
angles = [random.uniform(0, 2 * math.pi) for _ in range(6)]

alpha1_val, beta1_val, gamma1_val = angles[:3]
alpha2_val, beta2_val, gamma2_val = angles[3:]

A1_val = projected_A.subs({alpha: alpha1_val, beta: beta1_val, gamma: gamma1_val})
B1_val = projected_B.subs({alpha: alpha1_val, beta: beta1_val, gamma: gamma1_val})
C1_val = projected_C.subs({alpha: alpha1_val, beta: beta1_val, gamma: gamma1_val})
D1_val = projected_D.subs({alpha: alpha1_val, beta: beta1_val, gamma: gamma1_val})

A2_val = projected_A.subs({alpha: alpha2_val, beta: beta2_val, gamma: gamma2_val})
B2_val = projected_B.subs({alpha: alpha2_val, beta: beta2_val, gamma: gamma2_val})
C2_val = projected_C.subs({alpha: alpha2_val, beta: beta2_val, gamma: gamma2_val})
D2_val = projected_D.subs({alpha: alpha2_val, beta: beta2_val, gamma: gamma2_val})

print(type(B1_val))
print(B1_val.shape)

# 确保 A1_val, B1_val, C1_val, D1_val 是 (x, y) 的数值
A1_coords = (A1_val[0], A1_val[1])  # (x, y)
B1_coords = (B1_val[0], B1_val[1])  
C1_coords = (C1_val[0], C1_val[1])  
D1_coords = (D1_val[0], D1_val[1])  

A2_coords = (A2_val[0], A2_val[1])  
B2_coords = (B2_val[0], B2_val[1])  
C2_coords = (C2_val[0], C2_val[1])  
D2_coords = (D2_val[0], D2_val[1])  

# 定义两个 square，使用 AB, BC, CD, DA 的顺序构造
square1 = Polygon([A1_coords, B1_coords, C1_coords, D1_coords])  
square2 = Polygon([A2_coords, B2_coords, C2_coords, D2_coords])

# # 计算两个多边形的交集
intersection = square1.intersection(square2)

if not intersection.is_empty:
    intersection_coords = list(intersection.exterior.coords)
    print(f"交叠区域的边界顶点: {intersection_coords}")

    # 获取交叠区域的质心
    intersection_centroid = intersection.centroid
    print(f"交叠区域的质心坐标: {intersection_centroid}")

else:
    print("没有交叠区域")


# def calculate_eccentricity_inertia(polygon):

#     # 获取多边形的顶点
#     coords = np.array(polygon.exterior.coords)

#     # 计算质心平移
#     centroid = polygon.centroid.coords[0]
#     translated_coords = coords - centroid

#     # 计算惯性矩的分量 (Ix, Iy, Ixy)
#     Ix = np.sum((translated_coords[:, 1] ** 2))  # y^2 分量
#     Iy = np.sum((translated_coords[:, 0] ** 2))  # x^2 分量
#     Ixy = -np.sum(translated_coords[:, 0] * translated_coords[:, 1])  # xy 分量

#     inertia_matrix = np.array([[Ix, Ixy], [Ixy, Iy]])

#     # 计算惯性矩矩阵的特征值（主惯性矩）
#     eigenvalues, _ = np.linalg.eig(inertia_matrix)

#     # 计算二阶偏心度（eccentricity = sqrt(1 - min(eigenvalue) / max(eigenvalue)))
#     max_eigenvalue = np.max(eigenvalues)
#     min_eigenvalue = np.min(eigenvalues)
    
#     # 计算二阶偏心度
#     eccentricity = np.sqrt(1 - min_eigenvalue / max_eigenvalue)

#     return eccentricity


# def calculate_eccentricity(polygon, n=2):
#     # 获取多边形的顶点
#     coords = np.array(polygon.exterior.coords)
    
#     # 计算质心平移
#     centroid = polygon.centroid.coords[0]
#     translated_coords = coords - centroid

#     # 初始化各项的累加
#     sum_rn_cosnphi = 0
#     sum_rn_sinnphi = 0
#     sum_rn = 0

#     # 遍历多边形的所有顶点
#     for coord in translated_coords:
#         x, y = coord
#         r = np.sqrt(x**2 + y**2)  # 距离
#         phi = atan2(y, x)  # 极角

#         # 计算 r^n * cos(n*phi) 和 r^n * sin(n*phi)
#         sum_rn_cosnphi += (r**n) * cos(n * phi)
#         sum_rn_sinnphi += (r**n) * sin(n * phi)
#         sum_rn += r**n  # r^n 的和

#     # 计算偏心度
#     numerator = np.sqrt(sum_rn_cosnphi**2 + sum_rn_sinnphi**2)  # 分子部分
#     denominator = sum_rn  # 分母部分

#     epsilon_n = numerator / denominator if denominator != 0 else 0

#     return epsilon_n