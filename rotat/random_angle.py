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


# 生成6个随机角度(0到2π之间的弧度)，
angles = [random.uniform(0, 2 * math.pi) for _ in range(6)]

alpha1_val, beta1_val, gamma1_val = angles[:3]
alpha2_val, beta2_val, gamma2_val = angles[3:]

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
