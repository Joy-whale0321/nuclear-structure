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


