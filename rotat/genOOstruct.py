import numpy as np
import math

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

    return np.array(positions)i
