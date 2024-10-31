import matplotlib.pyplot as plt
import os

def plot_nucleons(x_coords1, y_coords1, x_coords2, y_coords2, x_coords_part, y_coords_part, radius_nucleon=0.85,radius_part=0.80, output_filename="nucleon_plot.pdf"):
    # 创建输出目录
    output_dir = "/mnt/e/git-repo/nuclear-structure/output/display"
    os.makedirs(output_dir, exist_ok=True)
    
    fig, ax = plt.subplots()
    
    # 绘制核子位置
    for (x1, y1) in zip(x_coords1, y_coords1):
        circle1 = plt.Circle((x1, y1), radius_nucleon, edgecolor='green', linewidth=0.01, facecolor='none', zorder=2)
        ax.add_patch(circle1)
        ax.text(x1, y1, 'G', color='green', ha='center', va='center')  # 标注绿色圆的坐标
    
    for (x2, y2) in zip(x_coords2, y_coords2):        
        circle2 = plt.Circle((x2, y2), radius_nucleon, edgecolor='blue', linewidth=0.01, facecolor='none', zorder=2)
        ax.add_patch(circle2)
        ax.text(x2, y2, 'B', color='blue', ha='center', va='center')  # 标注蓝色圆的坐标

    for (x_part, y_part) in zip(x_coords_part, y_coords_part):        
        circle_part = plt.Circle((x_part, y_part), radius_part, edgecolor='red', linewidth=0.01, facecolor='red', zorder=1)
        ax.add_patch(circle_part)
    
    # 设置边界
    all_x_coords = list(x_coords1) + list(x_coords2) + list(x_coords_part)
    all_y_coords = list(y_coords1) + list(y_coords2) + list(y_coords_part)
    x_min, x_max = -10, 10
    y_min, y_max = -10, 10

    ax.set_aspect('equal', 'box')
    ax.set_xlim(x_min, x_max)  
    ax.set_ylim(y_min, y_max)  
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Nucleon Positions')
    plt.grid(True)
    
    # 保存图片
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Figure saved to {output_path}")



