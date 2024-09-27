import numpy as np
from docx import Document

def extract_xyz_from_docx(textpage=1):
    """
    从 docx 文件中提取指定页（段）的 x, y, z 坐标，并返回一个 20x3 的 numpy 数组。
    默认读取第一页（第1个段落群）。
    """

    doc = Document('/mnt/c/Users/12896/Desktop/cluster/cluster-Ne20-EQMD.docx')

    xyz_coords = []

    start_reading = False
    line_count = 0  # 记录行数

    # 根据 textpage 参数来选择从哪个部分开始读取
    current_page = 0  # 当前处理的页面或段落群编号

    # 遍历文档的所有段落
    for paragraph in doc.paragraphs:
        line = paragraph.text.strip()

        # 每次遇到特定的关键字 "X", "z", "pz" 作为新的部分的开头
        if "X" in line and "z" in line and "pz" in line:
            current_page += 1
            # 继续寻找我们指定的页面（段落群）
            if current_page < textpage:
                continue
            start_reading = True
            line_count = 0  # 重置行计数
            continue
        
        # 如果已经进入指定的 textpage 页面，提取前3列
        if start_reading and line:
            # 只读取第2行到第21行
            if 1 <= line_count <= 20:
                values = line.split()
                if len(values) >= 3:  # 确保至少有3个值
                    x, y, z = values[:3]
                    xyz_coords.append([float(x), float(y), float(z)])  # 转换为浮点数并存储

            line_count += 1

        if line_count > 20:
            break

    return np.array(xyz_coords)

