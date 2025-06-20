import os
import shutil
import re

# 源文件夹路径
source_dir = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【1】德康\德康-武乐种场"
# 目标文件夹路径(与源文件夹相同)
target_dir = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【1】德康\德康-武乐种场"

# 用于存储找到的文件对
file_pairs = []

# 遍历源文件夹及其子文件夹
for root, dirs, files in os.walk(source_dir):
    map_files = [f for f in files if f == "plink.map"]
    ped_files = [f for f in files if f == "plink_ForwardStrand.ped"]
    
    # 如果在当前文件夹中找到了两个文件
    if map_files and ped_files:
        map_file = os.path.join(root, map_files[0])
        ped_file = os.path.join(root, ped_files[0])
        
        # 将文件对添加到列表中
        file_pairs.append((map_file, ped_file))

# 复制并重命名文件
for i, (map_file, ped_file) in enumerate(file_pairs, 1):
    # 构建新的文件名
    new_map_name = f"{i}.map"
    new_ped_name = f"{i}.ped"
    
    # 构建目标文件的完整路径
    target_map = os.path.join(target_dir, new_map_name)
    target_ped = os.path.join(target_dir, new_ped_name)
    
    # 复制文件
    shutil.copy2(map_file, target_map)
    shutil.copy2(ped_file, target_ped)
    
    print(f"已复制并重命名第{i}对文件:")
    print(f"  {map_file} -> {target_map}")
    print(f"  {ped_file} -> {target_ped}")

print(f"总共找到并处理了{len(file_pairs)}对文件。")