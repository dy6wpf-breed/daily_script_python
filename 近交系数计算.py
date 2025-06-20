import pandas as pd
import numpy as np
from tqdm import tqdm

def calculate_inbreeding(pedigree):
    inbreeding = {}
    calculated_relationships = {}
    
    # 预处理：创建一个字典，用于快速查找父母
    parent_dict = dict(zip(pedigree['个体号'], zip(pedigree['父号'], pedigree['母号'])))
    
    for individual in tqdm(pedigree['个体号'], desc="计算近交系数"):
        inbreeding[individual] = calculate_individual_inbreeding(individual, parent_dict, inbreeding, calculated_relationships)
    return inbreeding

def calculate_individual_inbreeding(individual, parent_dict, inbreeding, calculated_relationships):
    if individual in inbreeding:
        return inbreeding[individual]
    
    if individual not in parent_dict:
        inbreeding[individual] = 0
        return 0
    
    sire, dam = parent_dict[individual]
    if pd.isna(sire) or pd.isna(dam):
        inbreeding[individual] = 0
        return 0
    
    f_sire = calculate_individual_inbreeding(sire, parent_dict, inbreeding, calculated_relationships)
    f_dam = calculate_individual_inbreeding(dam, parent_dict, inbreeding, calculated_relationships)
    
    common_ancestors = find_common_ancestors(sire, dam, parent_dict)
    
    f = 0
    for ancestor in common_ancestors:
        p_sire = calculate_relationship(ancestor, sire, parent_dict, calculated_relationships)
        p_dam = calculate_relationship(ancestor, dam, parent_dict, calculated_relationships)
        f_ancestor = inbreeding.get(ancestor, 0)
        f += 0.5 * p_sire * p_dam * (1 + f_ancestor)
    
    inbreeding[individual] = f
    return f

def find_common_ancestors(individual1, individual2, parent_dict):
    ancestors1 = get_ancestors(individual1, parent_dict)
    ancestors2 = get_ancestors(individual2, parent_dict)
    return list(set(ancestors1) & set(ancestors2))

def get_ancestors(individual, parent_dict):
    ancestors = set()
    queue = [individual]
    while queue:
        current = queue.pop(0)
        if current in parent_dict:
            sire, dam = parent_dict[current]
            if pd.notna(sire):
                ancestors.add(sire)
                queue.append(sire)
            if pd.notna(dam):
                ancestors.add(dam)
                queue.append(dam)
    return list(ancestors)

def calculate_relationship(ancestor, descendant, parent_dict, calculated_relationships):
    if ancestor == descendant:
        return 1
    
    key = (ancestor, descendant)
    if key in calculated_relationships:
        return calculated_relationships[key]
    
    if descendant not in parent_dict:
        calculated_relationships[key] = 0
        return 0
    
    sire, dam = parent_dict[descendant]
    if pd.isna(sire) or pd.isna(dam):
        calculated_relationships[key] = 0
        return 0
    
    r_sire = calculate_relationship(ancestor, sire, parent_dict, calculated_relationships)
    r_dam = calculate_relationship(ancestor, dam, parent_dict, calculated_relationships)
    
    r = 0.5 * (r_sire + r_dam)
    calculated_relationships[key] = r
    return r

# 读取Excel文件
file_path = r"E:\desk\系谱文件.xlsx"
pedigree = pd.read_excel(file_path)

print("开始计算近交系数...")

# 计算近交系数
inbreeding_coefficients = calculate_inbreeding(pedigree)

# 计算平均近交系数
average_inbreeding = np.mean(list(inbreeding_coefficients.values()))

print(f"\n平均近交系数: {average_inbreeding}")

# 输出每个个体的近交系数
print("\n各个体的近交系数:")
for individual, coefficient in inbreeding_coefficients.items():
    print(f"个体 {individual} 的近交系数: {coefficient}")
    # ... existing code ...

# 输出每个个体的近交系数
print("\n各个体的近交系数:")
for individual, coefficient in inbreeding_coefficients.items():
    print(f"个体 {individual} 的近交系数: {coefficient}")

# 将结果保存到Excel文件
output_df = pd.DataFrame(list(inbreeding_coefficients.items()), columns=['个体号', '近交系数'])
output_file_path = r"E:\desk\近交系数结果.xlsx"
output_df.to_excel(output_file_path, index=False)
print(f"\n近交系数结果已保存到 {output_file_path}")