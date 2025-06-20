import pandas as pd
import os

def classify_pig(row):
    dd = row['DD']
    ll = row['LL']
    yy = row['YY']
    
    # 纯种杜洛克
    if dd > 0.90:
        return "纯种杜洛克"
    
    # 纯种长白
    if ll >= 0.85 and yy < 0.11 and dd < 0.07:
        return "纯种长白"
    
    # 纯种大白
    if yy >= 0.85 and ll < 0.11 and dd < 0.07:
        return "纯种大白"
    
    # 二元种猪
    if ll + yy >= 0.93 and min(ll, yy) >= 0.15:
        return "二元种猪"
    
    # 三元猪
    if dd > 0.10 and ll > 0.15 and yy > 0.15:
        return "三元猪"
    
    # 其他
    return "其他"

# 读取 CSV 文件
input_path = r"Y:\新版中芯一号\血统分析\KPSNY132025018C_KPS202503057K_new(1)\admixture.3.csv"
df = pd.read_csv(input_path)

# 应用分类函数
df['品种类别'] = df.apply(classify_pig, axis=1)

# 准备输出路径
output_dir = os.path.dirname(input_path)
output_path = os.path.join(output_dir, "pig_classification_results.csv")

# 输出结果到新的CSV文件
df.to_csv(output_path, index=False, encoding='utf-8-sig')

print(f"分类结果已保存到: {output_path}")

# 打印统计信息
print("\n各品种数量统计:")
print(df['品种类别'].value_counts())