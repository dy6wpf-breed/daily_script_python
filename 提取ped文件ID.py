# 提取PED文件第二列ID的Python脚本
input_file = r"Y:\种质资源数据管理平台2024\9-自贡德康\德康芯片数据释放-KPSNY032024002C-250516\merged_part5.ped"
output_file = r"Y:\种质资源数据管理平台2024\9-自贡德康\德康芯片数据释放-KPSNY032024002C-250516\ids.txt"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        columns = line.strip().split()
        if len(columns) >= 2:  # 确保每行至少有两列
            outfile.write(columns[1] + '\n')  # 写入第二列的ID

print("ID提取完成！")