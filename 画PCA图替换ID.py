# 定义文件路径
input_file_path = r"Z:\临时下载文件\大好河山芯片数据\大好河山芯片数据\pca_output.eigenvec"
replace_file_path = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【1】德康\德康互助场-挪威LL\【4】KPSNY032024002C_KPS202412117_德康-互猪种猪场_朱康平&罗毅&陈鹏宇&李肖晗_113个-20241230\KPSNY032024002C_朱康平、罗毅、陈鹏宇、李肖晗_德康-互猪种猪场_113个_中芯一号PLUS_report_20241230\3.Analysis\td.txt"
output_file_path = r"Z:\临时下载文件\大好河山芯片数据\大好河山芯片数据\pca_output_replaced.eigenvec"

# 读取替换表
replace_dict = {}
with open(replace_file_path, "r") as f:
    for line in f:
        old_id, new_id = line.strip().split()
        replace_dict[old_id] = new_id

# 打开输入文件并创建输出文件
with open(input_file_path, "r") as input_file, open(output_file_path, "w") as output_file:
    for line in input_file:
        columns = line.strip().split()
        fam_id = columns[0]
        ind_id = columns[1]

        # 如果 ind_id 在替换表中，进行替换
        if ind_id in replace_dict:
            columns[1] = replace_dict[ind_id]

        # 写入替换后的行
        output_file.write(" ".join(columns) + "\n")

print("ID替换完成，结果已保存到 pca_output_replaced.eigenvec")
