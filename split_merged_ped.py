import os

def split_ped_by_ids(merged_ped_file, id_files_dir, output_dir):
    """
    根据提供的ID列表文件分割merged.ped文件。

    :param merged_ped_file: merged.ped文件的路径。
    :param id_files_dir: 包含Sample_ID列表的txt文件所在的目录。
    :param output_dir: 分割后的ped文件输出目录。
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"创建输出目录: {output_dir}")

    # 预先读取merged.ped文件的所有行到内存，如果文件非常巨大，可以考虑逐行读取和处理
    print(f"正在读取主文件: {merged_ped_file}")
    with open(merged_ped_file, 'r') as f_merged:
        merged_lines = f_merged.readlines()
    print(f"主文件读取完成，共 {len(merged_lines)} 行。")

    for i in range(1, 6):  # 假设有5个ID文件
        id_file_name = f"1-样本基因型对应关系表{i}_Sample_IDs.txt"
        id_file_path = os.path.join(id_files_dir, id_file_name)
        
        output_ped_filename = f"merged_part{i}.ped"
        output_ped_path = os.path.join(output_dir, output_ped_filename)

        if not os.path.exists(id_file_path):
            print(f"ID文件 {id_file_path} 不存在，跳过。")
            continue

        print(f"正在处理ID文件: {id_file_path}")
        with open(id_file_path, 'r') as f_ids:
            # 读取ID，并去除每行末尾的换行符
            sample_ids = set(line.strip() for line in f_ids)
        
        print(f"从 {id_file_name} 中读取了 {len(sample_ids)} 个Sample_ID。")
        
        lines_written = 0
        with open(output_ped_path, 'w') as f_out:
            for line in merged_lines:
                parts = line.strip().split() # PED文件通常是空格或制表符分隔
                if len(parts) > 1:
                    ped_sample_id = parts[1] # 第二列是样本ID
                    if ped_sample_id in sample_ids:
                        f_out.write(line)
                        lines_written += 1
        
        if lines_written > 0:
            print(f"成功将 {lines_written} 行数据写入到: {output_ped_path}")
        else:
            print(f"没有找到与 {id_file_name} 匹配的数据行，未创建或内容为空: {output_ped_path}")

if __name__ == "__main__":
    base_dir = r"y:\种质资源数据管理平台2024\9-自贡德康\德康芯片数据释放-KPSNY032024002C-250516"
    merged_ped_file_path = os.path.join(base_dir, "merged.ped")
    # ID文件和merged.ped在同一个目录下
    id_files_directory = base_dir 
    # 输出文件也放在这个目录下，可以根据需要修改
    output_directory = base_dir 

    split_ped_by_ids(merged_ped_file_path, id_files_directory, output_directory)
    print("所有文件处理完成。")