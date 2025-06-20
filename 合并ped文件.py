import os

def merge_ped_files(output_file='merged.ped'):
    # 获取当前文件夹中的所有PED文件，排除输出文件本身
    ped_files = [f for f in os.listdir() if f.endswith('.ped') and f != output_file]
    
    if not ped_files:
        print("当前文件夹中没有找到PED文件")
        return
    
    print(f"找到的PED文件: {ped_files}")
    
    # 打开输出文件
    with open(output_file, 'w') as outfile:
        # 遍历每个PED文件并合并内容
        for i, ped_file in enumerate(ped_files):
            print(f"正在处理: {ped_file}")
            with open(ped_file, 'r') as infile:
                content = infile.read().rstrip()  # 移除末尾的换行符
                outfile.write(content)
                # 只有不是最后一个文件时才添加换行符
                if i < len(ped_files) - 1:
                    outfile.write('\n')
    
    print(f"成功合并 {len(ped_files)} 个PED文件到 {output_file}")

if __name__ == "__main__":
    merge_ped_files()