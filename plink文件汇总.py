import os
import shutil
import glob

def find_and_copy_plink_files(root_folder, output_folder, start_number=1):
    """
    遍历指定文件夹下的3.Analysis文件夹，找到plink.map和plink_ForwardStrand.ped文件
    复制并重命名为n.map、n.ped、(n+1).map、(n+1).ped等，从start_number开始
    
    参数:
        root_folder: 根文件夹路径
        output_folder: 输出文件夹路径
        start_number: 开始的编号，默认为1
    """
    print(f"正在搜索文件夹: {root_folder}")
    
    # 确保输出文件夹存在
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"创建输出文件夹: {output_folder}")
    
    # 查找所有的plink.map和plink_ForwardStrand.ped文件
    found_file_pairs = []
    for dirpath, dirnames, filenames in os.walk(root_folder):
        map_file_path = os.path.join(dirpath, "plink.map")
        ped_file_path = os.path.join(dirpath, "plink_ForwardStrand.ped")
        
        if os.path.isfile(map_file_path) and os.path.isfile(ped_file_path):
            found_file_pairs.append((map_file_path, ped_file_path))
    
    print(f"找到 {len(found_file_pairs)} 对plink文件")
    
    # 计数器，用于重命名文件，从指定的起始编号开始
    counter = start_number
    
    # 处理每个找到的文件对
    for map_file, ped_file in found_file_pairs:
        print(f"处理文件对: {map_file} 和 {ped_file}")
        
        # 创建新的文件名
        new_map_file = os.path.join(output_folder, f"{counter}.map")
        new_ped_file = os.path.join(output_folder, f"{counter}.ped")
        
        # 复制并重命名文件
        shutil.copy2(map_file, new_map_file)
        shutil.copy2(ped_file, new_ped_file)
        
        print(f"已复制并重命名为: {new_map_file} 和 {new_ped_file}")
        
        # 增加计数器
        counter += 1
    if not found_file_pairs:
        print("未找到任何plink文件对。")

def main():
    # 获取当前脚本所在目录作为根目录
    root_folder = os.path.dirname(os.path.abspath(__file__))
    
    # 询问用户是否要使用其他目录
    use_current_dir = input(f"是否使用当前目录 ({root_folder}) 作为搜索根目录? (y/n): ").strip().lower()
    
    if use_current_dir != 'y':
        root_folder = input("请输入要搜索的根目录路径: ").strip()
        # 检查输入的路径是否存在
        if not os.path.exists(root_folder):
            print(f"错误: 目录 '{root_folder}' 不存在!")
            return
    
    # 设置输出文件夹
    output_folder = os.path.join(root_folder, "汇总的Plink文件")
    
    # 询问用户是否要使用其他输出目录
    use_default_output = input(f"是否使用默认输出目录 ({output_folder})? (y/n): ").strip().lower()
    
    if use_default_output != 'y':
        output_folder = input("请输入输出目录路径: ").strip()
    
    # 询问用户从几号开始命名
    start_number_str = input("请输入起始编号 (默认为1): ").strip()
    start_number = 1  # 默认值
    
    # 验证输入是否为有效的数字
    if start_number_str:
        try:
            start_number = int(start_number_str)
            if start_number < 1:
                print("警告: 起始编号必须大于等于1，将使用默认值1")
                start_number = 1
        except ValueError:
            print("警告: 输入的不是有效的数字，将使用默认值1")
    
    # 执行文件查找和复制
    find_and_copy_plink_files(root_folder, output_folder, start_number)
    
    print(f"处理完成! 所有文件已复制到: {output_folder}")
    print(f"文件编号从 {start_number} 开始")

if __name__ == "__main__":
    main()