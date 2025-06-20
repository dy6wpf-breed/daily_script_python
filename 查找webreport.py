import os
import shutil
import time

def find_and_rename_web_report_folders(root_dir, target_dir):
    # 确保目标目录存在
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    found = False
    # 遍历根目录下的所有文件夹
    for root, dirs, files in os.walk(root_dir):
        print(f"正在遍历目录: {root}")  # 添加调试信息，查看当前遍历的目录
        for dir_name in dirs:
            if dir_name.endswith("Web_report") and '.' in dir_name:
                parts = dir_name.split('.')
                if parts[0].isdigit():
                    found = True
                    web_report_path = os.path.join(root, dir_name)
                    # 获取文件夹的创建时间
                    creation_time = os.path.getctime(web_report_path)
                    # 将创建时间转换为可读的日期格式
                    creation_date = time.strftime("%Y-%m-%d", time.localtime(creation_time))
                    # 生成新的文件夹名称
                    new_folder_name = f"{dir_name}_{creation_date}"
                    new_folder_path = os.path.join(target_dir, new_folder_name)

                    # 复制文件夹到目标目录
                    try:
                        shutil.copytree(web_report_path, new_folder_path)
                        print(f"已将 {web_report_path} 复制到 {new_folder_path}")
                    except Exception as e:
                        print(f"复制 {web_report_path} 时出错: {e}")

    if not found:
        print("未找到符合条件的文件夹。")

if __name__ == "__main__":
    # 要遍历的根目录
    root_directory = r"G:\项目数据备份\芯片数据\中芯一号\GS分析\【9】自贡德康_尖山场"
    # 汇总的目标目录
    target_directory = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【1】德康\德康-4个场数据\webreport汇总"

    find_and_rename_web_report_folders(root_directory, target_directory)