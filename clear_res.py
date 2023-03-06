# import os
import shutil
from coordProcess import *

def clear_res():
    # 设置目录路径
    dir_path = args["res_file"]

    # 遍历目录并删除所有文件和文件夹
    for filename in os.listdir(dir_path):
        file_path = os.path.join(dir_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                # 如果是文件或符号链接，直接删除
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                # 如果是文件夹，递归删除
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")
