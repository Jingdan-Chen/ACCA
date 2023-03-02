import re
import os

def de_comment(my_string)->str:

    # 定义正则表达式
    pattern = r'^(.*?)(?=#)'

    # 使用正则表达式匹配
    match = re.search(pattern, my_string)

    # 获取匹配结果
    if match:
        result = match.group(1)
        return(result)  # 输出 "hello world "
    else:
        return my_string

def read_config(filename,value_type=str,name_idx=0,value_idx=1)->dict:
    res = {}
    with open(filename,"r",encoding="utf-8") as f:
        cont = f.readlines()
    for item in cont:
        item = de_comment(item.strip())
        if len(item)==0:
            continue
        for obj in ["=",",","，"]: #judge segment
            if obj in item:
                seg = obj
                break
            seg = " "
        item_ = item.split(seg)
        name = item_[name_idx].strip()
        value = value_type(item_[value_idx]) if value_type!=str else item_[1].strip()
        res[name] = value
    return res

args = read_config("./config.txt")
covR = read_config(args["colvR_file"],value_type=float)
mol_path = args["mol_path"]

def is_bond_exist(atom1, atom2, length,thresh=1.25,Covalence_Dict=covR):

    """ 给定长度和两端原子种类，判断该距离下两原子是否成键 """

    judge = length/(Covalence_Dict[atom1] + Covalence_Dict[atom2])
    if judge > thresh:
        return False
    else:
        return True
    
def find_files(dir_path,name): # Doing!!!!!!!!
    res = []
    # 遍历目录中的文件
    for root, dirs, files in os.walk(dir_path):
        for file in files:
            if file.find(name) != -1:
                # return符合条件的文件名
                #res.append(os.path.join(root, file))
                res.append(file)
    return res

class geom():
    def __init__(self,name,energy,file) -> None:
        pass



# Reading energy.txt involving filename and energy
file_list = read_config(mol_path+"/energy.txt",name_idx=1,value_idx=2)
del file_list["Name"]

all_geom = {key:geom()} #Doing!!!!!!!!!!
