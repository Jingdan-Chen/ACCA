import re
import os
import numpy as np
from scipy.spatial import distance_matrix
from scipy.sparse import csr_matrix
import networkx as nx
# import matplotlib.pyplot as plt
import pandas as pd

pt_dict = {1: 'H',
 2: 'He',
 3: 'Li',
 4: 'Be',
 5: 'B',
 6: 'C',
 7: 'N',
 8: 'O',
 9: 'F',
 10: 'Ne',
 11: 'Na',
 12: 'Mg',
 13: 'Al',
 14: 'Si',
 15: 'P',
 16: 'S',
 17: 'Cl',
 18: 'Ar',
 19: 'K',
 20: 'Ca',
 21: 'Sc',
 22: 'Ti',
 23: 'V',
 24: 'Cr',
 25: 'Mn',
 26: 'Fe',
 27: 'Co',
 28: 'Ni',
 29: 'Cu',
 30: 'Zn',
 31: 'Ga',
 32: 'Ge',
 33: 'As',
 34: 'Se',
 35: 'Br',
 36: 'Kr',
 37: 'Rb',
 38: 'Sr',
 39: 'Y',
 40: 'Zr',
 41: 'Nb',
 42: 'Mo',
 43: 'Tc',
 44: 'Ru',
 45: 'Rh',
 46: 'Pd',
 47: 'Ag',
 48: 'Cd',
 49: 'In',
 50: 'Sn',
 51: 'Sb',
 52: 'Te',
 53: 'I',
 54: 'Xe',
 55: 'Cs',
 56: 'Ba',
 57: 'La',
 58: 'Ce',
 59: 'Pr',
 60: 'Nd',
 61: 'Pm',
 62: 'Sm',
 63: 'Eu',
 64: 'Gd',
 65: 'Tb',
 66: 'Dy',
 67: 'Ho',
 68: 'Er',
 69: 'Tm',
 70: 'Yb',
 71: 'Lu',
 72: 'Hf',
 73: 'Ta',
 74: 'W',
 75: 'Re',
 76: 'Os',
 77: 'Ir',
 78: 'Pt',
 79: 'Au',
 80: 'Hg',
 81: 'Tl',
 82: 'Pb',
 83: 'Bi',
 84: 'Po',
 85: 'At',
 86: 'Rn',
 87: 'Fr',
 88: 'Ra',
 89: 'Ac',
 90: 'Th',
 91: 'Pa',
 92: 'U',
 93: 'Np',
 94: 'Pu',
 95: 'Am',
 96: 'Cm',
 97: 'Bk',
 98: 'Cf',
 99: 'Es',
 100: 'Fm',
 101: 'Md',
 102: 'No',
 103: 'Lr',
 104: 'Rf',
 105: 'Db',
 106: 'Sg',
 107: 'Bh',
 108: 'Hs'}

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


def is_bond_exist(atom1, atom2, length,thresh=float(args["bonding_thresh"]),Covalence_Dict=covR):

    """ 给定长度和两端原子种类，判断该距离下两原子是否成键 """

    judge = length/(Covalence_Dict[atom1] + Covalence_Dict[atom2])
    if judge > thresh:
        return False
    else:
        return True
    

def extract(file,beginl=0): #extract and process xyz and SCF energy information in IRC out
# cont::list/["line_1","line_2"..."line_x"] the content of file to be processedd
# beginl::int/ begin line of the target part
# endl::int/ end line of the targer part
# num::int/ The index of the irc structure to be processed
    with open(file,"r") as f:
        cont = f.readlines()
    endl = len(cont)

    # scf_lines = [line.strip("\n") for line in cont if 'SCF Done' in line]
    # if scf_lines:
    #     spe = scf_lines[-1].split()[4]
    # else:
    #     spe = 'N/A'
    # Gaussian output args
    geom_segl = "---------------------------------------------------------------------"
    geom_begin = "Number     Number       Type             X           Y           Z"
    
    # Args define
    raw_xyz = []
    processed_xyz = []
    scf = None
    geom_flag = 0
    idx = beginl

    # detect the line of information (xyz and scf)
    while idx<endl:
        # print(idx)
        temp = cont[idx].strip("\n").strip()
        if geom_flag==0: # in front of xyz
            if temp==geom_begin:
                geom_flag = 1 #recieving
                idx+=1
        elif geom_flag==1: # in xyz
            if temp!=geom_segl:
                raw_xyz.append(temp)
            else:
                geom_flag=2
        else: # behind xyz
            if temp[0:8]=="SCF Done":
                scf = float(temp.split()[4])
                break
        idx+=1

    # process xyz
    ele_id_lis = [] # record the order of the elements (written in numbers in out file)
    for item in raw_xyz:
        temp = item.split()
        ele = pt_dict[int(temp[1])]
        ele_id_lis.append(ele)
        # new_lis = [ele]+temp[3:]+["\n"]
        # processed_xyz.append("    ".join(new_lis))
    # processed_xyz = [str(len(raw_xyz))+"\n","IRC:"+str(num)+"\n"]+processed_xyz

# raw_xyz::list/ raw xyz information in out file
# processed_xyz::list/ xyz information in the formal .xyz style
# ele_id_lis::list/ elements with order of the molecule
    return raw_xyz,scf,ele_id_lis

def xyz2mat(xyz): # convert xyz information to numpy matrix
# xyz: formal .xyz; result: numpy.array

    temp = list(map(lambda x:list(map(float,x.split()[-3:])),xyz[:]))
    result = np.array([l for l in temp])
    return result

def compute_internal(xyz_lis,args,numpy_flag=True): # compute internal arguments for a given xyz(numpy based)
# xyz_lis::list/ xyz information in formal .xyz style(or numpy format)
# args::list/[(arg1),(arg2)] bond args to be calculated 
# numpy_flag::bool/ False:input formal .xyz style; True:numpy format

    res = []
    for item in args:
        if len(item)==2: # Bond Distance
            if numpy_flag:
                vec1 = xyz_lis[item[0]-1]
                vec2 = xyz_lis[item[1]-1]
            else:
                new_xyzlis = xyz_lis[1:]
                xyz1 = new_xyzlis[item[0]].split()[-3:]
                xyz2 = new_xyzlis[item[1]].split()[-3:]
                vec1 = np.array(list(map(float,xyz1)))
                vec2 = np.array(list(map(float,xyz2)))
            dist = np.linalg.norm(vec1-vec2)
            res.append(dist)

        elif len(item)==3: # Bond Angle
            if numpy_flag:
                vec1 = xyz_lis[item[0]-1]-xyz_lis[item[1]-1]
                vec2 = xyz_lis[item[2]-1]-xyz_lis[item[1]-1]
            else:
                xyz1 = new_xyzlis[item[0]].split()[-3:]
                xyz2 = new_xyzlis[item[1]].split()[-3:]
                xyz3 = new_xyzlis[item[2]].split()[-3:]
                vec1 = np.array(list(map(float,xyz1)))-np.array(list(map(float,xyz2)))
                vec2 = np.array(list(map(float,xyz3)))-np.array(list(map(float,xyz2)))
            angle = np.arccos(vec1.dot(vec2)/(np.linalg.norm(vec1) * np.linalg.norm(vec2)))/np.pi*180
            res.append(angle)

        elif len(item)==4: # Dihedral Angle
            if numpy_flag:
                vec1 = xyz_lis[item[0]-1]-xyz_lis[item[1]-1]
                vec2 = xyz_lis[item[1]-1]-xyz_lis[item[2]-1]
                vec3 = xyz_lis[item[3]-1]-xyz_lis[item[2]-1]
            else:
                xyz1 = new_xyzlis[item[0]].split()[-3:]
                xyz2 = new_xyzlis[item[1]].split()[-3:]
                xyz3 = new_xyzlis[item[2]].split()[-3:]
                xyz4 = new_xyzlis[item[3]].split()[-3:]
                vec1 = -np.array(list(map(float,xyz2)))+np.array(list(map(float,xyz1))) # 2->1
                vec2 = -np.array(list(map(float,xyz2)))+np.array(list(map(float,xyz3))) # 2->3
                vec3 = -np.array(list(map(float,xyz3)))+np.array(list(map(float,xyz4))) # 3->4
            cross1 = np.cross(vec1,vec2)
            cross2 = np.cross(-vec2,vec3)
            cross3 = np.cross(vec1,vec3)
            dihe = np.sign(cross3.dot(vec2))*np.arccos(cross1.dot(cross2)/(np.linalg.norm(cross1) * \
                np.linalg.norm(cross2)))/np.pi*180
            res.append(dihe)

        else: #  input error?
            res.append(np.nan)

# information of the bond args in corresponding order
    return res

def is_IC_same(obj1,obj2):
    try:
        length = len(obj1)
    except:
        print("ERROR ADDING! Wrong obj length!")
        return
    assert len(obj1)==len(obj2), "IC compare length wrong"
    obj1 = list(obj1)
    obj2 = list(obj2)
    if length==2 or length==3:
        res = True if obj1==obj2 or obj1==obj2[::-1] else False
    elif length==4:
        res = True if obj1[1:3]==obj2[1:3] or obj1[1:3]==obj2[2:0:-1] else False
    return res

def get_geom_id(geomx,IClis):
    res = [geomx.energy]
    for item in IClis:
        res.append(geomx.IClis[item])
    return res

def redun_add(aim_set,string):
    obj = eval(string)
    for item in obj:
        for exist in aim_set: # 这里exist和item都是从1数起的内坐标，类型为tuple
            if is_IC_same(item,exist):
                break
        else:
            aim_set.add(item)
    return aim_set

class geom():
    def __init__(self,file,name=None) -> None:
        self.file = file
        self.name = name if name!=None else file.split(".")[0]
        self.xyz2graph()
        self.rotable_dict = self.find_rotable()
        # self.redun_add(args["redundantIC"])

        return
    
    def xyz2graph(self):
        [temp,self.energy,self.ele_lis] = extract(args["mol_path"]+"/"+self.file) 
        # temp is raw_xyz in g16 output; self.energy in a.u.; self.ele_lis -> list of elements
        self.xyzM = xyz2mat(temp)

        # 使用 distance_matrix 函数计算距离矩阵
        self.dist_mat = distance_matrix(self.xyzM, self.xyzM)

        # 将距离矩阵转换为稀疏矩阵
        self.sparse_mat = csr_matrix(self.dist_mat)

        self.G = nx.Graph()
        # 遍历稀疏矩阵中的每个元素，并根据阈值添加节点和边
        for i, j, d in zip(*self.sparse_mat.nonzero(), self.sparse_mat.data):
            if is_bond_exist(self.ele_lis[i],self.ele_lis[j],d):
                self.G.add_edge(i, j, weight=d)
        
        degree_dict = dict(nx.degree(self.G))
        for i in degree_dict.keys():
            if degree_dict[i]==0:
                self.fix_bond(i)

    def fix_bond(self,atom1):
        lis_conv = np.array([covR[atom1]+covR[atom2] for atom2 in range(len(self.G.nodes))])
        lis_dist = np.array(self.dist_mat[atom1])
        judge_lis = lis_dist/lis_conv
        val = 1000
        obj_atom = None
        for item in range(len(judge_lis)):
            if judge_lis[item]>0 and judge_lis[item] < val:
                val = judge_lis[item]
                obj_atom = item
        assert obj_atom not in [None,atom1], "atom fix wrong"
        self.G.add_edge(atom1,obj_atom,weight=val)
        return

    def get_neighbor(self,atom,exclude={}): # exclude:list/tuple
        try: 
            _ = len(exclude) #其实是类型检查
            pass
        except:
            exclude = {exclude}
        res = list(set(nx.neighbors(self.G, atom))-set(exclude))
        return res

    def find_rotable(self):
        rotable_dict = dict()
        edge_lis = list(self.G.edges)
        for item in edge_lis:
            [atom1,atom2] = item
            num_nei1 = len(self.get_neighbor(atom1,exclude=atom2))
            num_nei2 = len(self.get_neighbor(atom2,exclude=atom1))
            if num_nei1==0 or num_nei2==0:
                continue
            elif num_nei1<num_nei2:
                nei1 = self.get_neighbor(atom1,exclude=atom2)
                nei2 = self.get_neighbor(atom2,exclude={atom1,nei1[0]})
            else:
                nei2 = self.get_neighbor(atom2,exclude=atom1)
                nei1 = self.get_neighbor(atom1,exclude={atom2,nei2[0]})
            dihe_rep = (nei1[0]+1,atom1+1,atom2+1,nei2[0]+1)
            rotable_dict[(atom1,atom2)] = [dihe_rep,
                                            compute_internal(self.xyzM,[dihe_rep])[0]]
            # (atom1,atom2) starts from 0 (consistent with self.G); 
            # (nei1[0]+1,atom1+1,atom2+1,nei2[0]+1) starts from 1 (consistent with chem)
        return rotable_dict

    # def redun_add(self,string):
    #     obj = eval(string)
    #     for item in obj:
    #         for exist in self.IClis.keys(): # 这里exist和item都是从1数起的内坐标，类型为tuple
    #             if is_IC_same(item,exist):
    #                 break
    #         else:
    #             self.IClis[item] = compute_internal(self.xyzM,item)[0]
    #     return
    
    def generate_IClis(self,IClis):
        self.IClis = dict()
        for item in IClis:
            self.IClis[item] = compute_internal(self.xyzM,[item])[0]
        return
