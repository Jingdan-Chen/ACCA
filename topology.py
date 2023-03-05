# -*- coding=utf-8 -*-
import numpy as np

# Cited from Dalton Trans., 2008, 2832.
Covalence_Dict = {'H': 0.31, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.73, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Na': 1.66,
                  'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Mn': 1.39, 'Co': 1.26,
                  'Ni': 1.24, 'Cu': 1.32, 'Br': 1.20}

def is_bond_exist(atom1, atom2, length,thresh=1.25):

    """ 给定长度和两端原子种类，判断该距离下两原子是否成键 """

    judge = length/(Covalence_Dict[atom1] + Covalence_Dict[atom2])
    if judge > thresh:
        return False
    else:
        return True

def find_small_index(atom1, atom2):

    if atom1.get_index() > atom2.get_index():
        return [atom2, atom1]
    elif atom1.get_index() < atom2.get_index():
        return [atom1, atom2]
    else:
        raise ValueError("Two same atoms")

class Atom():

    def __init__(self, idx, name, loc):

        self.index = idx
        self.name = name
        self.loc = loc

    def get_loc(self):

        return self.loc

    def get_index(self):

        return self.index

    def get_name(self):

        return self.name

    def set_name(self, newname):

        self.name = newname
        return None

    def __repr__(self):

        return "Atom idx {} name {}".format(self.index, self.name)

class Bond():

    """ 表示两原子间的距离以及成键状况 """

    def __init__(self, atomA: Atom, atomB: Atom):

        """ 初始化对象，确定对象两端点和距离值 """

        self.terminal = find_small_index(atomA, atomB)
        self.length = np.linalg.norm(atomA.get_loc() - atomB.get_loc())

    def get_length(self):

        return self.length

    def get_atom(self):

        return self.terminal

    def another(self, atom):

        for i in self.terminal:
            if i.get_name() == atom:
                getidx = self.terminal.index(i)
        return self.terminal[1-getidx]

    def __repr__(self):

        return "Bond terminal: {} & {}; length: {} A.\n".format(self.terminal[0], self.terminal[1], self.length)

class Xyzlocation():

    """ 本对象用于存储一个分子结构的xyz信息 """

    def __init__(self, string):

        """ 初始化对象：从一段xyz格式文件中读取分子结构（针对不同应用场景建议重载该函数） """

        self.atomlist = []              # 存储xyz坐标信息：（原子序号，原子名）：xyz坐标
        self.bond_matrix = []           # 存储原子间距离信息：初次调用get_bond_matrix方法即可
        self.set_topo = False
        location_list = string.split('\n')          # xyz文件中每行一个原子坐标，因此根据分行符分割出每个原子的坐标
        for i in location_list:
            if i != "":
                atom_index = location_list.index(i)     # 第i行为第i个原子的坐标，从0计数
                atom_location = i.split()               # 分割每行：第一个元素为原子名，后面为坐标信息
                atom_name = atom_location[0]
                atom_coord_0 = atom_location[1:4]       # 初始坐标信息，元素类型为str
                atom_coord = np.around(list(map(lambda x: float(x), atom_coord_0)), decimals=8)
                # 将列表转换为numpy数组类型，元素为float类型，保留小数点后8位
                atom = Atom(idx=atom_index, name=atom_name, loc=atom_coord)
                self.atomlist.append(atom)              # 存储atom对象

    def __getitem__(self, index):

        """ 重写序号下标方法 """

        return self.atomlist[index]

    def get_bond_matrix(self):

        """ 获取xyz文件中的所有原子之间的距离 """

        for i in self.atomlist:             # 对所有原子
            atombind = []
            for j in self.atomlist:         # 对单个原子，检索所有原子
                if i.get_index() < j.get_index():       # 根据编号确保所有距离只计算一遍
                    bond = Bond(i, j)                   # 建立Bond对象
                    atombind.append(bond)               # 建立根据原子区分的子列表，方便后续查找
            if len(atombind) != 0:
                self.bond_matrix.append(atombind)           # 存储所有子列表

    def find_bond(self, atomA, atomB):

        """ 查找两个原子之间的Bond对象 """

        handle = find_small_index(atomA, atomB)                 # 对输入的原子进行排序
        atom_part = self.atomlist.index(handle[0])              # 选择较小序号的原子
        bond_part = self.bond_matrix[atom_part]                 # 直接找到对应bondmatrix的子列表
        search = 0                                              # 搜索参数，如果没有搜索到则报错
        for i in bond_part:
            if i.get_atom()[1] == handle[1]:                    # 匹配条件
                return i
            else:
                search += 1
        if search == len(bond_part):
            raise ValueError("cannot find the bond")

    def translate_xyz(self, vector):
        for i in self.locdict:
            self.locdict[i] += vector

class Topology():

    """ 该对象表示分子拓扑结构 """

    def __init__(self, xyzdata: Xyzlocation):

        """ 初始化对象：self.atom: 体系中的原子列表，相当于图论中的顶
                      self.bond: 体系中的成键列表，相当于图论中的边
                      atomlist: 输入的原子列表
                      bond_matrix: 输入的原子距离矩阵"""

        self.atom = xyzdata.atomlist
        self.bond_matrix = xyzdata.bond_matrix
        self.bond = []
        self.path = []

    def set_bond(self):

        """ 从bondmatrix中找出成键的距离，将对应的原子对保存下来 """

        for idx in range(len(self.atom)-1):
            i = self.atom[idx]
            bondlist = self.bond_matrix[idx]
            for j in bondlist:
                another = j.another(i.get_name()).get_name()
                if is_bond_exist(i.get_name(), another, j.length):
                    self.bond.append(j.get_atom())

    def get_carbon_topo(self, C_bond, found):

        """ 判断碳原子周围的成键情况"""

        carbonsearch = 0
        for i in [1, 2, 3]:
            C_category = 'C-{}'.format(i)
            C_link = 0
            C_topo = []
            printidx = 1
            for j in C_bond:
                if j.another('C').get_name() == 'C':
                    another = 'C-1'
                else:
                    another = j.another('C').get_name()
                if is_bond_exist(C_category, another, j.length):
                    C_link += 1
                    printidx += 1
                    C_topo.append(j.get_atom())
                else:
                    printidx += 1
            C_link += found
            if C_link == 5-i:
                return (C_topo, i)
            else:
                carbonsearch += 1
        if carbonsearch == 3:
            raise ValueError("cannot determine the hybrid type of carbon")

    def find_topo(self, atom):

        found = []
        for i in self.bond:
            if atom in i:
                found.append(i)
        return found

    def get_link_atoms(self, atom):

        found = self.find_topo(atom)
        link = []
        for i in found:
            idx = i.index(atom)
            link.append(i[1-idx])
        return link

    def getCH3(self):

        CH3list = []
        for i in self.atom:
            if i.get_name() == "C":
                link = self.get_link_atoms(i)
                Hnum = 0
                for j in link:
                    if j.get_name() == 'H':
                        Hnum += 1
                if Hnum == 3:
                    CH3list.append(i)
        return CH3list

    def stepstep(self, previous, now, step):

        step.append(now)
        link = self.get_link_atoms(now)
        link.remove(previous)
        if len(link) == 0:
            print(self.path)
            print(step)
            self.path.append(step.copy())
            print(self.path)
            step.remove(now)
        else:
            for i in link:
                if i in step:
                    self.path.append(step.copy())
                else:
                    print(self.path)
                    self.stepstep(now, i, step)

    def stepstart(self, start):

        link = self.get_link_atoms(self.atom[start])
        for j in link:
            step = [self.atom[start]]
            self.stepstep(self.atom[start], j, step)
        return self.path

    def get_ring(self, start):

        totalring = []
        def stepring(previous, now, step):
            step.append(now)
            link = self.get_link_atoms(now)
            link.remove(previous)
            link1 = []
            for i in link:
                if i.get_name() != 'H':
                    link1.append(i)
            if len(link1) == 0:
                #print("get terminal")
                #print(step)
                step.remove(now)
            else:
                for i in link1:
                    if i in step:
                        #print("get ring")
                        #print(step)
                        ringstart = step.index(i)
                        ring = set(step[ringstart:])
                        if ring in totalring:
                            continue
                        else:
                            totalring.append(ring)
                    else:
                        #print("further go")
                        #print(step)
                        stepring(now, i, step)
                step.remove(now)
        link = self.get_link_atoms(self.atom[start])
        link1 = []
        for j in link:
            if j.get_name() != 'H':
                link1.append(j)
        for j in link:
            step = [self.atom[start]]
            stepring(self.atom[start], j, step)
        return totalring

    def get_ring_atom(self):

        ringatom = []
        for i in self.get_ring(0):
            for j in i:
                if j not in ringatom:
                    ringatom.append(j)
        print(ringatom)
        return ringatom

    def get_dihedral(self, ring):
        totaldihedral = []
        for i in self.atom:
            def step4(previous, now, step):
                if len(step) < 3:
                    step.append(now)
                    link = self.get_link_atoms(now)
                    link.remove(previous)
                    link1 = []
                    for i in link:
                        if i.get_name() != 'H':
                            link1.append(i)
                    if len(link1) == 0:
                        #print("get terminal")
                        #print(step)
                        step.remove(now)
                    else:
                        for i in link1:
                            step4(now, i, step)
                        step.remove(now)
                elif len(step) == 3:
                    step.append(now)
                    totaldihedral.append(step.copy())
                    step.remove(now)
            link = self.get_link_atoms(i)
            link1 = []
            for j in link:
                if j.get_name() != 'H':
                    link1.append(j)
            for j in link:
                step = [i]
                step4(i, j, step)
        final = []
        for k in totaldihedral:
            num = 0
            for l in k:
                if l in ring:
                    num += 1
            if num < 3:
                final.append(k)
        return final
            

if __name__ == "__main__":
    filepath = r"D:\111 g16w_files\21 Cobalt-R complex\CoBrCF3Hx2\CoBrCF3Hx2.xyz"
    with open(filepath, 'r') as f:
        s0 = f.readlines()
        s1 = s0[2:]
        s2 = ""
        for i in s1:
            s2 += i
        xyz = Xyzlocation(s2)
        xyz.get_bond_matrix()
        topo = Topology(xyz)
        topo.set_bond()
        ring = topo.get_ring_atom()
        a = topo.get_dihedral(ring)
        print(len(a))
        op = ""
        for i in a:
            for j in i:
                op += str(j.get_index()) + ' -> '
            op += '\n'
        print(op)

