""" 本文档用来处理xyz文件信息，包括读写和数据操作"""

import os
import numpy as np

class xyzread():

    """ 本对象是用来根据文件路径读取并缓存xyz信息的 """

    def __init__(self, path, flag):

        """ 对象初始化：确定文件路径和读取方法 """

        self.path = path
        self.flag = flag  # 0: crest构象搜索得到的conformer.xyz文件
        self.framedict = {}

    def get_frame(self):

        """ 根据flag选择合适的读取和存储数据的方法"""

        assert self.path, "Filepath has not been set."  # 如果对象没有合理初始化，则报错
        fileread = open(self.path, 'r')                 # 打开文件
        if self.flag == 0:                              # flag=0指数据从crest读取
            self.read_crest_xyz(fileread)
            fileread.close()                            # 关闭文件释放内存
        elif self.flag == 1:                            # flag=1指数据从progdyn的traj读取
            self.read_prog_xyz()
            fileread.close()
        else:
            fileread.close()
            raise ValueError("flag out of range")       # 没有flag不能决定接下来的方法，所以报错

    def read_crest_xyz(self, fileread, **kwargs):

        """ 本方法是get_frame方法下的方法，用于读取crest文件"""

        options = kwargs.keys()         # 确定方法各项参数，方法会反复调用而第一次调用操作与之后不同，通过参数区分
        if "frameidx" in options:       # 如果输入frameidx参数，则不是第一次调用
            if "atomnum" in options:    # 此时xyz第一行已经被读取了，atomnum不能直接读取而需要通过参数输入
                frameidx = kwargs["frameidx"]
                atomnum = kwargs["atomnum"]
            else:                       # 如果没有输入atomnum，则循环不能建立只能报错
                raise ValueError("Unable to define the total number of atoms")
        else:                           # 如果没有输入frameidx参数，则是第一次调用
            fileread.seek(0)            # 将文件指针调至文件开头
            frameidx = 0                # 设定frameidx默认值为0
            atomnum = int(fileread.readline())       # 读取文件第一行，确定atomnum的值
        fileread.readline()                 # 读取文件第二行，根据xyz格式，从第三行开始才是有用信息
        atomidx = 0
        xyzdict = {}                        # 初始化每一帧结构对应的词典
        while atomidx < int(atomnum):       # 当atomidx小于atomnum，表明原子还没有读取完
            xyzline = self.read_xyz_line(fileread.readline())   # 调用方法得到原子坐标的tuple
            xyzdict[(atomidx, xyzline[0])] = xyzline[1]         # 保存进xyzdict
            atomidx += 1                                        # 迭代参数+1
        self.framedict[frameidx] = xyzdict      # 读取完原子坐标后保存xyzdict
        while fileread.readline() != "":        # 向下读取文件，还没有读取完
            self.read_crest_xyz(fileread, frameidx=frameidx+1, atomnum=atomnum)
            # 迭代方法，frameidx+1表明读取下一帧，并传递atomnum

    def read_xyz_line(self, string):

        """ 对xyz文本进行读取，得到坐标信息以字典形式保存"""

        loclist = string.split()
        atom = loclist[0]
        location = list(map(lambda x: float(x), loclist[1:]))
        return (atom, location)

    def writeout_crest_xyz(self, frame, filepath):

        """ 选择crest文件的某一帧结构写入其它文件 """

        assert self.flag == 0, "The file is not read from a crest log."
        assert len(self.framedict) != 0, "The data is empty, no file read."
        assert os.path.exists(filepath), "The filepath does not exist."
        filename = "crest_frame_{}.xyz".format(frame)
        filewrite = open(filepath+'\\'+filename, 'a')
        datawrite = self.framedict[frame]
        passage = str(len(datawrite))+'\n'
        passage += "XYZ write from CREST frame {} atomatically\n".format(frame)
        for i in datawrite:
            line = " " + i[-1]
            location = datawrite[i]
            for j in location:
                if j >= 0:
                    line += "      "+str(j)
                else:
                    line += "     "+str(j)
            line += '\n'
            passage += line
        filewrite.write(passage)
        filewrite.close()

    def get_frame_atoms(self, frame=0):

        """ 本方法用于获取frame下的原子种类和排序 """

        if len(self.framedict) == 0:            # 必须确保xyz坐标已经被读取，否则报错
            raise ValueError("No frame read.")
        else:
            framedata = self.framedict[frame]       # 取第一个frame
            atomlist = []                       # 初始化原子列表
            for i in framedata:
                atomlist.insert(i[0], i[1])     # 按照编号写入原子名
            return atomlist

    def check_frame_atoms(self):

        """ 本方法用于检查各帧的atomlist是否正确"""

        if len(self.framedict) != 0:                        # 只能在读取信息之后该方法才能执行
            atomlist = self.get_frame_atoms(0)              # 选择其中一帧获取atomlist
            framelist = self.framedict.keys().remove(0)     # 减少计算次数，把获取atomlist的这帧除去
            for i in framelist:
                if self.get_frame_atoms(i) != atomlist:     # 如果两帧atomlist不一致，则获取信息过程出错
                    return False
            return atomlist                                 # 所有循环均未返回值表明所有帧atomlist一致，返回atomlist
        else:
            return False                                    # 没有读取信息，无法得到atomlist


    def compare_crest_frames(self, frame1, frame2):
        frameA_dict = self.framedict[frame1]
        frameB_dict = self.framedict[frame2]
        frame_diff = {}
        frameA_atoms = frameA_dict.keys()
        frameB_atoms = frameB_dict.keys()
        for i in range(len(frameA_atoms)):
            if frameA_atoms[i] != frameB_atoms[i]:
                raise ValueError("The frames in one crest conformer file has different atoms.")
            else:
                LocAtomA = np.array(frameA_dict(frameA_atoms[i]))
                LocAtomB = np.array(frameB_dict(frameB_atoms[i]))
                LocAtomDiff = LocAtomB - LocAtomA
                frame_diff[frameA_atoms[i]] = LocAtomDiff
        return frame_diff

    def align_crest_frames(self, frame1, frame2, atomlist):

        """ 本方法可以按照atomlist提供的原子对frame2的坐标按照frame1的方向进行对齐"""

        atomlist = self.check_frame_atoms()
        if atomlist == False:
            raise ValueError("Cannot get the list of atoms: no frame read or wrong frame read. ")
        if len(atomlist) == 1:          # atomlist只有一个原子，即frame2中指定原子的坐标要与frame1中的对齐
            if atomlist[0] in atomlist:   # 需要确定指定原子的序号是否相符
                pass


if __name__ == "__main__":
    p = r"D:\111 g16w_files\21 Cobalt-R complex\CoBrCF3Hx2\crest_conformers.xyz"
    a = xyzread(p, 0)
    a.get_frame()
    a.writeout_crest_xyz(1, r"D:\111 g16w_files\21 Cobalt-R complex\CoBrCF3Hx2")
