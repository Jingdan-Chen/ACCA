# Automatic *Crest* Conformation Analysis

## Workflow

1. Run Crest + Collect All spe Data(√)
2. Collect xyz coordinates and calculate internal coords （√）
3. Molecule Graph building 
   1. **scipy 生成稀疏矩阵+networkx生成分子图**
   2. 检查是否还有未成键的原子（配位数=0），取"最近的原子成一根键"
4. Read in internal coordinates
   1. 检查每根键两边的原子，如果有一个原子配位数=1，则此键旋转不存在，统计需要旋转的键（M个，约为2N吧）
   2. 对于M根键，原则上都可以找到一个二面角确定其旋转方向
   3. 录入冗余内坐标，确定其是否与已有二面角同构，若无则加入
   4. 读取每个结构的内坐标信息，将其与能量一起输出到csv
5. Do a recursion using internal args, and defining their importance
   1. 看看每个二面角的方差，方差太小的舍去该特征
   2. 用LASSO模型（或者是SISSO）做一个回归分析，确定哪些指标的重要性很大，把他们拿出来
   3. 用这些特征的id进行聚类分析（密度聚类）
6. Manual Judge

## Docs Clarification

### dirs and configs:

data/: 包含一些物理化学常数（如范德华半径），目前是从rdkit.chem直接导出的

molecules/: 需要分析的对象，应包含每个构象的单点的.gjf输入文件和输出文件

result/: 分析的结果都会放在这里

config.txt: 设置文件，超参数可调

### run_clear.py: 直接运行一次clear_res.py

​	clear_res.py: 用于清空result/文件夹

### main.py:主程序

​	run_coordP: 运行coordP.py，提取分子信息，构建分子图，输出结果到csv文件

​	coordP.py: 用于定义实现从输出文件中取出几何/能量，构建分子图的类

​	cluster_sort.py: 用于实现特征筛选和分子的聚类，输出结果到csv文件和cluster.out文件，基于分类结果创建目录归纳分子

