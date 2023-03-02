# Automatic *Crest* Conformation Analysis

1. Run Crest + Collect All spe Data
2. Collect xyz coordinates and calculate internal coords
   1. 利用lz的数据结构，从xyz信息构建分子图结构（只考虑重原子和选定H原子）
      1. 每个原子都与N-1个原子进行判定，从而进行成键行为，成键的thresh由参数决定，统计每个原子的配位数；
      2. 检查是否还有未成键的原子（配位数=0），取"最近的原子成一根键"
      3. 检查每根键两边的原子，如果有一个原子配位数=1，则此键旋转不存在，统计需要旋转的键（M个，约为2N吧）
      4. 对于M根键，原则上都可以找到一个二面角确定其旋转方向
      5. 录入冗余内坐标，确定其是否与已有二面角同构，若无则加入
      6. 读取每个结构的内坐标信息，将其与能量一起输出到csv
3. Do a recursion using internal args, and defining their importance
   1. 看看每个二面角的方差，方差太小的舍去该特征
   2. 用LASSO模型（或者是SISSO）做一个回归分析，确定哪些指标的重要性很大，把他们拿出来
   3. 用“区分度定义”进行简单的结构聚类分析？（聚为一类的结构在此重要特征的超空间上相近）
   4. 同时也输出一个经典的二面角结构聚类结果
4. Manual Judge