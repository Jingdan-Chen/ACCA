# import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
from scipy.stats import f_oneway
from sklearn.linear_model import Lasso
from coordProcess import *
from sklearn.cluster import DBSCAN
import shutil

def cluster_sort(all_geom,scf_wrong):
    worksheet = pd.read_csv(args["res_file"]+"/mol_id.csv")

    name_l = worksheet.iloc[:,0]
    y = np.array(worksheet.iloc[:,1]).reshape(-1)
    X = np.array(worksheet.iloc[:,2:])
    X_shape_ori = X.shape
    norm_y = (y-np.min(y))/(np.max(y)-np.min(y)) # 归一化
    X_var = np.var(np.abs(X),axis=0) #方差计算

    # 能量筛选文件
    energy_window = float(args["spe_window"])
    ener_out_geom = list(name_l[(y-np.min(y))*627.5>energy_window])

    # 二面角方差筛选
    threshold = float(args["dihe_var_thre"])
    f_values, p_values = f_oneway(*X.T)
    selected_features = np.where(X_var > threshold)[0] 
    X_var_pick = X[:,selected_features]


    use_lasso = eval(args['use_lasso'])
    if use_lasso:
    # LASSO筛选特征
        alpha = float(args["lasso_alpha"])  # L1正则化强度参数,可能需要变化
        lasso = Lasso(alpha=alpha, 
                    max_iter=int(args["lasso_max_iter"]))
        lasso.fit(np.cos(X_var_pick), norm_y)
        feature_importance = np.abs(lasso.coef_)
        res_feature = selected_features[feature_importance>0]
    else:
        res_feature = selected_features

    new_sheet = worksheet.iloc[:,[0,1]+list(res_feature+2)]
    new_sheet.to_csv(args["res_file"]+"/selected_molid.csv")


    # 密度聚类分析
    X_var_pick = X[:,res_feature]
    cluster_X = (X_var_pick-np.min(X_var_pick,axis=0))/(
        np.max(X_var_pick,axis=0)-np.min(X_var_pick,axis=0))

    # 创建一个DBSCAN聚类模型
    model = DBSCAN(eps=float(args["dbscan_eps"]),
                    min_samples=1) # eps参数可能需要调整

    # 拟合模型并预测簇标签
    labels = model.fit_predict(cluster_X)

    # 将结构进行归类标记
    clustered_dict = dict()
    for i in range(len(labels)):
        if labels[i] not in clustered_dict:
            clustered_dict[labels[i]] = list()
        clustered_dict[labels[i]].append(name_l[i])

    # 输出与lasso分类相关的信息，并依据分类结果创建文件夹
    out_info = ["Input Geom: {}; Original IC: {}; Selected IC: {}\n".format(
        X_shape_ori[0],X_shape_ori[1],len(res_feature)),"\n"]
    
    os.mkdir(args["res_file"]+"/best")
    for i in clustered_dict.keys():
        out_info.append("Type {}:\n".format(i))
        temp_glis = sorted(clustered_dict[i],key=lambda a:float(all_geom[a].energy))
        out_info.append(";".join(temp_glis)+"\n")
        out_info.append("\n")
        folder_name = args["res_file"]+"/type{}".format(i)
        os.mkdir(folder_name)
        for file in temp_glis:
            source_file = args["mol_path"]+"/"+file+".gjf"
            destination_folder = folder_name
            shutil.copy(source_file, destination_folder)

        # 将 "source_file.txt" 复制到 "destination_folder" 中
        file_best = temp_glis[0]
        if file_best not in ener_out_geom:
            source_file = args["mol_path"]+"/"+file_best+".gjf"
            destination_folder = args["res_file"]+"/best"
            shutil.copy(source_file, destination_folder)
        

    out_info.append("\n")

    # 输出energy out的geometry
    out_info.append("Energy Out Files:\n")
    out_info.append(";".join(ener_out_geom)+"\n")
    out_info.append("\n")

    # scf_wrong 信息
    out_info.append("SCF wrong files: "+";".join(scf_wrong)+"\n")
    # 输出结束语
    out_info.append("Thank you for using ACCA!\n")

    with open(args["res_file"]+"/cluster.out","w") as f:
        f.writelines(out_info)

