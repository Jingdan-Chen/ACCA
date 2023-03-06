# 从这边把类和函数导过来
from coordProcess import *

def run_coordP():
    # read in all geom and energies
        # 获取指定目录下的所有文件
    molfiles = os.listdir(args['mol_path'])
    all_geom=dict()
    scf_wrong = []
    for file in molfiles:
        [basename,suffix] = file.split(".")
        if suffix == args['mol_suffix']:
            obj = geom(file=file)
            if obj.energy is None:
                scf_wrong.append(obj.name)
                continue
            all_geom[basename] = obj # 录入所有分子的结构和能量

    # 获取此批结构的重要id信息，对每个结构生成id，并进行整合
    all_IClis = set()
    geom_lis = list(all_geom.values())
    for  i in range(len(geom_lis)):
        all_IClis = all_IClis |\
        {list(geom_lis[i].rotable_dict.values())[j][0] for j in range(len(geom_lis[i].rotable_dict))}

    all_IClis = redun_add(all_IClis,args["redundantIC"])

    for item in all_geom.values():
        item.generate_IClis(all_IClis)

    def id_to_csv(all_geom):
        temp = list(all_geom.values())[0]
        col_idx = ["spe(a.u.)"]+list(map(str,temp.IClis.keys()))
        row_idx = list(all_geom.keys())
        res_book = pd.DataFrame([get_geom_id(list(all_geom.values())[i],all_IClis) for i in range(len(all_geom))],
                    columns=col_idx,index=row_idx)
        res_book.to_csv(args["res_file"]+"/mol_id.csv")
        
    # 在当前目录下生成 mol_id.csv 文件，作为汇总结果
    id_to_csv(all_geom)

    return all_geom,scf_wrong # return allgeom to cluster_sort
