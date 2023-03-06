if __name__=='__main__':
    print("Initializing process!")
    from coordProcess import *
    from clear_res import *

    if eval(args["begin_clear"]):
        print("clear all the results!")
        clear_res()

    print("running coord reading and molecular graph generation!")
    from run_coordP import *
    [all_geom,scf_wrong] = run_coordP()

    print("begin clustering the geoms")
    from cluster_sort import *
    cluster_sort(all_geom,scf_wrong)
    print("process terminated!")
    