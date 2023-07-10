# Automatic *Crest* Conformation Analysis

## Purpose

Grimme's crest software is a handful tool for conformation generation, but sometime it generates too many structures and many of their structure very similar. Here I propose a method to process generated "crest_conformer.xyz" to re-sampling valuable conformations in chemical space based on dihedral angle

## Workflow

1. Run Crest + Collect All spe Data(√)
2. Collect xyz coordinates and calculate internal coords （√）
3. Molecule Graph building 
   1. **scipy generates sparse matric + networkx generates molecule graph representation**
   2. check if there is any isolated atom (coordination number, CN=0)，bond it to the closest atom
4. Read in internal coordinates
   1. check the two atoms connected by a bond, if one of them owns CN=1, the bond is not rotatable, count the rotatable bonds（M~~2N）
   2. for M bonds, specify the 
   3. read in redundant internal coord, incorporate those which are not included (juding equivalence)
   4. read the dihedral angle of each geom, print to csv together with energy
5. Do a recursion using internal args, and defining their importance
   1. exclude dihedral angle with low variance (argument with low value)
   2. Use LASSO to evaluate the value of dihedral angle (those with high correlation with energy is valuable)
   3. print the sorted dihedral angles
   4. conduct density clustering to pick up those representative structures in chemical space


## Docs Clarification

### dirs and configs:

data/: physical/chemical contents(e.g. vdw radius), directly exported from rdkit.chem

molecules/: files to be analysed(.gjf/.out of Gaussian calculation containing single point energy)

result/: result of analyze

config.txt: configure file, editable

### run_clear.py: run clear_res.py

​	clear_res.py: empty result/

### main.py:

​	run_coordP: run coordP.py, extract molecular information, construct molecule graph, print to csv

​	coordP.py:define classes and functions used to read geom/energy and construct molecular graph

​	cluster_sort.py: sort character and cluster, print results to csv file and cluster.out file. Make new directories and re-allocate and move files based on clustering result 

