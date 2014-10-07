""" This script load "COMBAT_allECL_CG_TESRA_Expression_9-26-14.csv"
"""

import numpy as np
import csv
import scipy.io
from sys import platform as _platform
from load_gene_expression import load_gene_expression

if _platform == 'darwin':
    basepath = "/Users/changyale/dataset/ECL_CG_TESRA/"
elif _platform == 'linux2' or _platform == 'linux':
    basepath = "/home/changyale/dataset/ECL_CG_TESRA/"
else:
    basepath = "basepath Error!"

# file containing ECL, CG, TESRA
filename_1 = "COMBAT_allECL_CG_TESRA_Expression_9-26-14.csv"

# file containing 1439 probe names, which are the overlap between ECL and
# STRING network
filename_2 = "sub_gene_probe_name.csv"

# load gene expression dataset
tmp = load_gene_expression(basepath+filename_1)

# extract ECL dataset
data_ecl = tmp[0][:,0:232]
case_id_ecl = tmp[1][0:232]

# extract COPDGene dataset
data_cg = tmp[0][:,232:368]
case_id_cg = tmp[1][232:368]

# extract tesra dataset
data_tesra = tmp[0][:,368:615]
case_id_tesra = tmp[1][368:615]

# extract probe names
probe_name_3096 = tmp[2]
data_ecl_cg_tesra = tmp[0]
case_id_ecl_cg_tesra = tmp[1]

# save dataset, case_ids and probe names into .mat format
scipy.io.savemat(basepath+"data_ecl_cg_tesra.mat",\
        mdict={'data_ecl_cg_tesra':data_ecl_cg_tesra})

# save case id
file_case_id = open(basepath+"case_id_ecl_cg_tesra.csv","wb")
wr = csv.writer(file_case_id)
for i in range(len(case_id_ecl_cg_tesra)):
    wr.writerow([case_id_ecl_cg_tesra[i]])
file_case_id.close()

# save probe name
file_probe_name = open(basepath+"probe_name_3096.csv","wb")
wr = csv.writer(file_probe_name)
for i in range(len(probe_name_3096)):
    wr.writerow([probe_name_3096[i]])
file_probe_name.close()

# Extract data corresponding to 1439 probes
probe_name_1439 = []
csvfile = open(basepath+filename_2,"rU")
csvreader = csv.reader(csvfile)
for row in csvreader:
    probe_name_1439 = probe_name_1439 + row
csvfile.close()

idx_probe_1439 = []
for i in range(len(probe_name_1439)):
    idx_probe_1439.append(probe_name_3096.index(probe_name_1439[i]))
data_ecl_cg_tesra_1439 = data_ecl_cg_tesra[idx_probe_1439,:]

# save this dataset in a .mat file
scipy.io.savemat(basepath+"data_ecl_cg_tesra_1439.mat",\
        mdict={"data_ecl_cg_tesra_1439":data_ecl_cg_tesra_1439})


