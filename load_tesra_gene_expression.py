import csv
import numpy as np
from sys import platform as _platform
import scipy.io

if _platform == 'darwin':
    basepath = "/Users/changyale/dataset/TESRA/"
elif _platform == 'linux2' or _platform == 'linux':
    basepath = "/home/changyale/dataset/TESRA/"
else:
    basepath = "basepath Error!"

# TESRA gene expression dataset: 3096 probes, 247 patients
filename_1 = "TESRA_Expression_for_Yale_9-2-14.csv"

# Probe names of 1439 overlappling probes between TESRA and COPD datasets
filename_2 = "sub_gene_probe_name.csv"

# Load gene expression values
file_csv = open(basepath+filename_1,"rb")
csvreader = csv.reader(file_csv,delimiter=',')
lines = [line for line in csvreader]
file_csv.close()

# probe names
probe_name = lines[0][1:len(lines[0])]

# case ids
case_id = []

# gene expression values
data_exp = []
for i in range(1,len(lines)):
    case_id.append(lines[i][0])
    data_exp.append(lines[i][1:len(lines[i])])

# store expression data in genes x patients format
data_exp = np.array(data_exp).astype(np.float).T

# save gene expression data in a .csv file
np.savetxt("data_gene_expression_tesra.csv",data_exp,delimiter=",")

# save gene expression data in a .mat file
scipy.io.savemat("data_gene_expression_tesra.mat",mdict={'data_exp_tesra':data_exp})

# save patient id
file_case_id = open("case_id_tesra.csv","wb")
wr = csv.writer(file_case_id)
for i in range(len(case_id)):
    wr.writerow([case_id[i]])
file_case_id.close()

# save probe name
file_probe_name = open("probe_name_tesra.csv","wb")
wr = csv.writer(file_probe_name)
for i in range(len(probe_name)):
    wr.writerow([probe_name[i]])
file_probe_name.close()

############################## filename_2 ##################################
probe_name_overlap = []
csvfile = open(basepath+filename_2,"rU")
csvreader = csv.reader(csvfile)
for row in csvreader:
    probe_name_overlap = probe_name_overlap + row

# extract 1439 genes from the original 3096 dataset
idx_probe_overlap = []
for i in range(len(probe_name_overlap)):
    idx_probe_overlap.append(probe_name.index(probe_name_overlap[i]))

data_tesra_overlap_with_ecl = data_exp[idx_probe_overlap,:]

# Dump this into a csv file
np.savetxt("data_tesra_overlap_with_ecl.csv",data_tesra_overlap_with_ecl,delimiter=",")

# save this dataset in a mat file
scipy.io.savemat("data_tesra_overlap_with_ecl.mat",\
        mdict={'data_tesra_overlap_with_ecl':data_tesra_overlap_with_ecl})


