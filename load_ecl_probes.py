import csv
import numpy as np
from sys import platform as _platform
import scipy.io

if _platform == 'darwin':
    basepath = "/Users/changyale/dataset/nbs_gene_expression/"
elif _platform == 'linux2' or _platform == 'linux':
    basepath = "/home/changyale/dataset/nbs_gene_expression/"
else:
    basepath = "basepath Error!"

# Gene expression dataset: 3096 probes, 232 patients
filename_1 = "ECL_244_Blood_05_NBS_4-19-14.csv"

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
np.savetxt("data_gene_expression_ecl.csv",data_exp,delimiter=",")

# save gene expression data in a .mat file
scipy.io.savemat("data_gene_expression_ecl.mat",mdict={'data_exp_ecl':data_exp})

# save patient id
# save patient id
file_case_id = open("case_id_ecl.csv","wb")
wr = csv.writer(file_case_id)
for i in range(len(case_id)):
    wr.writerow([case_id[i]])
file_case_id.close()

# save probe name
file_probe_name = open("probe_name_ecl.csv","wb")
wr = csv.writer(file_probe_name)
for i in range(len(probe_name)):
    wr.writerow([probe_name[i]])
file_probe_name.close()

