import csv
import numpy as np
from sys import platform as _platform
import scipy.io

if _platform == 'darwin':
    basepath = '/Users/changyale/dataset/copd_gene_expression/'
elif _platform == 'linux2' or _platform == 'linux':
    basepath = '/home/changyale/dataset/copd_gene_expression/'
else:
    basepath = 'basepath Error!'

# COPDGene expression dataset with 12531 genes
filename_1 = "COPDGENE_PBMC_GEO_PROCESSED.txt"

# Deprecated: this file is not correct
# COPDGene expression dataset with 1487 genes(overlap with ECLIPSE)
filename_2 = "COPDGene_ECLIPSEoverlap_Probes_for_Yale.csv"

# COPDGene expression dataset with 3096 probes
filename_3 = "COPDGene_Full_Expression_For_Yale_7.22.14.csv"

# Probe names of 1439 overlapping probes between ECL and COPD datasets
filename_4 = "sub_gene_probe_name.csv"

# Probe names in ECL dataset
filename_5 = "probe_name_ecl.csv"

############################# filename_1 #####################################
# COPDGene Gene Expression Data
text_file = open(basepath+filename_1,"r")
lines = text_file.read().split('\n')

# Extract patient id
tmp = lines[0].split('\t')
case_id = tmp[1:len(tmp)]

# Extract probe_id, probe names and expression data
probe_id = []
probe_name = []
data_exp = []
for i in range(1,len(lines)-1):
    tmp = lines[i].split('\t')
    probe_id.append(float(tmp[0]))
    probe_name.append(tmp[1])
    data_exp.append(tmp[2:len(tmp)])

# store data in an array
data_exp = np.array(data_exp).astype(np.float)

# Dump expression data in a csv file
np.savetxt("data_gene_expression_copd.csv",data_exp,delimiter=",")

# save gene expression data in a .mat file
scipy.io.savemat('data_gene_expression_copd.mat',mdict={'data_exp':data_exp})

# save patient id
file_case_id = open("case_id_copd.csv","wb")
wr = csv.writer(file_case_id)
for i in range(len(case_id)):
    wr.writerow([case_id[i]])
file_case_id.close()

# save probe name
file_probe_name = open("probe_name_copd.csv","wb")
wr = csv.writer(file_probe_name)
for i in range(len(probe_name)):
    wr.writerow([probe_name[i]])
file_probe_name.close()

############################# filename_2 #####################################
# read probe names from csv file to a list
probe_name_1500 = []
csvfile = open(basepath+filename_2,"rU")
csvreader = csv.reader(csvfile)
for row in csvreader:
    probe_name_1500 = probe_name_1500 + row

# extract subset of expression data from the full set
# first find indices of probes in the original dataset
idx_probe_1500 = []
for i in range(len(probe_name_1500)):
    idx_probe_1500.append(probe_name.index(probe_name_1500[i]))
# extract subset
data_exp_1500 = data_exp[idx_probe_1500,:]

# Dump this in a csv file
np.savetxt("data_gene_expression_copd_1500.csv",data_exp_1500,delimiter=",")

# save the same data in a .mat file
scipy.io.savemat("data_gene_expression_copd_1500.mat",mdict={'data_exp_1500':data_exp_1500})

############################# filename_3 #####################################
csvfile = open(basepath+filename_3,"rU")
csvreader = csv.reader(csvfile)
lines = [line for line in csvreader]
csvfile.close()

# patient ids
case_id_3000 = lines[0][1:len(lines[0])]

# extract probe id and expression data
probe_name_3000 = []
data_exp_3000 = []
for i in range(1,len(lines)):
    probe_name_3000.append(lines[i][0])
    data_exp_3000.append(lines[i][1:len(lines[i])])

# store expression data in an array
data_exp_3000 = np.array(data_exp_3000).astype(np.float)

# reorder rows according to probe names in ECL dataset
data_exp_3000_reorder = np.zeros(data_exp_3000.shape)

# load ECL probe names
file_1 = open(filename_5,"rb")
lines = file_1.read().split('\r\n')
probe_ecl = lines[0:len(lines)-1]
file_1.close()

assert len(probe_ecl) == data_exp_3000.shape[0]
for i in range(len(probe_ecl)):
    idx = probe_name_3000.index(probe_ecl[i])
    data_exp_3000_reorder[i][:] = data_exp_3000[idx][:]

# save this data in a mat file
scipy.io.savemat("data_gene_expression_copd_3000.mat",mdict={'data_exp_3000':data_exp_3000})
scipy.io.savemat("data_gene_expression_copd_3000_reorder.mat",\
        mdict={'data_exp_3000_reorder':data_exp_3000_reorder})

# save patient id to a csv file
file_case_id = open("case_id_copd_3000.csv","wb")
wr = csv.writer(file_case_id)
for i in range(len(case_id_3000)):
    wr.writerow([case_id_3000[i]])
file_case_id.close()

# save probe names to a csv file
file_probe_name = open("probe_name_copd_3000.csv","wb")
wr = csv.writer(file_probe_name)
for i in range(len(probe_name_3000)):
    wr.writerow([probe_name_3000[i]])
file_probe_name.close()

########################## filename_4 #########################################
probe_name_overlap = []
csvfile = open(basepath+filename_4,"rU")
csvreader = csv.reader(csvfile)
for row in csvreader:
    probe_name_overlap = probe_name_overlap + row

# extract subset of expression data from the 3096 dataset, first find indices
# of probes in the 3096 COPDGene dataset
idx_probe_overlap = []
flag_overlap = True
for i in range(len(probe_name_overlap)):
    if probe_name_overlap[i] in probe_name_3000:
        idx_probe_overlap.append(probe_name_3000.index(probe_name_overlap[i]))
    else:
        flag_overlap = False
        break

assert flag_overlap == True

data_copd_overlap_with_ecl = data_exp_3000[idx_probe_overlap,:]

# Dump this in a csv file
np.savetxt("data_copd_overlap_with_ecl.csv",data_copd_overlap_with_ecl,delimiter=",")

# save this dataset in a mat file
scipy.io.savemat("data_copd_overlap_with_ecl.mat",\
        mdict={'data_copd_overlap_with_ecl':data_copd_overlap_with_ecl})


