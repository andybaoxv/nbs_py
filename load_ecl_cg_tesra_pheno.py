""" This script load "All_ECL_CG_TESRA_Phenotypes_11-1-14.csv"
"""

import numpy as np
import csv
import scipy.io
from sys import platform as _platform

if _platform == 'darwin':
    basepath = "/Users/changyale/dataset/ECL_CG_TESRA/"
elif _platform == 'linux2' or _platform == 'linux':
    basepath = '/home/changyale/dataset/ECL_CG_TESRA/'
else:
    basepath = 'basepath Error!'

# file containing phenotype data for patients in three cohorts: ECLIPSE,
# COPDGene and TESRA
filename_1 = "All_ECL_CG_TESRA_Phenotypes_11-1-14.csv"

# Load csv file
# case ids for patients in three cohorts
case_id = []
# phenotype feature names
feature_name = []
file_pheno = open(basepath + filename_1)
csvreader = csv.reader(file_pheno)
lines = [line for line in csvreader]
file_pheno.close()

# array of phenotype data
data_pheno = np.zeros((len(lines)-1,len(lines[0])-1))

# masked array for missing values
data_mask = np.ones(data_pheno.shape)

# first line in the file is phenotype feature names
feature_name = lines[0][1:len(lines[0])]

for i in range(1,len(lines)):
    
    # first column is case id
    case_id.append(lines[i][0])
    
    # 'NA' represents missing values
    for j in range(1,len(lines[i])):
        if lines[i][j] == 'NA':
            data_pheno[i-1,j-1] = -999
            data_mask[i-1,j-1] = 0
        elif lines[i][j] == 'ECLIPSE':
            data_pheno[i-1,j-1] = 1
        elif lines[i][j] == 'COPDGene':
            data_pheno[i-1,j-1] = 2
        elif lines[i][j] == 'TESRA':
            data_pheno[i-1,j-1] = 3
        else:
            data_pheno[i-1,j-1] = float(lines[i][j])

