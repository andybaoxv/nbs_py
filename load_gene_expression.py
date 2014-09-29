""" This function load data from a .csv format gene expression dataset
"""

import numpy as np
import csv

def load_gene_expression(file_name):
    """
    Parameters
    ----------
    file_name: string
        file_name of gene expression dataset
    
    Returns
    -------
    data_exp: array
        data matrix containing gene expression data
    
    case_id: case id of patients
    
    probe_name: probe names
    """
    file_csv = open(file_name,"rb");
    csvreader = csv.reader(file_csv,delimiter=",")
    lines = [line for line in csvreader]
    file_csv.close()

    # Probe names
    probe_name = lines[0][1:len(lines[0])]

    # case ids
    case_id = []

    # gene expression values
    data_exp = []
    for i in range(1,len(lines)):
        case_id.append(lines[i][0])
        data_exp.append(lines[i][1:len(lines[i])])
    data_exp = np.array(data_exp).astype(np.float).T

    return (data_exp,case_id,probe_name)


