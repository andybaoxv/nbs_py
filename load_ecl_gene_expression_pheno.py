import csv
import numpy as np
from sys import platform as _platform
from sklearn.metrics import normalized_mutual_info_score as nmi
import matplotlib.pyplot as plt

#filename_1 = "ECL_244_Blood_05_NBS_3-23-14.csv"
#filename_2 = "ECL_NBS_ProbeGeneEntrezEnsemblMap_3-23-14.csv"
#filename_3 = "ECL_244_Blood_CellAdj_05_NBS_3-23-14.csv"
#filename_4 = "ECL_NBS_ProbeGeneMap_3-23-14.csv"

# Updates on April 23, 2014
if _platform == 'darwin':
    basepath = "/Users/changyale/dataset/nbs_gene_expression/"
elif _platform == 'linux2' or _platform == 'linux':
    basepath = "/home/changyale/dataset/nbs_gene_expression/"
else:
    basepath = "basepath Error!"

# Gene expression dataset
filename_1 = "ECL_244_Blood_05_NBS_4-19-14.csv"

# Gene ID Mapping dataset
filename_2 = "ECL_NBS_ProbeGeneEntrezEnsemblMap_4-19-14.csv"

# Phenotype dataset
filename_3 = "ECL_244_Blood_PhenoforYale.csv"

# 20 probes
filename_4 = "copd_associated_probes_age_packs_sex_cellcounts-4-19-14.csv"

## Load gene expression values
filename = basepath+filename_1
file_csv = open(filename,"rb")
csvreader = csv.reader(file_csv,delimiter=',')
lines = [line for line in csvreader]
file_csv.close()

# probe names in Gene expression dataset
probe_exp = lines[0][1:len(lines[0])]

# patient id
patient_id = []

# gene expression values
gene_exp = []
for i in range(1,len(lines)):
    patient_id.append(lines[i][0])
    gene_exp.append(lines[i][1:len(lines[i])])

gene_exp = np.array(gene_exp).astype(np.float)

## Load gene IDs mappings
file_name = basepath+filename_2
file_csv = open(file_name,"rb")
csvreader = csv.reader(file_csv,delimiter=',')
lines = [line for line in csvreader]
file_csv.close()

# Gene ID types
type_geneid = lines[0]
entrez = []
hgnc = []
probe = []
ensembl = []
for i in range(1,len(lines)):
    entrez.append(lines[i][0])
    hgnc.append(lines[i][1])
    probe.append(lines[i][2])
    ensembl.append(lines[i][3])

# Load a subset of probes
probe_subset = []
entrez_subset = []
hgnc_subset = []
ensembl_subset = []
file_name = basepath+filename_4
file_csv = open(file_name,"rb")
csvreader = csv.reader(file_csv,delimiter=',')
lines = [line for line in csvreader]
file_csv.close()
for i in range(len(lines)):
    if lines[i][0] in probe:
        probe_subset.append(lines[i][0])
        tmp = probe.index(lines[i][0])
        entrez_subset.append(entrez[tmp])
        hgnc_subset.append(hgnc[tmp])
        ensembl_subset.append(ensembl[tmp])
    else:
        print "NOT in the list: ",lines[i][0]
gene_exp_subset = np.zeros((gene_exp.shape[0],len(lines)))
n_probes_exist = 0
for j in range(len(lines)):
    if lines[j][0] in probe_exp:
        tmp = probe_exp.index(lines[j][0])
        gene_exp_subset[:,n_probes_exist] = gene_exp[:,tmp]
        n_probes_exist += 1
gene_exp_subset = gene_exp_subset[:,0:n_probes_exist]

file_csv = open("gene_exp_subset.csv","wb")
csvwriter = csv.writer(file_csv)
for i in range(gene_exp_subset.shape[0]):
    csvwriter.writerow(list(gene_exp_subset[i,:]))
file_csv.close()


# Use entrez ID as key and the list of probes as values
map_gene_probe = {}
for i in range(len(probe_exp)):
    if probe_exp[i] in probe:
        idx = probe.index(probe_exp[i])
        if ensembl[idx] != 'NA'and entrez[idx] !='NA' and hgnc[idx] != 'NA':
            if entrez[idx] not in map_gene_probe.keys():
                map_gene_probe[entrez[idx]] = [probe[idx]]
            else:
                map_gene_probe[entrez[idx]].append(probe[idx])

# Choose a single probe from a set of probes for genes that have multiple
# probes, resulting 1-to-1 mapping between entrez,probe,hgnc,ensembl
entrez_sel = []
probe_sel = []
hgnc_sel = []
ensembl_sel = []

for i in map_gene_probe.keys():
    # There's only one probe for this gene
    if len(map_gene_probe[i]) == 1:
        entrez_sel.append(i)
        probe_sel.append(map_gene_probe[i][0])
        idx = probe.index(map_gene_probe[i][0])
        hgnc_sel.append(hgnc[idx])
        ensembl_sel.append(ensembl[idx])
    # There're multiple probes for one gene
    else:
        entrez_sel.append(i)
        # Choose the probe that gives highest average expression value 
        tmp_avg = []
        for j in map_gene_probe[i]:
            tmp_avg.append(np.mean(gene_exp[:,probe_exp.index(j)]))
        tmp = map_gene_probe[i][tmp_avg.index(max(tmp_avg))]
        probe_sel.append(tmp)
        idx = probe.index(tmp)
        hgnc_sel.append(hgnc[idx])
        ensembl_sel.append(ensembl[idx])

# construct resulting gene expression data matrix
file_csv = open("eclipse_3.csv","wb")
csvwriter = csv.writer(file_csv)
csvwriter.writerow(['Probe','HGNC','Ensembl','Entrez']+patient_id)
for i in range(len(probe_sel)):
    tmp = list(gene_exp[:,probe_exp.index(probe_sel[i])])
    csvwriter.writerow([probe_sel[i]]+[hgnc_sel[i]]+[ensembl_sel[i]]+\
            [entrez_sel[i]]+tmp)
file_csv.close()

# Load Phenotype data
file_name = basepath+filename_3
file_csv = open(file_name,"rb")
csvreader = csv.reader(file_csv,delimiter=',')
lines = [line for line in csvreader]
file_csv.close()

# feature names in phenotype dataset
feature_name = lines[0][1:len(lines[0])]

# patient id from phenotype data, need to be matched with that in gene
# expression data
patient_id_pheno = []

# phenotype data
data_pheno = []

# subject_group
subject_group = []

# gold_stage
gold_stage = []

# Add label to patients according to case/control, label case as 1 and control
# as 0
patient_label = []

for i in range(1,len(lines)):
    patient_id_pheno.append(lines[i][0])
    assert lines[i][0] == patient_id[i-1]
    data_pheno.append(lines[i][1:-2])
    subject_group.append(lines[i][-2])
    if lines[i][-2] == 'COPD Subjects':
        patient_label.append(1)
    elif lines[i][-2] == 'Smoker Controls':
        patient_label.append(0)
    else:
        pass
    if lines[i][-1] == 'Stage I':
        gold_stage.append(1)
    elif lines[i][-1] == 'Stage II':
        gold_stage.append(2)
    elif lines[i][-1] == 'Stage III':
        gold_stage.append(3)
    elif lines[i][-1] == 'Stage IV':
        gold_stage.append(4)
    elif lines[i][-1] == 'NA':
        gold_stage.append(0)
    else:
        print "gold_stage error"
assert len(patient_label) == len(patient_id)

# There're missing values in phenotype data
#data_pheno = np.array(data_pheno).astype(np.float)

# Store patient label in .csv file for future use
file_csv = open("patient_label.csv","wb")
csvwriter = csv.writer(file_csv)
csvwriter.writerow(patient_label)
file_csv.close()

# Store gold_stage label in .csv file for future use
file_csv = open("gold_stage.csv","wb")
csvwriter = csv.writer(file_csv)
csvwriter.writerow(gold_stage)
file_csv.close()

# Missing values in phenotype data
feature = []
for i in range(len(data_pheno[0])):
    feature.append([])
for i in range(len(data_pheno)):
    for j in range(len(data_pheno[i])):
        feature[j].append(data_pheno[i][j])

data = np.zeros((len(data_pheno),len(data_pheno[0])))
assert len(data_pheno[0]) == len(feature)

cnt = 0
for j in range(data.shape[1]):
    tmp_sum = 0.
    tmp_num = 0.
    for i in range(len(feature[j])):
        if feature[j][i] != 'NA':
            tmp_sum += float(feature[j][i])
            tmp_num += 1.
    tmp_imputation = tmp_sum/tmp_num
    for i in range(len(feature[j])):
        if feature[j][i] != 'NA':
            data[i,j] = float(feature[j][i])
        else:
            data[i,j] = tmp_imputation
            cnt += 1
#print cnt*1./(data.shape[0]*data.shape[1])

# Store patient phenotype data in csv file
file_csv = open("data_pheno.csv","wb")
csvwriter = csv.writer(file_csv)
for i in range(data.shape[0]):
    csvwriter.writerow(list(data[i,:]))
file_csv.close()

label_pred_4 = np.loadtxt('label_pred_4.csv')
label_nbs_lf4 = np.loadtxt('label_nbs_lf4.csv')

nmi_1 = []
nmi_2 = []
nmi_3 = []
nmi_4 = []

for j in range(data.shape[1]):
    nmi_1.append(nmi(data[:,j],label_pred_4))
    nmi_2.append(nmi(data[:,j],label_nbs_lf4))
    nmi_3.append(nmi(data[:,j],patient_label))
    nmi_4.append(nmi(data[:,j],gold_stage))

index = np.arange(1,13)
bar_width = 0.2
labels = ['BD.FEV1','oxygen','ExacTrunc','BD.FEV.FVC','FracVol.950U',\
        'Lowest.15.','Emphysema','Neutrophils','Lymphocytes',\
        'Monocytes','Eosinophils','Basophils']

bar_1 = plt.bar(index,nmi_1,bar_width,color='b',label='Normalization+NMF')
bar_2 = plt.bar(index+bar_width,nmi_2,bar_width,color='r',label='NMF')
bar_3 = plt.bar(index+bar_width*2,nmi_3,bar_width,color='g',label='Case/Control')
bar_4 = plt.bar(index+bar_width*3,nmi_4,bar_width,color='y',label='Gold Stage')
plt.xlabel('Phenotype Features')
#plt.ylabel('NMI between Features and Labels')
plt.title('NMI between Phenotype Features and Labels')
plt.xticks(index+bar_width*2,labels)
plt.legend()
plt.tight_layout()
plt.show()
