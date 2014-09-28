import numpy as np
import csv

# Extract patient ID from ECLIPSE dataset
file_eclipse = open("/Users/changyale/dataset/eclipse/"+\
        "ECLIPSE_121_Blood_Raw_Top_ForNBS_1.25.14.csv","rb")
reader = csv.reader(file_eclipse)
lines = [line for line in reader]
file_eclipse.close()
patient_id = lines[0][4:len(lines[0])]

method = "NBS_with_CC"

# Read labels from csv file
file_labels = open("/Users/changyale/matlab/nbs_release_v0.2_ECLIPSE/"+\
        method+"_labels.csv","rb")
reader = csv.reader(file_labels)
lines = [line for line in reader]
patient_label = []
for i in range(len(lines)):
    patient_label.append(lines[i][0])

assert len(patient_id) == len(patient_label)

# Write patient id and label together to csv file
file_id_label = open("/Users/changyale/matlab/nbs_release_v0.2_ECLIPSE/"+\
        method+"_ID_labels.csv","wb")
csvwriter = csv.writer(file_id_label)
for i in range(len(patient_id)):
    csvwriter.writerow([patient_id[i],patient_label[i]])
file_id_label.close
