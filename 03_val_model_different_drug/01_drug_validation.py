from sklearn.externals import joblib
from sklearn import preprocessing
from scipy.stats import spearmanr
import numpy as np
import pandas as pd
import math
import os


# load different models
model_list = os.listdir("..//..//02_drug_model//prot_model//")
model_list = sorted(model_list)

# load drug info data (for standardization)
all_df = pd.read_csv("..//..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
scaler = preprocessing.StandardScaler().fit(all_properties)

# load drug info per proteins
drug_file_list = os.listdir("..//..//02_drug_model//prot_model//")
drug_info_dic = {}

for drug_file in drug_file_list:
    prot_name = drug_file.replace("_model_MACCS_properties", "")
    df = pd.read_csv("..//..//02_drug_model//drug_per_prot//" + prot_name + ".txt", sep="\t", header=None, lineterminator="\n")
    MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
    properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
    Ki = np.array(list(df[df.columns[3]]))
    # data standardization
    properties = scaler.transform(properties)
    MACCS_properties = np.array(np.append(MACCS, properties, axis=1))
    x_test = MACCS_properties
    y_test = Ki
    if len(y_test) >= 10:
        drug_info_dic[prot_name] = [x_test, y_test]

counter = 1
for model in model_list:
    print(counter)
    validation_drug_prot = model.replace("_model_MACCS_properties", "")

    # use different models except validation_model
    with open(".//result//" + validation_drug_prot + "_model_pred_other_drug.txt", "w") as o1:
        print("Model_used", "Drug_Number", "Self_SCC", "Self_Z", "Val_Drug", "Val_Drug_Number", "SCC_of_Val_real", "Val_Z", "Val_Z_Bigger", sep="\t", file=o1)

        SVR_linear = joblib.load("..//..//02_drug_model//prot_model//" + model)
        self_result = list(SVR_linear.predict(drug_info_dic[validation_drug_prot][0]))
        self_scc, self_pval = spearmanr(self_result, drug_info_dic[validation_drug_prot][1])
        self_arctanh = np.arctanh(self_scc)
        Number_of_Drugs_to_Build = len(drug_info_dic[validation_drug_prot][1])
        self_z = self_arctanh*(math.pow(((Number_of_Drugs_to_Build-3)/1.06), 0.5))

        for drug in drug_info_dic:
            val_number = len(drug_info_dic[drug][1])
            pred = list(SVR_linear.predict(drug_info_dic[drug][0]))
            val_scc, val_pval = spearmanr(pred, drug_info_dic[drug][1])
            val_arctanh = np.arctanh(val_scc)
            val_z = val_arctanh*(math.pow(((val_number-3)/1.06), 0.5))

            bigger = ("yes" if val_z > self_z else "no")

            print(validation_drug_prot, Number_of_Drugs_to_Build, round(self_scc, 4), self_z, drug, val_number, round(val_scc, 4), val_z, bigger, sep="\t", file=o1)
    
    counter += 1