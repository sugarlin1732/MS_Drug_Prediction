from sklearn.externals import joblib
from sklearn import preprocessing
from scipy.stats import spearmanr
import numpy as np
import pandas as pd
import os


# load different models
model_list = os.listdir("..//..//02_drug_model//prot_model//")
model_list = sorted(model_list)


counter = 1
for model in model_list:
    print(counter)
    validation_drug_prot = model.replace("_model_MACCS_properties", "")

    # load drug info data (for standardization)
    all_df = pd.read_csv("..//..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
    all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
    scaler = preprocessing.StandardScaler().fit(all_properties)

    # load drugs
    df = pd.read_csv("..//..//02_drug_model//drug_per_prot//" + validation_drug_prot + ".txt", sep="\t", header=None, lineterminator="\n")
    MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
    properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
    Ki = np.array(list(df[df.columns[3]]))

    # data standardization
    properties = scaler.transform(properties)
    MACCS_properties = np.array(np.append(MACCS, properties, axis=1))

    x_test = MACCS_properties
    y_test = Ki

    # use different models except validation_model

    with open(".//result//real_" + validation_drug_prot + "_drugs_w_all_model_scc.txt", "w") as o1:
        print("drug numbers =", len(y_test), file=o1)
        print("Prot_model_compared", "SCC_of_predicted_ki", "p_value", sep="\t", file=o1)

        SVR_linear = joblib.load("..//..//02_drug_model//prot_model//" + model)
        result = SVR_linear.predict(x_test)
        validation_predict = list(result)
        self_scc, self_pval = spearmanr(Ki, validation_predict)
        print(model, self_scc, self_pval, sep="\t", file=o1)

        for other_model in model_list:
            if validation_drug_prot in other_model:
                pass
            else:
                SVR_linear = joblib.load("..//..//02_drug_model//prot_model//" + other_model)
                result = SVR_linear.predict(x_test)
                y_pred = list(result)
                scc, pval = spearmanr(Ki, y_pred)
                print(other_model.replace("_model_MACCS_properties", ""), scc, pval, sep="\t", file=o1)

    counter += 1
