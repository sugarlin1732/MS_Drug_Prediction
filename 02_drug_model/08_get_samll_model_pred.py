import numpy as np
import pandas as pd
import os
from sklearn.svm import SVR
from sklearn import preprocessing
from joblib import dump, load
from scipy.stats import spearmanr
import random
from sklearn.model_selection import train_test_split
import math
import statistics
import sys
np.set_printoptions(threshold=sys.maxsize)


# load drug number
prot_list = os.listdir(".//drug_per_prot//")


# load drug info data (for standardization)
all_df = pd.read_csv("..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
scaler = preprocessing.StandardScaler().fit(all_properties)


counter = 1
with open(".//small_model_data.txt", "w") as o1:
    print("protein", "drug_count", "100%_SCC", "100%pval", "avg_30%SCC", "median_30%_SCC", "Avg_pval_30%", "Std_of_30%_distribution", sep="\t", file=o1)

    for i in prot_list:
        print(counter)
        with open(".//drug_per_prot//" + i, "r") as file1:
            drug_count = len(file1.readlines())
            if (drug_count < 10) and (drug_count > 2):
                protein = i.replace("_HUMAN.txt", "")
                print(protein, drug_count, sep="\t", end="\t", file=o1)

                svr_all = load(".//prot_all_model//" + protein + "_HUMAN_model_MACCS_properties")
                # loading protein_drug data
                df = pd.read_csv(".//drug_per_prot//" + i, sep="\t", header=None, lineterminator="\n")
                MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
                properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
                Ki = np.array(list(df[df.columns[3]]))

                # data standardization
                properties = scaler.transform(properties)
                MACCS_properties = np.array(np.append(MACCS, properties, axis=1))
                x_test = MACCS_properties
                y_test = Ki

                # predict
                result = svr_all.predict(x_test)
                y_pred = list(result)
                result_all = spearmanr(y_test, y_pred)
                print(result_all[0], result_all[1], sep="\t", end="\t", file=o1)
                #############################################################################

                train_rho_result = list()
                train_pvalue = list()

                rho_result = list()
                pvalue = list()

                for i in range(0, 1000):
                    print(i)
                    
                    MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
                    properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))

                    Ki = np.array(list(df[df.columns[3]]))

                    # set standardization scale
                    properties = scaler.transform(properties)
                    feature = np.array(np.append(MACCS, properties, axis=1))

                    # shuffle data
                    temp = list(zip(feature, Ki))
                    random.shuffle(temp)
                    feature, Ki = zip(*temp)
                    feature = np.asarray(feature)
                    Ki = np.asarray(Ki)

                    # data split
                    x_train, x_test, y_train, y_test = train_test_split(feature, Ki, test_size=.3)

                    # setting model parameters
                    svr_linear = SVR(kernel="linear")

                    # training
                    svr_linear.fit(x_train, y_train)

                    #testing result
                    result = svr_linear.predict(x_test)
                    y_pred = list(result)

                    #training result
                    training_result = svr_linear.predict(x_train)
                    train_predict = list(training_result)

                    # record training result
                    train_rho, train_pval = spearmanr(y_train, train_predict)

                    train_rho_result.append(float(train_rho))
                    train_pvalue.append(float(train_pval))

                    # record testing result
                    rho, pval = spearmanr(y_test, y_pred)
                    rho_result.append(float(rho))
                    pvalue.append(float(pval))
                    
                train_rho_result = np.array(train_rho_result)
                train_pvalue = np.array(train_pvalue)
                rho_result = np.array(rho_result)
                pvalue = np.array(pvalue)

                print(statistics.mean(rho_result), statistics.median(rho_result), statistics.mean(pvalue), np.std(rho_result), sep="\t", file=o1)

        counter += 1
