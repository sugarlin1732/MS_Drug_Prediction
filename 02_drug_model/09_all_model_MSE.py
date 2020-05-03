import sys
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
from scipy.stats import spearmanr
from sklearn import preprocessing
from joblib import dump, load
from sklearn.metrics import mean_squared_error
import random
import os
import numpy as np
import pandas as pd
import math
import statistics
np.set_printoptions(threshold=sys.maxsize)


# load drug info data (for standardization)
all_df = pd.read_csv("..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
scaler = preprocessing.StandardScaler().fit(all_properties)

with open(".//all_model_MSE.txt", "w") as o2:
    print("model", "drug_number", "average_70%_MSE", "average_30%_MSE", "100%_MSE", sep="\t", file=o2)
    model_list = os.listdir(".//prot_all_model//")

    for model in model_list:
        prot = model.replace("_model_MACCS_properties", "")
        # loading data
        df = pd.read_csv(".//drug_per_prot//" + prot + ".txt", sep="\t", header=None, lineterminator="\n")
        MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
        properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
        Ki = np.array(list(df[df.columns[3]]))

        if len(Ki) < 2:
            continue

        # set standardization scale
        properties = scaler.transform(properties)
        feature = np.array(np.append(MACCS, properties, axis=1))


        train_rho_result = list()
        train_pvalue = list()

        rho_result = list()
        pvalue = list()

        mean_train_MSE_list = list()
        mean_test_MSE_list = list()


        for i in range(0, 1000):
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
            mean_train_MSE_list.append(mean_squared_error(y_train, train_predict))

            # record testing result
            rho, pval = spearmanr(y_test, y_pred)

            rho_result.append(float(rho))
            pvalue.append(float(pval))
            mean_test_MSE_list.append(mean_squared_error(y_test, y_pred))

        with open(".//all_model_MSE//" + prot + "_SVR_MSE_0.7_0.3_1000loop.txt", "w") as o1:

            print("min/max training spearman rho results", min(train_rho_result), max(train_rho_result), file=o1)
            print("min/max training spearman p value", min(train_pvalue), max(train_pvalue), file=o1)
            print("average training spearman results", statistics.mean(train_rho_result), file=o1)
            print("mdeian training spearman results", statistics.median(train_rho_result), file=o1)
            print("average training MSE", statistics.mean(mean_train_MSE_list), file=o1)


            print("SCC", "MSE", sep="\t", file=o1)
            for i in range(len(train_rho_result)):
                print(train_rho_result[i], mean_train_MSE_list[i], file=o1)

            print("min/max testing spearman rho results", min(rho_result), max(rho_result), file=o1)
            print("min/max testing spearman p value", min(pvalue), max(pvalue), file=o1)
            print("average testing spearman results", statistics.mean(rho_result), file=o1)
            print("mdeian testing spearman results", statistics.median(rho_result), file=o1)
            print("average testing MSE", statistics.mean(mean_test_MSE_list), file=o1)



            print("SCC", "MSE", sep="\t", file=o1)
            for j in range(len(rho_result)):
                print(rho_result[j], mean_test_MSE_list[j], file=o1)
        
        
        # load 100% model
        svr_all = load(".//prot_all_model//" + model)
        result_all = list(svr_all.predict(feature))
        MSE_all = mean_squared_error(Ki, result_all)
        print(prot.replace("_HUMAN", ""), str(len(Ki)), statistics.mean(mean_train_MSE_list), statistics.mean(mean_test_MSE_list), MSE_all, sep="\t", file=o2)
