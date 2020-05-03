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

# list drug data of proteins
prot_list = os.listdir(".//drug_per_prot//")

# load drug data of proteins
for protein in prot_list:
    with open(".//drug_per_prot//" + protein, "r") as file1:
        if len(file1.readlines()) < 3:            
            continue
        # make protein directory if not exist
        prot_name = protein.replace("_HUMAN.txt", "")
        save_path = ".//all_model_70_data//"+ prot_name
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        # load training data
        df = pd.read_csv(".//drug_per_prot//" + protein, sep="\t", header=None, lineterminator="\n")
        MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
        properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
        Ki = np.array(list(df[df.columns[3]]))

        # train 70% model for 1,000 times
        for loop in range(0, 1000):
            
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
            svr = SVR(kernel="linear")

            # training
            svr.fit(x_train, y_train)

            #testing result
            result = svr.predict(x_test)
            y_pred = list(result)

            #training result
            training_result = svr.predict(x_train)
            train_predict = list(training_result)

            # write data
            with open(save_path + "//" + str(loop+1) + ".txt", "w") as o1:
                print("experimental_-logKi", "predicted_-logKi", sep="\t", file=o1)
                print("training result========================================================", file=o1)
                for i in range(len(y_train)):
                    print(y_train[i], train_predict[i], sep="\t", file=o1)
                print("testing result========================================================", file=o1)
                for j in range(len(y_test)):
                    print(y_test[j], y_pred[j], sep="\t", file=o1)
        