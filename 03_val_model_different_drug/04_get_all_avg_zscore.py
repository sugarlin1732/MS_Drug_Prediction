import os
from statistics import mean
model_list = os.listdir(".//result//")

counter = 1
with open("all_mean_zscore.txt", "w") as o1:
    for model in model_list:
        print(counter)
        with open(".//result//" + model, "r") as file1:
            prot = model.replace("_HUMAN_model_pred_other_drug.txt", "")
            f1 = file1.readlines()[1:]
            zscore_list = []
            for line in f1:
                data = line.split("\t")
                if data[0] == data[4]:
                    continue
                else:
                    zscore = data[-2]
                    zscore_list.append(float(zscore))
            print(prot, mean(zscore_list), sep="\t", file=o1)
        counter += 1