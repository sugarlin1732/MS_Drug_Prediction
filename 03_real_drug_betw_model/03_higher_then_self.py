import os

result_list = os.listdir(".//result//")


with open("self_predict_lower_other_models.txt", "w") as o1:
    print("Drug_of_Protein", "Drug_number", "Self_pred_scc", "Model_Use", "SCC", sep="\t", file=o1)
    for result in result_list:
        prota = result.replace("real_", "").replace("_drugs_w_all_model_scc.txt", "")
        with open(".//result//" + result) as file1:
            f = file1.readlines()
            durg_number = f[0].replace("\n", "").split(" ")[-1]
            pascc = float(f[2].split("\t")[1])

            for line in f[3:]:
                data = line.split("\t")
                scc = float(data[1])
                if scc > pascc:
                    protb = data[0]
                    print(prota, durg_number, pascc, protb, scc,  sep="\t", file=o1)
