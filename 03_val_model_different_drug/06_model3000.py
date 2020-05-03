model_dic = {}


with open("all_result.txt", "r") as file1:
    for line in file1.readlines()[1:]:
        data = line.split("\t")
        model, drug_number, val_z = data[0], data[1], data[7]
        try:
            model_dic[model + "," + drug_number].append(float(val_z))
        except:
            model_dic[model + "," + drug_number] = [float(val_z)]

for key in model_dic:
    model_dic[key] = (sum(1 for i in model_dic[key] if i > 5))/len(model_dic[key])


with open("model_drug_number_N_z3percent.txt", "w") as o1:
    print("percent of predicting other drugs zscore > 3", file=o1)
    print("model_used", "drug_number", "percent_zscore>3", sep="\t", file=o1)
    for key in model_dic:
        print(key.split(",")[0], (key.split(",")[1]), model_dic[key], sep="\t", file=o1)