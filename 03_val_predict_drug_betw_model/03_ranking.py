import os


result_list = os.listdir(".//result//")



with open("ranking_scc_pair.txt", "w") as o1:
    print("Drug_of_Protein", "Drug_number", "Model_Use", "Drugs_used_to_build", "SCC", sep = "\t", file=o1)

    scc_pair = {}
    bulid_number_dic = {}
    for result in result_list:
        prot_name = result.replace("predicted_", "").replace("_drugs_w_all_model_scc.txt", "")
        with open("..//..//02_drug_model//drug_per_prot//" + prot_name + ".txt", "r") as file2:
            f2 = file2.readlines()
            build_drug_number = len(f2)
            bulid_number_dic[prot_name] = str(build_drug_number)

    for result in result_list:
        with open(".//result//" + result) as file1:
            prota = result.replace("predicted_", "").replace("_drugs_w_all_model_scc.txt", "")
            f = file1.readlines()
            drug_number = f[0].replace("\n", "").split(" ")[-1]
            for line in f[2:]:
                data = line.split("\t")
                protb = data[0]
                scc = data[1]
                scc_pair[prota + "," + protb + "," + drug_number] = [float(scc)]
    sort_scc_pair = list(sorted(scc_pair.items(), key=lambda d: d[1], reverse=True))

    for pair in sort_scc_pair:
        data = pair[0].split(",")
        print(data[0], data[2], data[1], bulid_number_dic[data[1]], pair[1][0], sep="\t", file=o1)
