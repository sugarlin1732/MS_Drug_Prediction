
prot_dic = {}
drug_dic = {}

with open("06_all_data_processed.txt", "r") as f1:
    for line in f1.readlines():
        data = line.split("\t")
        drug = data[0]
        prot = data[1]
        if prot in prot_dic:
            prot_dic[prot].append(drug)
        else:
            prot_dic[prot] = [drug]
        if drug in drug_dic:
            drug_dic[drug].append(prot)
        else:
            drug_dic[drug] = [prot]

with open("08_drug_number_per_prot.txt", "w") as o1:
    for key in prot_dic:
        print(key, len(prot_dic[key]), sep="\t", file=o1)

with open("08_prot_number_per_drug.txt", "w") as o2:
    for key in drug_dic:
        print(key, len(drug_dic[key]), sep="\t", file=o2)