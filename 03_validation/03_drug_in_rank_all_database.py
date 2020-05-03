import numpy as np
from scipy.stats import fisher_exact

model_prot = "AA2AR_HUMAN"

# load drug id
with open("..//01_raw_data//05_all_drug_info.txt", "r") as file1:
    drug_id = [line.split("\t")[0] for line in file1]

# load predicted all data ki
with open(model_prot + "_model_all_database_drug.txt", "r") as file2:
    drug_ki = [line.replace("\n", "") for line in file2.readlines()[1:]]

drug_dic = {}
for i in range(len(drug_id)):
    drug_dic[int(drug_id[i])] = float(drug_ki[i])

sort_list = list(map(list, sorted(drug_dic.items(), key=lambda d: d[1], reverse=True)))

# load protein drug id
with open("..//02_drug_model//drug_per_prot//" + model_prot + ".txt", "r") as file3:
    prot_drug_id = [int(line.split("\t")[0]) for line in file3]


# computing ranking
count_list = [0]*10
lowest = sort_list[-1][1]
highest = sort_list[0][1]
range_list = []
for i in range(0, 11):
    range_list.append(highest - ((highest-lowest)/10)*i)

print(highest, lowest)
print(range_list)

"""
for j in range(len(sort_list)):
    if sort_list[j][0] in prot_drug_id:
        # count_list[int((highest-sort_list[j][1])/((highest-lowest)/10))-1] += 1

print(count_list)
"""
"""
# computing ranking
count_list = [0]*10
for i in range(len(sort_list)):
    if sort_list[i][0] in prot_drug_id:
        percent = int(((i+1)/len(sort_list))*10)
        count_list[percent] += 1


# compute fisher exact test
total_drug_number = len(drug_id)
total_model_prot_number = len(prot_drug_id)


ten_percent = total_drug_number // 10
ninety_percent = total_drug_number - ten_percent


with open("distribution_" + model_prot + "_in_all_" + model_prot + "_model.txt", "w") as o1:
    print("Model_drug_number", "odds_ratio", "p_value", sep="\t", file=o1)
    for i in count_list:
        oddsratio, pvalue = fisher_exact([[i, ten_percent-i], [total_model_prot_number - i, ninety_percent - total_model_prot_number + i]], alternative="greater")
        print(i, oddsratio, pvalue, sep="\t", file=o1)
"""