import numpy as np
from statistics import mean
from scipy.stats import fisher_exact

# load protein name
with open("..//02_drug_model//02_protein_order.txt", "r") as file1:
    prot_name = np.array([line.replace("\n", "") for line in file1])

# load ppi network
# protein A as key
# [protein B, score] as value
ppi_dic = {}
counter = 1
with open(".//InBio_Map_core_2016_09_12//core.psimitab", "r") as file3: 
    print("load ppi network")       
    for line in file3:
        print(counter)
        data = line.replace("\n", "").split("\t")
        prot_A_name = data[2].split("|")[0].replace("uniprotkb:", "")
        prot_B_name = data[3].split("|")[0].replace("uniprotkb:", "")
        method = data[6]
        score = max([float(i) for i in data[-2].replace("-", "").split("|") if i.replace(".", "").isdigit()])
        
        # ppi score threshold
        if "experimental interaction detection" in method:
            if score >= 0:
                value = [prot_B_name, score]

                if prot_A_name == prot_B_name:
                    counter += 1
                    continue
                if prot_A_name in prot_name:
                    try:
                        ppi_dic[prot_A_name].append(value)
                    except:
                        ppi_dic[prot_A_name] = [value]
        counter += 1

# load protein proportionality rho corr matrix
rho_matrix = np.loadtxt("..//04_prot_corr//03_proportionality_rho//prot_rho_modified_matrix.txt", delimiter=",")
corr_cutoff = list(set(np.round_(rho_matrix, 4).flatten()))


# calculate odds ratio
odds_ratio_matrix = []

print("calculating")

calculate_count = 1
for cutoff in corr_cutoff[1:]:
    print(calculate_count)
    temp_odds_ratio_list = []

    with open(".//rho_cutoff_odds_ratio//" + str(cutoff) + "_rho_cutoff_odds_ratio.txt", "w") as o1:

        for protein_A in range(len(rho_matrix)):
            member_ppi = int()
            member_n_ppi = int()
            n_member_ppi = int()
            n_member_n_ppi = int()

            member_name = prot_name[np.where(rho_matrix[protein_A] <= cutoff)[0]]
            member_name = np.delete(member_name, np.argwhere(member_name == prot_name[protein_A]))
            for protein_B in member_name:
                if protein_B in ppi_dic:
                    member_ppi += 1
                else:
                    member_n_ppi += 1

            n_member_name = prot_name[np.where(rho_matrix[protein_A] > cutoff)[0]]
            for protein_B in n_member_name:
                if protein_B in ppi_dic:
                    n_member_ppi += 1
                else:
                    n_member_n_ppi += 1
            oddsratio, pvalue = fisher_exact([[member_ppi, member_n_ppi], [n_member_ppi, n_member_n_ppi]])
            if (str(oddsratio) != "inf") and (str(oddsratio) != "nan"):
                temp_odds_ratio_list.append(oddsratio)
                print(oddsratio, pvalue, file=o1)
            
    
    try:
        odds_ratio_matrix.append(mean(temp_odds_ratio_list))
    except:
        odds_ratio_matrix.append("nan")
    calculate_count += 1


with open("rho_cutoff_pip_odds_ratio.txt", "w") as o2:
    print("cutoff", "average_odds_ratio", sep="\t", file=o2)
    for i in range(len(corr_cutoff)):
        print(corr_cutoff[i], odds_ratio_matrix[i], sep="\t", file=o2)
print("END")
