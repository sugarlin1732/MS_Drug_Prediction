import numpy as np
from statistics import mean
from scipy.stats import fisher_exact

# load protein name
with open("..//02_drug_model//02_protein_order.txt", "r") as file1:
    prot_name = np.array([line.replace("\n", "") for line in file1])


# load ppi network
# protein A as key
# [protein B] as value
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

                if prot_A_name == prot_B_name:
                    counter += 1
                    continue
                if prot_A_name in prot_name:
                    try:
                        ppi_dic[prot_A_name].append(prot_B_name)
                    except:
                        ppi_dic[prot_A_name] = [prot_B_name]
        counter += 1


# load protein proportionality rho corr matrix
p2p_rho = np.loadtxt("..//04_prot_corr//03_proportionality_rho//prot2prot_proportionality.txt", delimiter="\t")
rho_list = p2p_rho[:, 2]
corr_cutoff = list(set(np.round_(rho_list, 4).flatten()))
corr_cutoff.sort()


protein_pair = []
for i in p2p_rho[:, 0:2]:
    pair = [int(x) for x in i]
    protein_pair.append(prot_name[pair].tolist())
protein_pair = np.array(protein_pair)

####################################
with open("prot_pair_rho_ppi.txt", "w") as o3:
    print("proteint_A", "proteint_B", "proportionality", "PPI", sep="\t", file=o3)
    
    for i in range(len(protein_pair)):
        print(protein_pair[i, 0], protein_pair[i, 1], rho_list[i], sep="\t", end="\t", file=o3)
        if protein_pair[i, 0] in ppi_dic:
            if protein_pair[i, 1] in ppi_dic[protein_pair[i, 0]]:
                print("1", file=o3)
            else:
                print("0", file=o3)
        else:
            print("nan", file=o3)
####################################

"""
# calculate odds ratio
odds_ratio_matrix = []
cutoff_count = 1
for cutoff in corr_cutoff:
    print("cutoff ", str(cutoff_count))
    pair_ppi = int()
    pair_n_ppi = int()
    n_pair_ppi = int()
    n_pair_n_ppi = int()

    for prot_pair in protein_pair[np.where(rho_list >= cutoff)[0]]:
        if prot_pair[0] in ppi_dic:
            if prot_pair[1] in ppi_dic[prot_pair[0]]:
                pair_ppi += 1
            else:
                pair_n_ppi += 1

    for prot_pair in protein_pair[np.where(rho_list < cutoff)[0]]:
        if prot_pair[0] in ppi_dic:
            if prot_pair[1] in ppi_dic[prot_pair[0]]:
                n_pair_ppi += 1
            else:
                n_pair_n_ppi += 1

    oddsratio, pvalue = fisher_exact([[pair_ppi, pair_n_ppi], [n_pair_ppi, n_pair_n_ppi]])
    odds_ratio_matrix.append([oddsratio, pvalue])

    with open("rho_cutoff_odds_ratio//" + str(cutoff) + ".txt", "w") as o1:
        print("#pair_ppi", str(pair_ppi), sep="\t", file=o1)
        print("#pair_n_ppi", str(pair_n_ppi), sep="\t", file=o1)
        print("#n_pair_ppi", str(n_pair_ppi), sep="\t", file=o1)
        print("#n_pair_n_ppi", str(n_pair_n_ppi), sep="\t", file=o1)
        print("oddsratio", str(oddsratio), sep="\t", file=o1)
        print("pvalue", str(pvalue), sep="\t", file=o1)

    cutoff_count += 1
with open("rho_cutoff_pip_odds_ratio.txt", "w") as o2:
    print("cutoff", "odds_ratio", "p_value", sep="\t", file=o2)
    for i in range(len(corr_cutoff)):
        print(corr_cutoff[i], odds_ratio_matrix[i][0], odds_ratio_matrix[i][1], sep="\t", file=o2)
pri
####################################





####################################
nt("END")
"""