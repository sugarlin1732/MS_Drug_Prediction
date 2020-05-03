from statistics import mean
from scipy.stats import spearmanr

prot_pair_list = []
coef_list = []
partner_jaccard = []

# jaccard
def jaccard(pa, pb):
    try:
        pa_partner = ppi_dic[pa]
        pb_partner = ppi_dic[pb]
        unions = len(set(pa_partner + pb_partner))

        intersec = 0
        for i in pa_partner:
            if i in pb_partner:
                intersec += 1
        return(1. * intersec/unions)
    except:
        return("nan")


# load feature coef
with open(".//model_to_model_coeff_scc.txt", "r") as file1:
    f1 = file1.readlines()[1:]
    for line in f1:
        data = line.replace("\n", "").split("\t")
        coef = data[-1]

        prot_pair_list.append(data[:2])
        coef_list.append(coef)

protein_lsit = list(set(j for i in prot_pair_list for j in i))


# load ppi network
# protein A as key
# [protein B] as value
ppi_dic = {}
with open("..//..//05_ppi//InBio_Map_core_2016_09_12//core.psimitab", "r") as file3: 
    print("load ppi network")       
    for line in file3:
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
                if prot_A_name in protein_lsit:
                    try:
                        ppi_dic[prot_A_name].append(prot_B_name)
                    except:
                        ppi_dic[prot_A_name] = [prot_B_name]


# count jaccard
for pair in prot_pair_list:
    jac = jaccard(pair[0], pair[1])
    partner_jaccard.append(jac)

# write
with open("scc_of_coefSCC_partner_jaccard.txt", "w") as o1:
    print("SCC between feature_coef_SCC and protein_partner_jaccard = ", spearmanr(coef_list, partner_jaccard), sep="\n", file=o1)
    print("model_1", "model_2", "SCC_of_coed", "prot_partner_jaccard", sep="\t", file=o1)
    for i in range(len(prot_pair_list)):
        print(prot_pair_list[i][0], prot_pair_list[i][1], coef_list[i], partner_jaccard[i], sep="\t", file=o1)
        
