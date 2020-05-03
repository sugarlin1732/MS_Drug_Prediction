from statistics import mean
from scipy.stats import spearmanr

prot_pair = []
prop_list = []
jaccard_list = []

# load prop
with open("prot_pair_rho_ppi.txt", "r") as file1:
    f1 = file1.readlines()[1:]
    for line in f1:
        data = line.replace("\n", "").split("\t")
        prot_pair.append(data[0:2])
        prop_list.append(float(data[2]))

protein_lsit = list(set(j for i in prot_pair for j in i))


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


# load ppi network
# protein A as key
# [protein B] as value
ppi_dic = {}
with open(".//InBio_Map_core_2016_09_12//core.psimitab", "r") as file3: 
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
for pair in prot_pair:
    jaccard_list.append(jaccard(pair[0], pair[1]))


# write file
with open("prop_ppi_partner_scc.txt", "w") as o1:
    print("scc between proportionality and ppi_partner_jaccard = ", file=o1)
    print(spearmanr(prop_list, jaccard_list), file=o1)
    print("protein A", "protein B", "proportionality", "PPI_partner_jaccard", sep="\t", file=o1)
    for i in range(len(prot_pair)):
        print(prot_pair[i][0], prot_pair[i][1], prop_list[i], jaccard_list[i], sep="\t", file=o1)