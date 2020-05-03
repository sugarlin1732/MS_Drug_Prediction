from statistics import mean
from scipy.stats import spearmanr

protein_ebc_dic = {}
with open("prot_pair_data_of_betweeness_network.txt", "r") as file1:
    for line in file1.readlines()[1:]:
        data = line.split("\t")
        protein_ebc_dic[data[0] + "," + data[1]] = data[3]



# protein_pair = []
# ebc_lsit = []
# PPIpartner_jaccard = []
# with open("prop_ppi_partner_scc.txt", "r") as file2:
#     for line in file2.readlines()[3:]:
#         data = line.replace("\n", "").split("\t")
#         if (data[0] + "," + data[1]) in protein_ebc_dic:
#             protein_pair.append(data[0] + "," + data[1])
#             ebc_lsit.append(float(protein_ebc_dic[data[0] + "," + data[1]]))
#             PPIpartner_jaccard.append(data[3])
#         elif (data[1] + "," + data[0]) in protein_ebc_dic:
#             protein_pair.append(data[0] + "," + data[1])
#             ebc_lsit.append(float(protein_ebc_dic[data[1] + "," + data[0]]))
#             PPIpartner_jaccard.append(float(data[3]))
            

# with open("prot_ebc_PPIpartner_SCC.txt", "w") as o1:
#     print("scc between eBC and ppi_partner_jaccard =", file=o1)
#     print(spearmanr(ebc_lsit, PPIpartner_jaccard), file=o1)
#     print("protein A", "proteinB", "eBC", "PPI_partner_jaccard", file=o1)
#     for i in range(len(protein_pair)):
#         pa, pb = protein_pair[i].split(",")
#         print(pa, pb, ebc_lsit[i], PPIpartner_jaccard[i], sep="\t", file=o1)

protein_pair = []
ebc_lsit = []
PPIpartner_jaccard = []
with open("prop_ppi_partner_scc.txt", "r") as file2:
    for line in file2.readlines()[3:]:
        data = line.replace("\n", "").split("\t")
        if float(data[-1]) > 0:
            if (data[0] + "," + data[1]) in protein_ebc_dic:
                protein_pair.append(data[0] + "," + data[1])
                ebc_lsit.append(float(protein_ebc_dic[data[0] + "," + data[1]]))
                PPIpartner_jaccard.append(data[3])
            elif (data[1] + "," + data[0]) in protein_ebc_dic:
                protein_pair.append(data[0] + "," + data[1])
                ebc_lsit.append(float(protein_ebc_dic[data[1] + "," + data[0]]))
                PPIpartner_jaccard.append(float(data[3]))


with open("prot_ebc_PPIpartner_SCC_1.txt", "w") as o1:
    print("scc between eBC and ppi_partner_jaccard =", file=o1)
    print(spearmanr(ebc_lsit, PPIpartner_jaccard), file=o1)
    print("protein A", "proteinB", "eBC", "PPI_partner_jaccard", file=o1)
    for i in range(len(protein_pair)):
        pa, pb = protein_pair[i].split(",")
        print(pa, pb, ebc_lsit[i], PPIpartner_jaccard[i], sep="\t", file=o1)