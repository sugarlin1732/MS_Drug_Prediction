import networkx as nx
from statistics import mean
import os

"""
# load GO term
GO_dic = {}
file_list = os.listdir("..//06_go//prot_raw_BP_GO")
for files in file_list:
    with open("..//06_go//prot_raw_BP_GO//" + files, "r") as gofile:
        go_terms = set(line.replace("\n", "") for line in gofile)
        GO_dic[files.replace(".txt", "")] = go_terms

"""
# load protein order
with open("..//02_drug_model//02_protein_order.txt", "r") as file2:
    protein_lsit = [line.replace("\n", "") for line in file2]

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
for key in ppi_dic:
    ppi_dic[key] = set(ppi_dic[key])


# jaccard
def jaccard(path_list):
    try:
        PPI_path_list = [ppi_dic[key] for key in path_list]
    except:
        return 0
    intersections = len(set.intersection(*PPI_path_list))
    union_list = []
    
    if len(PPI_path_list) > 2:
        for i in range(0, len(PPI_path_list)-1):
            for j in range(i+1, len(PPI_path_list)):
                union_list.append(set.intersection(set(PPI_path_list[i]), set(PPI_path_list[j])))
        unions = len(set().union(*union_list))
    else:
        unions = len(set().union(set(PPI_path_list[0]), set(PPI_path_list[1])))
    try:
        return 1. * intersections / unions
    except:
        return 0



# load network
G = nx.read_graphml("..//08_betweenness_network//edge_betweenness_name_network_loop20.graphml")
edges = list(map(list, G.edges()))



# load shortest path
most_lenth = 0
path_dic = {}
with open("..//08_betweenness_network//shortest_ebc_path_length_betweenness_network.txt", "r") as file1:
    f1 = file1.readlines()[1:]
    for line in f1:
        data = line.replace("\n", "").split("\t")
        pair = protein_lsit[int(data[0])] + "," + protein_lsit[int(data[1])]
        path = [protein_lsit[int(i)] for i in data[3].split(",")]
        path_dic[pair] = path
        if len(path) > most_lenth:
            most_lenth = len(path)
        


# compute all path
final_result = {}
total_count = 1
for lenth in range(2, most_lenth+1):
    print(total_count)
    total_sim = []
    with open(".//ppi_partner_by_steps//" + str(lenth-1) + "step_PPIpart_result_EBC_network.txt", "w") as o1:
        for key in path_dic:
            if len(path_dic[key]) == lenth:
                for i in path_dic[key]:
                    print(i, end="\t", file=o1)
                PPIpart_sim = jaccard(path_dic[key])
                print(PPIpart_sim, file=o1)
                total_sim.append(PPIpart_sim)
        print("average PPIpart_sim =", mean(total_sim), file=o1)
    final_result[lenth-1] = mean(total_sim)
    total_count += 1
                
with open("step_number_PPIpart.txt", "w") as o2:
    print("step_number", "PPIpart", sep="\t", file=o2)
    for key in final_result:
        print(key, final_result[key], sep="\t", file=o2)



