import networkx as nx
from statistics import mean
import os


# load GO term
GO_dic = {}
file_list = os.listdir(".//prot_raw_BP_GO")
for files in file_list:
    with open(".//prot_raw_BP_GO//" + files, "r") as gofile:
        go_terms = set(line.replace("\n", "") for line in gofile)
        GO_dic[files.replace(".txt", "")] = go_terms



# jaccard
def jaccard(path_list):
    GO_path_list = [GO_dic[key] for key in path_list]
    
    intersections = len(set.intersection(*GO_path_list))
    union_list = []
    
    if len(GO_path_list) > 2:
        for i in range(0, len(GO_path_list)-1):
            for j in range(i+1, len(GO_path_list)):
                union_list.append(set.intersection(set(GO_path_list[i]), set(GO_path_list[j])))
        unions = len(set().union(*union_list))
    else:
        unions = len(set().union(set(GO_path_list[0]), set(GO_path_list[1])))
    return 1. * intersections / unions
    


# load network
G = nx.read_graphml("..//08_betweenness_network//edge_betweenness_name_network_loop20.graphml")
edges = list(map(list, G.edges()))

# load protein order
with open("..//02_drug_model//02_protein_order.txt", "r") as file2:
    prot_name = [line.replace("\n", "") for line in file2]


# load shortest path
most_lenth = 0
path_dic = {}
with open("..//08_betweenness_network//shortest_ebc_path_length_betweenness_network.txt", "r") as file1:
    f1 = file1.readlines()[1:]
    for line in f1:
        data = line.replace("\n", "").split("\t")
        pair = prot_name[int(data[0])] + "," + prot_name[int(data[1])]
        path = [prot_name[int(i)] for i in data[3].split(",")]
        path_dic[pair] = path
        if len(path) > most_lenth:
            most_lenth = len(path)
        



# compute all path
total_count = 1

with open("compare_direct_step_GOsim.txt", "w") as o2:
    print("start", "end", "n_step", "n_step_GOsim", "direct_GOsim", sep="\t", file=o2)
    for lenth in range(2, most_lenth+1):
        print(total_count)        
        for key in path_dic:
            if len(path_dic[key]) == lenth:
                GO_sim = jaccard(path_dic[key])
                pa = key.split(",")[0]
                pb = key.split(",")[1]
                
                print(pa, pb, str(lenth-1), GO_sim, jaccard([pa, pb]), sep="\t", file=o2)
        
        total_count += 1



