import networkx as nx
from scipy.stats import spearmanr

# calculate shortest path with betweenness
G0 = nx.read_graphml("edge_betweenness_network_loop20.graphml")

# make proportionality weight
G1 = G0.copy()
edge_list = list(G1.edges())
rho_dic = {}

with open("02_protein_order.txt", "r") as file2:
    prot_name = [line.replace("\n", "") for line in file2]

with open("prot_pair_data.txt", "r") as file1:
    f1 = file1.readlines()[1:]
    for line in f1:
        data = line.split("\t")
        pa = data[0]
        pb = data[1]
        rho = data[2]        
        rho_dic[str(prot_name.index(pa)) + "," + str(prot_name.index(pb))] = rho

for edge in edge_list:
    key1 = str(edge[0]) + "," + str(edge[1])
    key2 = str(edge[1]) + "," + str(edge[0])
    if (key1 in rho_dic):
        G1.add_edge(str(edge[0]), str(edge[1]), weight=float(rho_dic[key1]))        
    elif (key2 in rho_dic):
        G1.add_edge(str(edge[0]), str(edge[1]), weight=float(rho_dic[key2]))
        


# write file
b_path_list = []
with open("shortest_ebc_path_length_betweenness_network.txt", "w") as o1:
    print("Start node", "End node", "Betweenness weight length", "Shortest path", sep="\t", file=o1)
    for i in range(0, len(G0.nodes())-1):
        for j in range(i+1, len(G0.nodes())):
            b_shortest = nx.shortest_path_length(G0, source=str(i), target=str(j), weight= "weight")
            b_full_lenth = nx.shortest_path_length(G0, source=str(i), target=str(j))
            b_path = nx.shortest_path(G0, source=str(i), target=str(j), weight= "weight")
            print(i, j, b_full_lenth-b_shortest, ",".join(map(str,b_path)), sep="\t", file=o1)

r_path_list = []
with open("shortest_rho_path_length_betweenness_network.txt", "w") as o2:
    print("Start node", "End node", "Proportionality weight length", "Shortest path", sep="\t", file=o2)
    for i in range(0, len(G1.nodes())-1):
        for j in range(i+1, len(G1.nodes())):
            r_shortest = nx.shortest_path_length(G1, source=str(i), target=str(j), weight= "weight")
            r_full_lenth = nx.shortest_path_length(G1, source=str(i), target=str(j))
            r_path = nx.shortest_path(G1, source=str(i), target=str(j), weight= "weight")
            print(i, j, r_full_lenth-r_shortest, ",".join(map(str,r_path)), sep="\t", file=o2)

##################
for i in range(len(b_path_list)):
    if b_path_list[i] != r_path_list:
        print(i)