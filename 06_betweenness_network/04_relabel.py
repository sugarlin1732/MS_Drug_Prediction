import networkx as nx

G0 = nx.read_graphml("edge_betweenness_network_loop20.graphml")

with open("02_protein_order.txt", "r") as file2:
    prot_name = [line.replace("\n", "") for line in file2]

mapping = {}
node = list(G0.nodes())

for i in range(len(node)):
    mapping[node[i]] = prot_name[int(node[i])].replace("_HUMAN", "")

G0 = nx.relabel_nodes(G0, mapping)

nx.write_graphml(G0, "edge_betweenness_name_network_loop20.graphml")
