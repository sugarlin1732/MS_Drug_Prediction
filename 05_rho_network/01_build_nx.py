import networkx as nx
import numpy as np


rho_corr = np.loadtxt("bbb.txt", delimiter = ",")


G = nx.Graph()

for i in range(len(rho_corr)):
    for j in range(len(rho_corr)):
        if i == j:
            continue
        G.add_edge(i, j, weight = rho_corr[i, j])

print(G.number_of_nodes())
print(G.number_of_edges())


nx.write_graphml(G, "rho_corr_network.graphml")