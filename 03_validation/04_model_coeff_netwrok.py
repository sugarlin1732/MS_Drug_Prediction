import networkx as nx
import numpy as np

data = []
with open("model_to_model_coeff_scc.txt", "r") as file1:
    for line in file1.readlines()[1:]:
        aaa = line.replace("\n", "").replace("_HUMAN", "").split("\t")
        if float(aaa[2]) > 0.5:
            data.append([aaa[0], aaa[1], float(aaa[2])])

G = nx.Graph()

for i in data:
    G.add_edge(i[0], i[1], weight = i[2])

nx.write_graphml(G, "model_coeff_network.graphml")