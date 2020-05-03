import networkx as nx
from sklearn import preprocessing
import numpy as np

goal_loop = 1

# originalgraph
G0 = nx.read_graphml("rho_corr_network.graphml")


EB = nx.edge_betweenness_centrality(G0, k=None, normalized=True, weight="weight", seed=None)
betweenness_list = [EB[key] for key in EB]
betweenness_list = np.array(betweenness_list).reshape(1, -1)
betweenness_list = preprocessing.normalize(betweenness_list)
EB_old = betweenness_list[0]


#  loop
loop = 1
"""
while True:
                
    G = nx.Graph()
    index = 0
    for pair in EB:
        if EB_old[index] == 0:
            G.add_node(pair[0])
            G.add_node(pair[1])
            index += 1
            continue
        else:
            G.add_edge(pair[0], pair[1], weight = EB_old[index])
            index += 1

    print(EB_old)
    EB = nx.edge_betweenness_centrality(G, k=None, normalized=True, weight="weight", seed=None)
    betweenness_list = [EB[key] for key in EB]
    betweenness_list = np.array(betweenness_list).reshape(1, -1)
    betweenness_list = preprocessing.normalize(betweenness_list)

    print(betweenness_list[0])

    MSE = 1
    if len(EB_old) == len(betweenness_list[0]):
        mse = np.square(np.subtract(EB_old,betweenness_list[0])).mean() 
        print("Loop", str(loop), ", MSE =", mse)
        MSE = mse        
        
    print("Loop", str(loop) ,"number of edges =", G.number_of_edges())
    print()
    EB_old = betweenness_list[0]    

    if loop == (goal_loop+1):
        nx.write_graphml(G, "edge_betweenness_network_loop" + str(goal_loop) + ".graphml")            
        break

    # if MSE < 0.05:
    #     print("Converge")
    #     print("MSE =", MSE)
    #     print("Number of edges =", G.number_of_edges())
    #     print("Loop =", str(loop))
    #     #nx.write_graphml(G, "edge_betweenness_network.graphml")
    #     break

    loop += 1
"""


# 1 loop

G = nx.Graph()
index = 0
for pair in EB:
    if EB_old[index] == 0:
        G.add_node(pair[0])
        G.add_node(pair[1])
        index += 1
        continue
    else:
        G.add_edge(pair[0], pair[1], weight = EB_old[index])
        index += 1

nx.write_graphml(G, "edge_betweenness_network_loop1.graphml")            



print("END")