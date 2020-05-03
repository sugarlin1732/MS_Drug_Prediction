import networkx as nx


# load prot name
with open("..//02_drug_model//02_protein_order.txt", "r") as f1:
    prot_name = [line.replace("\n", "") for line in f1.readlines()]


# load betweenness network
G = nx.read_graphml("edge_betweenness_network_loop20.graphml")
edge_data = list(G.edges.data())

betweenness_dic = {}
for pair in edge_data:
    pair_name = prot_name[int(pair[0])] + "," + prot_name[int(pair[1])]
    betweenness_dic[pair_name] = pair[2]["weight"]


# load proportionality PPI and domain_sin info, and GO semsim
# ex: {'5HT1A_HUMAN,5HT1B_HUMAN': '1,1.0,10.719368456811244'}
info = {}

with open("..//07_pfam_domain//prot_pair_rho_ppi_domainS.txt", "r") as f2:
    file2 = f2.readlines()[1:-1]
    with open("..//06_go//prot_pair_rho_semsim_biological_process_Resnik_max.txt", "r") as f3:
        file3 = f3.readlines()[2:]
        for i in range(len(file2)):
            data = file2[i].replace("\n", "").split("\t")
            prot_pair = data[0] + "," + data[1]
            GOsemsin = file3[i].replace("\n", "").split("\t")[-1]
            info[prot_pair] = data[2] + "," + data[3] + "," + data[4] + "," + GOsemsin


# write file
with open("prot_pair_data_of_betweeness_network.txt", "w") as o1:
    print("ProteinA", "ProteinB", "Proportionality", "Betweenness", "PPI", "Domain_sim", "GO_semsim", sep="\t", file=o1)
    for key in betweenness_dic:
        prot_A = key.split(",")[0]
        prot_B = key.split(",")[1]
        if key in info:   
            rho = info[key].split(",")[0]         
            ppi = info[key].split(",")[1]
            dom_sim = info[key].split(",")[2]
            GO_ss = info[key].split(",")[3]
            print(prot_A, prot_B, rho, betweenness_dic[key], ppi, dom_sim, GO_ss, sep="\t", file=o1)

        else:
            info_key = prot_B + "," +prot_A
            rho = info[info_key].split(",")[0]         
            ppi = info[info_key].split(",")[1]
            dom_sim = info[info_key].split(",")[2]
            GO_ss = info[info_key].split(",")[3]
            print(prot_A, prot_B, rho, betweenness_dic[key], ppi, dom_sim, GO_ss, sep="\t", file=o1)

