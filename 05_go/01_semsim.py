import fastsemsim
import numpy as np
from scipy.stats import spearmanr


# SS settings
ontology = fastsemsim.load_ontology(ontology_type = "GeneOntology")
ac = fastsemsim.load_ac(ontology, species="human")
semsim_type="obj"
semsim_measure="Cosine"
"""
Resnik
Lin
Jiang-Conrath
SimGIC
SimUI
SimIC
SimRel
Dice
SimTO
SimNTO
Jaccard
Czekanowski-Dice
Cosine
GSESAME
SimICND
SimICNP
"""
mixing_strategy="max"
"""
Average
BMA (Best Match Average)
max
"""

GO_type = "molecular_function"
"""
molecular_function
cellular_component
biological_process
"""

# Initializing semantic similarity
ss = fastsemsim.init_semsim(ontology = ontology, ac = ac, semsim_type = semsim_type, semsim_measure = semsim_measure, mixing_strategy = mixing_strategy)


# load protein name
with open(".//uniprot_to_uniprot_entry.txt", "r") as file1:
    prot_name = np.array([line.replace("\n", "").split("\t")[1] for line in file1])


# load protein proportionality rho corr matrix
p2p_rho = np.loadtxt("prot2prot_bbb.txt", delimiter="\t")
rho_list = p2p_rho[:, 2]

protein_pair = []
for i in p2p_rho[:, 0:2]:
    pair = [int(x) for x in i]
    protein_pair.append(prot_name[pair].tolist())
protein_pair = np.array(protein_pair)



# Calculating SS
ss_list = []
pair_count = 1
for pairs in protein_pair:
    print(pair_count)
    res = ss.SemSim(pairs[0], pairs[1], root=GO_type)
    ss_list.append(res)
    pair_count += 1
ss_list = np.array(ss_list)



# calculate SCC
ss_for_scc = []
rho_for_scc = []

for i in range(len(ss_list)):
    if type(ss_list[i]) == type(1.1): # check if ss value is Nonetype
        ss_for_scc.append(ss_list[i])
        rho_for_scc.append(rho_list[i])

SCC = spearmanr(rho_for_scc, ss_for_scc)
"""
####################################
with open("prot_pair_rho_semsim_" + GO_type + "_" + semsim_measure + "_" + mixing_strategy + ".txt", "w") as o3:
    print(SCC, file=o3)
    print("proteint_A", "proteint_B", "proportionality", "semsim", sep="\t", file=o3)
    
    for i in range(len(protein_pair)):
        print(protein_pair[i, 0], protein_pair[i, 1], rho_list[i], ss_list[i], sep="\t", file=o3)
        
####################################
"""
