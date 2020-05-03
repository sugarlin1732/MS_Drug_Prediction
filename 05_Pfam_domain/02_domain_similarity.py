from scipy import spatial
from scipy.stats import spearmanr
import numpy as np


# load domain file
prot_dom_dic = {}
with open("prot_pfam_domain.txt", "r") as file1:
    for line in file1.readlines()[1:]:
        data = line.replace("\n", "").split("\t")
        prot_dom_dic[data[0]] = data[1].split(",")

prop_list = []
consine_list = []

counter = 1
# load all data
with open("prot_pair_rho_ppi_domainS.txt", "w") as o1:
    print("proteint_A", "proteint_B", "proportionality", "PPI", "domain_sim", sep="\t", file=o1)
    with open("prot_pair_rho_ppi.txt", "r") as file2:        
        for line in file2.readlines()[1:]:
            #print(counter)
            data = line.replace("\n", "").split("\t")
            pa, pb, rho, ppi = data
            prop_list.append(rho)
            pa_dom, pb_dom = prot_dom_dic[pa], prot_dom_dic[pb]

            com_list = []
            for i in pa_dom:
                if i not in com_list:
                    com_list.append(i)
            for j in pb_dom:
                if j not in com_list:
                    com_list.append(j)

            pa_bit = []
            pb_bit = []
            for k in com_list:
                if k in pa_dom:
                    pa_bit.append(1)
                else:
                    pa_bit.append(0)
                if k in pb_dom:
                    pb_bit.append(1)
                else:
                    pb_bit.append(0)
            consine_s = 1 - spatial.distance.cosine(pa_bit, pb_bit)
            consine_list.append(consine_s)
            print(pa, pb, rho, ppi, consine_s, sep="\t", file=o1)
            #counter += 1

    print(spearmanr(prop_list, consine_list), file=o1)