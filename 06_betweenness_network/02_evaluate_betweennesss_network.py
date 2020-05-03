# compare prot pairs with betweenness or not
# proportionaily, domain_sim, GO_semsim use wilcoxon and effect size (odds ratio)
# ppi use fisher test (odds ratio)
import numpy as np
from scipy import stats
import statistics

# cohen's d
#correct if the population S.D. is expected to be equal for the two groups.
def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)



# load betweenness data
betweenness_rho_list = []
betweenness_ppi_list = []
betweenness_domain_sim_list = []
betweenness_GO_semsim_list = []
betweenness_prot_pair_list = []
with open("prot_pair_data_of_betweeness_network.txt", "r") as file1:
    for line in file1.readlines()[1:]:
        data = line.replace("\n", "").split("\t")
        betweenness_rho_list.append(float(data[2]))
        betweenness_ppi_list.append(float(data[4]))
        betweenness_domain_sim_list.append(float(data[5]))
        betweenness_GO_semsim_list.append(float(data[6]))
        betweenness_prot_pair_list.append(data[0] + "," + data[1])
        betweenness_prot_pair_list.append(data[1] + "," + data[0])


# load non_betweenness data
NB_rho_list = []
NB_ppi_list = []
NB_domain_sim_list = []
NB_GO_semsim_list = []

count = 0
with open("prot_pair_data.txt", "r") as file2:
    for line in file2.readlines()[1:]:
        data = line.replace("\n", "").split("\t")
        prot_pair = data[0] + "," + data[1]
        if prot_pair in betweenness_prot_pair_list:
            count +=1
            continue
        else:
            NB_rho_list.append(float(data[2]))
            NB_ppi_list.append(float(data[3]))
            NB_domain_sim_list.append(float(data[4]))
            NB_GO_semsim_list.append(float(data[5]))


# PPI & betweenness fisher exact test
b_p = 0
b_np = 0
nb_p = 0
nb_np = 0
for i in betweenness_ppi_list:
    if i == 1:
        b_p += 1
    else:
        b_np += 1
for j in NB_ppi_list:
    if j == 1:
        nb_p += 1
    else:
        nb_np += 1


with open("result.txt", "w") as o1:
    print("Fisher exact test of PPI =", stats.fisher_exact([[b_p, nb_p], [b_np, nb_np]]), sep="\t", file=o1)
    
    # print wilcoxon
    print("proportionality wilcoxon =", stats.mannwhitneyu(betweenness_rho_list, NB_rho_list), file=o1)
    print("domain_sim wilcoxon =", stats.mannwhitneyu(betweenness_domain_sim_list, NB_domain_sim_list), file=o1)
    print("GO_semsim wilcoxon =", stats.mannwhitneyu(betweenness_GO_semsim_list, NB_GO_semsim_list), file=o1)

    # print cohens' d
    print("proportionality cohens' d =", cohen_d(betweenness_rho_list, NB_rho_list), file=o1)
    print("domain_sim cohens' d =", cohen_d(betweenness_domain_sim_list, NB_domain_sim_list), file=o1)
    print("GO_semsim cohens' d =", cohen_d(betweenness_GO_semsim_list, NB_GO_semsim_list), file=o1)


