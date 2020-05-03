import numpy as np

# load protein predicted -logKi matrix
ki_matrix = np.loadtxt(".//prot_SCC_modified_matrix.txt", delimiter=",")


# generate SCC matrix
scc_raw_matrix = []

counter = 1
# calculate SCC
for i in range(len(ki_matrix)):
    print(counter)
    scc_temp = []
    for j in range(len(ki_matrix)):
        if i == j:
            scc_temp.append(1)
        else:            
            scc_temp.append(np.corrcoef(ki_matrix[i], ki_matrix[j])[0, 1])
    scc_raw_matrix.append(scc_temp)
    counter += 1

scc_raw_matrix = np.array(scc_raw_matrix)
scc_modified_matrix = 0.5 + 0.5*scc_raw_matrix
scc_modified_matrix[np.where(scc_modified_matrix == 1)] = 0


np.savetxt("prot_SCC_raw_matrix.txt", scc_raw_matrix, fmt = "%.5e", delimiter=",")
np.savetxt("prot_SCC_modified_matrix.txt", scc_modified_matrix, fmt = "%.5e", delimiter=",")