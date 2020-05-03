import numpy as np
from scipy.stats.mstats import gmean


# load file
print("loading file...")
matrix = np.loadtxt("..//02_zscore_cutoff_2_matrix//protein_ki_zscore_2_matrix.txt", delimiter=",")
matrix = matrix + 1


# tranpose to fit gmean
matrix = np.transpose(matrix)


# clr(x)
print("generating clr matrix...")
new_transpose_matrix = []
for i in matrix:
    gx = gmean(i)
    new_transpose_matrix.append(list(map(lambda x: np.log(x/gx), i)))


# transpose back
new_matrix = np.transpose(new_transpose_matrix)



# rho
counteri = 1

print("calcuating rho correlation...")
corr_matrix = []
for i in new_matrix:
    rho_temp = []
    counterj = 2
    for j in new_matrix:
        print("protein " + str(counteri) + " to protein " + str(counterj))
        rho = 1 - (np.var(i-j)/(np.var(i) + np.var(j)))
        rho_temp.append(rho)
        counterj += 1
    corr_matrix.append(rho_temp)
    counteri += 1

corr_matrix = np.array(corr_matrix)
rho_modified_matrix = 0.5 + 0.5*corr_matrix
rho_modified_matrix[np.where(rho_modified_matrix == 1)] = 0

print("writing file...")

#np.savetxt(".//prot_rho_raw_matrix.txt", corr_matrix, fmt = "%.5e", delimiter=",")
#np.savetxt(".//prot_rho_modified_matrix.txt", rho_modified_matrix, fmt = "%.5e", delimiter=",")

np.savetxt(".//aaa.txt", corr_matrix, fmt = "%.5e", delimiter=",")
np.savetxt(".//bbb.txt", rho_modified_matrix, fmt = "%.5e", delimiter=",")

print("done...")