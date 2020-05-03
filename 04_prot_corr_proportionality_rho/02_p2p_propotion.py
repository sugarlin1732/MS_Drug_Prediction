import numpy as np

# read rho corr matrix
rho_matrix = np.loadtxt("bbb.txt", delimiter=",")

"""
# read protein uniprot name
uniprot_name = []
with open("uniprot_name.txt", "r") as file1:
    for line in file1:
        uniprot_name.append(line.replace("\n", ""))
"""



counter = 1
# name & corr data
with open("prot2prot_bbb.txt", "w") as o1:
    for i in range(len(rho_matrix)-1):
        for j in range(i+1, len(rho_matrix)):
            print(i, end="\t", file=o1)
            print(j, end="\t", file=o1)
            print((rho_matrix[i, j]), file=o1)