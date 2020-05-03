import numpy as np
import pandas as pd
from scipy import stats


# load predict -ligKi matrix
ki_matrix = np.loadtxt("..//..//02_drug_model//03_predicted_ki_matrix.txt", delimiter=",")

# generate zscore matrix
z_score_matrix = []

# calculate zscore
counter = 1
for prot in ki_matrix:
    print(counter)
    z_score = stats.zscore(prot)
    z_score[z_score < 2] = 0
    z_score_matrix.append(z_score)

    counter += 1


z_score_matrix = np.array(z_score_matrix)
np.savetxt(".//protein_ki_zscore_2_matrix.txt", z_score_matrix, fmt = "%.5e", delimiter=",")




