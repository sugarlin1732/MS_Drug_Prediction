from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["savefig.dpi"] = 300 #图片像素
plt.rcParams["figure.dpi"] = 300 #分辨率



rho_go_pair_list = []

with open("prot_pair_rho_semsim_biological_process_Resnik_max.txt", "r") as file1:
    for line in file1.readlines()[2:]:
        rho_go_pair_list.append(list(map(float, line.replace("\n", "").split("\t")[2:])))

rho_go_pair_list = sorted(rho_go_pair_list,key=lambda l:l[0], reverse=False)

low_pair = []
high_pair = []
for i in rho_go_pair_list:
    if i[0] < 0.5:
        low_pair.append(i)
    else:
        high_pair.append(i)

low_pair = np.array(low_pair)
high_pair = np.array(high_pair)

low_rho, low_ss = np.hsplit(low_pair, 2)
high_rho, high_ss = np.hsplit(high_pair, 2)

ss_list = [list(low_ss), list(high_ss)]

fig, ax = plt.subplots()
ax.boxplot(ss_list)

plt.show()
