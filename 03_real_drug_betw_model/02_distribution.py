import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import os 

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["savefig.dpi"] = 300 #图片像素
plt.rcParams["figure.dpi"] = 300 #分辨率


file_list = os.listdir(".//result//")

counter = 1
for files in file_list:
    # load data files
    with open(".//result//" + files) as file1:
        print(counter)
        f1 = file1.readlines()
        drug_number = f1[0].replace("\n", "").split(" ")[-1]
        scc_lsit = [float(line.split("\t")[1]) for line in f1[3:]]
        self_pred = str(round(float(f1[2].split("\t")[1]), 2))
        
        prot_drug = files.replace("real_", "").replace("_HUMAN_drugs_w_all_model_scc.txt", "")

        # image parameter
        fig, ax = plt.subplots(figsize=(20,10))
        ax.hist(scc_lsit, bins=50, normed=0, facecolor="blue", edgecolor="black", alpha=0.7)
        


        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%d"))
        ax.tick_params(direction="in", labelsize=15)
        xstart, xend = ax.get_xlim()
        ystart, yend = ax.get_ylim()
        ax.xaxis.set_ticks(np.arange(xstart, xend, 0.1))
        ax.set_xlabel("SCC", fontsize=20)
        ax.set_ylabel("Count",fontsize=20)
        ax.set_title("SCC_of_" + prot_drug + "_drug_predicted_other_models",fontsize=20)
        ax.text(xend*0.43, yend*0.87, "Number of Drugs = " + drug_number + "\nself_pred_scc = " + self_pred, fontsize=25)

        fig.tight_layout()
        fig.savefig(".//distribution_fig//SCC_of_real_" + prot_drug + "_drug_predicted_other_models.png")
        plt.clf()
        counter += 1