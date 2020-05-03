import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os
import statistics

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["savefig.dpi"] = 300 #图片像素
plt.rcParams["figure.dpi"] = 300 #分辨率


model_list = os.listdir(".//result//")
zscore_dic = {}
median_dic = {}

# load zscore of model to drugs
for model in model_list:
    with open(".//result//" + model, "r") as file1:
        prot = model.replace("_HUMAN_model_pred_other_drug.txt", "")
        f1 = file1.readlines()[1:]
        zscore_list = []
        for line in f1:
            zscore = line.split("\t")[-2]
            zscore_list.append(float(zscore))
        zscore_dic[prot] = zscore_list
        median_dic[prot] = [statistics.median(zscore_list)]

# sort and pick top 30 median
new_median = list(sorted(median_dic.items(), key=lambda d: d[1], reverse=True))[:30]
# new_median = list(sorted(median_dic.items(), key=lambda d: d[1], reverse=True))[:2]
top30_prot = []
top30_zscore = []


for i in new_median:
    top30_prot.append(i[0])
    top30_zscore.append(zscore_dic[i[0]])

# boxplot
data = top30_zscore
fig1, ax1 = plt.subplots(figsize=(20,10))

ax1.set_title("Z-score of top-30-median Models")
ax1.boxplot(data, labels=top30_prot)

plt.savefig("zscore_boxplot.png")
