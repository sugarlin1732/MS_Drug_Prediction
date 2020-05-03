from scipy.stats import wilcoxon
from scipy import stats
import statistics
import numpy as np

# cohen's d
def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)



total_1 = []
total_not1 = []
step_dic = {}
with open("compare_direct_step_GOsim.txt", "r") as file1:
    f1 = file1.readlines()[1:]
    for line in f1:
        data = line.replace("\n", "").split("\t")
        step = data[2]
        if step == "1":
            continue
        else:
            total_1.append(float(data[4]))
            total_not1.append(float(data[3]))
            try:
                step_dic[int(step)][0].append(float(data[4]))
                step_dic[int(step)][1].append(float(data[3]))
            except:
                step_dic[int(step)] = [[float(data[4])], [float(data[3])]]
        


with open("compare_1_nstep_wilcoxon.txt", "w") as o1:
    print("1_step VS not_1_step (" + str(len(total_1)) + " pairs)", file=o1)
    print(wilcoxon(total_1, total_not1), file=o1)
    print("Cohen's d =", cohen_d(total_1, total_not1), file=o1)
    print(file=o1)
    for key in step_dic:
        print("1_step VS " + str(key) + "_step (" + str(len(step_dic[key][0])) + " pairs)", file=o1)
        print(wilcoxon(step_dic[key][0], step_dic[key][1]), file=o1)
        print("Cohen's d =", cohen_d(step_dic[key][0], step_dic[key][1]), file=o1)
        print(file=o1)
