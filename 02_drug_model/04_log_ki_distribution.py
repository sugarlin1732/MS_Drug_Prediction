import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["savefig.dpi"] = 300 #图片像素
plt.rcParams["figure.dpi"] = 300 #分辨率



# load protein predicted -logKi matrix
ki_matrix = np.loadtxt("03_predicted_ki_matrix.txt", delimiter=",").ravel()


# image parameter
plt.hist(ki_matrix, bins=100, normed=0, facecolor="blue", edgecolor="black", alpha=0.7)
#plt.yticks(np.linspace(0, 100000, 6))
plt.xlim((-20, 10))
plt.xticks(np.linspace(-20, 10, 7))
plt.xlabel("predicted_-logKi")
plt.ylabel("Count")
plt.title("predict_-logKi&count_distribution")

# save image
plt.tight_layout()
plt.savefig(".//predict_logKi_count.png")
