import math

with open("05_all_data.txt", "r") as file1:
    with open("06_all_data_processed.txt", "w") as o1:
        for line in file1:
            data = line.split("\t")
            for i in range(len(data)):
                if i == 3:
                    print(-math.log(float(data[3]), 10), end="\t", file=o1)
                    continue
                if i == 5:
                    print(data[i], end="", file=o1)
                else:
                    print(data[i], end="\t", file=o1)



