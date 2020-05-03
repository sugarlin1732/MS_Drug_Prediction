import numpy as np

dic = {}
counter = 1
with open("05_all_data.txt", "r") as file1:
    with open("05_all_drug_info.txt", "w") as o1:
        for line in file1:
            print(counter)
            data = line.replace("\n", "").split("\t")
            cid = data[0]
            inichi = data[2]
            fpt = data[4]
            properties = data[5]

            if cid in dic:
                continue
            else:
                dic[cid] = cid + "\t" + inichi + "\t" + fpt + "\t" + properties
            counter += 1
            

        for key in dic:
            print(dic[key], file=o1)
