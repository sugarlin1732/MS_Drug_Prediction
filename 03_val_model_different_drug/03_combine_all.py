import os
import math

file_list = os.listdir(".//result//")

with open("all_result.txt", "w") as o1:
    print("Model_used", "Drug_Number","Self_SCC","Self_Z","Val_Drug","Val_Drug_Number","SCC_of_Val_real","Val_Z","Val_Z_Bigger","log_zratio", sep="\t", file=o1) 
    for files in file_list:
        with open(".//result//" + files, "r") as file1:
            for line in file1.readlines()[1:]:
                data = line.split("\t")
                try:
                    print(line.replace("\n", ""), math.log(float(data[7])/float(data[3]), 10), sep="\t", file=o1)
                except:
                    pass