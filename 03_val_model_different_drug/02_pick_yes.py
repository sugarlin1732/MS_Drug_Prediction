import os

file_list = os.listdir(".//result//")

with open("model_different_drugs.txt", "w") as o1:
    print("Model_used", "Drug_Number", "Self_SCC", "Self_Z", "Val_Drug", "Val_Drug_Number", "SCC_of_Val_real", "Val_Z", "Val_Z_Bigger", sep = "\t", file=o1)
    for files in file_list:
        with open(".//result//" + files, "r") as file1:
            for line in file1.readlines()[1:]:
                if "yes" in line.split("\t")[-1]:
                    print(line, end="", file=o1)
                