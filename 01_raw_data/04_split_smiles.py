import os

os.mkdir("SMILES")

# load data
with open("try.txt", "r") as file1:
    for line in file1:
        cid = line.replace("\n", "").split(",")[0]
        smiles = line.replace("\n", "").split(",")[-1]
        with open(".//SMILES//" + cid + ".smi", "w") as o1:
            print(smiles, file=o1)