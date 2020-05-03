import numpy as np
import pubchempy as pcp

counter = 1
# load data
with open("022_no_duplicate.sdf", "r") as f1:
    with open("031_data_with_SMILES.txt", "w") as o1:
        for line in f1:
            print(counter)
            cid = line.split(",")[0]
            
            try:
                comp = pcp.Compound.from_cid(cid)
                smiles = comp.isomeric_smiles
                print(line.replace("\n", ""), end=",", file=o1)
                print(smiles, file=o1)
            except:
                print("error")
                continue
            
            counter += 1

        