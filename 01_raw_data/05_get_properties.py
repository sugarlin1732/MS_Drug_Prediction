import numpy as np
import pandas as pd
import os
from chembl_webresource_client.new_client import new_client

os.mkdir("Properties")
# properties that are chosen
choose = ["acd_logd", "acd_logp", "alogp", "aromatic_rings", "full_mwt",
          "hba", "hba_lipinski", "hbd", "hbd_lipinski", "heavy_atoms", "mw_freebase", "mw_monoisotopic", "num_lipinski_ro5_violations", "num_ro5_violations", "psa", "qed_weighted", "rtb"]

# load data
drug_inchi = []
with open("022_no_duplicate.sdf", "r") as file1:
    for line in file1:
        inchi = line.split(",")[2]
        drug_inchi.append(inchi)

drug_inchi = list(set(drug_inchi))

counter = 1
error = 0
for drug in drug_inchi:
    print("counter " + str(counter))
    print("errror " + str(error))
    try:
        properties_list = []        
        molecule = new_client.molecule
        for properties in choose:
            properties_list.append(molecule.get(drug)["molecule_properties"][properties])
        with open(".//Properties//" + drug + ".txt", "w") as o1:
            for i in properties_list:
                print(float(i), end=",", file=o1)
    except:
        print("error")
        error += 1
        pass
    counter += 1

print("END")
