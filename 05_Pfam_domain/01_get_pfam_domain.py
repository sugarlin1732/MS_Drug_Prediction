import requests
import time
import numpy as np

requests.adapters.DEFAULT_RETRIES = 5 # retry times
session = requests.Session()
session.keep_alive = False # remove extra connection

# load protein name
with open("..//02_drug_model//02_protein_order.txt", "r") as file1:
    prot_name_list = [line.replace("\n", "") for line in file1]

# crawler
counter = 1
with open("prot_pfam_domain.txt", "w") as o1:
    print("protein_name", "pfam_domain", sep="\t", file=o1)
    for protein in prot_name_list:
        print(counter)
        url = "https://pfam.xfam.org/protein/" + protein

        data = session.get(url).text.split("<tbody>")[2].split('<tr class="odd">')
        domain_line = [i.replace("            ", "").split("\n")[2] for i in data[1:]]
        domain_data = [j.split(">")[2].replace("</a", "") for j in domain_line if '<td><span class="inactive">n/a</span></td>' not in j]
        
        print(protein, end="\t", file=o1)
        for k in range(len(domain_data)):
            if k == (len(domain_data) - 1):
                print(domain_data[k], file=o1)
            else:
                print(domain_data[k], end=",", file=o1)
        counter += 1
