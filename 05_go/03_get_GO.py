from bioservices import QuickGO

s = QuickGO()

# load gene product id
uniprot = []
ID = []
with open("uniprot_to_uniprot_entry.txt", "r") as file1:
    for line in file1.readlines():
        data = line.replace("\n", "").split("\t")
        uniprot.append(data[0])
        ID.append(data[1])



# get go term id
for i in range(len(ID)):
    print(ID[i])
    result = s.Annotation(limit=-1, aspect="P", geneProductId=ID[i], taxonId="9606")
    with open(".//prot_raw_BP_GO//" + uniprot[i] + ".txt", "w") as o1:
        GO_list = []
        for key in result:
            if key == "results":
                for i in result[key]:
                    if i["goId"] not in GO_list:
                        GO_list.append(i["goId"])
        for i in GO_list:
            print(i, file=o1)

