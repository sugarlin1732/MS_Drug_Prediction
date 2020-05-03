

p1 = input("protein1")
p2 = input("protein2")

with open("..//02_drug_model//drug_per_prot//" + p1 + ".txt", "r") as file1:
    drug1 = sorted([int(line.split("\t")[0]) for line in file1])

with open("..//02_drug_model//drug_per_prot//" + p2 + ".txt", "r") as file2:
    drug2 = sorted([int(line.split("\t")[0]) for line in file2])


with open(".//drugs_of_" + p1 + "_" + p2 + ".txt", "w") as o1:
    print(p1, len(drug1), p2, len(drug2), sep="\t", file=o1)
    for i in drug1:
        if i not in drug2:
            print(i, "", sep="\t", file=o1)
    for j in drug2:
        if j not in drug1:
            print("", j, sep="\t", file=o1)

    for k in list(set.intersection(set(drug1), set(drug2))):
        print(k, k, sep="\t", file=o1)

