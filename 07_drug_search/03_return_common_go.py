p1 = "DRD2_HUMAN"
p2 = "BCL2_HUMAN"
p3 = "BCL2_HUMAN"

with open(".//DRD2_HUMAN.txt", "r") as file1:
    go1 = set(line.replace("\n", "") for line in file1)

with open(".//BCL2_HUMAN.txt", "r") as file2:
    go2 = set(line.replace("\n", "") for line in file2)

with open(".//BCL2_HUMAN.txt", "r") as file3:
    go3 = set(line.replace("\n", "") for line in file3)


int12 = list(set.intersection(go1, go2))
int23 = list(set.intersection(go2, go3))
int123 = list(set.intersection(go1, go2, go3))

with open(".//go-basic.obo", "r") as file4:
    f4 = file4.readlines()[25:545598]
    string = ""
    for i in f4:
        if ("[Term]" in i) or ("id: GO:" in i) or ("name" in i) or ("alt_id" in i) or ("namespace" in i) or ("def: " in i):
            string += i
    data_list = string.split("[Term]")

    with open("DRD2_BCL2_BCL2_shared_GO.txt", "w") as o2:
        for term in data_list:            
            if term == "":
                pass
            else:
                term_data = term.split("GO:")
                go_list = []
                for j in term_data[1:]:
                    termid = "GO:" + j.split("\n")[0]
                    go_list.append(termid)

                for goterm in go_list:
                    if goterm in int123:
                        print(term, end="", file=o2)




"""
with open("shared_goterm.txt", "w") as o1:
    for i in int12:
        print(i, end="\t", file=o1)
    print(file=o1)
    for j in int23:
        print(j, end="\t", file=o1)
    print(file=o1)
    for k in int123:
        print(k, end="\t", file=o1)
"""