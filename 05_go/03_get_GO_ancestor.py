

# load gene product id
uniprot = []
ID = []
with open("uniprot_to_uniprot_entry.txt", "r") as file1:
    for line in file1.readlines():
        data = line.replace("\n", "").split("\t")
        uniprot.append(data[0])
        ID.append(data[1])


# make ID and term dic
GO_dict = {}
with open("goa_human.gaf", "r") as file2:
    temp_dic = {}
    for line in file2:
        if "!" in line:
            continue
        else:
            name = line.split("\t")[1]
            term = line.split("\t")[4]
        
        if name not in temp_dic:
            temp_dic[name] = []
        if term in temp_dic[name]:
            pass
        else:
            temp_dic[name].append(term)
    
    for key in temp_dic:
        if key in ID:
            GO_dict[key] = temp_dic[key]




# get ancestor
ancestor_dic = {}
with open("go-basic.obo", "r") as file3:
    all_data = file3.readlines()[25:545599]
    result = []
    tmp = []
    
    switch = 0
    for line in all_data:
        if line == "[Term]":
            result.append(tmp)
            tmp = []
        elif "id: " in line:
            tmp.append(line.replace("\n", "").split(" ")[-1])
            continue
        elif "namespace: " in line:
            if "biological_process" in line:
                switch = 1
            else:
                switch = 0
            continue
        elif "GO:" in line:
            go = "GO:" + line.replace("\n", "").split("GO:")[1].split(" !")[0]
            tmp.append(go)
        elif line == "\n":
            if switch == 1:
                result.append(tmp)
                tmp = []
            else:
                tmp = []
                switch = 0

    for i in result:
        ancestor_dic[i[0]] = i[1:]


def add_term(term_list, new_list):
    for i in term_list:
        if i not in new_list:
            new_list.append(i)
        try:
            ancestor_list = ancestor_dic[i]
            add_term(ancestor_list, new_list)
        except:
            pass



for key in GO_dict:
    new_list = []
    new_list = GO_dict[key]
    add_term(GO_dict[key], new_list)
    GO_dict[key] = new_list


for i in range(len(ID)):
    with open(".//prot_raw_BP_GO//" + uniprot[i] + ".txt", "w") as o1:
        for j in GO_dict[ID[i]]:
            print(j, file=o1)
