
dic = {}
with open("E:\\Research\\total_new\\01_raw_data\\06_all_data_processed.txt", "r") as file1:
    for line in file1:
        data = line.replace("\n", "").split("\t")
        target = data[1]
        
        try:
            dic[target].append(data)
        
        except:
            dic[target] = [data]
    

for key in dic:
    with open("drug_per_prot//" + key + ".txt", "w") as o1:
        for value in dic[key]:
            for i in range(len(value)):
                if value[i] == value[-1]:
                    print(value[i], file=o1)
                else:
                    print(value[i], end="\t", file=o1)