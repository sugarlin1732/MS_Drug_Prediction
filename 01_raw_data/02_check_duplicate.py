# load data
dic = {}
with open("02_processed_data.sdf", "r") as f1:
    for line in f1:
        line_list = line.replace("\n", "").split(",")
        key = line_list[0] + "," + line_list[1] + "," + line_list[2]
        value = float(line_list[-1])

        if key not in dic:
            dic[key] = value
        else:
            if dic[key] < value:
                dic[key] = value
            else:
                pass

# write data
with open("022_no_duplicate.sdf", "w") as o1:
    for key in dic:
        print(key, end=",", file=o1)
        print(dic[key], file=o1)
        