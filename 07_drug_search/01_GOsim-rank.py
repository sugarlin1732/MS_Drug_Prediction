with open("compare_direct_step_GOsim.txt", "r") as file1:
    f1 = file1.readlines()[1:30720]

    with open("within_2steps_GOsim_rank.txt", "w") as o1:
        print("start", "end", "n_step", "n_step_GOsim", "direct_GOsim", sep="\t", file=o1)
        info_dic = {}
        for line in f1:
            data = line.replace("\n", "").split("\t")
            nsim = data[3]
            osim = data[4]
            if osim > nsim:
                continue
            info_dic[data[0] + "," + data[1] + "," + data[2] + "," + data[4]] = data[3]

        rank_list = sorted(info_dic.items(), key=lambda d: d[1], reverse=True)

        for i in rank_list:
            data = i[0].split(",")
            print(data[0], data[1], int(data[2]), float(i[1]), float(data[3]), sep="\t", file=o1)