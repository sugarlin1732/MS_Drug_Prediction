
# load prot name
with open("..//02_drug_model//02_protein_order.txt", "r") as f1:
    prot_name = [line.replace("\n", "") for line in f1.readlines()]


# load proportionality PPI and domain_sin info, and GO semsim
info = {}

with open("..//07_pfam_domain//prot_pair_rho_ppi_domainS.txt", "r") as f2:
    file2 = f2.readlines()[1:-1]
    with open("..//06_go//prot_pair_rho_semsim_biological_process_Resnik_max.txt", "r") as f3:
        file3 = f3.readlines()[2:]
        for i in range(len(file2)):
            data = file2[i].replace("\n", "").split("\t")
            prot_pair = data[0] + "," + data[1]
            GOsemsin = file3[i].replace("\n", "").split("\t")[-1]
            info[prot_pair] = data[2] + "," + data[3] + "," + data[4] + "," + GOsemsin


# write file
with open("prot_pair_data.txt", "w") as o1:
    print("ProteinA", "ProteinB", "proportionality", "PPI", "Domain_sim", "GO_semsim", sep="\t", file=o1)
    for key in info:
        prot_name = key.split(",")
        rho = info[key].split(",")[0]
        prot_pair = info[key].split(",")[1]
        dom_sim = info[key].split(",")[2]
        go_sim = info[key].split(",")[3]
        print(prot_name[0], prot_name[1], rho, prot_pair, dom_sim, go_sim, sep = "\t", file=o1)
    