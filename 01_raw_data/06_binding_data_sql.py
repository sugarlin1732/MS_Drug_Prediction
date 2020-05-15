import MySQLdb

# connect mySQL
db = MySQLdb.connect(db="ms_thesis", host="localhost", user="root", passwd="1123581321")
cur = db.cursor()

# create table
code_list = []
create_table = "create table if not exists drug_processed_data (\
    PubChem_Cid int NOT NULL, \
    Target_Name varchar(15) NOT NULL, \
    InChI_Key varchar(35) NOT NULL, \
    minus_logKi float(25, 20) NOT NULL, \
    MACCS_hex varchar(64) NOT NULL \
    );"
code_list.append(create_table)
# create 17 property columns
for i in range(17):
    code_list.append("alter table drug_processed_data add column property_" + str(i+1) \
        + " float(10, 4) NOT NULL;")

# execute sql code: create table
for i in code_list:
    cur.execute(i)
    results = cur.fetchall()

    # print result
    for record in results:
        col1 = record[0]
        col2 = record[1]
        print("%s, %s" % (col1, col2))
code_list = []


# load drug data
with open(".//06_all_data_processed.txt", "r") as file1:
    data_tuple = ""
    for line in file1.readlines():
        data = line.replace("\n", "").split("\t")
        cid = str(data[0])
        target = str(data[1])
        inchi = str(data[2])
        ki = str(data[3])
        MACCS = data[4].replace("[", "").replace("]", "").replace(", ", "")
        MACCS = "%0*X" % ((len(MACCS) + 3) // 4, int(MACCS, 2))
        properties = [str(i) for i in data[-1].replace("[", "").replace("]", "").split(", ")]
        tmp = cid + ",'" + target + "','" + inchi + "'," + ki + ",'" + MACCS + "',"
        for i in properties:
            tmp += i + ","
        tmp = tmp.rstrip(",")
        
        # insert into sql
        code = "insert into drug_processed_data values (" + tmp + ");"
        cur.execute(code)



db.commit()
db.close()