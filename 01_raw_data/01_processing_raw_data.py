def printCheck(tmp_data):
    ki_check = tmp_data['Ki_nM'].strip() != '' and not '>' in tmp_data['Ki_nM'] and not '<' in tmp_data['Ki_nM'] and float(tmp_data['Ki_nM']) > 0
    uniprot_check = '_HUMAN' in tmp_data['UniProt_SwissProt_Entry_Name_of_Target_Chain']
    cid_check = tmp_data['PubChem_CID'].strip() != ''
    ichi_check = tmp_data['Ligand_InChI_Key'].strip() != ''
    
    if ki_check and uniprot_check and cid_check and ichi_check:
        with open("processed_data.sdf", "a") as o1:
            print(tmp_data['PubChem_CID'].strip(), tmp_data['UniProt_SwissProt_Entry_Name_of_Target_Chain'].strip(),\
                tmp_data['Ligand_InChI_Key'].strip(), tmp_data['Ki_nM'].strip(), sep=',', file=o1)

with open('BindingDB_All_3D.sdf', 'r') as f:
    tmp_data = {'Structure':''}
    subkey = 'Structure'
    for line in f:
        if '$$$$' in line:
            printCheck(tmp_data)
            tmp_data = {'Structure':''}
            subkey = 'Structure'
        elif '> <' in line:
            subkey = line.strip().replace('> <','').replace('>','').replace('(','').replace(')','').replace(' ','_')
            tmp_data[subkey] = ''
        else:
            tmp_data[subkey] += line
            
         
            
