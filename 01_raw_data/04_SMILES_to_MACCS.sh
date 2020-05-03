
cd /home/sugarlin1732/work/total_new/01_raw_data
mkdir fpt_MACCS
cd /home/sugarlin1732/work/total_new/01_raw_data/SMILES

for SMILES_file in /home/sugarlin1732/work/total_new/01_raw_data/SMILES/*
    do
        filename=$(basename "$SMILES_file" .smi)
        input="$filename".smi
        output="$filename".fpt

        babel -ismi $input -ofpt $output -xfMACCS
        mv $output /home/sugarlin1732/work/total_new/01_raw_data/fpt_MACCS
    done
exit 0