#!/bin/sh

# Initialize an empty array
pdb_list=()

# Read the file line by line
while IFS= read -r line; do
    # Add the line (pdb) to the array
    pdb_list+=("$line")
done < "pdb_list.txt"

# Specify the output directory
output_dir="`pwd`/pdb_files"

# Ensure the output directory exists
mkdir -p "$output_dir"

# Loop through the PDB IDs
for pdb in "${pdb_list[@]}"; do
    # Define the output file path
    output_file="$output_dir/$pdb/$pdb.pdb"

    # Check if the file already exists and is not empty
    if [ -s "$output_file" ]; then
        echo "PDB file $pdb.pdb already exists and is not empty. Skipping download."
    else
        # Download the PDB file to the specified output directory
        wget "http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=$pdb" -O "$output_file"
    fi
done