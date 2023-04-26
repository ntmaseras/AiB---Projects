#!/bin/bash

quicktree_path="/home/nonatorruella/miniconda3/bin/quicktree"
rapidnj_path="/home/nonatorruella/miniconda3/bin/rapidnj"

# Create a CSV file to store the execution times
echo "file,QuickTreeNJ,RapidNJ" > execution_times.csv

for file in *.phy; do
   
    # Run QuickTree and record the execution time
    quicktree_output="${file}.quicktree"
    quicktree_start_time=$(date +%s.%N)
    "${quicktree_path}" -in m "${file}" > "${quicktree_output}"
    quicktree_end_time=$(date +%s.%N)
    quicktree_time=$(echo "${quicktree_end_time}-${quicktree_start_time}" | bc)

    # Run RapidNJ and record the execution time
    rapidnj_output="${file}.rapidnj"
    rapidnj_start_time=$(date +%s.%N)
    "${rapidnj_path}" "${file}" > "${rapidnj_output}"
    rapidnj_end_time=$(date +%s.%N)
    rapidnj_time=$(echo "${rapidnj_end_time}-${rapidnj_start_time}" | bc)
    
    # Write the execution times to the CSV file
    echo "${file},${quicktree_time},${rapidnj_time}" >> execution_times.csv
  

done

