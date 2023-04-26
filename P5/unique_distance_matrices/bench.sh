#!/bin/bash

quicktree_path="/home/nonatorruella/miniconda3/bin/quicktree"
rapidnj_path="/home/nonatorruella/miniconda3/bin/rapidnj"

# Create a CSV file to store the execution times
echo "file,QuickTreeNJ,RapidNJ" > execution_times.csv

for file in *.phy; do
    #quicktree_output="${file}.quicktree"
    rapidnj_output="${file}.rapidnj"

    # Run QuickTree and record the execution time
    quicktree_output="${file}.quicktree"
    quicktree_start_time=$(date +%s.%N)
    "${quicktree_path}" -in m "${file}" > "${quicktree_output}"
    quicktree_end_time=$(date +%s.%N)
    quicktree_time=$(echo "${quicktree_end_time}-${quicktree_start_time}" | bc)

    #quicktree_time="$(/usr/bin/time -f "%e" "${quicktree_path}" -in phylip "${file}" > "${quicktree_output}" 2>&1)"


    # Run RapidNJ and record the execution time
    rapidnj_start_time=$(date +%s.%N)
    "${rapidnj_path}" "${file}" > "${rapidnj_output}"
    rapidnj_end_time=$(date +%s.%N)
    rapidnj_time=$(echo "${rapidnj_end_time}-${rapidnj_start_time}" | bc)
    #quicktree_time = 0
    # Write the execution times to the CSV file
    echo "${file},${quicktree_time},${rapidnj_time}" >> execution_times.csv
  

done

