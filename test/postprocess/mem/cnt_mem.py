import os
import re
import sys
import pandas as pd

def find_max_memory_in_file(file_path, method="BSX"):
    """Find the maximum value of `max_memory` in the file."""
    max_memory = 0
    if (method != "GUP"):
        pattern = re.compile(r"\s*max_memory:\s*([0-9]+)\s*KB")
    else:
        pattern = re.compile(r"\tMaximum resident set size \(kbytes\):\s*([0-9]+)")

    with open(file_path, 'r') as file:
        for line in file:
            match = pattern.search(line)
            if match:
                current_memory = int(match.group(1))
                if current_memory > max_memory:
                    max_memory = current_memory
                    print(f"Found max_memory: {current_memory} KB in {file_path}")

    return max_memory

def find_max_memory_in_directory(directory, method):
    """Recursively traverse the directory to find the maximum max_memory value in all files."""
    max_memory = 0

    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            try:
                # Find the maximum max_memory value of the current file.
                current_max = find_max_memory_in_file(file_path, method)
                if current_max > max_memory:
                    max_memory = current_max
            except Exception as e:
                print(f"Error reading file {file_path}: {e}")
    
    return max_memory

if __name__ == "__main__":
    results = pd.DataFrame(columns=['method', 'dataset', 'mem'])
    methods = ["BSX", "VEQ", "GUP", "KSS", "BICE", "RM"]
    datasets = ["citeseer", "dblp", "youtube", "HPRD", "human", "maayan-figeys", "twitch", "web-Stanford", "wordnet-words", "YeastS"]
    dataspace = "/var/lib/docker/subgraph/output/bsx/mem/"

    idx = 0
    for method in methods:
        for dataset in datasets:
            datapath = os.path.join(dataspace, method, dataset)
            
            # Find the maximum max_memory value for all files in this directory.
            max_memory = find_max_memory_in_directory(datapath, method)

            results.loc[idx] = [method, dataset, max_memory]
            idx += 1
            print(f"The maximum max_memory value for {method} in {dataset} is: {max_memory} KB")

    results.to_csv('mem.csv', index=False)
