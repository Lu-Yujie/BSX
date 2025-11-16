import os
import csv
import sys

if __name__ == '__main__':
    """
    Write all experimental results under a single data-graph folder into one file,
    and add corresponding information.
    Input is the data folder path; experiment data under the path are organized as Qxx/xx.csv.
    """
    folder_path = sys.argv[1]
    output_file = os.path.join(folder_path, 'result.csv')
    # Get all subfolders that start with 'Q' in the folder
    subfolders = [name for name in os.listdir(folder_path)
                  if os.path.isdir(os.path.join(folder_path, name)) and name.startswith('Q')]

    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        # Example of a more detailed header (commented out)
        writer.writerow(["size_", "id", "unsolved", "filter_", "order", "engine", "total_time", "result_num"])

        # Iterate through each subfolder
        for subfolder in subfolders:
            subfolder_path = os.path.join(folder_path, subfolder)
            csv_files = [file for file in os.listdir(subfolder_path) if file.endswith('.csv')]
            print(csv_files)

            # Iterate through each .csv file
            for csv_file in csv_files:
                csv_file_path = os.path.join(subfolder_path, csv_file)

                # Read the .csv file content
                with open(csv_file_path, 'r') as file:
                    reader = csv.reader(file)
                    data_started = False
                    graph_id = ""

                    # Iterate through each row
                    for row in reader:
                        if len(row) <= 0:
                            continue
                        # Skip lines that start with Data, load, end, or Query
                        if row[0].startswith('Data') or row[0].startswith('load') or row[0].startswith('end') \
                           or row[0].startswith('Query'):
                            continue
                        else:
                            # Extract id part from the first column and insert folder size prefix
                            row[0] = row[0].split('.')[0].split('-')[-1]
                            row.insert(0, subfolder[1:])
                            writer.writerow(row)
