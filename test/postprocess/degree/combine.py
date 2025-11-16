import os
import csv
import sys

if __name__ == '__main__':
  folder_path = sys.argv[1]
  output_file = folder_path + '/result.csv'
  dFolders = [name for name in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, name)) and name.startswith('D')]

  with open(output_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["size_", "degree", "id", "unsolved", "filter_", "order", "engine", "total_time", "result_num"])

    for dFolder in dFolders:
      dFolder_path = os.path.join(folder_path, dFolder)
      qFolders = [name for name in os.listdir(dFolder_path) if os.path.isdir(os.path.join(dFolder_path, name)) and name.startswith('Q')]
      for qFolder in qFolders:
        qFolder_path = os.path.join(dFolder_path, qFolder)
        csv_files = [file for file in os.listdir(qFolder_path) if file.endswith('.csv')]
        print(csv_files)

        for csv_file in csv_files:
          csv_file_path = os.path.join(qFolder_path, csv_file)

          with open(csv_file_path, 'r') as file:
            reader = csv.reader(file)
            data_started = False
            graph_id = ""

            for row in reader:
              if len(row) <= 0:
                continue
              if row[0].startswith('Data') or row[0].startswith('load') or row[0].startswith('end') \
                   or row[0].startswith('Query'):
                continue
              else:
                row[0] = row[0].split('.')[0]
                row.insert(0,dFolder[1:])
                row.insert(0,qFolder[1:])
                writer.writerow(row)
