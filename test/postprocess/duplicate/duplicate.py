from datasketch import MinHash, MinHashLSH
import sys
import numpy as np
import os


# Function to read sets from a file
def read_sets_from_file(file_path):
  with open(file_path, 'r', encoding='utf8') as f:
    data = f.read()
  set_data = data.split('------')
  sets = [set(s.split()) for s in set_data if s.strip()]
  return sets

def process_query(query_path, qsize):
  threshold = 0.98
  results = []
  for file_idx in range(qsize-1):
    # Read sets from the file
    file_path = query_path+'/'+str(file_idx)+'.txt'
    sets = read_sets_from_file(file_path)
    if len(sets) == 0:
      continue
    if len(sets) > 100000:
      return []

    # Create MinHash objects
    minhashes = []
    for s in sets:
      m = MinHash(num_perm=128)
      for d in s:
        m.update(d.encode('utf8'))
      minhashes.append(m)

    # Create LSH index
    lshs = []
    lsh = MinHashLSH(threshold=threshold, num_perm=128)
    lsh.insert("m0", minhashes[0])
    lshs.append(lsh)

    for idx, minhash in enumerate(minhashes[1:], start=1):  # start=2 to name as "m2", "m3", etc.
      lshs_num = len(lshs)
      build_lsh = True
      for i in range(lshs_num):
        result = lshs[i].query(minhash)
        if len(result) > 0:
          # if use the first added minHash as representative, comment next line
          lshs[i].insert(f"m{idx}", minhash)  # cumulate minHash to each lsh
          build_lsh = False
          break
      if build_lsh:
        lsh = MinHashLSH(threshold=threshold, num_perm=128)
        lsh.insert(f"m{idx}", minhash)
        lshs.append(lsh)
    results.append([len(sets), len(lshs), len(lshs)/len(sets)])

  return results

if __name__ == '__main__':
  data_path = str(sys.argv[1])
  graph     = str(sys.argv[2])
  method    = str(sys.argv[3])
  i         = int(sys.argv[4])

  files_path = os.path.join(data_path, graph, method)

  qsize = 10
  query_path = os.path.join(files_path, f"Q{qsize}-{i+1}")
  print(query_path)
  results = process_query(query_path, 10)
  if len(results) == 0:
    exit()

  array = np.array(results)
  np.savetxt(query_path+'/array.txt', array, fmt='%.6f')
