import numpy as np
import os

num_files = 100
qsize = 10
method = "bsx" # "bs1"
graph = "citeseer" # "Figeys" "YeastS"
src_path = "/var/lib/docker/subgraph/output/bsx/duplicate/"+graph+"/"+ method+"/"
file_paths = []

for i in range(100):
    file_paths.append(src_path + "Q" + str(qsize) + "-" + str(i + 1) + "/"+"array.txt")

data_num = [0 for _ in range(qsize)]
data_sum = [[0, 0, 0] for _ in range(qsize)]
print(data_sum)
data_mean = [[0, 0, 0] for _ in range(qsize)]
for file_path in file_paths:
    if not os.path.exists(file_path):
      continue
    data = np.loadtxt(file_path)
    if data.shape == (3,):
      data = data.reshape(1,3)
    for i in range(data.shape[0]):
      for j in range(3):
        data_sum[i][j] += data[i][j]
      data_num[i] += 1

print(data_sum)
print(data_num)

for i in range(qsize):
  if (data_num[i] == 0):
    continue
  for j in range(2):
    data_mean[i][j] = data_sum[i][j] / data_num[i]
  data_mean[i][2] = 1 - (data_sum[i][2] / data_num[i])

# 输出结果
np.savetxt(method+".txt", data_mean, fmt='%.6f')
