# Batch-Backtracking Search

Subgraph matching is a fundamental problem in graph analysis. Recently, many algorithms have been developed, often using classic backtracking search. This traditional **b**acktracking **s**earch matches **1** vertex at a time, denoted as B $\circledS$ 1, which can lead to redundant computations due to overlapping search spaces. To address this problem, we propose a novel _batch-**b**acktracking **s**earch_ framework that enables matching a set of data vertices **$X$**, denoted as B $\circledS X$, in each backtracking step. B $\circledS X$ models the search space as a ``search box'', allowing for flexible search space exploration and significantly minimizing the overlap between search space. It effectively selects batches to cluster data vertices with similar search spaces. For each search box, we introduce a refinement method to filter out unpromising candidate mappings. Furthermore, we propose a homomorphism termination to break the backtracking process as early as possible and an efficient embedding enumeration method to list all embeddings within the search box simultaneously. Extensive experiments on real-world graphs demonstrate that B $\circledS X$ significantly outperforms existing state-of-the-art algorithms, achieving a speedup of up to 9 orders of magnitude under the EPS metric.

## Compile

Under the root directory of the project, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```

## Correctness Verification

We provide 200 test cases along with the corresponding test script. The usage is as follows:

```bash
python test.py ../build/matching/BS
```

If all 200 cases pass correctly, the following text will be displayed:

```bash
{your_method_name} engine passed the correctness check.
```

Additionally, we provide a script(_check_result.py_) for comparing results between different methods as an auxiliary tool for verifying the correctness of the code. Please refer to the script comments for specific usage instructions.

## Execute

After compiling the source code, you can find the binary file 'BS' under the 'build/matching' directory.
Execute the binary with the following command ./BS -d data_graph -q query_graph
-filter filter_technique -order order_technique -engine engine_technique -num max_number_of_embeddings -time_limit max_execute_time,
in which -d specifies the input of the data graphs and -q specifies the input of the query graphs.
The -filter parameter gives the filtering method, the -order specifies the ordering method, and the -engine
sets the enumeration method. The -num parameter sets the maximum number of embeddings that you would like to find. The time_time constrains the maximum execution time. If the number of embeddings enumerated reaches the limit or all the results have been found or the time limit is reached, then the program will terminate.
Set -num as 'MAX' to find all results.

Example (Use the filtering method of CFL and order method of GraphQL to generate the candidate vertex sets and the matching order respectively.
Enumerate results with the set-intersection based local candidate computation method):

```zsh
./BS -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter CFL -order GQL -engine BS1 -num MAX
```

## Input

Both the input query graph and data graph are vertex-labeled.
Each graph starts with 't N M' where N is the number of vertices and M is the number of edges. A vertex and an edge are formatted
as 'v VertexID LabelId Degree' and 'e VertexId VertexId' respectively. Note that we require that the vertex
id is started from 0 and the range is [0,N - 1] where V is the vertex set. The following
is an input sample. You can also find sample data sets and query sets under the test folder.

Example:

```zsh
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1
e 0 2
e 1 2
e 1 3
e 2 4
e 3 4
```

## Techniques Supported

The filtering methods that generate candidate vertex sets.

|Parameter of Command Line (-filter) | Description |
| :-----------------------------------: | :-------------: |
|NLF| the neighborhood label frequency filter |
|CFL| the filtering method of CFL|
|DPiso| the filtering method of DP-iso |

The ordering methods that generate matching order.

|Parameter of Command Line (-order) | Description |
| :-----------------------------------: | :-------------: |
|GQL| the ordering method of GraphQL |

The enumeration methods that find all results.

|Parameter of Command Line (-engine) | Description |
| :-----------------------------------: | :-------------: |
|BS1| Naive-Backtracking Search |
|BSX| Batch-Backtracking Search |

## Generalizability

BSX supports different tightness ratios for refinement and varying similarity ratios for candidates within a batch. In the _supplement_ directory, there are two patches that provide support for these extensions.

## Experiment Datasets

We have placed all the datasets used for testing in the paper at this link: [dataset_BSX]([https://github.com/Lu-Yujie/BSX_dataset](https://anonymous.4open.science/r/BSX_dataset-F0E2)).

```bash
tar -xJf dataset_BSX.tar.xz
```
