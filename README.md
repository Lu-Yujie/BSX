# Batch-Backtracking Search

Subgraph matching is a fundamental problem in graph analysis. Recently, many algorithms have been developed, often based on classic backtracking search. Traditional **b**acktracking **s**earch matches **1** vertex at a time (denoted B $\circledS$ 1), which can cause redundant computations due to overlapping search spaces. To address this, we propose a novel *batch-**b**acktracking **s**earch* framework that allows matching a set of data vertices **$X$** (denoted B $\circledS X$) in each backtracking step. B $\circledS X$ models the search space as a “search box,” enabling flexible exploration of the search space and substantially reducing overlap between search subspaces. It selects batches to group data vertices that have similar search spaces. For each search box, we introduce a refinement procedure to filter out unpromising candidate mappings. In addition, we propose a homomorphism termination to stop the backtracking process as early as possible and an efficient embedding enumeration method to list all embeddings within a search box simultaneously. Extensive experiments on real-world graphs show that B $\circledS X$ significantly outperforms existing state-of-the-art algorithms, achieving speedups up to nine orders of magnitude under the EPS metric.

## Compile

From the project root directory, run the following commands to compile the source code. The compiled executable is `build/bin/BS`.

```zsh
mkdir build
cd build
cmake ..
make
```

## Run

To improve experimental efficiency, we implemented multithreaded experiment scripts in C++ and placed them under `test/run`. Before releasing the code, I tested every script on our server; they all run correctly there. I kept our file paths in the scripts—if you want to reproduce our results, you only need to replace the paths.

### Preliminary

* First run `bash test/run/create_out_dir.sh [output_path]` to create the output directory.
* Then edit the experiment configuration files under `test/conf` to choose the filters, order, and engine for the current run. Their usage is described in `test/conf/readme.md`.
* Note that our code integrates VEQ, RM, KSS, BS1, and BSX into a single framework. The other two comparison methods in the paper, GUP and BICE, are not supported by this framework; you should download their original open-source code and run them separately, then convert their outputs to the format our code expects. To help reproduction, we provide processed experimental results in `test/draw` in CSV format.

### Overall (Fig. 10, 14, 16, Table 3)

We provide a relatively simple parallel script `overall.sh` to compute the overall experiment results used for Fig. 10, 13, 14, and 16. All subsequent `run` scripts are derived from this script. For clarity, parameter descriptions are included in the script as comments.
To increase parallelism, we implemented the experiment scripts in C++; this was a matter of familiarity, and users may port them to other languages (bash/python/…) if preferred.

For peak memory consumption experiments, enable *ANALYZE_MEMORY* in `configuration/config.h`, then rerun the experiments. Use `test/postprocess/mem/cnt_mem.py` to calculate the final results.

**conf**:

* `excluded.txt`: remove *BSX* from *engine*, remove *CLF*, *DPiso*, *NLF* from *filter*.
* `wanted.txt`: keep *RMRMRM* and *DPisoDPisoKSS* (each on separate lines).
* `src.txt`: keep *VEQ*.

For the Peak Memory Consumption experiment, keep only the corresponding combinations in `conf`:

* VEQ: keep *VEQ* in `src.txt`.
* BSX: clear `wanted.txt` and `src.txt`, and remove *CFL* from `filter` and *BSX* from `engine` in `excluded.txt`.
* KSS: clear `src.txt`, keep *DPisoDPisoKSS* in `wanted.txt`.
* RM: clear `src.txt`, keep *RMRMRM* in `wanted.txt`.
* BICE: we embedded the same statistics functions into its source code and then re-ran it.
* GUP: since its implementation language differs, we used the `time` tool to run its open-source code and then converted the results to the output structure used by our framework.

### Global (Fig. 11)

Use `global.cpp`. Based on **Overall**, we change `query_num` to 8; query IDs are 1–8.

**conf** is the same as **Overall**, except remove *RMRMRM* from `wanted.txt`.

### Degree / Density (Fig. 12)

Use `degree.cpp`. On top of **Overall**, add a degree parameter with values in the range 1–5. We removed query size 10 because for query size 10 it is hard to find subgraphs with large degree variation.

**conf** is the same as **Overall**.

### Label Size (Fig. 13)

Use `ratio.cpp`. Replace the `label_size` parameter on top of **Overall**. The tested values are: 10, 50, 100, 200, $10^3$, $10^4$, $10^5$.

**conf**:

* For `label_size <= 1000`, use the same configuration as **Overall**.
* For `label_size > 10000`, remove VEQ from `src.txt`.

### Data Graph Size (Fig. 15)

Use `scale.cpp`. This experiment uses EvoGraph to scale up real graphs: `"human"`, `"HPRD"`, `"citeseer"`, `"YeastS"`, `"Figeys"`, and adds a `scale` variable.

**conf** is the same as **Overall**.

### Similarity & Refinement (Fig. 17)

This group of experiments requires merging two patches from the master branch into the `multiple` branch first (you will need to resolve merge conflicts). `similar.cpp` corresponds to the similarity experiments and adds a `similar_ratio` parameter with values `{0, 25, 50, 75, 100}`. `refine.cpp` corresponds to refinement experiments and adds a `refinement_ratio` parameter with values `{0, 25, 50, 75, 100}`. After specifying the output path, also modify the macro definitions in `configuration/config.h` to match the parameters.

For convenience, you can switch to the `similar` and `refinement` branches of the code directly—these are the versions used for our experiments.

**conf**:

* `excluded.txt`: remove *BSX* from *engine*, remove *CLF* from *filter*.
* `wanted.txt`: clear.
* `src.txt`: clear.

### Overlap Ratio (Table 4)

* First, use the master branch code. Enable **ANALYZE_DUPLICATE** in `configuration/config.h` (remember to recompile after switching branches).
* Run the experiments using `run/duplicate.sh`. The script will print the path where results are stored.
* After generating experimental results, compile and run the executable built from `postprocess/duplicate/process_dup.cpp` to process the data and compute numbers of **buckets** and **boxes**.
* Finally, run `postprocess/duplicate/avg.py` to compute average results and obtain the numbers shown in the paper.

## Postprocess

After producing experiment outputs, use our postprocessing scripts under `test/postprocess` to structure the results into formatted files. The data-processing scripts correspond to the run scripts (by name) and are straightforward to execute; each script contains comments describing required parameters. Typically, each experiment setting includes two basic scripts:

* `combine.py`: merges results from different query sizes into a single CSV.
* `combine_all.sh`: calls `combine.py` to produce the merged file.
* `avg.py`: averages the merged experimental data to produce plotting inputs.

Different experiment settings may require small, obvious adjustments (e.g., adding specific paths or adjusting column names).

### Peak Memory Consumption

Use `test/postprocess/mem/cnt_mem.py` to aggregate memory-use results. This script will search for the relevant memory output files and produce the summarized data.

## Plots

We provide the code to draw all figures shown in the paper.

### Scripts for the paper’s figures

|            Script            |                            Figure                           |
| :--------------------------: | :---------------------------------------------------------: |
|       `eps/draw_eps.py`      |        Fig. 10. Relative EPS ratio (VEQ as baseline).       |
|    `global/draw_global.py`   |     Fig. 11. QPT of different methods with output limit.    |
|    `degree/draw_degree.py`   | Fig. 12. Performance under different query graph densities. |
|     `ratio/draw_ratio.py`    |          Fig. 13. Effect of data graph label size.          |
|     `qsize/draw_qsize.py`    |             Fig. 14. Effect of query graph size.            |
|     `scale/draw_scale.py`    |             Fig. 15. Effect of data graph size.             |
| `filter/draw_filter_base.py` |        Fig. 16. Relative EPS Ratio (DPiso-BⓈ𝑋 = 1).        |
|  `flexible/draw_flexible.py` |   Fig. 17. Performance of BⓈ𝑋 with different parameters.   |

Notes:

* Only minor modifications to the data-processing scripts are typically needed for different experiment settings, such as adding specific paths or adjusting column names.
* The flexible experiments were conducted after applying the two patches included in the supplement folder.
