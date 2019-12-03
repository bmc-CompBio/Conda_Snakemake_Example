# Conda Snakemake Example

Demonstration of Conda package management and Snakemake pipeline execution on HPC cluster.
This pipeline will should also run on your workstation.

Outline: One sample ChIP-Seq of MSL2 in Drosophila cell lines. Fetch data from SRA. Align the reads, call peaks, create a bigwig file for browsing and generate an HTML report from a Rmd script.

### Prerequisites

In your cluster project directory install Miniconda as described here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

```
bash Miniconda3-latest-Linux-x86_64.sh
```

### Installing

Clone this repository in your cluster project directory

```
git clone https://github.com/bmc-CompBio/Conda_Snakemake_Example.git
```

Create the conda environment and activate it

```
conda env create --name chipseq --file environment.yaml
conda activate chipseq
```

Install a github R package:
with the chipseq environment activated run R

```
R
```

within R execute

```
library(devtools)
install_github("musikutiv/tsTools")
```

and quit the R session

```
q()
```

## Running the pipeline

Check the pipeline by running
```
snakemake -np
snakemake --dag | dot -Tsvg > dag.svg
```

Run the pipeline

```
snakemake -j 999 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} --mem {cluster.mem}"
```
