# PopPIPE: Population analysis PIPEline 🛠🧬

[![Build and push Docker image](https://github.com/johnlees/PopPIPE/workflows/Build%20and%20push%20Docker%20image/badge.svg?branch=master)](https://hub.docker.com/repository/docker/poppunk/poppipe)

Downstream analysis of [PopPUNK](https://www.poppunk.net/) results. Produces subclusters and visualisations of all strains.

Further documentation can be found in the [PopPUNK docs](https://poppunk.readthedocs.io/en/latest/subclustering.html).

## Pipeline description

The pipeline consists of the following steps:
- Split files into their [PopPUNK](https://www.poppunk.net/) strains.
- Use [pp-sketchlib](https://github.com/bacpop/pp-sketchlib) to calculate core and accessory distances within each strain.
- Use core distances and [rapidnj](https://birc.au.dk/software/rapidnj/) to make a neighbour-joining tree.
- (lineage_clust mode) Generate clusters from core distances with lineage clustering in PopPUNK.
- Use [ska](https://github.com/simonrharris/SKA) to generate within-strain alignments.
- Use [IQ-TREE](http://www.iqtree.org/) to generate an ML phylogeny using this alignment, and the NJ tree as a starting point.
- Use [fastbaps](https://github.com/gtonkinhill/fastbaps) to generate subclusters which are partitions of the phylogeny.
- Create an overall visualisation with both core and accessory distances, as in PopPUNK. The final tree consists of refining the NJ tree by grafting the maximum likelihood trees for subclusters to their matching nodes.

### Example pipeline DAG

![pipeline dag](snakemake_dag.png)

`ska index` steps (one per sample, and a dependency of all `ska align` steps) have been omitted for simplicity.

## Installation

The supported method is to use conda, which is most easily accessed by first
installing [miniconda](https://conda.io/miniconda.html). PopPIPE simply depends
upon snakemake and pandas:
```
conda install snakemake pandas
```

Other dependencies will be automatically installed by conda the first time
you run the pipeline. You can also install them yourself and omit the `--use-conda`
directive to snakemake:
```
conda create -n poppipe --file=environment.yml
```

If the package cannot be found you will need to add the necessary channels:
```
conda config --add channels r
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Running inside a container

An alternative, if you are having trouble with the above, is to use the PopPIPE docker
container. If you are comfortable running commands inside docker containers and mounting
your external files, the whole pipeline is in the container available by running:
```
docker pull poppunk/poppipe:latest
```

You can also follow the above process and make a local clone of snakemake, and replace
`--use-conda` with `--use-singularity`, which will automatically pull this container, and run
each step inside it.

Use `--singularity-args` if you need to bind directories.
## Usage

1. Modify `config.yml` as appropriate.
2. Run `snakemake --cores <n_cores> --use-conda`.

In particular, check the three `poppunk_` arguments, which should be set to the
full path of the `--r-files` argument, strain clusters .csv file and `.h5` database file,
from your PopPUNK run.

On a cluster or the cloud, you can use snakemake's built-in `--cluster` argument:
```
snakemake --cluster qsub -j 16 --use-conda
```
See the [snakemake docs](https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html)
for more information on your cluster/cloud provider.

### Creating visualisations
The default target is the first in the Snakefile: `cluster_summary`. This
produces subclusters but no visualisations. To continue the run forward from here
use the `microreact` target:
```
snakemake --use-conda make_microreact
```

This will create a phylogeny, embedding and format your strains and their subclusters
for [microreact](https://microreact.org/) and save these files to the output. The phylogeny
and clusters will be sent to microreact, and a link to your page will be output to the terminal
and saved in `output/microreact_url.txt`.

**NB** From 2021-10-27 Microreact requires an API key for the final step to work. See the
[microreact docs](https://docs.microreact.org/api/access-tokens) for instructions on how to generate one for your account.

## Config file

### PopPIPE configuration

* `script_location`: The `scripts/` directory, if not running from the root of this repository
* `poppunk_db`: The PopPUNK HDF5 database file, without the `.h5` suffix.
* `poppunk_clusters`: The PopPUNK cluster CSV file, usually `poppunk_db/poppunk_db_clusters.csv`.
* `poppunk_rfile`: The `--rfile` used with PopPUNK, which lists sample names and files, one per line, tab separated.
* `min_cluster_size`: The minimum size of a cluster to run the analysis on (recommended at least 6).

### SKA configuration

* `fastq_qual`: With read input, the `-q` option, which ignores k-mers with bases below this score.
* `fastq_cov`: With read input, the `-c` option, which sets a minimum k-mer count.

### IQ-TREE configuration

* `enabled`: Set to `false` to turn off ML tree generation, and use the NJ tree throughout.
* `mode`: Set to `full` to run with the specified model, set to `fast` to run using `--fast` (like fasttree).
* `model`: A string for the `-m` parameter describing the model. Adding `+ASC` is recommended.

### fastbaps configuration

* `levels`: Number of levels of recursive subclustering.
* `script`: Location of the `run_fastbaps` script. Find by running `system.file("run_fastbaps", package = "fastbaps")` in R.

### mandrake configuration

* `perplexity`: Perplexity parameter for t-SNE (between 5 and 50).
* `knn`: Number of nearest neighbours (at least two).
* `maxIter`: Iterations in the optimisation (at least 10000, default 100000).

### Microreact configuration

* `name`: Title of the Microreact to produce
* `website`: Website link to give in Microreach
* `email`: Contact email to list in Microreact
* `api_token`: The API token from your Microreact account

## Updating a run with results from poppunk_assign

You can use the helper script `poppipe_assign.py` to help you re-run after query assignment. For example, if you assigned to a database with:

```
poppunk_assign --db listeria_rlist --query qlist.txt --output listeria_qlist
```

Give the same arguments, and your config file to `python poppipe_assign.py`:

```
python poppipe_assign.py --db listeria_rlist  --query qlist.txt --output listeria_qlist --config config.yml
```

This will generate combined input files, and a new config file. Run snakemake again:

```
snakemake --configfile configv2rfmbr9.yml
```

Here, the first snakemake pipeline ran on four strains consisting of 83 samples:
```
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 4
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	cluster_summary
	4	fastbaps
	4	generate_nj
	4	iq_tree
	4	ska_align
	58	ska_index
	4	sketchlib_dists
	4	split_strains
	83
```

The second one, with the new queries, had one more strain, and nineteen new samples to index:
```
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	cluster_summary
	5	fastbaps
	5	generate_nj
	5	iq_tree
	5	ska_align
	19	ska_index
	5	sketchlib_dists
	5	split_strains
	50
```

**NB**: This will re-run all downstream steps for each strain, other than the `ska index` steps. If you have a small number of strains being changed this is likely to be inefficient. If you would like us to support this type of analysis please get in touch.
