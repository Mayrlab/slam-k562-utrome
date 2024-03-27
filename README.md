## slam-k562-utrome
This repository provides a pipeline for reprocessing the SLAM-seq dataset from [Wu et al., 2019](https://doi.org/10.7554/eLife.45396) at 3'UTR isoform resolution using [the human UTRome annotation](https://github.com/Mayrlab/hcl-utrome) developed by Fansler et al. 2023.

## Reproducing the Pipeline
### Cloning
This repository can be cloned with:

```bash
git clone https://github.com/Mayrlab/slam-k562-utrome.git
```

### Prerequisite Software
This requires Conda/Mamba and Snakemake. If you do not already have a Conda installation, we strongly recommend [Miniforge](https://github.com/conda-forge/miniforge#miniforge). We provide a YAML for the Snakemake version we used (v6.8), but feel free to use later versions.

### Configuration
Three configuration options in `config.yaml` should be adjusted by the user prior to running:

  - `tmp_dir`: temporary directory for scratch
  - `genome_fa`: FASTA file of hg38 genome
  - `utrome_gtf`: GTF file of human UTRome (obtainable through [`scUTRquant`](https://github.com/Mayrlab/scUTRquant/tree/e793dcbce24adad832bca42901b11cd9bd66ebe9/extdata/targets/utrome_hg38_v1))

### Running
The full pipeline can be executed with simply

```bash
snakemake --use-conda
```

We encourage HPC users to configure [a Snakemake profile](https://github.com/Snakemake-Profiles) and use this via a `--profile` argument.
