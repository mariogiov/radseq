radseq
======

Scripts and wrappers for radtags code to deal with SciLifeLab directory structures.

## Snakemake pipeline installation
create a new conda environment with python3 and install snakemake:

```
conda create -n snakemake python=3.4.2 pip
source activate snakemake
pip install snakemake
```

clone mario's radseq repo and find the Snakefile:

```
cd && git clone https://github.com/mariogiov/radseq.git
```

copy the Snakefile to whichever directory you're processing:

```
cp ~/radseq/snakemake/Snakefile <wherever>
```

## usage
just type 'snakemake' in the same directory as the Snakefile:

```
snakemake
```
