# Aim of the pipeline

1. intersect pairwise genome alignment results in UCSC chain format (obtained for instance using cactus: https://github.com/ComparativeGenomicsToolkit/cactus) for different species, two being typically refs and others outgroups, in order to obtain common sections
2. get ancestral genome for refs using outgroups
3. analyse pattern of mutations in regions of interest, here NIEBs.

# Dependencies

- snakemake
- bedtools
- R
- Python 3
- cmake
- compiler compatible with C++ 17

# Usage

Clone the git repo and then issue the snakemake command. You *must* specify number of cores to use through the `-nx` (with `x` being the quantity of cores) or `--cores x`. The pipeline has only like 8 steps that can be launched in parallel so no need to specify more (they won't be used).

You can get to a particular step of the pipeline by calling its name or its output : for instance :

- if you want the bed with all the common regions between outgroups, you need to call `snakemake -n1 data/intersected.bed`
- for the ancestral genomes just call the full pipeline.
