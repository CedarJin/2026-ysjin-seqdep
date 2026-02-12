# 2026-ysjin-seqdep
Microbiome metagenome and metatranscriptome sequencing depth

## Test workflow using SRA data
### Set up sratoolkit

Please refer to https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit to download and install sratoolkit. 


### Download and Downsample the SRA data to 10M, 20M, 30M, 40M, 50M

Use the following code to download SRA data:
```bash
uv run ./download_test_metagenome.sh # alternatively, use the following snakemake workflow. 
```

Snakemake workflow:
```bash
module load seqtk
uv run snakemake -s testdata_download_downsample.smk -n # -n : Dry run, just test if DAG can be build.
uv run snakemake -s testdata_download_downsample.smk -j 1 # -j : specify the # of cores used in the workflow. This workflow is memory intensive and time-consuming. You can submit a sbatch job and see if 4 cores and 64G mem are faster.
```
Submit a batch job:
```bash
sbatch run_testdata_downsample.sbatch
squeue -u $USER # see the progress
scancel <jobid> # cancel if needed # or use scancel -u $USER
```
Once downloading and downsampling are finished, please head several files to double check the read pairs, see if the # are correct and if they are matched.


