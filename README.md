## Lampropeltis getula genomes

### Sean Harrington


<br>

This repository contains code for whole genome analyses of Lampropeltis getula, primary to identify selection across the genome. This is all a follow up to RAD data analysis performed in [Harrington and Burbrink](https://onlinelibrary.wiley.com/doi/full/10.1111/jbi.14536) looking at lineage structure across this complex.

The vast majority of these are slurm scripts to be run on the Beartooth cluster at the University of Wyoming. You will need to adapt these to your own computational architecture if you wish to run them.

* Side note to myself - some paths in scripts might be slightly broken since reorganizing this directory into subdirectories.




<br>
<br>

### Raw data quality assessment and trimming

* `qc_trim` directory

1. `fastqc_array.slurm` runs fastqc to look at quality of the raw reads

2. `trim_array.slurm` runs Trimmomatic to trim reads for quality and to remove adapters

3. `fastqc_array_bbduck_trims.slurm` runs fastqc on reads after trimming. We used the adapter file from bbduk, which is where the name of this file comes from.


<br>

### Mapping, duplicate removal

* `mapping` directory

1. `BWA_to_ratsnake.slurm` runs BWA to map the trimmed reads to the Pantherophis genome from HERE.  **DROP IN LINK TO RATSNAKE GENOME**

2. `idx_ratmapped_bam.slurm` indexes each of the .bam files created by the previous script.

3. `flagstat_ratmapped.slurm` runs `samtools flagstat` on each bam file to get some mapping stats.


4. `picard_rmd_ratmapped.slurm` runs `picard MarkDuplicates` to identify and remove duplicate reads.


<br>

### Variant calling and filtering

* `var_call` directory

1. `var_call_ratmapped_all.slurm` calls variants from the bam files with duplicates removed. This is done for each scaffold of the reference genome separately here as a slurm job array.

2. `sort_ratmapped_all_vcfs.slurm` sorts each of the resulting vcf files generated in the previous step.

3. `combine_ratmapAll_vcf.slurm` combines the individual scaffold vcf files into a single vcf of the whole genome. 

4. `get_vcf_stats_ratmapAll.slurm` gets some statistics about the vcf file of all scaffolds combined. `r_vcf_stats_ratmapAllFilt.R` will then create plots of those stats.

5. Filter vcf files: `filter_vcf_ratmapAll.slurm` filters the combined vcf of all scaffolds. `filter_each_vcf_ratmapAll.slurm` filters each of the scaffold vcf files using the same filters.

6. `get_vcf_stats_ratmap_filtered.slurm` calculates vcf stats of the filtered whole genome vcf and then runs `r_vcf_stats_cmdline.R` (in main `scripts` directory) to make plots.


- **add details of filtering parameters**

<br>

### Genome scans: use pixy to get Fst, dxy, and pi across the genome Pixy requires a vcf file with all sites to run, need to generate this.

* `gen_scan` directory

1. `var_call_allsites.slurm` calls variants from the bam files with duplicates removed, generating a vcf with all sites, including invariants. This is done for each scaffold of the reference genome separately here as a slurm job array.

2. `sort_allsites_NC_vcfs.slurm` sorts each vcf. After running this script, delete unsorted vcf files manually. I removed this from the script in case part of it fails.

3. `filter_each_vcf_allsites.slurm` filters all of the vcf files and tabix index each.

4. `Combine_NC_vcf_get_stats_allsites.slurm` combines the pre-filtered and post-filtered vcfs from NC scaffolds and gets some statistics about the vcf files and then runs `r_vcf_stats_cmdline.R` to make plots.

5. `Pixy.slurm` runs pixy as a job array across all NC scaffolds. `Pixy100_combined.slurm` does the same, but in 100 kb windows. `Pixy100_combined.slurm` runs pixy on the vcf containing the whole genome in a single file.

- `pixy_popfiles` directory contains populations files for pixy

6. `plot_pixy_cmd.R` plots pixy for scaffolds from the command line. `plot_pixy_combined_chroms.R` plots the combined genome pixy output.

For each scaffold, ran this interactively:

```
salloc -A getpop -t 0-03:00 --mem=1G --cpus-per-task=1
module load gcc/12.2.0 r/4.2.2
# ensure tidyverse is installed first
cd /project/getpop/scripts
Rscript plot_pixy_cmd.R /project/getpop/pixy_out100 /project/getpop/pixy_plots
```


<br>



### R analyses

* `r_analyses` directory

1. LFMM. Use LFMM in the LEA R package to find loci possibly under environmental selection. File processing: ran into issues with reading in my vcf files using vcf2lfmm. So now, first convert to `.ped` using plink. E.g.,

```
cd /project/getpop/vcf_rat_map_all
module load miniconda3
conda activate plink # my plink conda env
plink --vcf filtered_ratmapAll_NConly.vcf.gz --recode --out /project/getpop/lfmm/filtered_ratmapAll_NConly --allow-extra-chr
```

extract environmental data using `extract_enviro_data.R`






<br>
<br>
<br>


## tests/development:

- `testpixy.txt` code for testing pixy out interactively

Installing pixy required extra messing around in conda:

```{bash}
# standard pixy install
conda create -n pixy -y
conda activate pixy
conda install --yes -c conda-forge pixy python=3.8 -y
conda install --yes -c bioconda htslib -y

# install numpy and samtools with versions that work with pixy
conda install numpy=1.21 -y
conda install -c bioconda samtools=1.9 --force-reinstall -y
```

I believe that this is now fixed after I posted this on Github as an issue.


## To do
- Check that I have stats of the filtered vcf - pretty sure I do
- double check that the vcf is the way I want it
	- remove any individuals with high missing data?
- Pixy plots - get the pi for each pop in the all chromo plot


## Notes to me:
`filter_each_vcf_allsites_TEST.slurm` tests to make sure that including min alleles as zero doesn't get rid of things I want to keep by removing that argument altogether - it seems that it does not matter, file sizes are identical.

# plotting pixy output:

```
salloc -A getpop -t 0-03:00 --mem=1G --cpus-per-task=1
module load gcc/12.2.0 r/4.2.2
# start R and install tidyverse if necessary, then quit
cd /project/getpop/scripts
Rscript plot_pixy_cmd.R /project/getpop/pixy_out100 /project/getpop/pixy_plots
```





- need to tabix each file before pixy: format is: tabix <vcf_file> # samtools needs to be loaded up 


- what is the most appropriate window size for this??
- do we want to run things just on NC chromosomes? Or on all, including the really small ones

- allele surfing??    Pereira P., Teixeira J., Velo-Antón G. 2018. Allele surfing shaped the genetic structure of the European pond turtle via colonization and population expansion across the Iberian Peninsula from Africa. J. Biogeogr. 45:2202–2215.




