## Lampropeltis getula genomes

### Sean Harrington


<br>

This repository contains code for whole genome analyses of Lampropeltis getula, primary to identify selection across the genome. This is all a follow up to RAD data analysis performed in [Harrington and Burbrink](https://onlinelibrary.wiley.com/doi/full/10.1111/jbi.14536) looking at lineage structure across this complex.

The vast majority of these are slurm scripts to be run on the Beartooth cluster at the University of Wyoming. You will need to adapt these to your own computational architecture if you wish to run them.


<br>
<br>

### Raw data quality assessment and trimming

1. `fastqc_array.slurm` runs fastqc to look at quality of the raw reads

2. `trim_array.slurm` runs Trimmomatic to trim reads for quality and to remove adapters

3. `fastqc_array_bbduck_trims.slurm` runs fastqc on reads after trimming. We used the adapter file from bbduk, which is where the name of this file comes from.


<br>

### Mapping, duplicate removal

1. `BWA_to_ratsnake.slurm` runs BWA to map the trimmed reads to the Pantherophis genome from HERE.  **DROP IN LINK TO RATSNAKE GENOME**

2. `idx_ratmapped_bam.slurm` indexes each of the .bam files created by the previous script.

3. `flagstat_ratmapped.slurm` runs `samtools flagstat` on each bam file to get some mapping stats.


4. `picard_rmd_ratmapped.slurm` runs `picard MarkDuplicates` to identify and remove duplicate reads.


<br>

### Variant calling and filtering

1. `var_call_ratmapped_all.slurm` calls variants from the bam files with duplicates removed. This is done for each scaffold of the reference genome separately here as a slurm job array.

2. `sort_ratmapped_all_vcfs.slurm` sorts each of the resulting vcf files generated in the previous step.

3. `combine_ratmapAll_vcf.slurm` combines the individual scaffold vcf files into a single vcf of the whole genome. 

4. `get_vcf_stats_ratmapAll.slurm` gets some statistics about the vcf file of all scaffolds combined. `r_vcf_stats_ratmapAllFilt.R` will then create plots of those stats.

5. Filter vcf files: `filter_vcf_ratmapAll.slurm` filters the combined vcf of all scaffolds. `filter_each_vcf_ratmapAll.slurm` filters each of the scaffold vcf files using the same filters.

6. `get_vcf_stats_ratmap_filtered.slurm` calculates vcf stats of the filtered whole genome vcf and then runs `r_vcf_stats_cmdline.R` to make plots.


- **add details of filtering parameters**

<br>

### Pixy - get Fst, dxy, and pi across the genome Pixy requires a vcf file with all sites to run, need to generate this.

1. `var_call_allsites.slurm` calls variants from the bam files with duplicates removed, generating a vcf with all sites, including invariants. This is done for each scaffold of the reference genome separately here as a slurm job array.

2. `sort_allsites_NC_vcfs.slurm` sorts each vcf. After running this script, delete unsorted vcf files manually. I removed this from the script in case part of it fails.

3. `filter_each_vcf_allsites.slurm` filters all of the vcf files and tabix index each.

4. `Combine_NC_vcf_get_stats_allsites.slurm` combines the pre-filtered and post-filtered vcfs from NC scaffolds and gets some statistics about the vcf files and then runs `r_vcf_stats_cmdline.R` to make plots.


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




## To do
- Check that I have stats of the filtered vcf - pretty sure I do
- double check that the vcf is the way I want it
	- remove any individuals with high missing data?
- plot pixy - modify R code on the website


## Notes to me:
`filter_each_vcf_allsites_TEST.slurm` tests to make sure that including min alleles as zero doesn't get rid of things I want to keep by removing that argument altogether - it seems that it does not matter, file sizes are identical.





