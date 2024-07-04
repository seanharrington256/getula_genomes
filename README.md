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



### Subsample reads for MVZ137799 - started at ~80x coverage, to get down to ~30x, subsample 180000000 reads

```
# some directory management
cd /project/getpop/trim_out_bbduk_adapters
mv trim_read_MVZ137799 trim_read_MVZ137799_not_downsampled
mkdir trim_read_MVZ137799

cd /project/getpop/trim_out_bbduk_adapters/trim_read_MVZ137799_not_downsampled

# FASTQ R1
seqtk sample -s 123 paired-MVZ137799_S64_L003_R1_001.fastq.gz 160000000 > ../trim_read_MVZ137799/paired-MVZ137799_S64_L003_R1_001.fastq
gzip ../trim_read_MVZ137799/paired-MVZ137799_S64_L003_R1_001.fastq

# FASTQ R2
seqtk sample -s 123 paired-MVZ137799_S64_L003_R2_001.fastq.gz 160000000 > ../trim_read_MVZ137799/paired-MVZ137799_S64_L003_R2_001.fastq
gzip ../trim_read_MVZ137799/paired-MVZ137799_S64_L003_R2_001.fastq

mv /project/getpop/trim_out_bbduk_adapters/trim_read_MVZ137799_not_downsampled /project/getpop/
```

<br>


### Mapping, duplicate removal

* `mapping` directory

1. `BWA_to_ratsnake.slurm` runs BWA to map the trimmed reads to the Pantherophis genome from HERE. Also Indexes the resulting bam files. **DROP IN LINK TO RATSNAKE GENOME**

2. `flagstat_ratmapped.slurm` runs `samtools flagstat` on each bam file to get some mapping stats.

3. `picard_rmd_ratmapped.slurm` runs `picard MarkDuplicates` to identify and remove duplicate reads.




<br>

### Variant calling and filtering

**Make sure to bgzip all files, otherwise tabix fails for pixy**

```
bgzip -c your_file.vcf > your_file.vcf.gz
```

* `var_call` directory

1. `var_call_ratmapped_all.slurm` calls variants from the bam files with duplicates removed. This is done for each scaffold of the reference genome separately here as a slurm job array.


**1 RERUNNING NOW^**


2. `sort_ratmapped_all_vcfs.slurm` sorts each of the resulting vcf files generated in the previous step.

3. `combine_ratmapAll_vcf.slurm` combines the individual scaffold vcf files into a single vcf of the whole genome. 

4. `get_vcf_stats_ratmapAll.slurm` gets some statistics about the vcf file of all scaffolds combined. `r_vcf_stats_ratmapAllFilt.R` will then create plots of those stats.

5. Filter vcf files: `filter_vcf_ratmapAll.slurm` filters the combined vcf of all scaffolds. `filter_each_vcf_ratmapAll.slurm` filters each of the scaffold vcf files using the same filters.

6. `get_vcf_stats_ratmap_filtered.slurm` calculates vcf stats of the filtered whole genome vcf and then runs `r_vcf_stats_cmdline.R` (in main `scripts` directory) to make plots.


- **add details of filtering parameters** May want to use qual 20?

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
plink --vcf filtered_ratmapAll_NConly.vcf.gz --recode12 --out /project/getpop/lfmm/filtered_ratmapAll_NConly --allow-extra-chr
```

extract environmental data using `extract_enviro_data.R`

* Testing:

to see how long LFMM will take, copied the `filtered_sort_rat_map_all_NC_045554.1_RagTag.vcf.gz` file into lfmm directory, filter to complete data with vcftools:

## note that this includes the wrong recode without 12

```
cd /project/getpop/lfmm
module load vcftools
vcftools --gzvcf filtered_sort_rat_map_all_NC_045554.1_RagTag.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out filtered_SNPs_no_missing_data

# then convert to plink
module load miniconda3
conda activate plink # my plink conda env
plink --vcf filtered_SNPs_no_missing_data.recode.vcf --recode --out /project/getpop/lfmm/no_missingNC_045554 --allow-extra-chr
```

or for the file with missing data:

```
cd /project/getpop/vcf_rat_map_all/filt_scafs_vcfs
plink --vcf filtered_sort_rat_map_all_NC_045557.1_RagTag.vcf.gz --recode12 --out /project/getpop/lfmm/rat_map_all_NC_045557 --allow-extra-chr
```

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

## may want to re-run variant calling, etc. with MVZ with huge num reads subsampled
- filter to biallelic only in filtering? thought I had this done, but LEA finds SNPs with a single state


- need to tabix each file before pixy: format is: tabix <vcf_file> # samtools needs to be loaded up 


- what is the most appropriate window size for this??
- do we want to run things just on NC chromosomes? Or on all, including the really small ones

- allele surfing??    Pereira P., Teixeira J., Velo-Antón G. 2018. Allele surfing shaped the genetic structure of the European pond turtle via colonization and population expansion across the Iberian Peninsula from Africa. J. Biogeogr. 45:2202–2215.


<br>
<br>
<br>
<br>
<br>



<br>


## Coverage issues

A bunch of individuals are like 8x. This is not ideal. Try out some different filtering on a single allsites scaffold:

```
cd /project/getpop/vcf_allsites
salloc -A inbreh -t 0-05:00 --mem=48G --cpus-per-task=1
module load gcc vcftools
```
- missing 20%
	all sites	
	`vcftools --gzvcf  sort_NC_045557.1_RagTag.vcf.gz --minDP 10 --max-missing 0.8 --max-alleles 2`
	kept 770,350 out of a possible 10637025 Sites
	Biallelic	
	`vcftools --gzvcf  sort_NC_045557.1_RagTag.vcf.gz --minDP 10 --max-missing 0.8 --max-alleles 2 --min-alleles 2`
	kept 101,750 out of a possible 10637025 Sites
	
	13% variable sites

- missing 20% but with depth cutoff of 8 instead of 10
	minDP 8  !!!!! a lower depth cutoff can keep WAY more sites
	all sites
	`vcftools --gzvcf  sort_NC_045557.1_RagTag.vcf.gz --minDP 8 --max-missing 0.8 --max-alleles 2`
	kept 3,433,603 out of a possible 10637025 Sites
	
	biallelic
	`vcftools --gzvcf  sort_NC_045557.1_RagTag.vcf.gz --minDP 8 --max-missing 0.8 --max-alleles 2 --min-alleles 2`
	kept 330,651 out of a possible 10637025 Sites	
	
	9% variable sites
		
		
- missing 20% but with depth cutoff of 7???
	minDP 7
	`vcftools --gzvcf  sort_NC_045557.1_RagTag.vcf.gz --minDP 7 --max-missing 0.8 --max-alleles 2`
	kept 5,569,858 out of a possible 10,637,025 Sites

	`vcftools --gzvcf  sort_NC_045557.1_RagTag.vcf.gz --minDP 7 --max-missing 0.8 --max-alleles 2 --min-alleles 2`
	kept 517,336 out of a possible 10637025 Sites	
	
	9.2% variable sites
	


- missing 30%, 10x
	all sites
	`vcftools --gzvcf  sort_NC_045557.1_RagTag.vcf.gz --minDP 10 --max-missing 0.7 --max-alleles 2`
	kept 2,357,569 out of a possible 10637025 Sites

	Biallelic
	`vcftools --gzvcf  sort_NC_045557.1_RagTag.vcf.gz --minDP 10 --max-missing 0.7 --max-alleles 2 --min-alleles 2`
	kept 243,355 out of a possible 10637025 Sites

	10% variable sites


- missing 30%, 8x
	all sites
	`vcftools --gzvcf  sort_NC_045557.1_RagTag.vcf.gz --minDP 8 --max-missing 0.7 --max-alleles 2`
	kept 6,204,385 out of a possible 10637025 Sites

	Biallelic
	`vcftools --gzvcf  sort_NC_045557.1_RagTag.vcf.gz --minDP 8 --max-missing 0.7 --max-alleles 2 --min-alleles 2`
	kept 578,425 out of a possible 10637025 Sites	

	9.3% variable sites











## Go back and map the reads to the Arizona genome and see if coverage was better


## Directory Arizona_mapping

### Mapping, duplicate removal

* `mapping` directory

1. `BWA_to_Arizona.slurm` runs BWA to map the trimmed reads to the *Arizona elegans* genome from HERE.  **DROP IN LINK TO GENOME**

2. `idx_arizona_map_bam.slurm` indexes each of the .bam files created by the previous script.

3. `flagstat_arizona_map.slurm` runs `samtools flagstat` on each bam file to get some mapping stats.

4. `picard_rmd_arizona_map.slurm` runs `picard MarkDuplicates` to identify and remove duplicate reads.


### Variant calling and filtering

* `var_call` directory

1. `var_call_arizona_mapAllSites.slurm` calls variants from the bam files with duplicates removed. This is done for each scaffold of the reference genome separately here as a slurm job array. Calls allsites to be used for pixy or filtered to only variant sites as desired.

2. `sort_arimappedAllSites_vcfs.slurm` sorts each of the resulting vcf files generated in the previous step.

3. `combine_arimapAllSites_vcf.slurm` combines the individual scaffold vcf files into a single vcf of the whole genome. 

4. `get_vcf_stats_arimapAllSites.slurm` gets some stats from the vcf file. I deleted some of the output from this that is on individual sites - files are huge

5. `filt_vcf_pixy.slurm` filters the vcf for use in Pixy (allsites vcf) based on coverage, missing data 20%, and quality. Removes indels. **RUNNING**

6. `filt_vcf_snps.slurm` filters the vcf down to just SNPs based on coverage, missing data 20%, and quality. Removes indels. **RUNNING**

7. `filt_vcf_snps_2_MAC.slurm` filters the vcf from 6 down to only sites with a minor allele count at least 2 for population structure stuff

8. `filt_vcf_snps_2MAC_thin.slurm` thins the vcf from 7 down to 1 snp per 10 kb window

9. `filt_vcf_snps_2MAC_thin1k.slurm` thins the vcf from 7 down to 1 snp per 1 kb window



- just to test out var calling without rmduplicates - very slightly higher coverage, but probably not good:

`var_call_arizona_mapAllSites_NORMD.slurm` calls variants on the bam files that have not had markduplicates run

`sort_arimappedAllSites_vcfs_NORMD.slurm` sorts the vcf files

`combine_arimapAllSites_vcf_NORMD.slurm` combines those into whole genome

`get_vcf_stats_arimapAllSites_NORMD.slurm` gets idepth for the vcf file only.


* I deleted all of the NORMD called variants - seems obvious that this is worse
* also deleted the bam files mapped to arizona without dups removed



### Arizona pixy analyses

* Arimap_pixy

1. `Pixy10_combined.slurm` runs pixy on the vcf containing the whole genome in a single file.

- `pixy_popfiles` directory contains populations files for pixy
- make sure vcf is tabiz indexed

2. `plot_pixy_combined_chroms.R` plots the combined genome pixy output.



### Arizona LFMM analyses

Convert vcf file to ped, have had trouble converting vcf to SNMF/LFMM

```
cd /project/getpop/lfmm

# then convert to plink
module load miniconda3
conda activate plink # my plink conda env


VCF_FILE=/project/getpop/vcf_allsites_Arizona_map/filt_getula_arimap_SNPs_MAC2_1k_THIN_full_genome.vcf.gz
plink --vcf $VCF_FILE --recode12 --out arimap_1k_genome --allow-extra-chr
```

Have been running LFMM stuff in main scripts for for it - **imputation of genoptypes for LFMM is running now**



<br>
<br>
<br>
<br>
<br>
<br>

Is it possible/desirable to filter the invariants separately from the variants? E.g., lower cov for the invariants?
no, doesn't seem like.

Run a quick check by calculating depth at the variant sites (biallelic):


```
cd /project/getpop/vcf_allsites_Arizona_map
salloc -A inbreh -t 0-05:00 --mem=200G --cpus-per-task=1
module load gcc/12.2.0 vcftools/0.1.16

OUT=/project/getpop/vcf_stats/vcf_stats_allsites_Arizona_map
OUTFILE=vcf_stats_allsites_Arizona_map_SNPS


vcftools --gzvcf getula_arimap_AllSites_full_genome.vcf.gz --min-alleles 2 --max-alleles 2 --depth --out $OUT/$OUTFILE
```

<br>



# Ok, plans on the remapped data:

Definitely going to want to drop out 

- FHSM7738 (~3x coverage)

Probably drop

- FHSM13623 (~4x)
- Mayybe FHSM11300 (~5x)

- add in the 6x cov filter and then look an imiss before dropping anybody





# For just WCH, slam out some stuff - don't super worry about removing individuals & stuff

Use Arizona mapping, the filtering in there doesn't remove inds


