The following pipeline describes the process I use to analyze publicly available RNAseq data, which is required by FungiExpresZ:

### Data Download: 
I utilize the prefetch and fastq-dump commands to download and decompress the RNAseq datasets. These commands allow me to obtain the raw reads needed for further analysis.

### Mapping to the Genome: 
To align the raw reads with the reference genome, I employ the hisat2 software. This step enables me to map the reads to their corresponding genomic locations, facilitating downstream analysis.

### Quantification of Gene Expression:
Once the reads have been mapped, I utilize stringtie to calculate the Fragments Per Kilobase of transcript per Million mapped reads (FPKM) value. This metric provides an estimation of gene expression levels by normalizing for transcript length and library size.

### Expression Matrix Construction: 
In order to generate the expression matrix required by FungiExpresZ, I have developed a custom script. This script takes the output from the previous steps and constructs the matrix, which serves as input for subsequent analyses.

### Softwares / Tools
The following softwares or tools are required:
(1) SRA Toolkit
(2) HISAT2
(3) StringTie
(4) Samtools
(5) Snakemake

Please make sure to install these tools before running snakemake.

## 1. Collect RNAseq data access number
updating ......

## 2. Download genome fasta file and gff file
For most fungi species, the genome information file could be downloaded from Fungidb

To create a HISAT2 index, you can use the following command:
```
hisat2-build species_genome.fasta species_hisat2_index
```

## 3.Download and process RNAseq data with snakemake
Edit the Snakefile to change the file path within it.

Run a dry run to check for any errors:
```
snakemake -s FungiExpresZ_single_end.py --dryrun
snakemake -s FungiExpresZ_paired_end.py --dryrun
```

If there are no errors, run the Snakefile using the following commands:
```
snakemake -s FungiExpresZ_single_end.py --latency-wait 10 --keep-going --rerun-incomplete
snakemake -s FungiExpresZ_paired_end.py --latency-wait 10 --keep-going --rerun-incomplete
```



## 4. Selecting file according to mapping rate and construct an expression matrix

updating ......
