The following pipeline describes the process I use to analyze publicly available RNAseq data, which is required by FungiExpresZ:

    Data Download: I utilize the prefetch and fastq-dump commands to download and decompress the RNAseq datasets. These commands allow me to obtain the raw reads needed for further analysis.

    Mapping to the Genome: To align the raw reads with the reference genome, I employ the hisat2 tool. This step enables me to map the reads to their corresponding genomic locations, facilitating downstream analysis.

    Quantification of Gene Expression: Once the reads have been mapped, I utilize stringtie to calculate the Fragments Per Kilobase of transcript per Million mapped reads (FPKM) value. This metric provides an estimation of gene expression levels by normalizing for transcript length and library size.

    Expression Matrix Construction: In order to generate the expression matrix required by FungiExpresZ, I have developed a custom script. This script takes the output from the previous steps and constructs the matrix, which serves as input for subsequent analyses.

## 1. collect RNAseq dataset
updating ......

## 2.download, mapping and calculate FPKM

Here is the snakemake file FungiExpresZ_double_end.py
The file is uploaded to script folder
what you need to do is change the four config parameters
```
# config
hisat2_index = '/home/genome_info/Calbicans/Calbicans_hisat2_index'
gff_file = '/home/genome_info/Calbicans/FungiDB-66_CalbicansSC5314.gff'
output_path = '/home/FungiExpresZ/Candida_albicans/test/output'
sra_file_path = '/home/FungiExpresZ/Candida_albicans/test/test.txt'


# read_sra_number
def read_sra_numbers():
    with open(sra_file_path) as f:
        return [line.strip() for line in f]

sra_numbers = read_sra_numbers()


# main
rule all:
    input:
        expand(os.path.join(output_path, "{sra_number}.xls"), sra_number=sra_numbers)

rule download:
    output:
        xls = os.path.join(output_path, "{sra_number}.xls"),
        stat = os.path.join(output_path, "{sra_number}_alignmentstat")
    params:
        bam = os.path.join(output_path, "{sra_number}.bam"),
        Index = hisat2_index,
        gff = gff_file
    shell:
        """
        prefetch {wildcards.sra_number} -O ./
        mv ./{wildcards.sra_number}/{wildcards.sra_number}.sra ./{wildcards.sra_number}.sra
        rm -r {wildcards.sra_number}
        fasterq-dump ./{wildcards.sra_number}.sra --outdir ./
        hisat2 -p 8 -x {params.Index} -1 {wildcards.sra_number}_1.fastq -2 {wildcards.sra_number}_2.fastq | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o {params.bam}
        stringtie {params.bam} -e -G {params.gff} > {output.xls}
        samtools flagstat {params.bam} > {output.stat}
        rm -f {wildcards.sra_number}.sra
        rm -f {wildcards.sra_number}_1.fastq
        rm -f {wildcards.sra_number}_2.fastq
        """

```

## 3. Selecting file according to mapping rate and construct an expression matrix

updating ......
