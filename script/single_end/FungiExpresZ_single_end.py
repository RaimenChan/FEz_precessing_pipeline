# config
hisat2_index = '/home/genome_info/Calbicans/Calbicans_hisat2_index'
gff_file = '/home/genome_info/Calbicans/FungiDB-66_CalbicansSC5314.gff'
output_path = '/home/FungiExpresZ/Candida_albicans/single_end/output'
sra_file_path = '/home/FungiExpresZ/Candida_albicans/single_end/Calbicans_RNAseq_run_single_end.txt'


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
        hisat2 -p 8 -x {params.Index} -U {wildcards.sra_number}.fastq | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o {params.bam}
        stringtie {params.bam} -e -G {params.gff} > {output.xls}
        samtools flagstat {params.bam} > {output.stat}
        rm -f {wildcards.sra_number}.sra
        rm -f {wildcards.sra_number}.fastq
        """