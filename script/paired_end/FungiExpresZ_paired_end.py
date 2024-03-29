# config
hisat2_index = '~/genome_info/Anidulans/Anidulans_hisat2_index'
gff_file = '~/genome_info/Anidulans/FungiDB-65_AnidulansFGSCA4.gff'
output_path = '~/FungiExpresZ/Aspergillus_nidulans/paired_end/2/output'
sra_file_path = '~/FungiExpresZ/Aspergillus_nidulans/paired_end/2/Aspergillus_nidulans_RNAseq-paired_end_2.txt'


# read_sra_number
def read_sra_numbers():
    with open(sra_file_path) as f:
        return [line.strip() for line in f]

sra_numbers = read_sra_numbers()


# main
rule all:
    input:
        expand(os.path.join(output_path,"expression", "{sra_number}.xls"), sra_number=sra_numbers)

rule download:
    output:
        xls = os.path.join(output_path, "expression", "{sra_number}.xls"),
        stat = os.path.join(output_path, "alignment", "{sra_number}_alignmentstat")
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
        rm -f {params.bam}
        """