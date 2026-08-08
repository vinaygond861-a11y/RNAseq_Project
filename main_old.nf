nextflow.enable.dsl=2

process FASTQC {

    publishDir "results/fastqc", mode: 'copy'

    input:
    path reads

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc $reads
    """
}


process FASTP {

    publishDir "trimmed", mode: 'copy'

    input:
    path reads

    output:
    path "*.fastq"

    script:
    """
    fastp -i ${reads} -o trimmed_${reads.baseName}.fastq
    """
}

process HISAT2_ALIGN {

    publishDir "results/alignment", mode: 'copy'

    input:
    path reads

    output:
    path "*.sam"

    script:
    """
    hisat2 \
      -x /home/vinay_bioinformatic/project/refrence/hisat2_index/ecoli \
      -U ${reads} \
      -S ${reads.baseName}.sam
    """
}

process SAMTOOLS {

    publishDir "results/bam", mode: 'copy'

    input:
    path sam

    output:
    path "*.sorted.bam"
    

    script:
    """
    samtools view -bS ${sam} | samtools sort -o ${sam.baseName}.sorted.bam
    samtools index ${sam.baseName}.sorted.bam
    """
}
process FEATURECOUNTS {

    publishDir "results/counts", mode: 'copy'

    input:
    path bam

    output:
    path "gene_counts.txt"

    script:
    """
    featureCounts \
    -F GFF \
    -t CDS \
    -g locus_tag \
    -a /home/vinay_bioinformatic/project/refrence/GCF_000005845.2_ASM584v2_genomic.gff \
    -o gene_counts.txt \
    ${bam}
    """
}

workflow {

    reads = Channel.fromPath("fastq/*.fastq")

    FASTQC(reads)

    trimmed = FASTP(reads)

    sam_files = HISAT2_ALIGN(trimmed)

bam_files = SAMTOOLS(sam_files)

FEATURECOUNTS(bam_files)

}
