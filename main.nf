nextflow.enable.dsl=2

params.reads = "fastq/*_{1,2}.fastq"

process FASTQC {

    publishDir "results/fastqc", mode: "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.html")
    path("*_fastqc.zip")

    script:
    """
    fastqc ${reads[0]} ${reads[1]}
    """
}

process FASTP {

    publishDir "results/fastp", mode: "copy"
    publishDir "trimmed", mode: "copy", pattern: "*.fastq"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("trimmed_${sample_id}_1.fastq"),
          path("trimmed_${sample_id}_2.fastq")

    script:
    """
    fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o trimmed_${sample_id}_1.fastq \
        -O trimmed_${sample_id}_2.fastq \
        -h ${sample_id}_fastp.html \
        -j ${sample_id}_fastp.json
    """
}

process HISAT2_ALIGN {

    publishDir "results/alignment", mode: "copy"

    input:
    tuple val(sample_id),
          path(read1),
          path(read2)

    output:
    tuple val(sample_id),
          path("${sample_id}.sam")

    script:
    """
    hisat2 \
        -x /home/vinay_bioinformatic/project/refrence/hisat2_index/ecoli \
        -1 ${read1} \
        -2 ${read2} \
        -S ${sample_id}.sam
    """
}
process SAMTOOLS {

    publishDir "results/bam", mode: "copy"

    input:
    tuple val(sample_id), path(sam)

    output:
tuple val(sample_id),
      path("${sample_id}.sorted.bam"),
      path("${sample_id}.sorted.bam.bai")

    script:
    """
    samtools view -bS ${sam} | samtools sort -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}
process FEATURECOUNTS {

    publishDir "results/counts", mode: "copy"

    input:
    tuple val(sample_id),
          path(bam),
          path(bai)

    output:
    path("${sample_id}_gene_counts.txt")

    script:
    """
    featureCounts \
-p \
-B \
-C \
-F GFF \
-t CDS \
-g locus_tag \
-a /home/vinay_bioinformatic/project/refrence/GCF_000005845.2_ASM584v2_genomic.gff \
-o ${sample_id}_gene_counts.txt \
${bam}
    """
}
workflow {

    reads = Channel.fromFilePairs(params.reads, checkIfExists: true)

    FASTQC(reads)

    trimmed_reads = FASTP(reads)

    sam_files = HISAT2_ALIGN(trimmed_reads)

    bam_files = SAMTOOLS(sam_files)

    FEATURECOUNTS(bam_files)
}
