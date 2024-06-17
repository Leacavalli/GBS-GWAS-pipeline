process Trimmomatic {
  conda "java"

  publishDir "${params.path_nextflow_dir}/1.TRIMMED_READS" , overwrite: true

  input:
  val base

  output:
  val base
  path "${base}_trimmed_1.fq.gz"
  path "${base}_trimmed_2.fq.gz"

  """
  java -jar /n/holylfs05/LABS/hanage_lab/Lab/holyscratch01/lcavalli/Software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog ${base}_1.log ${params.path_nextflow_dir}/0.RAW_READS/${base}_1.fastq.gz ${params.path_nextflow_dir}/0.RAW_READS/${base}_2.fastq.gz ${base}_trimmed_1.fq.gz ${base}_unpaired_1.fq.gz ${base}_trimmed_2.fq.gz ${base}_unpaired_2.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > ${base}_trim.out
  rm ${base}_1.log ${base}_unpaired_1.fq.gz ${base}_unpaired_2.fq.gz ${base}_trim.out
  """
}


process Kraken {
  publishDir "${params.path_nextflow_dir}/2.KRAKEN_REPORTS" , overwrite: true

  input:
  val base
  val trimmed_read1
  val trimmed_read2

  output:
  val base
  val trimmed_read1
  val trimmed_read2
  path "${base}_kraken.out"
  path "${base}_report.txt.gz"

  script:
  """
  kraken2 --db /n/holylfs05/LABS/hanage_lab/Lab/holyscratch01/lcavalli/Software/kraken2/standard_db --report ${base}_report.txt --paired ${trimmed_read1} ${trimmed_read2} > ${base}_kraken.out
  gzip ${base}_report.txt
  """
}

process Kraken_tools {
  conda "biopython"

  publishDir "${params.path_nextflow_dir}/3.FILTERED_READS", overwrite: true

  input:
  val base
  val trimmed_read1
  val trimmed_read2
  val kraken_output
  val kraken_report

  output:
  val base
  path "${base}_filtered_1.fq.gz"
  path "${base}_filtered_2.fq.gz"

  script:
  """
  extract_kraken_reads.py -k ${kraken_output} -s1 ${trimmed_read1} -s2 ${trimmed_read2} -o ${base}_filtered_1.fq -o2 ${base}_filtered_2.fq --fastq-output --taxid 1301
  gzip ${base}_filtered_1.fq
  gzip ${base}_filtered_2.fq
  """
}

process FastQC {
  conda "fastqc"

  input:
  val base
  val filtered_read1
  val filtered_read2

  output:
  val base
  val filtered_read1
  val filtered_read2

  script:
  """
  fastqc ${filtered_read1} -o ${params.path_nextflow_dir}/4.FASTQC
  fastqc ${filtered_read2} -o ${params.path_nextflow_dir}/4.FASTQC
  """
}

process Serotype_MLST {
  conda "srst2"

  publishDir "${params.path_nextflow_dir}/5.SEROTYPE_MLST", overwrite: true

  input:
  val base
  val filtered_read1
  val filtered_read2

  output:
  val base
  val filtered_read1
  val filtered_read2
  path "${base}__genes__GBS-SBG__results.txt"
  path "${base}__mlst__Streptococcus_agalactiae__results.txt"

  script:
  """
  srst2 --input_pe ${params.path_nextflow_dir}/0.RAW_READS/${base}_{1,2}.fastq.gz --output ${base} --log --gene_db ${params.path_nextflow_dir}/Files/GBS-SBG.fasta --mlst_db ${params.path_nextflow_dir}/Files/Streptococcus_agalactiae.fasta --mlst_definitions ${params.path_nextflow_dir}/Files/profiles_csv --mlst_delimiter '_'
  """
}

process SNIPPY {
  conda "snippy"

  input:
  val base
  val filtered_read1
  val filtered_read2
  val serotypes_out
  val mlst_out

  output:
  val base
  val filtered_read1
  val filtered_read2

  script:
  """
  snippy --outdir ${params.path_nextflow_dir}/6.SNIPPY/${base} --R1 ${filtered_read1} --R2 ${filtered_read2} --ref ${params.path_nextflow_dir}/Files/AP018935.1.fa --cpus 16
  """
}

process ASSEMBLY {
  conda "unicycler"

  publishDir "${params.path_nextflow_dir}/7.ASSEMBLIES", overwrite: true

  input:
  val base
  val filtered_read1
  val filtered_read2


  output:
  val base
  path "${base}.fasta"

  script:
  """
  unicycler -1 ${filtered_read1}  -2 ${filtered_read2} -o ${base}.unicycler.out
  mv ${base}.unicycler.out/assembly.fasta ${base}.fasta
  rm -r ${base}.unicycler.out
  """
}


process ANNOTATION {
  conda "prokka"

  publishDir "${params.path_nextflow_dir}/8.ANNOTATION", overwrite: true

  input:
  val base
  val fasta_file

  output:
  val base
  val fasta_file
  path "${base}.gff"

  script:
  """
  prokka ${fasta_file} --outdir ${base}_prokka  --prefix ${base}
  mv ${base}_prokka/*.gff .
  rm -r ${base}_prokka
  """
}




workflow {
  Channel.of(params.sample)
    | Trimmomatic
    | Kraken
    | Kraken_tools
    | FastQC
    | Serotype_MLST
    | SNIPPY
    | ASSEMBLY
    | ANNOTATION
}
