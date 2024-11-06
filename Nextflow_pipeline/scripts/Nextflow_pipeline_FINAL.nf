params.OUTGROUP = false
params.main = false
params.sub1 = false
params.Fasttree = false
params.Mash = false
params.SC = false
params.sub2 = false

process Make_architecture {
  input:
  val base

  output:
  val base

  """
  mkdir -p ${params.path_nextflow_dir}/1.TRIMMED_READS
  mkdir -p ${params.path_nextflow_dir}/2.FASTQC
  mkdir -p ${params.path_nextflow_dir}/3.SEROTYPE_MLST
  mkdir -p ${params.path_nextflow_dir}/4.SNIPPY
  mkdir -p ${params.path_nextflow_dir}/5.ASSEMBLIES
  mkdir -p ${params.path_nextflow_dir}/6.ANNOTATION
  mkdir -p ${params.path_nextflow_dir}/7.FastANI
  mkdir -p ${params.path_nextflow_dir}/8.QUAST
  """
}

process Serotype_MLST {
  conda "srst2"

  publishDir "${params.path_nextflow_dir}/3.SEROTYPE_MLST", overwrite: true

  input:
  val base

  output:
  path "${base}__genes__GBS-SBG__results.txt"
  path "${base}__mlst__Streptococcus_agalactiae__results.txt"

  script:
  """
  srst2 --input_pe ${params.path_nextflow_dir}/0.RAW_READS/${base}_{1,2}.fastq.gz --output ${base} --log --gene_db ${params.path_nextflow_dir}/Files/GBS-SBG.fasta --mlst_db ${params.path_nextflow_dir}/Files/Streptococcus_agalactiae.fasta --mlst_definitions ${params.path_nextflow_dir}/Files/profiles_csv --mlst_delimiter '_'
  """
}

process Trimmomatic {
  conda "openjdk=11"

  publishDir "${params.path_nextflow_dir}/1.TRIMMED_READS" , overwrite: true

  input:
  val base

  output:
  val base
  path "${base}_trimmed_1.fq.gz"
  path "${base}_trimmed_2.fq.gz"

  """
  java -jar ${params.path_nextflow_dir}/Files/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog ${base}_1.log ${params.path_nextflow_dir}/0.RAW_READS/${base}_1.fastq.gz ${params.path_nextflow_dir}/0.RAW_READS/${base}_2.fastq.gz ${base}_trimmed_1.fq.gz ${base}_unpaired_1.fq.gz ${base}_trimmed_2.fq.gz ${base}_unpaired_2.fq.gz ILLUMINACLIP:${params.path_nextflow_dir}/Files/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:50 > ${base}_trim.out
  """
}

process FastQC {
  conda "fastqc"

  input:
  val base
  val trimmed_read1
  val trimmed_read2

  output:
  val base

  script:
  """
  fastqc ${trimmed_read1} --extract -o ${params.path_nextflow_dir}/2.FASTQC
  fastqc ${trimmed_read2} --extract -o ${params.path_nextflow_dir}/2.FASTQC
  """
}

process SNIPPY {
  conda "${params.path_nextflow_dir}/Files/Snippy_environment.yml"

  input:
  val base
  val trimmed_read1
  val trimmed_read2

  output:
  val base

  script:
  """
  snippy --outdir ${params.path_nextflow_dir}/4.SNIPPY/${base} --R1 ${trimmed_read1} --R2 ${trimmed_read2} --ref ${params.path_nextflow_dir}/Files/AP018935.1.fa
  """
}

process ASSEMBLY {
  conda "${params.path_conda_envs}/unicycler"

  publishDir "${params.path_nextflow_dir}/5.ASSEMBLIES", overwrite: true

  input:
  val base
  val trimmed_read1
  val trimmed_read2

  output:
  val base
  path "${base}.fasta"

  script:
  """
  unicycler -1 ${trimmed_read1}  -2 ${trimmed_read2} -o ${base}.unicycler.out
  mv ${base}.unicycler.out/assembly.fasta ${base}.fasta
  """
}


process ANNOTATION {
  conda "prokka"

  publishDir "${params.path_nextflow_dir}/6.ANNOTATION", overwrite: true

  input:
  val base
  val fasta_file

  output:
  path "${base}.gff"

  script:
  """
  prokka ${fasta_file} --outdir ${base}_prokka  --prefix ${base}
  mv ${base}_prokka/*.gff .
  """
}

process FastANI {
  conda "fastani"

  publishDir "${params.path_nextflow_dir}/7.FastANI", overwrite: true

  input:
  val base
  val fasta_file

  output:
  path "${base}.fastani.out"

  script:
  """
  fastANI -q ${fasta_file} -r ${params.path_nextflow_dir}/Files/AP018935.1.fa -o ${base}.fastani.out
  """
}

process QUAST {
  conda "quast anaconda::joblib"

  publishDir "${params.path_nextflow_dir}/8.QUAST", overwrite: true

  input:
  val base
  val fasta_file

  output:
  val base
  path "${base}.quast.out"

  script:
  """
  quast -o ${base}.quast.out ${fasta_file}
  """
}


process Filter_out_QC {
  input:
  val base
  val sero_MLST_data
  val FASTQC_data
  val SNIPPY_data
  val ANNOTATION_data
  val FastANI_data

  output:
  val base

  script:
  """
  # Create necessary directories if they don't already exist
    mkdir -p ${params.path_nextflow_dir}/4.SNIPPY/Filter_out_QC
    mkdir -p ${params.path_nextflow_dir}/5.ASSEMBLIES/Filter_out_QC
    mkdir -p ${params.path_nextflow_dir}/6.ANNOTATION/Filter_out_QC

  # Extract information from the QUAST report
    n50=\$(grep -w "N50" ${params.path_nextflow_dir}/8.QUAST/${base}.quast.out/report.txt | awk '{print \$2}')
    gc=\$(grep -w "GC (%)" ${params.path_nextflow_dir}/8.QUAST/${base}.quast.out/report.txt | awk '{print \$3}')
    total_length=\$(grep -w "Total length" ${params.path_nextflow_dir}/8.QUAST/${base}.quast.out/report.txt | grep -v "(>= " | awk '{print \$3}')
    contigs=\$(grep "^# contigs " ${params.path_nextflow_dir}/8.QUAST/${base}.quast.out/report.txt | tail -1 | awk '{print \$3}')
    ANI=\$(cat ${params.path_nextflow_dir}/7.FastANI/${base}.fastani.out | awk '{print \$3}')
    Adapter_content_1=\$(grep "Adapter Content" ${params.path_nextflow_dir}/2.FASTQC/${base}_trimmed_1_fastqc/summary.txt | awk '{print \$1}')
    Adapter_content_2=\$(grep "Adapter Content" ${params.path_nextflow_dir}/2.FASTQC/${base}_trimmed_2_fastqc/summary.txt | awk '{print \$1}')

  # Check the conditions and move files if necessary
    if (( \$(echo "\$ANI < 95" | bc -l) )) || [ "\$n50" -lt 30000 ] || [ "\$contigs" -gt 250 ] || (( \$(echo "\$gc <= 30" | bc -l) )) || (( \$(echo "\$gc >= 40" | bc -l) )) || [ "\$total_length" -lt 1700000 ] || [ "\$total_length" -gt 2400000 ] || [[ "\$Adapter_content_1" == "FAIL" || "\$Adapter_content_2" == "FAIL" ]]; then
      # Check if the SNIPPY directory exists, and move it if it does
      if [ -d "${params.path_nextflow_dir}/4.SNIPPY/${base}" ]; then
        mv ${params.path_nextflow_dir}/4.SNIPPY/${base} ${params.path_nextflow_dir}/4.SNIPPY/Filter_out_QC
      fi
      # Check if the assembly fasta file exists, and move it if it does
      if [ -f "${params.path_nextflow_dir}/5.ASSEMBLIES/${base}.fasta" ]; then
      mv ${params.path_nextflow_dir}/5.ASSEMBLIES/${base}.fasta ${params.path_nextflow_dir}/5.ASSEMBLIES/Filter_out_QC
      fi
      # Check if the annotation gff file exists, and move it if it does
      if [ -f "${params.path_nextflow_dir}/6.ANNOTATION/${base}.gff" ]; then
      mv ${params.path_nextflow_dir}/6.ANNOTATION/${base}.gff ${params.path_nextflow_dir}/6.ANNOTATION/Filter_out_QC
      fi
      # Remove Sample from GWAS phenotype file
      grep -v "^${base}\s" ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt >> ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt
    fi
  """
}


process POPPUNK {
  cpus 8

  conda "poppunk==2.6.0"

  publishDir "${params.path_nextflow_dir}/9.POPPUNK", overwrite: true

  input:
  val base

  output:
  val base

  script:
  """
  # Prepare input file
  mkdir -p ${params.path_nextflow_dir}/9.POPPUNK
  for i in ${params.path_nextflow_dir}/5.ASSEMBLIES/*fasta
  do
  sample=\$(basename \$i | cut -d '.' -f1-1)
  echo -e "\$sample\t\$i"
  done > ${params.path_nextflow_dir}/9.POPPUNK/qfile.txt

  # Run
  # 1. Create sketch database
  poppunk --create-db --output ${params.path_nextflow_dir}/9.POPPUNK/GBS_GWAS_db --r-files ${params.path_nextflow_dir}/9.POPPUNK/qfile.txt --threads 8
  # 2. Run QC on the database
  poppunk --qc-db --ref-db ${params.path_nextflow_dir}/9.POPPUNK/GBS_GWAS_db --length-range 1500000 2500000  --max-zero-dist 1 --threads 8
  # 3. Fit a model, check cluster (core+accessory)
  poppunk --fit-model dbscan --ref-db ${params.path_nextflow_dir}/9.POPPUNK/GBS_GWAS_db --threads 8
  # 4. Fit a model, check cluster (core only)
  poppunk --fit-model refine --ref-db ${params.path_nextflow_dir}/9.POPPUNK/GBS_GWAS_db --threads 8 --indiv-refine core
  # Keep only the samples with classified SC for pyseer
  tail -n +2 ${params.path_nextflow_dir}/9.POPPUNK/GBS_GWAS_db/GBS_GWAS_db_core_clusters.csv | cut -d ',' -f1 | sort > sorted_file1.csv
  cut -f1 ${params.path_nextflow_dir}/Files/phenotypes.txt | sort > sorted_file2.csv
  comm -3 sorted_file1.csv  sorted_file2.csv | sed 's/\t//g' > removed_samples.txt
  grep -v -F -f removed_samples.txt ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt > ${params.path_nextflow_dir}/Files/phenotypes_filtered_SC.txt
  """
}

process ASSIGN_CC {
input:
val base

output:
val base

script:
  """
  # Create Output Directory
  mkdir -p ${params.path_nextflow_dir}/10.ASSIGN_CC

  # Append all results from mlst to a single file
  awk '(NR == 1) || (FNR > 1)' ${params.path_nextflow_dir}/3.SEROTYPE_MLST/*__mlst__Streptococcus_agalactiae__results.txt > ${params.path_nextflow_dir}/3.SEROTYPE_MLST/combined_results.txt
  cut -f1,2 ${params.path_nextflow_dir}/3.SEROTYPE_MLST/combined_results.txt | tail -n +2 | sed 's/\t/,/g' | sort -t',' -k2,2 | sed  's/\\*//g'  > ${params.path_nextflow_dir}/10.ASSIGN_CC/STs.txt

  # Define the file paths
  cut -d',' -f1,9 ${params.path_nextflow_dir}/Files/GBS_CC_profiles.csv | tail -n +2  | sort -t',' -k1,1 > ${params.path_nextflow_dir}/10.ASSIGN_CC/GBS_ST_CC_ref.txt

  # Process the files and join them based on the 'ST' column
  echo "ST,Sample,clonal_complex" > ${params.path_nextflow_dir}/10.ASSIGN_CC/combined_results_with_cc.txt
  join -t',' -1 2 -2 1 -a 1 -e 'NA' -o 1.2,1.1,2.2 ${params.path_nextflow_dir}/10.ASSIGN_CC/STs.txt ${params.path_nextflow_dir}/10.ASSIGN_CC/GBS_ST_CC_ref.txt -a1  >> ${params.path_nextflow_dir}/10.ASSIGN_CC/combined_results_with_cc.txt
  """
}

process SNIPPY_MULTI {
conda "${params.path_nextflow_dir}/Files/Snippy_environment.yml"

  publishDir "${params.path_nextflow_dir}/11.SNIPPY_MULTI", overwrite: true

  input:
  val base

  output:
  path "core.aln"
  path "core.full.aln"
  path "core.tab"
  path "core.vcf"
  path "core.txt"
  path "core.ref.fa"
  path "clean.full.aln"

  script:
  """
  snippy-core --ref ${params.path_nextflow_dir}/Files/AP018935.1.fa ${params.path_nextflow_dir}/4.SNIPPY/*RR*
  snippy-clean_full_aln core.full.aln > clean.full.aln
  """
}


process SNIPPY_MULTI_OUTGROUP {
conda "${params.path_nextflow_dir}/Files/Snippy_environment.yml"

  publishDir "${params.path_nextflow_dir}/11.SNIPPY_MULTI", overwrite: true

  input:
  val base

  output:
  path "core.aln"
  path "core.full.aln"
  path "core.tab"
  path "core.vcf"
  path "core.txt"
  path "core.ref.fa"
  path "clean.full.aln"

  script:
  """
  snippy-core --ref ${params.path_nextflow_dir}/Files/AP018935.1.fa --prefix core_outgroup ${params.path_nextflow_dir}/4.SNIPPY/* ${params.path_nextflow_dir}/4.SNIPPY/Outgroup
  snippy-clean_full_aln core_outgroup.full.aln > clean_outgroup.full.aln
  """
}

process RAxML_MLsearch {
  cpus 50

  input:
  val snippy_out

  output:
  val snippy_out

  script:
  """
  mkdir -p ${params.path_nextflow_dir}/12.Phylogeny
  mkdir -p ${params.path_nextflow_dir}/12.Phylogeny/12.1.RAxML
  # Generate 100 ML trees on distinct starting trees and output the best likelihood tree
  ${params.path_nextflow_dir}/Files/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 50 -m GTRGAMMA -p 12345 -# 100 -s ${params.path_nextflow_dir}/11.SNIPPY_MULTI/clean.full.aln -w ${params.path_nextflow_dir}/12.Phylogeny/12.1.RAxML -n T1
  """
}

process RAxML_bootstrap {
  cpus 50

  input:
  val snippy_out

  output:
  val snippy_out

  script:
  """
  mkdir -p ${params.path_nextflow_dir}/12.Phylogeny
  mkdir -p ${params.path_nextflow_dir}/12.Phylogeny/12.1.RAxML
  # Generate 250 bootstrap tree to infer statistical support of the branches.
  ${params.path_nextflow_dir}/Files/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 50 -m GTRGAMMA -p 12345 -b 12345 -# 250 -s ${params.path_nextflow_dir}/11.SNIPPY_MULTI/clean.full.aln -w ${params.path_nextflow_dir}/12.Phylogeny/12.1.RAxML -n T2
  """
}

process RAxML_bipartition {
  cpus 50

  input:
  val RAxML_MLsearch_out
  val RAxML_bootstrap_out

  output:
  val RAxML_MLsearch_out

  script:
  """
 # Draw bipartitions on the best ML tree
  ${params.path_nextflow_dir}/Files/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 50 -m GTRCAT -p 12345 -f b -t ${params.path_nextflow_dir}/12.Phylogeny/12.1.RAxML/RAxML_bestTree.T1 -z ${params.path_nextflow_dir}/12.Phylogeny/12.1.RAxML/RAxML_bootstrap.T2 -w ${params.path_nextflow_dir}/12.Phylogeny/12.1.RAxML -n T3
    """
}




process RAxML_OUTGROUP {
  cpus 50

  input:
  val snippy_out

  output:
  val snippy_out

  script:
  """
  mkdir -p ${params.path_nextflow_dir}/12.Phylogeny
  mkdir -p ${params.path_nextflow_dir}/12.Phylogeny/12.3.RAxML_outgroup
  # Generate 100 ML trees on distinct starting trees and output the best likelihood tree
  ${params.path_nextflow_dir}/Files/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 50 -m GTRGAMMA -p 12345 -# 100 -s ${params.path_nextflow_dir}/11.SNIPPY_MULTI/clean_outgroup.full.aln -w ${params.path_nextflow_dir}/12.Phylogeny/12.3.RAxML_outgroup -n T1 -o Outgroup
  # Generate 250 bootstrap tree to infer statistical support of the branches.
  ${params.path_nextflow_dir}/Files/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 50 -m GTRGAMMA -p 12345 -b 12345 -# 250 -s ${params.path_nextflow_dir}/11.SNIPPY_MULTI/clean_outgroup.full.aln -w ${params.path_nextflow_dir}/12.Phylogeny/12.3.RAxML_outgroup -n T2 -o Outgroup
  # Draw bipartitions on the best ML tree
  ${params.path_nextflow_dir}/Files/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 50 -m GTRCAT -p 12345 -f b -t ${params.path_nextflow_dir}/12.Phylogeny/12.3.RAxML_outgroup/RAxML_bestTree.T1 -z ${params.path_nextflow_dir}/12.Phylogeny/12.3.RAxML_outgroup/RAxML_bootstrap.T2 -w ${params.path_nextflow_dir}/12.Phylogeny/12.3.RAxML_outgroup -n T3
  """
}

process FASTTREE {
  input:
  val snippy_out

  output:
  val snippy_out

  script:
  """
  mkdir -p ${params.path_nextflow_dir}/12.Phylogeny
  mkdir -p ${params.path_nextflow_dir}/12.Phylogeny/12.2.FastTree
  FastTree -gtr -nt ${params.path_nextflow_dir}/11.SNIPPY_MULTI/clean.full.aln > ${params.path_nextflow_dir}/12.Phylogeny/12.2.FastTree/GBS_GWAS_FastTree_phylo.tre
  """
}

process PANAROO {
  cpus 8

  conda "${params.path_nextflow_dir}/Files/Panaroo_environment.yml"

  input:
  val base

  output:
  val base

  script:
  """
  mkdir -p ${params.path_nextflow_dir}/13.Pangenome_analysis
  mkdir -p ${params.path_nextflow_dir}/13.Pangenome_analysis/13.1.Panaroo
  panaroo -t 8 -i ${params.path_nextflow_dir}/6.ANNOTATION/*.gff -o ${params.path_nextflow_dir}/13.Pangenome_analysis/13.1.Panaroo --clean-mode strict --alignment pan --aligner mafft
  """
}

process PANAROO_CLARC {
  conda "${params.path_conda_envs}/clarc_env"

  input:
  val base

  output:
  val base

  script:
  """
  mkdir -p ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC
  mkdir -p ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC/data
  mkdir -p ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC/Output
  # Prepare input Files
  cp ${params.path_nextflow_dir}/13.Pangenome_analysis/13.1.Panaroo/gene_presence_absence_roary.csv ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC/data
  cp ${params.path_nextflow_dir}/13.Pangenome_analysis/13.1.Panaroo/pan_genome_reference.fa ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC/data
  cp ${params.path_nextflow_dir}/13.Pangenome_analysis/13.1.Panaroo/gene_data.csv ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC/data
  for i in ${params.path_nextflow_dir}/5.ASSEMBLIES/*fasta
  do
  sample=\$(basename \$i | cut -d '.' -f1-1)
  echo "\$sample"
  done > ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC/data/needed_sample_names.txt

  # Run CLARC
  clarc --panaroo --input_dir ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC/data --output_dir ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC/Output
  """
}

process ROARY {
  conda "roary"

  input:
  val base

  output:
  val base

  script:
  """
  mkdir -p ${params.path_nextflow_dir}/13.Pangenome_analysis
  roary -e --mafft -p 8 ${params.path_nextflow_dir}/6.ANNOTATION/*.gff -f ${params.path_nextflow_dir}/13.Pangenome_analysis/13.2.ROARY
  """
}

process ROARY_CLARC {
  conda "${params.path_conda_envs}/clarc_env"

  input:
  val base

  output:
  val base

  script:
  """
  mkdir -p ${params.path_nextflow_dir}/13.Pangenome_analysis/13.4.ROARY_CLARC
  mkdir -p ${params.path_nextflow_dir}/13.Pangenome_analysis/13.4.ROARY_CLARC/data
  mkdir -p ${params.path_nextflow_dir}/13.Pangenome_analysis/13.4.ROARY_CLARC/Output
  # Prepare input Files
  cp ${params.path_nextflow_dir}/13.Pangenome_analysis/13.2.ROARY/gene_presence_absence.csv ${params.path_nextflow_dir}/13.Pangenome_analysis/13.4.ROARY_CLARC/data
  cp ${params.path_nextflow_dir}/13.Pangenome_analysis/13.2.ROARY/pan_genome_reference.fa ${params.path_nextflow_dir}/13.Pangenome_analysis/13.4.ROARY_CLARC/data
  for i in ${params.path_nextflow_dir}/5.ASSEMBLIES/*fasta
  do
  sample=\$(basename \$i | cut -d '.' -f1-1)
  echo "\$sample"
  done > ${params.path_nextflow_dir}/13.Pangenome_analysis/13.4.ROARY_CLARC/data/needed_sample_names.txt

  # Run CLARC
  clarc --input_dir ${params.path_nextflow_dir}/13.Pangenome_analysis/13.4.ROARY_CLARC/data --output_dir ${params.path_nextflow_dir}/13.Pangenome_analysis/13.4.ROARY_CLARC/Output
  """
}

process UNITIG {
  cpus 8

  conda "unitig-counter"

  input:
  val base

  output:
  val base

  script:
  """
  # Prepare input file
  for i in ${params.path_nextflow_dir}/5.ASSEMBLIES/*fasta
  do
  sample=\$(basename \$i | cut -d '.' -f1-1)
  echo -e "\$sample\t\$i"
  done > qfile.txt

  # Run Unitig
  unitig-counter -strains qfile.txt -output ${params.path_nextflow_dir}/14.UNITIGS -nb-cores 8
  cdbg-ops extend --graph ${params.path_nextflow_dir}/14.UNITIGS/graph --unitigs ${params.path_nextflow_dir}/14.UNITIGS/unitigs.txt > ${params.path_nextflow_dir}/14.UNITIGS/extended.txt
  """
}

process RAxML_dist {
  conda "pyseer"

  input:
  val input1
  val input2
  val input3

  output:
  val input1

  script:
  """
  # create output folders
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis

  # create phylogeny tsv file
  python ${params.path_nextflow_dir}/Files/pyseer/scripts/phylogeny_distance.py --lmm ${params.path_nextflow_dir}/12.Phylogeny/12.1.RAxML/RAxML_bestTree.T1 > ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis/RAxML_phylogeny_K.tsv
  """
}

process SNPGWAS_RAxML {
  cpus 8

  conda "pyseer"

  input:
  val base

  script:
  """
  pyseer --lmm --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --vcf ${params.path_nextflow_dir}/11.SNIPPY_MULTI/core.vcf --similarity ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis/RAxML_phylogeny_K.tsv --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis/SNPGWAS_RAxML.txt
  """
}

process DBGWAS_RAxML {
  cpus 8

  conda "pyseer"

  input:
  val base

  script:
  """
  pyseer --lmm --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --kmers ${params.path_nextflow_dir}/14.UNITIGS/unitigs.txt --uncompressed --similarity ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis/RAxML_phylogeny_K.tsv --output-patterns ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis/DBGWAS_RAxML_unitig_patterns.txt --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis/DBGWAS_RAxML.txt
  """
}

process PanGWAS_Panaroo_RAxML {
  cpus 8

  conda "pyseer"

  input:
  val base

  script:
  """
  pyseer --lmm --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --pres ${params.path_nextflow_dir}/13.Pangenome_analysis/13.1.Panaroo/gene_presence_absence.Rtab --similarity ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis/RAxML_phylogeny_K.tsv --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis/PanGWAS_Panaroo_RAxML.txt
  """
}


process Fasttree_dist {
  conda "pyseer"

  input:
  val input1
  val input2
  val input3

  output:
  val input1

  script:
  """
  # create output folders
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.1.Fasttree

  # create phylogeny tsv file
  python ${params.path_nextflow_dir}/Files/pyseer/scripts/phylogeny_distance.py --lmm ${params.path_nextflow_dir}/12.Phylogeny/12.2.FastTree/GBS_GWAS_FastTree_phylo.tre > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.1.Fasttree/FastTree_phylogeny_K.tsv
  """
}

process SNPGWAS_Fasttree {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  pyseer --lmm --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --vcf ${params.path_nextflow_dir}/11.SNIPPY_MULTI/core.vcf --similarity ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.1.Fasttree/FastTree_phylogeny_K.tsv --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.1.Fasttree/SNPGWAS_Fasttree.txt
  """
}

process DBGWAS_Fasttree {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  pyseer --lmm --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --kmers ${params.path_nextflow_dir}/14.UNITIGS/unitigs.txt --uncompressed --similarity ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.1.Fasttree/FastTree_phylogeny_K.tsv --output-patterns ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.1.Fasttree/DBGWAS_Fasttree_unitig_patterns.txt --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.1.Fasttree/DBGWAS_Fasttree.txt
  """
}

process PanGWAS_Panaroo_Fasttree {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  pyseer --lmm --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --pres ${params.path_nextflow_dir}/13.Pangenome_analysis/13.1.Panaroo/gene_presence_absence.Rtab --similarity ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.1.Fasttree/FastTree_phylogeny_K.tsv --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.1.Fasttree/PanGWAS_Panaroo_phylo_Fasttree.txt
  """
}

process Mash_dist {
  conda "pyseer"

  input:
  val input1
  val input2
  val input3

  output:
  val input1

  script:
  """
  # create output folders
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash

  # create distance file
  mash sketch -s 10000 -o ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash/mash_sketch ${params.path_nextflow_dir}/5.ASSEMBLIES/*.fasta
  mash dist ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash/mash_sketch.msh ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash/mash_sketch.msh  | square_mash > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash/mash.tsv
  """
}

process SNPGWAS_Mash {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  pyseer --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --vcf ${params.path_nextflow_dir}/11.SNIPPY_MULTI/core.vcf --distances ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash/mash.tsv --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash/SNPGWAS_mash.txt
  """
}

process DBGWAS_Mash {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  pyseer --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --kmers ${params.path_nextflow_dir}/14.UNITIGS/unitigs.txt --uncompressed --output-patterns ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash/DBGWAS_mash_unitig_patterns.txt --distances ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash/mash.tsv --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash/DBGWAS_mash.txt
  """
}

process PanGWAS_Panaroo_Mash {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  pyseer --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --pres ${params.path_nextflow_dir}/13.Pangenome_analysis/13.1.Panaroo/gene_presence_absence.Rtab --distances ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash/mash.tsv --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.2.Mash/PanGWAS_Panaroo_mash.txt
  """
}

process SC_dist {
  conda "pyseer"

  input:
  val input1
  val input2
  val input3

  output:
  val input1

  script:
  """
  # create output folders
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.3.SC

  # create distances file
  sed 's/,/\t/g' ${params.path_nextflow_dir}/9.POPPUNK/GBS_GWAS_db/GBS_GWAS_db_core_clusters.csv > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.3.SC/SC_clusters.txt
  """
}

process SNPGWAS_SC {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  pyseer --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered_SC.txt --vcf ${params.path_nextflow_dir}/11.SNIPPY_MULTI/core.vcf --covariates ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.3.SC/SC_clusters.txt --use-covariates 2 --no-distances --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.3.SC/SNPGWAS_SC.txt
  """
}

process DBGWAS_SC {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  pyseer --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered_SC.txt --kmers ${params.path_nextflow_dir}/14.UNITIGS/unitigs.txt --uncompressed --output-patterns ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.3.SC/DBGWAS_SC_unitig_patterns.txt --covariates ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.3.SC/SC_clusters.txt --use-covariates 2 --no-distances --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.3.SC/DBGWAS_SC.txt
  """
}

process PanGWAS_Panaroo_SC {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  pyseer --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered_SC.txt --pres ${params.path_nextflow_dir}/13.Pangenome_analysis/13.1.Panaroo/gene_presence_absence.Rtab --covariates ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.3.SC/SC_clusters.txt --use-covariates 2 --no-distances --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.2.Subanalysis_1/16.2.3.SC/PanGWAS_Panaroo_SC.txt
  """
}

process PanGWAS_Roary_RAxML {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  # create output folders
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer/15.3.Subanalysis_2

  # Run Pyseer
  pyseer --lmm --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --pres ${params.path_nextflow_dir}/13.Pangenome_analysis/13.2.ROARY/gene_presence_absence.Rtab --similarity ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis/RAxML_phylogeny_K.tsv --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.3.Subanalysis_2/PanGWAS_Roary_RAxML.txt
  """
}

process PanGWAS_Panaroo_CLARC_RAxML {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  # create output folders
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer/15.3.Subanalysis_2

  # Make gene presence/absence input (convert CLARC .csv output to Rtab)
  awk 'BEGIN {FS=OFS=","} NR==1 {\$1="Gene" \$1} 1' ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC/Output/clarc_results/clarc_condensed_presence_absence.csv | awk '
   {
       for (i=1; i<=NF; i++)  {
           a[NR,i] = \$i
       }
   }
   NF>p { p = NF }
   END {
       for (i=1; i<=p; i++) {
           for (j=1; j<=NR; j++) {
               printf "%s%s", a[j,i], (j==NR ? "\\n" : "\\t")
           }
       }
   }
   ' FS=',' OFS='\\t' > ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC/Output/clarc_results/clarc_condensed_presence_absence.Rtab


  # Run Pyseer
  pyseer --lmm --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --pres ${params.path_nextflow_dir}/13.Pangenome_analysis/13.3.Panaroo_CLARC/Output/clarc_results/clarc_condensed_presence_absence.Rtab --similarity ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis/RAxML_phylogeny_K.tsv --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.3.Subanalysis_2/PanGWAS_Panaroo_CLARC_RAxML.txt
  """
}

process PanGWAS_Roary_CLARC_RAxML {
  cpus 8

  conda "pyseer"

  input:
  val base


  script:
  """
  # create output folders
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer
  mkdir -p ${params.path_nextflow_dir}/15.Pyseer/15.3.Subanalysis_2

  # Make gene presence/absence input (convert CLARC .csv output to Rtab)
  awk 'BEGIN {FS=OFS=","} NR==1 {\$1="Gene" \$1} 1' ${params.path_nextflow_dir}/13.Pangenome_analysis/13.4.ROARY_CLARC/Output/clarc_results/clarc_condensed_presence_absence.csv | awk '
   {
       for (i=1; i<=NF; i++)  {
           a[NR,i] = \$i
       }
   }
   NF>p { p = NF }
   END {
       for (i=1; i<=p; i++) {
           for (j=1; j<=NR; j++) {
               printf "%s%s", a[j,i], (j==NR ? "\\n" : "\\t")
           }
       }
   }
   ' FS=',' OFS='\\t' > ${params.path_nextflow_dir}/13.Pangenome_analysis/13.4.ROARY_CLARC/Output/clarc_results/clarc_condensed_presence_absence.Rtab

  # Run Pyseer
  pyseer --lmm --phenotypes ${params.path_nextflow_dir}/Files/phenotypes_filtered.txt --pres ${params.path_nextflow_dir}/13.Pangenome_analysis/13.4.ROARY_CLARC/Output/clarc_results/clarc_condensed_presence_absence.Rtab --similarity ${params.path_nextflow_dir}/15.Pyseer/15.1.Main_analysis/RAxML_phylogeny_K.tsv --cpu 8 > ${params.path_nextflow_dir}/15.Pyseer/15.3.Subanalysis_2/PanGWAS_Roary_CLARC_RAxML.txt
  """
}

workflow flow1 {
take: data
main:
    // Step 1: Start with Make_architecture process
    Make_architecture(data)

    // Step 2: Serotype_MLST and Trimmomatic both depend on Make_architecture's output
    Serotype_MLST(Make_architecture.out)
    Trimmomatic(Make_architecture.out)

    // Step 3: FastQC, SNIPPY, and ASSEMBLY depend on Trimmomatic's output
    FastQC(Trimmomatic.out)
    SNIPPY(Trimmomatic.out)
    ASSEMBLY(Trimmomatic.out)

    // Step 4: ANNOTATION, FastANI, and QUAST depend on ASSEMBLY's output
    ANNOTATION(ASSEMBLY.out)
    FastANI(ASSEMBLY.out)
    QUAST(ASSEMBLY.out)

    emit:
      QUAST.out[0]
      Serotype_MLST.out[0]
      FastQC.out[0]
      SNIPPY.out[0]
      ANNOTATION.out[0]
      FastANI.out[0]
}


workflow flow2 {
take:
  quast_data
  sero_MLST_data
  FASTQC_data
  SNIPPY_data
  ANNOTATION_data
  FastANI_data
main:
    Filter_out_QC(quast_data, sero_MLST_data, FASTQC_data, SNIPPY_data, ANNOTATION_data, FastANI_data)
emit: Filter_out_QC.out
}

workflow flow3 {
    take: data
    main:
    // Step 1: Run all processes that take 'data' as input in parallel
    POPPUNK(data)
    ASSIGN_CC(data)
    SNIPPY_MULTI(data)
    PANAROO(data)
    ROARY(data)
    UNITIG(data)

    // Step 2: FASTTREE and RAxML both depend on the output of SNIPPY_MULTI
    FASTTREE(SNIPPY_MULTI.out[0])
    RAxML_MLsearch(SNIPPY_MULTI.out[0])
    RAxML_bootstrap(SNIPPY_MULTI.out[0])
    RAxML_bipartition(RAxML_MLsearch.out[0], RAxML_bootstrap.out[0])

    // If an outgroup is specified, run RAxML with outgroup rooting
    if ( params.OUTGROUP) {
    SNIPPY_MULTI_OUTGROUP(data)
    RAxML_OUTGROUP(SNIPPY_MULTI_OUTGROUP.out[0])
    }

    // Step 3: PANAROO_CLARC depends on PANAROO output, and ROARY_CLARC depends on ROARY output
    PANAROO_CLARC(PANAROO.out[0])
    ROARY_CLARC(ROARY.out[0])

    emit:
    RAxML_bipartition.out[0]
    UNITIG.out[0]
    PANAROO.out[0]
    SNIPPY_MULTI.out[0]
    POPPUNK.out[0]
    FASTTREE.out[0]
    PANAROO_CLARC.out[0]
    ROARY_CLARC.out[0]
}

workflow flow_main {
    take:
    raxml_data
    unitig_data
    panaroo_data

    main:
    RAxML_dist(raxml_data, unitig_data, panaroo_data)
    SNPGWAS_RAxML(RAxML_dist.out[0])
    DBGWAS_RAxML(RAxML_dist.out[0])
    PanGWAS_Panaroo_RAxML(RAxML_dist.out[0])
}

workflow flow_sub1_Fasttree {
    take:
    fasttree_data
    unitig_data
    panaroo_data

    main:
    Fasttree_dist(fasttree_data, unitig_data, panaroo_data)
    SNPGWAS_Fasttree(Fasttree_dist.out[0])
    DBGWAS_Fasttree(Fasttree_dist.out[0])
    PanGWAS_Panaroo_Fasttree(Fasttree_dist.out[0])
}

workflow flow_sub1_Mash {
    take:
    snippy_data
    unitig_data
    panaroo_data

    main:
    Mash_dist(snippy_data,unitig_data, panaroo_data)
    SNPGWAS_Mash(Mash_dist.out[0])
    DBGWAS_Mash(Mash_dist.out[0])
    PanGWAS_Panaroo_Mash(Mash_dist.out[0])
}

workflow flow_sub1_SC {
    take:
    poppunk_data
    unitig_data
    panaroo_data

    main:
    SC_dist(poppunk_data, unitig_data, panaroo_data)
    SNPGWAS_SC(SC_dist.out[0])
    DBGWAS_SC(SC_dist.out[0])
    PanGWAS_Panaroo_SC(SC_dist.out[0])
}

workflow flow_sub2 {
    take:
    raxml_data
    panaroo_clarc_data
    roary_clarc_data

    main:
    RAxML_dist(raxml_data , panaroo_clarc_data, roary_clarc_data)
    PanGWAS_Roary_RAxML(RAxML_dist.out[0])
    PanGWAS_Panaroo_CLARC_RAxML(RAxML_dist.out[0])
    PanGWAS_Roary_CLARC_RAxML(RAxML_dist.out[0])
}

workflow {
channel.fromFilePairs("${params.path_nextflow_dir}/0.RAW_READS/*_{1,2}.fastq.gz")
       .map { it[0] }
       | flow1

flow2(flow1.out)
       | collect
       | flow3

  if( params.main ) {
    flow_main(flow3.out[0], flow3.out[1], flow3.out[2])
    println "Running the Main Analysis: SNP-GWAS, DB-GWAS and Pan-GWAS, with the population structure controlled using a RAxML phylogeny and the accessory genome defined per Panaroo."
  } else {
    println "Skipping the Main Analysis (SNP-GWAS, DB-GWAS and Pan-GWAS, with the population structure controlled using a RAxML phylogeny and the accessory genome defined per Panaroo)."
  }
  if( params.sub1 ) {
    flow_sub1_Fasttree(flow3.out[5], flow3.out[1], flow3.out[2])
    flow_sub1_Mash(flow3.out[3], flow3.out[1], flow3.out[2])
    flow_sub1_SC(flow3.out[4], flow3.out[1], flow3.out[2])
    println "Running the Sub-Analysis #1: SNP-GWAS, DB-GWAS and Pan-GWAS, with the population structure controlled using a Fastree phylogeny, Mash distances and Sequence Cluster (SC) classifications."
  }
  if( params.Fasttree ) {
    flow_sub1_Fasttree(flow3.out[5], flow3.out[1], flow3.out[2])
    println "Running SNP-GWAS, DB-GWAS and Pan-GWAS, with the population structure controlled using a Fastree phylogeny and the accessory genome defined per Panaroo.."
  }
  if( params.Mash ) {
    flow_sub1_Mash(flow3.out[3],flow3.out[1], flow3.out[2])
    println "Running SNP-GWAS, DB-GWAS and Pan-GWAS, with the population structure controlled using Mash distances and the accessory genome defined per Panaroo.."
  }
  if( params.SC ) {
    flow_sub1_SC(flow3.out[4], flow3.out[1], flow3.out[2])
    println "Running SNP-GWAS, DB-GWAS and Pan-GWAS, with the population structure controlled using Sequence Cluster (SC) classifications and the accessory genome defined per Panaroo.."
  }
  if( params.sub2 ) {
    flow_sub2(flow3.out[0], flow3.out[6], flow3.out[7])
    println "Running the Sub-Analysis #2: Pan-GWAS, with the accessory genome defined per Roary, CLARC-redefined Panaroo genes, and CLARC-redefined Roary genes."
  }
}
