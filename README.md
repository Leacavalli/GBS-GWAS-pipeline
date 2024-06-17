# Pipeline Description
![alt text](https://github.com/Leacavalli/GBS-GWAS-pipeline/blob/main/Bioinfo_flowchart.png)

This Nextflow pipeline was built to integrate the following steps:
<br>   

| Steps  | Description  | Tool |
| ------------- | ------------- | ------------- |
| 1 | Trim sequencing primers  | [Trimmomatic](https://github.com/timflutre/trimmomatic)  |
| 2 | Identify contaminating reads  | [Kraken2](https://github.com/DerrickWood/kraken2) |
| 3 | Remove contaminating reads  | [KrakenTools](https://github.com/jenniferlu717/KrakenTools) |
| 4 | Reads Quality Control | [FastQC](https://github.com/s-andrews/FastQC)  |
| 5 | Determine the serotype and ST of isolates | [SRST2](https://github.com/katholt/srst2), with [GBS-SG](https://github.com/swainechen/GBS-SBG)  |
| 6 | Mapping to reference genome and variant calling | [Snippy](https://github.com/tseemann/snippy)  |
| 7 | Assemble reads  | [Unicycler](https://github.com/rrwick/Unicycler)  |
| 8 | Genome Annotation | [Prokka](https://github.com/tseemann/prokka)  |


# Pipeline Architecture

Nextflow_pipeline\
&nbsp;&nbsp;| --- scripts\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- nextflow.config\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Nextflow_pipeline.nf\
&nbsp;&nbsp;| --- Files\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- AP018935.1.fa\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- GBS-SBG.fasta\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Streptococcus_agalactiae.fasta\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- profiles_csv\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- phylogeny_distance.py\
&nbsp;&nbsp;|--- 0.RAW_READS\
&nbsp;&nbsp;|--- 1.TRIMMED_READS\
&nbsp;&nbsp;|--- 2.KRAKEN_REPORTS\
&nbsp;&nbsp;|--- 3.FILTERED_READS\
&nbsp;&nbsp;|--- 4.FASTQC\
&nbsp;&nbsp;|--- 5.SEROTYPE_MLST\
&nbsp;&nbsp;|--- 6.SNIPPY\
&nbsp;&nbsp;|--- 7.ASSEMBLIES\
&nbsp;&nbsp;|--- 8.ANNOTATION\
&nbsp;&nbsp;|--- 9.Post_Pipeline\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- 9.1.MultiQC\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- 9.2.Snippy_core\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- 9.3.Pangenome_analysis\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- 9.4.FastTree\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- 9.5.PopPunk\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- 9.6.Pyseer

# Running the Pipeline
The pipeline contains the steps that can be run for each sample individually. From read trimming and filtering to genome annotation.
## 1. Clone this Github
```
git clone https://github.com/Leacavalli/GBS-GWAS-pipeline.git
```
## 2. Download Raw reads from NCBI (Optional)
Note: SraAccList.txt contains the list of accessions you want
```
cd 0.RAW_READS
sbatch -p shared -t 1-00:00 --mem=100000 --wrap="prefetch --option-file SraAccList.txt"
for i in *RR*
do
  sbatch -p shared  -t 0-00:10 --mem=10000 --wrap="fasterq-dump --split-files $i/*sra"
done
```
## 3. Install Java v.11 for Nextflow
```
module load python
conda create -n java      # Only run once to set up; All following times, skip this
conda activate java
conda install openjdk=11  # Only run once to set up; All following times, skip this
java -version             # Check java version
```
## 4. Install Nextflow
```
# Download the executable package:
wget -qO- https://get.nextflow.io | bash   
# Make the binary executable on your system by running:  
chmod +x nextflow  
# Add executable to $PATH                           
nano  ~/.bashrc      
# add line at the end of file:
export PATH="<path_to_nextflow>:$PATH"
# exist file with ctlX and Y
source ~/.bashrc      
```
## 5. Run Nextflow
The following options are required to run the pipeline for each isolates:
```
Usage: nextflow run <Nextflow_script> [args...]
Arguments:
        -C                     <Nextflow configuration file (Available in scripts directory)>
        --path_nextflow_dir    <Path to your Nextflow Directory>
        --sample               <The Sample ID/Name>
```
For example, the command to run the pipeline for a single isolate is the following:  
```
nextflow run Nextflow_pipeline.nf -C nextflow.config --path_nextflow_dir <Path/to/your/Nextflow/Directory> --sample <Sample ID>
```
The following loop can be used  to submit jobs for all isolates:
```
for i in ../0.RAW_READS/*_1.fastq.gz
do
 base=$(echo $i | cut -d '/' -f3-3 | cut -d '_' -f1-1)
 sbatch -p shared -t 0-03:00 --mem=100000 --wrap="nextflow  -C nextflow.config run Nextflow_pipeline.nf--path_nextflow_dir <Path/to/your/Nextflow/Directory> --sample $base"
done
```

# Post-pipeline steps
Post-pipeline steps consists of bioinformatic steps that need to be run on the collection of pipeline outputs from all samples. Starting from merge QC reports from all samples and filtering bad quality controls, and ending with the GWAS.

| Steps  | Description  | Tool |
| ------------- | ------------- | ------------- |
| 1 | Combine QC results | [MultiQC](https://github.com/MultiQC/MultiQC)  |
| 2 | Produce a whole genome SNP alignment | [Snippy-core](https://github.com/tseemann/snippy)  |
| 3 | Pangenome Analysis | [Panaroo](https://github.com/gtonkinhill/panaroo)  |
| 4 | Call unitigs | [Unitig-counter](https://github.com/bacpop/unitig-counter) |
| 5 | Produce a phylogeny | [FastTree](http://www.microbesonline.org/fasttree/) |
| 6.1 | Conduct a GWAS using SNPs | [pyseer](https://github.com/mgalardini/pyseer) |
| 6.2 | Conduct a GWAS using unitigs (i.e. DBGWAS) | [pyseer](https://github.com/mgalardini/pyseer) |
| 6.3 | Conduct a GWAS using gene presence/absence (i.e. PanGWAS)| [pyseer](https://github.com/mgalardini/pyseer) |

```
path_nextflow_dir=$(pwd)
```
## 1. Genome quality control
Tool: MultiQC\
Output: One QC report for all samples in our dataset
```
# Install MultiQC
conda create -n multiqc
source activate multiqc
conda install multiqc

# Run MultiQC
cd $path_nextflow_dir/9.Post_Pipeline/9.1.MultiQC/
sbatch -p shared  -t 0-03:00 --mem=100000 --wrap="multiqc $path_nextflow_dir/4.FASTQC"
```
## 2. Merge Snippy results
Tool: snippy-core\
Output:  Whole genome SNP alignment
```
# Install Snippy
conda create -n multiqc
source activate snippy
conda install multiqc

# Run Snippy
cd $path_nextflow_dir/9.Post_Pipeline/9.2.Snippy_core
cp $path_nextflow_dir/6.SNIPPY/SRR5063682/ref.fa ./
sbatch -p shared -t 0-02:00 --mem=100000 --wrap="snippy-core --ref 'ref.fa' $path_nextflow_dir/6.SNIPPY/SRR*"
```
## 3. Pangenome analysis
Tool: panaroo\
Output: Gene Presence/Absence matrix for input into Pyseer
```
# Install Panaroo
conda create -n Panaroo python=3.9
source activate Panaroo
conda install -c conda-forge -c bioconda -c defaults 'panaroo>=1.3'

# Run Panaroo
cd $path_nextflow_dir/9.Post_Pipeline/9.3.Pangenome_analysis
sbatch -p shared -t 0-05:00 --mem=100000 --wrap="panaroo -i $path_nextflow_dir/8.ANNOTATION/*.gff -o $path_nextflow_dir/9.Post_Pipeline/9.3.Pangenome_analysis --clean-mode strict --alignment pan --aligner mafft"
```
## 4. Unitig calling
Tool: unitig-counter\
Output: Unitig list file for input into Pyseer
```
# Install unitig-counter
conda create -n unitig-counter
source activate unitig-counter
conda install unitig-counter

# Create input file
cd $path_nextflow_dir/9.Post_Pipeline/9.4.Unitig_counter
echo -e "ID\tPath" > strain_list.txt
for i in $path_nextflow_dir/7.ASSEMBLIES/*.fasta
do
  base=$(echo "$i" | cut -d '/' -f12-12 | cut -d '.' -f1)
  echo -e "$base\t$i"
done >> strain_list.txt

# Run unitig-counter
sbatch -p shared -c 8 -t 0-05:00 --mem=100000 --wrap="unitig-counter -strains strain_list.txt -output Unitig_counter -nb-cores 8"

# Extend short unitigs leftwards and rightwards by following the neightbouring nodes in the de Bruijn graph
cd Unitig_counter
sbatch -p shared -c 8 -t 0-05:00 --mem=100000 --wrap="cdbg-ops extend --graph graph --unitigs unitigs.txt > extended.txt"
```
## 5. Phylogeny - FastTree
Tool: FastTree
```
# install FastTree
curl -O http://www.microbesonline.org/fasttree/FastTree.c
gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm

# Run
cd $path_nextflow_dir/9.Post_Pipeline/9.5.FastTree
sbatch -p shared -t 0-02:00 --mem=100000 --wrap="FastTree -gtr -nt $path_nextflow_dir/9.Post_Pipeline/9.2.Snippy_core/core.full.aln > phylo.aln"
```
## 6. GWAS
Tool: pyseer
### Build similarity matrix from phylogeny
```
cd $path_nextflow_dir/9.Post_Pipeline/9.6.Pyseer
#installation
conda -n create pyseer
source activate pyseer
conda install pyseer

#create phylogeny tsv file
sbatch -p shared -t 0-02:00 --mem=100000 --wrap="python $path_nextflow_dir/Files/phylogeny_distance.py --lmm $path_nextflow_dir/9.Post_Pipeline/9.5.FastTree/phylo.aln > phylogeny_K.tsv"
```
### 6.1.SNP-based GWAS
```
pyseer --phenotypes phenotypes --vcf $path_nextflow_dir/9.Post_Pipeline/9.2.Snippy_core/core.vcf --lmm --similarity phylogeny_K.tsv > EOD_LOD_SNPs.txt
```
### 6.2.Pan-GWAS
```
pyseer --phenotypes phenotypes --pres $path_nextflow_dir/9.Post_Pipeline/9.3.Pangenome_analysis/gene_presence_absence.Rtab --lmm --similarity phylogeny_K.tsv > EOD_LOD_COGs.txt
```
### 6.3.Unigit-based GWAS (DBGWAS)
```
pyseer --lmm --phenotypes phenotypes --kmers $path_nextflow_dir/9.Post_Pipeline/9.4.Unitig_counter/Unitig_counter/unitigs.txt --similarity phylogeny_K.tsv --output-patterns unitig_patterns.txt --cpu 8 > phenotype_unitigs.txt
```
