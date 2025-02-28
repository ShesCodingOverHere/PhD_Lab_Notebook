# **All you need to know for basic bioinformatic analysis**

## **Command Line Basics Pipeline**

### **Introduction to the Command Line**
#### Goal: Understand what the command line is and how to access it.
- What is the command line?
- Difference between CLI (Command Line Interface) and GUI (Graphical User Interface).
- How to open the terminal:
  - **Windows**: Command Prompt (`cmd`), PowerShell, Windows Subsystem for Linux (WSL), Git Bash.
  - **Mac/Linux**: Terminal.

### **Navigation & File System**
#### Goal: Learn to move through directories and manipulate files.
🔹 **Commands to Learn**
- `pwd` → Print working directory.
- `ls` → List files in a directory.
- `cd` → Change directory (`cd ..` to move up, `cd /` for root).
- `mkdir <directory>` → Create a new directory.
- `rmdir <directory>` → Remove an empty directory.
- `touch <file>` (Linux/macOS) / `type nul > <file>` (Windows) → Create a new file.
- `rm <file>` / `del <file>` (Windows) → Delete a file.
- `mv <file> <destination>` → Move or rename a file.
- `cp <file> <destination>` → Copy a file.

🔹 **Exercises**
- Navigate to different directories using `cd`.
- Create a new directory and a text file inside it.
- Move and rename a file.
- Delete a file and a directory.

---

### **Viewing & Editing Files**
#### Goal: Learn how to read and edit text files in the command line.
🔹 **Commands to Learn**
- `cat <file>` → Display file contents.
- `less <file>` → View file with scrolling.
- `head <file>` / `tail <file>` → View first/last 10 lines.
- `nano <file>` / `vim <file>` / `notepad <file>` (Windows) → Edit a file.
- `echo "text" > file.txt` → Write text to a file.
- `echo "text" >> file.txt` → Append text to a file.

🔹 **Exercises**
- Create a text file and write some text into it.
- Use `cat`, `less`, `head`, and `tail` to view contents.
- Edit a file using `nano` or `vim`.

---

### **Permissions & Processes**
#### Goal: Understand how to manage file permissions and running processes.
🔹 **Commands to Learn**
- `chmod <permissions> <file>` → Change file permissions (Linux/macOS).
- `chown <user>:<group> <file>` → Change file owner (Linux/macOS).
- `ps` / `tasklist` (Windows) → View running processes.
- `kill <PID>` / `taskkill /PID <PID>` (Windows) → Kill a process.
- `top` / `htop` → Monitor system resources.

🔹 **Exercises**
- Change file permissions using `chmod`.
- View running processes using `ps` and `top`.

---

### **Searching & Redirection**
#### Goal: Learn how to filter and redirect command output.
🔹 **Commands to Learn**
- `grep "pattern" <file>` → Search for text in a file.
- `find <directory> -name "<file>"` → Search for files.
- `|` (Pipe) → Redirect output between commands.
- `>` (Redirect) → Save output to a file.
- `>>` (Append Redirect) → Append output to a file.
- `sort` → Sort output.
- `uniq` → Remove duplicate lines.

🔹 **Exercises**
- Search for a specific word inside a text file using `grep`.
- Redirect command output to a file using `>`.
- Combine `ls` and `sort` using pipes.

---

# **Basic Conda Pipeline**  

This simplified pipeline will help **beginners** learn **Conda**, focusing on essential commands for **environment management** and **package installation**.

---

## **Install Conda**
### **Install Miniconda**  
Miniconda is a lightweight version of Conda. Download from:  
🔗 [Miniconda Download](https://docs.conda.io/en/latest/miniconda.html)

**For Linux/Mac:**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

**For Windows:**  
Download and install **Miniconda** from the link above.

### **🔹 Step 2: Verify Installation**
```bash
conda --version
conda info
```

---

## **Conda Environment Basics**
### **Create an Environment**
```bash
conda create --name myenv python=3.9
```

### **Activate the Environment**
```bash
conda activate myenv
```

### **Deactivate the Environment**
```bash
conda deactivate
```

### **Remove an Environment**
```bash
conda remove --name myenv --all
```

---

## **Managing Packages**
### **Search for a Package**
```bash
conda search numpy
```

### **Install a Package**
```bash
conda install numpy
```

### **Update a Package**
```bash
conda update numpy
```

### **Remove a Package**
```bash
conda remove numpy
```

---

## **📜 Summary of Basic Conda Commands**
| **Task** | **Command** |
|----------|------------|
| Check Conda Version | `conda --version` |
| Create Environment | `conda create --name myenv python=3.9` |
| Activate Environment | `conda activate myenv` |
| Deactivate Environment | `conda deactivate` |
| Remove Environment | `conda remove --name myenv --all` |
| Search for Package | `conda search numpy` |
| Install Package | `conda install numpy` |
| Update Package | `conda update numpy` |
| Remove Package | `conda remove numpy` |
| Export Environment | `conda env export > environment.yml` |
| Recreate Environment | `conda env create -f environment.yml` |

---

# **AWK & Bioawk Pipeline**  

## **Introduction to AWK**
#### Goal: Understand the basics of AWK syntax and its role in text processing.
🔹 **Concepts to Learn**  
- What is AWK?  
- How AWK processes text line-by-line.  
- Basic AWK syntax:  
  ```bash
  awk '{print $0}' file.txt   # Print the entire line
  awk '{print $1, $2}' file.txt  # Print specific columns
  ```
- AWK execution methods:
  - One-liners (`awk 'pattern { action }' file.txt`).
  - Scripts (`awk -f script.awk file.txt`).

🔹 **Exercises**  
- Print the full content of a file.  
- Extract the first column of a tab-separated file.  
- Print the second and third column separated by a `-`.

---

### **AWK Data Selection & Filtering**  
#### Goal: Learn how to filter and process text efficiently.  
🔹 **Concepts to Learn**  
- Using `if` conditions in AWK:  
  ```bash
  awk '{ if ($3 > 50) print $0 }' file.txt
  ```
- Matching patterns using regular expressions:  
  ```bash
  awk '/pattern/ {print $0}' file.txt
  ```
- Filtering based on specific values:
  ```bash
  awk '$2 == "GeneX" {print $0}' genes.txt
  ```
  
🔹 **Exercises**  
- Print lines where column 3 is greater than 100.  
- Extract all lines where the second column contains "Drosophila".  
- Select lines where column 5 contains "ATG" or "TAA".

---

### **Processing Columns & Calculations in AWK**  
#### Goal: Perform mathematical operations and column manipulations.  
🔹 **Concepts to Learn**  
- Performing calculations:  
  ```bash
  awk '{print $1, $2, $3 * 2}' data.txt
  ```
- Calculating sums and averages:  
  ```bash
  awk '{sum+=$3} END {print "Total:", sum}' data.txt
  ```
- Counting occurrences of values:  
  ```bash
  awk '{count[$1]++} END {for (word in count) print word, count[word]}' file.txt
  ```

🔹 **Exercises**  
- Compute the sum of column 4.  
- Count the number of times each gene name appears in a dataset.  
- Calculate the average length of sequences from a FASTA-like file.

---

## **Introduction to Bioawk**
#### Goal: Learn how Bioawk extends AWK for bioinformatics.  
🔹 **Concepts to Learn**  
- What is Bioawk and how it differs from AWK?  
- Installing Bioawk:
  ```bash
  git clone https://github.com/lh3/bioawk.git
  cd bioawk && make
  ```
- Basic Bioawk commands:
  ```bash
  bioawk -c fastx '{print $name, length($seq)}' file.fasta
  ```
- Predefined data structures for:
  - **FASTA** (`-c fastx`)
  - **FASTQ** (`-c fastq`)
  - **GFF** (`-c gff`)
  - **SAM** (`-c sam`)

🔹 **Exercises**  
- Extract sequence lengths from a FASTA file.  
- Print the sequence name and first 10 bases.  
- Extract only protein-coding genes from a GFF file.

---

### **Advanced Bioawk - Genomic Data Processing**
#### Goal: Perform complex queries on biological datasets.  
🔹 **Concepts to Learn**  
- Extracting specific fields from GFF files:
  ```bash
  bioawk -c gff '$3 == "gene" {print $1, $4, $5, $9}' annotations.gff
  ```
- Filtering FASTQ files by quality scores:
  ```bash
  bioawk -c fastq '{if ($avgqual > 30) print $name}' file.fastq
  ```
- Parsing SAM files for mapped reads:
  ```bash
  bioawk -c sam '$3 != "*" {print $1, $3, $4}' alignments.sam
  ```

🔹 **Exercises**  
- Extract coding gene coordinates from a GFF file.  
- List all FASTQ reads with an average quality > 30.  
- Get the mapped read count from a SAM file.

---

### **Automating Bioinformatics Pipelines with AWK & Bioawk**
#### Goal: Combine AWK & Bioawk for batch processing.  
🔹 **Concepts to Learn**  
- Writing shell scripts with AWK/Bioawk:  
  ```bash
  # Extract sequence lengths from multiple FASTA files
  for file in *.fasta; do
      bioawk -c fastx '{print $name, length($seq)}' $file > ${file%.fasta}.lengths.txt
  done
  ```
- Integrating AWK into larger workflows:
  ```bash
  cat genes.txt | awk '$3 > 50 {print $1, $2}' | sort -k2,2n > filtered_genes.txt
  ```

🔹 **Exercises**  
- Automate extraction of gene lengths from multiple GFF files.  
- Write a script that filters short sequences from a set of FASTA files.  

---

## **🧬 Basic Genome Assembly Pipeline**
### **🛠 Tools Used**
1. **FastQC** - Quality check of raw reads.
2. **Fastp** - Trimming and filtering Illumina reads.
3. **Flye** - De novo assembly using long reads.
4. **Racon** - First round of polishing with long reads.
5. **Pilon** - Second round of polishing with short reads.
6. **BUSCO** - Genome assembly quality assessment.

---

## **📌 Step 1: Quality Control (FastQC)**
**Check the quality of raw reads before trimming.**
```bash
mkdir QC_reports
fastqc -o QC_reports raw_reads/*.fastq.gz
```

🔹 **Output:** Generates `.html` and `.zip` quality reports.

---

## **📌 Step 2: Read Trimming (Fastp)**
**Trim Illumina short reads for better polishing accuracy.**
```bash
fastp -i raw_reads/illumina_R1.fastq.gz -I raw_reads/illumina_R2.fastq.gz \
      -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz \
      --detect_adapter_for_pe --cut_right --json fastp_report.json --html fastp_report.html
```

🔹 **Output:** `trimmed_R1.fastq.gz` and `trimmed_R2.fastq.gz`

---

## **📌 Step 3: Genome Assembly (Flye)**
**Assemble genome using long reads (ONT/PacBio).**
```bash
flye --nano-hq raw_reads/long_reads.fastq.gz \
     --out-dir flye_assembly --genome-size 200m --threads 16
```

🔹 **Output:** `flye_assembly/assembly.fasta`

---

## **📌 Step 4: First Polishing Round (Racon)**
**Polish the assembly with long reads.**
```bash
minimap2 -ax map-ont flye_assembly/assembly.fasta raw_reads/long_reads.fastq.gz > long_reads.sam
racon -t 16 raw_reads/long_reads.fastq.gz long_reads.sam flye_assembly/assembly.fasta > racon_polished.fasta
```

🔹 **Output:** `racon_polished.fasta` (Improved assembly)

---

## **📌 Step 5: Second Polishing Round (Pilon)**
**Polish the genome with high-quality Illumina short reads.**
```bash
bwa index racon_polished.fasta
bwa mem -t 16 racon_polished.fasta trimmed_R1.fastq.gz trimmed_R2.fastq.gz > illumina.sam
samtools sort -o illumina.sorted.bam illumina.sam
samtools index illumina.sorted.bam

pilon --genome racon_polished.fasta --frags illumina.sorted.bam --output pilon_polished --threads 16
```

🔹 **Output:** `pilon_polished.fasta` (Final polished genome)

---

## **📌 Step 6: Quality Assessment (BUSCO)**
**Evaluate assembly completeness using the BUSCO database.**
```bash
busco -m genome -i pilon_polished.fasta -o busco_results -l drosophila_odb10
```

🔹 **Output:** `busco_results/short_summary.txt` (Completeness score)

---

### **🔁 Optional: Iterate for More Polishing**
- You can **repeat Racon and Pilon** multiple times for further refinement.

---

## **📜 Summary of Outputs**
| Step | Tool | Output |
|------|------|--------|
| Quality Control | FastQC | `QC_reports/*.html` |
| Read Trimming | Fastp | `trimmed_R1.fastq.gz`, `trimmed_R2.fastq.gz` |
| Assembly | Flye | `flye_assembly/assembly.fasta` |
| Polishing 1 | Racon | `racon_polished.fasta` |
| Polishing 2 | Pilon | `pilon_polished.fasta` |
| Quality Check | BUSCO | `busco_results/short_summary.txt` |

---

### 🚀 **Final Output:**  
The final **high-quality genome assembly** is `pilon_polished.fasta`, which has been polished using both **long and short reads** and validated using **BUSCO**.

### **QUAST: Genome Assembly Quality Assessment Step**  

**QUAST (Quality Assessment Tool for Genome Assemblies)** is used to evaluate genome assemblies by comparing them against a reference or in a reference-free mode.

---

## **Install QUAST**  
If you haven’t installed QUAST, use Conda:
```bash
conda install -c bioconda quast
```
Or install from source:
```bash
git clone https://github.com/ablab/quast.git
cd quast
python setup.py install
```

---

## **Run QUAST on an Assembly**
If you **have a reference genome**, use:
```bash
quast -o quast_output -r reference.fasta assembly.fasta
```

If you **don’t have a reference**, run QUAST in **de novo mode**:
```bash
quast -o quast_output assembly.fasta
```

### **🔹 Common Options**
- `-r reference.fasta` → Provide a reference genome (optional).
- `-o quast_output` → Specify an output directory.
- `--min-contig 500` → Set minimum contig length to 500 bp.
- `--threads 8` → Use 8 CPU threads for faster processing.

---

## **Interpret Results**
After running, open the **QUAST summary report**:

```bash
cat quast_output/report.txt
```

### **🔹 Key Metrics Explained**
| **Metric** | **Meaning** |
|------------|------------|
| **N50** | The contig length at which 50% of the genome is assembled. |
| **# contigs** | Total number of contigs in the assembly. |
| **Largest contig** | Length of the largest assembled contig. |
| **Total length** | Total base pairs in the assembly. |
| **GC content** | Percentage of GC bases. |
| **Misassemblies** | Structural errors compared to the reference. |

---

### **Full Code**

```
#!/bin/bash

# Set parameters
THREADS=16
GENOME_SIZE="200m"
BUSCO_DB="drosophila_odb10"

# Create output directories
mkdir -p QC_reports trimmed_reads assembly polishing busco_results

# Step 1: Quality Control (FastQC)
echo "Running FastQC..."
fastqc -o QC_reports raw_reads/*.fastq.gz

# Step 2: Read Trimming (Fastp)
echo "Trimming Illumina reads..."
fastp -i raw_reads/illumina_R1.fastq.gz -I raw_reads/illumina_R2.fastq.gz \
      -o trimmed_reads/trimmed_R1.fastq.gz -O trimmed_reads/trimmed_R2.fastq.gz \
      --detect_adapter_for_pe --cut_right --json QC_reports/fastp_report.json --html QC_reports/fastp_report.html

# Step 3: Genome Assembly (Flye)
echo "Running Flye for genome assembly..."
flye --nano-hq raw_reads/long_reads.fastq.gz \
     --out-dir assembly --genome-size $GENOME_SIZE --threads $THREADS

# Step 4: First Polishing Round (Racon)
echo "Running Racon polishing..."
minimap2 -ax map-ont assembly/assembly.fasta raw_reads/long_reads.fastq.gz > polishing/long_reads.sam
racon -t $THREADS raw_reads/long_reads.fastq.gz polishing/long_reads.sam assembly/assembly.fasta > polishing/racon_polished.fasta

# Step 5: Second Polishing Round (Pilon)
echo "Running Pilon polishing..."
bwa index polishing/racon_polished.fasta
bwa mem -t $THREADS polishing/racon_polished.fasta trimmed_reads/trimmed_R1.fastq.gz trimmed_reads/trimmed_R2.fastq.gz > polishing/illumina.sam
samtools sort -o polishing/illumina.sorted.bam polishing/illumina.sam
samtools index polishing/illumina.sorted.bam
pilon --genome polishing/racon_polished.fasta --frags polishing/illumina.sorted.bam --output polishing/pilon_polished --threads $THREADS

# Step 6: Quality Assessment (BUSCO)
echo "Running BUSCO for genome completeness check..."
busco -m genome -i polishing/pilon_polished.fasta -o busco_results -l $BUSCO_DB

# Step 7: QUAST
echo "Running QUAST..."
quast -o quast_output -r reference.fasta assembly.fasta

echo "Pipeline completed. Final genome assembly is in polishing/pilon_polished.fasta"
```

## **Annotation Pipeline**

```
#!/bin/bash

# Set parameters
THREADS=16
GENOME_FASTA="genome.fasta"
REPEAT_LIBRARY="repeats.lib"

# Create output directories
mkdir -p annotation repeatmasking helixer gff_processing

# Step 1: Build RepeatModeler Database
echo "Building RepeatModeler database..."
BuildDatabase -name repeat_db $GENOME_FASTA

# Step 2: Run RepeatModeler to identify repeat families
echo "Running RepeatModeler..."
RepeatModeler -database repeat_db -pa $THREADS -LTRStruct -dir annotation

# Step 3: Run RepeatMasker to mask repeats
echo "Running RepeatMasker..."
RepeatMasker -pa $THREADS -lib annotation/consensi.fa.classified -dir repeatmasking $GENOME_FASTA

# Step 4: Run Helixer for gene annotation
echo "Running Helixer..."
helixer predict --genome $GENOME_FASTA --output helixer/ --threads $THREADS

# Step 5: Process GFF using gffread
echo "Processing GFF annotations..."
gffread helixer/predictions.gff -T -o gff_processing/annotations.gtf

echo "Annotation pipeline completed. Results are in annotation, repeatmasking, helixer, and gff_processing directories."
```

### **Variant Calling Pipeline + R Code for Analysis**  

## **Read Quality Control (FastQC & Fastp)**
🔹 **Check raw reads quality**
```bash
mkdir QC_reports
fastqc -o QC_reports raw_reads/*.fastq.gz
```
🔹 **Trim low-quality reads**
```bash
fastp -i raw_reads/sample_R1.fastq.gz -I raw_reads/sample_R2.fastq.gz \
      -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz \
      --detect_adapter_for_pe --cut_right --json fastp_report.json --html fastp_report.html
```

---

## **Align Reads to Reference (BWA)**
🔹 **Index the reference genome**
```bash
bwa index reference.fasta
```
🔹 **Align reads to the reference**
```bash
bwa mem -t 16 reference.fasta trimmed_R1.fastq.gz trimmed_R2.fastq.gz > aligned.sam
```
🔹 **Convert, sort, and index BAM file**
```bash
samtools view -bS aligned.sam | samtools sort -o sorted.bam
samtools index sorted.bam
```

---

## **Mark Duplicates (GATK)**
```bash
gatk MarkDuplicates -I sorted.bam -O dedup.bam -M dedup_metrics.txt --REMOVE_DUPLICATES true
samtools index dedup.bam
```

---

## **Variant Calling**
### **🔹 Option 1: FreeBayes (haplotype-aware)**
```bash
freebayes -f reference.fasta dedup.bam > raw_variants.vcf
```
### **🔹 Option 2: GATK HaplotypeCaller**
```bash
gatk HaplotypeCaller -R reference.fasta -I dedup.bam -O raw_variants.vcf
```

---

## **Filter Variants**
🔹 **Filter out low-quality variants**
```bash
vcftools --vcf raw_variants.vcf --minQ 30 --recode --out filtered_variants
```

---

## **Variant Analysis in R**
### **🔹 Install & Load Required Packages**
```r
install.packages("vcfR")
install.packages("ggplot2")
library(vcfR)
library(ggplot2)
```

### **🔹 Read VCF File**
```r
vcf <- read.vcfR("filtered_variants.recode.vcf")
```

### **🔹 Summary of Variants**
```r
vcf_summary <- extract.info(vcf, "DP")
summary(vcf_summary)
```

### **🔹 Visualize Depth Distribution**
```r
ggplot(data.frame(Depth=vcf_summary), aes(x=Depth)) +
  geom_histogram(binwidth=5, fill="blue", color="black") +
  theme_minimal() +
  labs(title="Depth Distribution", x="Depth", y="Frequency")
```

### **🔹 Visualize Variant Types**
```r
variant_types <- table(getFIX(vcf)[, "ID"])
barplot(variant_types, las=2, col="steelblue", main="Variant Types")
```

---

## **📜 Summary of Outputs**
| Step | Tool | Output |
|------|------|--------|
| QC | FastQC | `QC_reports/*.html` |
| Trimming | Fastp | `trimmed_R1.fastq.gz, trimmed_R2.fastq.gz` |
| Alignment | BWA | `sorted.bam` |
| Deduplication | GATK | `dedup.bam` |
| Variant Calling | FreeBayes/GATK | `raw_variants.vcf` |
| Filtering | VCFtools | `filtered_variants.recode.vcf` |
| Analysis | R (vcfR) | Summary & Plots |

---

# **RNA-seq Pipeline Using Kallisto & R for Analysis**  
This pipeline covers **pseudo-alignment and quantification** with **Kallisto**, followed by **differential expression analysis** in **R** using `tximport`, `DESeq2`, and `ggplot2`.  

---

## **Quality Control (FastQC & Fastp)**
**Check read quality before processing:**
```bash
mkdir QC_reports
fastqc -o QC_reports raw_reads/*.fastq.gz
```

**Trim adapters and low-quality reads:**
```bash
fastp -i raw_reads/sample_R1.fastq.gz -I raw_reads/sample_R2.fastq.gz \
      -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz \
      --detect_adapter_for_pe --cut_right --json fastp_report.json --html fastp_report.html
```

---

## **Indexing the Transcriptome (Kallisto)**
Before quantification, build an index using a **transcriptome reference** (e.g., GENCODE, Ensembl).  
```bash
kallisto index -i transcriptome.idx reference_transcripts.fasta
```

---

## **Quantification with Kallisto**
Run Kallisto on paired-end RNA-seq data:
```bash
kallisto quant -i transcriptome.idx -o sample1_quant \
               -b 100 -t 8 \
               trimmed_R1.fastq.gz trimmed_R2.fastq.gz
```
🔹 `-b 100`: Bootstrap replicates for uncertainty estimation  
🔹 `-t 8`: Number of threads  

This will generate:
- `abundance.tsv`: Transcripts per million (TPM) & estimated counts  
- `abundance.h5`: HDF5 binary format for downstream analysis  

Repeat for all samples, organizing them into folders.

---

## **Import Kallisto Results into R**  
Use **tximport** to import results into **DESeq2** for differential expression analysis.

### **🔹 Install & Load Packages**
```r
install.packages("tximport")
install.packages("readr")
install.packages("DESeq2")
install.packages("ggplot2")

library(tximport)
library(readr)
library(DESeq2)
library(ggplot2)
```

---

## **Load & Prepare Data**  
Define file paths for multiple samples:
```r
sample_names <- c("sample1", "sample2", "sample3", "sample4")  
files <- file.path(sample_names, "abundance.tsv")  
names(files) <- sample_names  

txi <- tximport(files, type = "kallisto", txOut = TRUE)
```

---

## **Differential Expression Analysis (DESeq2)**
### **🔹 Load Metadata**
Create a metadata table with conditions:
```r
sample_info <- data.frame(
  row.names = sample_names,
  condition = c("control", "control", "treated", "treated")
)
```

### **🔹 Run DESeq2**
```r
dds <- DESeqDataSetFromTximport(txi, colData = sample_info, design = ~ condition)
dds <- DESeq(dds)
```

### **🔹 Extract Differentially Expressed Genes (DEGs)**
```r
res <- results(dds, contrast = c("condition", "treated", "control"))
res <- res[order(res$padj),]  # Sort by adjusted p-value
head(res)
```

---

## **Visualize Results**
### **🔹 Volcano Plot**
```r
res_df <- as.data.frame(res)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05)) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")
```

### **🔹 PCA Plot**
```r
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")
```

### **🔹 Heatmap of Top 50 Genes**
```r
library(pheatmap)
top50 <- head(order(res$padj), 50)
pheatmap(assay(vsd)[top50,], scale="row", cluster_cols=FALSE)
```

---

## **📜 Summary of Pipeline Steps**
| Step | Tool | Output |
|------|------|--------|
| QC | FastQC & Fastp | `QC_reports/*.html`, trimmed FASTQs |
| Indexing | Kallisto | `transcriptome.idx` |
| Quantification | Kallisto | `abundance.tsv` |
| Import & DE Analysis | R (tximport, DESeq2) | DEG tables, plots |
| Visualization | ggplot2, pheatmap | PCA, Volcano Plot, Heatmap |

---

## **Teaching Pipeline for scRNA-seq Analysis Using CellRanger & Seurat**  

This structured **teaching pipeline** will help learners process **single-cell RNA-seq (scRNA-seq) data** using **CellRanger for quantification** and **Seurat for downstream analysis**, including **quality control, normalization, clustering, and visualization**.

---

## **📌 Phase 1: Understanding the Workflow**
### **🔹 Goal:**
- Understand the **steps in scRNA-seq analysis**.
- Learn about **CellRanger** (for generating count matrices).
- Learn about **Seurat & Harmony** (for clustering & visualization).

### **🔹 Workflow Overview**
1. **CellRanger Preprocessing:**
   - Generate a **reference transcriptome**.
   - Run **CellRanger count** to generate a **gene-cell matrix**.
2. **Seurat Analysis:**
   - Load & create **Seurat objects**.
   - Perform **quality control & filtering**.
   - Normalize, scale, and find **highly variable genes**.
   - Perform **PCA, UMAP, and clustering**.
   - Use **Harmony for batch correction**.
   - Visualize clusters and **annotate cell types**.

---

## **📌 Phase 2: Preparing Data (CellRanger)**
### **🔹 Step 1: Install & Set Up CellRanger**
- Install **CellRanger** following official instructions:  
🔗 [CellRanger Download](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)

- Add **CellRanger to PATH** (Linux/Mac):
```bash
export PATH=/path/to/cellranger:$PATH
```
- Check installation:
```bash
cellranger --help
```

---

### **🔹 Step 2: Generate a Reference Transcriptome**
Before running **CellRanger count**, a reference transcriptome must be created:
```bash
cellranger mkref \
  --genome=dm6 \
  --fasta=/path/to/dm6.fa \
  --genes=/path/to/dm6.gtf \
  --nthreads=8 \
  --memgb=64
```
🔹 **Output:** `dm6/` reference directory

---

### **🔹 Step 3: Generate a Gene-Cell Count Matrix**
```bash
cellranger count \
  --id=sample_count \
  --transcriptome=/path/to/reference \
  --fastqs=/path/to/sample_fastqs \
  --sample=sample_name \
  --localcores=8 \
  --localmem=64
```
🔹 **Output:** `filtered_feature_bc_matrix/` directory (contains counts for Seurat)

---

## **📌 Phase 3: Seurat Analysis in R**
### **🔹 Step 4: Install Required Packages**
```r
install.packages("Seurat")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("BiocManager")
BiocManager::install("SeuratWrappers")
BiocManager::install("Nebulosa")
BiocManager::install("harmony")
```
```r
library(Seurat)
library(tidyverse)
library(ggplot2)
library(SeuratWrappers)
library(Nebulosa)
library(harmony)
```

---

### **🔹 Step 5: Load Data into Seurat**
```r
time8_data <- Read10X(data.dir = "/path/to/count_timepoint_8h/filtered_feature_bc_matrix")
seurat_8h <- CreateSeuratObject(time8_data, project = "8h")
```
Repeat for all time points (1d, 3d, 6d, 9d, 15d, 30d).

---

### **🔹 Step 6: Add Metadata (Replicates)**
```r
meta8 <- seurat_8h@meta.data
meta8$barcodes <- rownames(meta8)
bc_splits_8 <- data.frame(do.call('rbind', strsplit(as.character(meta8$barcodes), '-', fixed=TRUE)))
bc_splits_8$group <- ifelse(bc_splits_8$X2 == 1, "Rep1", "Rep2")
seurat_8h@meta.data$Replicate <- bc_splits_8$group
```
🔹 Repeat for all time points.

---

### **🔹 Step 7: Normalize Data**
```r
seurat_8h <- NormalizeData(seurat_8h)
```
🔹 Repeat for all time points.

---

### **🔹 Step 8: Merge All Time Points**
```r
merged_data <- merge(
  x = seurat_8h, 
  y = list(seurat_1d, seurat_3d, seurat_6d, seurat_9d, seurat_15d, seurat_30d),
  add.cell.ids = c("t_8h", "t_1d", "t_3d", "t_6d", "t_9d", "t_15d", "t_30d")
)
```

---

### **🔹 Step 9: Quality Control**
```r
merged_data$percent.mt <- PercentageFeatureSet(merged_data, pattern = '^mt:')
merged_data_filtered <- subset(
  merged_data, 
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25
)
```

---

### **🔹 Step 10: Identify Variable Genes**
```r
merged_data_filtered <- FindVariableFeatures(merged_data_filtered, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(merged_data_filtered), 10)
print(top10)
```

---

### **🔹 Step 11: Scale Data**
```r
merged_data_filtered <- ScaleData(merged_data_filtered)
```

---

### **🔹 Step 12: PCA & Clustering**
```r
merged_data_filtered <- RunPCA(merged_data_filtered, npcs = 50)
ElbowPlot(merged_data_filtered, ndims = 50)

merged_data_filtered <- RunUMAP(merged_data_filtered, dims = 1:30)
merged_data_filtered <- FindNeighbors(merged_data_filtered, dims = 1:30)
merged_data_filtered <- FindClusters(merged_data_filtered, resolution = 0.5)
```

---

### **🔹 Step 13: Batch Correction with Harmony**
```r
merged_data_harmony <- RunHarmony(merged_data_filtered, group.by.vars = 'orig.ident')
merged_data_harmony <- RunUMAP(merged_data_harmony, reduction = 'harmony', dims = 1:30)
merged_data_harmony <- FindClusters(merged_data_harmony, resolution = 0.5)
```

---

### **🔹 Step 14: Visualization**
```r
DimPlot(merged_data_harmony, reduction = 'umap', label = TRUE)
DimPlot(merged_data_harmony, reduction = 'umap', split.by = 'orig.ident')
```

---

### **🔹 Step 15: Find Cell Type Markers**
```r
all_markers <- FindAllMarkers(merged_data_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all_markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
```

---

### **🔹 Step 16: Annotate Clusters**
```r
new_cluster_ids <- c("Neuron", "Glia", "Microglia", "Astrocytes", "Oligodendrocytes")
names(new_cluster_ids) <- levels(merged_data_harmony)
merged_data_harmony <- RenameIdents(merged_data_harmony, new_cluster_ids)
DimPlot(merged_data_harmony, reduction = "umap", label = TRUE, pt.size = 0.5)
```

---

## **📜 Summary of Pipeline Steps**
| Step | Tool | Description |
|------|------|-------------|
| Preprocessing | CellRanger | Generates gene-cell count matrices |
| Load Data | Seurat | Create Seurat objects |
| QC | Seurat | Filter out low-quality cells |
| Normalization | Seurat | Normalize expression levels |
| Clustering | Seurat | Perform PCA, UMAP, and clustering |
| Batch Correction | Harmony | Remove batch effects |
| Visualization | Seurat | UMAP, FeaturePlots, and heatmaps |
| Annotation | Seurat | Assign cell types |

---

## **Manuals for each program**

- [Kallisto](https://pachterlab.github.io/kallisto/manual)
- [DESeq](https://bioconductor.org/packages/3.21/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [QUAST](https://quast.sourceforge.net/docs/manual.html)
- [Samtools and bcftools](https://www.htslib.org/doc/#manual-pages)
- [Bash](https://www.gnu.org/software/bash/manual/bash.html)
- [Awk](https://www.gnu.org/software/gawk/manual/gawk.html)
- [Bioawk](https://bioinformaticsworkbook.org/Appendix/Unix/bioawk-basics.html#gsc.tab=0)
- [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/index.html)
- [Fastp](https://github.com/OpenGene/fastp/blob/master/README.md)
- [Flye](https://gensoft.pasteur.fr/docs/Flye/2.9/USAGE.html)
- [Pilon and racon](https://denbi-nanopore-training-course.readthedocs.io/en/latest/artic/pilon/)
- [BUSCO](https://busco.ezlab.org/busco_userguide.html)
- [bwa](https://bio-bwa.sourceforge.net)
- [freebayes](https://open.bioqueue.org/home/knowledge/showKnowledge/sig/freebayes)
- [GATK](https://gatk.broadinstitute.org/hc/en-us/categories/360002310591)
- [minimap2](https://lh3.github.io/minimap2/minimap2.html)
- [Seurat official tutorials](https://satijalab.org/seurat/)
- [Scanpy tutorials](https://scanpy-tutorials.readthedocs.io/en/latest/)
- [Monocle3 trajectory inference](https://cole-trapnell-lab.github.io/monocle3/)
