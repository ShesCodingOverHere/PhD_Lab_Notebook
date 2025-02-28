# *D. affinis*: Genome line


## Assembly
```
#!/bin/bash
#SBATCH --job-name=hifiasm  # Job name
#SBATCH --partition=kucg      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL,BEGIN     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tconway@ku.edu   # Where to send mail
#SBATCH --ntasks=1          # Run on a single CPU
#SBATCH --mem=100gb           # Job memory request
#SBATCH --time=7-24:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=hifiasm_M%j.log  # Standard output and error log


module load conda
conda activate hifiasm

hifiasm -o /home/t043c581/scratch/data/DaffM.asm -t 32 /home/t043c581/scratch/data/D-affinis_M_F_HiFi.fastq.gz
```

## Remove bacterial contamination

```
blastn -task megablast -query foo.fa -remote -db nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -out foo.megablast
```

- When this is finished, go through the file and remove bacteria/virus using NCBI

### quast, busco, and other data I will undoubtedly be asked for...

```
# busco
nohup busco -i /hdd/Taylor/data/foo.fa -o /hdd/Taylor/data/BUSCO.foo -l /hdd/Taylor/data/diptera_odb10 -m genome --auto-lineage-euk -f &

# quast
python3 ../software/quast-5.2.0/quast.py /hdd/Taylor/data/foo.fa
```

| Species (sex)| Length| # of Scaffolds| N50 | BUSCO (complete)| BUSCO (single)| BUSCO (dup)|
|-------------:|:-----:|:-------------:|:---:|:---------------:|:-------------:|:----------:|
| _D. affinis_(male)| ~210Mb| 56| 16.8Mb| 99.1%| 98.0%| 1.1%|
| _D. affinis_(female)| ~200Mb| 57| 25.9Mb| 99.1%| 98.2%| 0.9%|

## Align fasta to reference via mummer to rename scaffolds locally

```
nucmer --maxgap=500 -mincluster=100 reference.fasta query.fasta
delta-filter -q -r out.delta > foo.filter
show-coords -B foo.filter > foo.coords
```

- Using the Rscript "chrom_mapping.R", check PID for each alignment and rename accordingly using:

```
sed 's/scaffold_to_be_renamed/rename_it_here/g' foo.fa >temp1
sed 's/scaffold_to_be_renamed/rename_it_here/g' temp1 >temp2
sed 's/scaffold_to_be_renamed/rename_it_here/g' temp2 >temp1

and so on...
```

## Coverage
### Method 1

Make .bam files with minimap and plot coverage using chromosome quotient method

#### Code:
```
# on the linux
# align female longreads to male reference
minimap2 -t 8 -ax map-hifi /hdd/Taylor/data/DaffM_hifi_nobact_.fa /hdd/Taylor/data/D-affinis_M_F_HiFi.fastq.gz |samtools view -bS > /hdd/Taylor/data/DalgF_hifi.bam

# align male longreads to male reference
minimap2 -t 8 -ax map-hifi /hdd/Taylor/data/DalgM_hifi_nobact_v2.fa /hdd/Taylor/data/D-algonquin_M_HiFi.fastq.gz |samtools view -bS > /hdd/Taylor/data/DaztM_hifi.bam

# sort bams and get coverage
samtools sort DalgM_hifi.bam >DalgM_hifi_sorted.bam
samtools coverage DalgM_hifi_sorted.bam >DalgM_hifi.cov
samtools sort DalgF_hifi.bam >DalgF_hifi_sorted.bam
samtools coverage DalgF_hifi_sorted.bam >DalgF_hifi.cov
```

Using RStudio:
```
# Load required libraries
library(ggplot2)
library(cowplot)
library(ggrepel)

# Load data
aztmale <- read.delim("/Users/conway/Desktop/CurrentWorkingDatasets/DalgM_hifi.cov", header = TRUE)
aztfem <- read.delim("/Users/conway/Desktop/CurrentWorkingDatasets/DalgF_hifi.cov", header = TRUE)

# Rename columns for clarity
colnames(algmale) <- c("scaffold", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")
colnames(algfem) <- c("scaffold", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")

# Normalize data
algmale$normalized <- algmale$numreads / sum(algmale$numreads)
algfem$normalized <- algfem$numreads / sum(algfem$numreads)

# Calculate scaffold sizes
algmale$size <- algmale$endpos - algmale$startpos
algfem$size <- algfem$endpos - algfem$startpos

# Compute log2 ratio of female/male normalized reads
algmale$log2fem_male <- log2(algfem$normalized / algmale$normalized)

# List of scaffolds to highlight
#highlight_scaffolds <- c("ptg000009l", "ptg000018l", "ptg000019l",
#                         "ptg000022l", "ptg000028l", "ptg000029l",
#                         "ptg000038l", "ptg000046l", "ptg000098l")

# Add a column to indicate if a scaffold should be highlighted
#algmale$highlight <- ifelse(algmale$scaffold %in% highlight_scaffolds, "yes", "no")

ggplot(algmale[algmale$size > 10000, ], aes(x = size, y = log2fem_male)) +
  # Add points with conditional coloring
  geom_point(aes(color = highlight), size = 3, alpha = 0.7) +
  # Use a log10 scale for the x-axis
  scale_x_continuous(trans = "log10", labels = scales::comma_format()) +
  # Add horizontal reference lines
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_hline(yintercept = log2(0.3), color = "red", linetype = "dashed") +
  geom_hline(yintercept = log2(2), color = "green", linetype = "dashed") +
  # Label only points below the red line
  geom_text_repel(aes(label = ifelse(log2fem_male < log2(0.3), scaffold, "")),
                  size = 5, 
                  box.padding = 0.4, 
                  point.padding = 0.4, 
                  max.overlaps = Inf) +
  # Customize color scale
  scale_color_manual(values = c("no" = "blue", "yes" = "blue"), guide = "none") +
  # Axis labels
  labs(
    x = "Scaffold Size (log10 scale)",
    y = "Log2(Female/Male) Normalized Reads",
    title = "Algonquin Normalized Coverage Using Chrom Quotient Concept"
  ) +
  # Set y-axis limits
  ylim(-7.5, 2.5) +
  # Apply clean theme
  theme_cowplot(font_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Center and bold title
    axis.title = element_text(size = 16),  # Larger axis titles
    axis.text = element_text(size = 14),   # Larger axis text
    panel.grid.major = element_line(color = "grey90", size = 0.5),  # Subtle grid lines
    legend.position = "none"               # Remove legend
  )
```

![alt text](image.png)

Points aligning around 1 should be X-linked, and points aligning around 0 should be autosomal. Anything under the red line is putative Y-linked. This is because females have 2 Xs when males have 1, and males have 1 Y while females have zero. You expect (when log2 tranformed) for the autosomal reads to cancel out and end up around 0.

### Method 2

Indexcov

#### Code:
```
goleft indexcov --directory ../indexcov_Dalg *.bam
```
![alt text](indexcov-indexcov-depth-ptg000009l.png) ![alt text](indexcov-indexcov-depth-ptg000018l.png) ![alt text](indexcov-indexcov-depth-ptg000019l.png) ![alt text](indexcov-indexcov-depth-ptg000022l.png) ![alt text](indexcov-indexcov-depth-ptg000028l.png) ![alt text](indexcov-indexcov-depth-ptg000029l.png) ![alt text](indexcov-indexcov-depth-ptg000038l.png) ![alt text](indexcov-indexcov-depth-ptg000046l.png) ![alt text](indexcov-indexcov-depth-ptg000098l.png) ![alt text](indexcov-indexcov-roc-ptg000009l.png) ![alt text](indexcov-indexcov-roc-ptg000018l.png) ![alt text](indexcov-indexcov-roc-ptg000019l.png) ![alt text](indexcov-indexcov-roc-ptg000022l.png) ![alt text](indexcov-indexcov-roc-ptg000028l.png) ![alt text](indexcov-indexcov-roc-ptg000029l.png) ![alt text](indexcov-indexcov-roc-ptg000038l.png) ![alt text](indexcov-indexcov-roc-ptg000046l.png) ![alt text](indexcov-indexcov-roc-ptg000098l.png)
## Putative Y scaffolds

### Confirmation of putative Y scaffolds
| Scaffold| Length| # of genes (unmasked)| # of genes (masked)| samtools coverage| indexcov| PCR|
|-------------:|:-----:|:-------------:|:---:|:---------------:|:-------------:|:----------:|
| | | | | | | |
| | | | | | | |
| | | | | | | |
| | | | | | | |
| | | | | | | |

## Locate and mask repeats with repeatmodelor and repeatmasker on the cluster

```
#!/bin/bash
#SBATCH --job-name=RMaffM  # Job name
#SBATCH --partition=kucg      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL,BEGIN     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tconway@ku.edu   # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=4-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=RMaffM_%j.log  # Standard output and error log

module load repeatmodeler
module load repeatmasker/4.0.9

#usage: sbatch RepeatMasker.args.job <fasta> <prefix>

cd $SCRATCH
mkdir RMaffinisM_pilon2

echo "STARTING"
cd RMaffinisM_pilon2
cp $HOME/$1 .

BuildDatabase -name $2 -engine ncbi $1

RepeatModeler -engine ncbi -pa 8 -database $2

RepeatMasker -pa 8 -gff -lib $2-families.fa -dir MaskerOutput$2 $1

echo done
```

- With this data, you can look at Y-linked repeat families.

## Annotate with helixer

- Go to https://www.plabipd.de/helixer_main.html
- Input fasta
- Change "Select Lineage-specific mode" to invertebrate
- Enter GFF label name and email address
- Submit job and wait
- grep gene foo.gff > genes.txt
- Import genes.txt into spreadsheet
- Convert gff to fasta using gffread (see below for code)
- blastx Y_transcripts.fa
- look up each gene on flybase and fill out spreadsheet

```
# gffread
gffread DaffM_unmasked_helixer.gff -g ../Daffinis/DaffM_hifi_nobact.fa -w DaffM_transcripts.fasta
```

## Renaming Transcripts
### Using Nilanjan's method

```
# gffread
gffread your_transcripts.gff -g genomic_reference.fasta -w your_transcripts.fasta​

# make a database
makeblastdb -in Dmel_translation_clean.fasta -dbtype prot -out dmel_protein_database

# Run blastx
blastx -query DaztM_masked_fixednames_transcripts.fa -db dmel_protein_database -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -num_threads 4 -out blast_azt_transcripts.txt

# retrieve protein sequences
cut -f2 blast_azt_transcripts.txt | sort | uniq > aztM_best_hit_proteins.list
seqtk subseq ../Dmel_translation_clean.fasta aztM_best_hit_proteins.list > aztM_best_hit_proteins.fa

makeblastdb -in DaztM_masked_fixednames_transcripts.fa -dbtype nucl -out DaztM_transcripts_db

tblastn -query aztM_best_hit_proteins.fa -db DaztM_transcripts_db -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -out aztM_blast_reciprocal.txt

awk '{print $1"\t"$2}' aztM_blast_reciprocal.txt > aztM_forward_hits.txt
awk '{print $2"\t"$1}' aztM_blast_reciprocal.txt > aztM_reciprocal_hits.txt
sort aztM_forward_hits.txt aztM_reciprocal_hits.txt | sed  's/-P[ABCDEFGHIJKLMNOPQRSTUVWXYZ]//g'  | uniq > reciprocal_best_hits.txt
awk '{if(a[$2]++){print $1"\t"$2"."a[$2]}else{print $0}}' reciprocal_best_hits.txt > aztM_RBH.txt
awk '{print $0 ".1"}' aztM_RBH.txt > aztM_RBH2.txt
awk '{sub(/\.1$/, "", $1); print}' RBH.txt  >RBH3.txt

gawk 'NR==FNR { mapping[$1] = $2; next } { for (key in mapping) gsub(key, mapping[key]) } 1' RBH2.txt ../DaztM_masked_helixer.gff > aztM_temp.gff

gawk 'NR==FNR { mapping[$1] = $2; next } { for (key in mapping) gsub(key, mapping[key]) } 1' RBH3.txt temp.gff > renamed.gff
```




## RNAseq Analysis
```
# Load required libraries
library(DESeq2)
library(ggplot2)
library(tidyr)
library(dplyr)

# Load counts data
counts <- read.table("/Users/conway/Desktop/CurrentWorkingDatasets/Daffinis/all.counts", 
                     header = TRUE)

# Set row names to the exact values in the first column
rownames(counts) <- counts[, 1]

# Remove the first column (now redundant in row names)
counts <- counts[, -1]

# Filter genes with low counts
keep <- rowSums(counts) > 20  # Keep genes with total counts > 20
cts <- round(counts[keep, ], digits = 0)

# Define experimental conditions
condition <- factor(c("XO", "XO", "XO", "XY", "XY", "XY"))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = DataFrame(condition), 
                              design = ~ condition)

# Normalize and filter
dds <- estimateSizeFactors(dds)
dds <- dds[rowSums(counts(dds, normalized = TRUE) > 20) >= 2, ]  # Filter low counts

# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Convert results to a data frame
df <- as.data.frame(res)
df$gene <- rownames(df)  # Add gene names as a column

# Add scaffold information by removing the last 8 characters from the gene column
df <- df %>%
  mutate(scaffold = substr(gene, 1, nchar(gene) - 9))  # Extract scaffold names

# Add another column `scaff` with the last 10 characters of the scaffold column
df <- df %>%
  mutate(scaff = substr(scaffold, nchar(scaffold) - 9, nchar(scaffold))) %>%  # Keep last 10 characters
  mutate(scaff = gsub("[_l]", "", scaff))  # Remove "_" or "l" from the name

# Define the putative Y scaffolds
putative_Y_scaffolds <- c("ptg000009", "ptg000018", "ptg000019", 
                          "ptg000022", "ptg000028", "ptg000029", 
                          "ptg000038", "ptg000046", "ptg000098")

# Add a new column `group` to label putative Y scaffolds
df <- df %>%
  mutate(group = ifelse(scaff %in% putative_Y_scaffolds, "Putative Y", "Other"))

# Create the boxplot using `group` for coloring and jittered points
ggplot(df, aes(x = reorder(scaff, log2FoldChange, FUN = median), 
               y = log2FoldChange, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +  # Create boxplot
  geom_jitter(aes(color = group, 
                  shape = ifelse(padj < 0.05, "Significant", "Not Significant")), 
              size = 1.5, width = 0.2, alpha = 0.6) +  # Add jittered points
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # Add reference line
  scale_fill_manual(values = c("Putative Y" = "deeppink", "Other" = "blue")) +  # Define colors for boxplots
  scale_color_manual(values = c("Putative Y" = "deeppink4", "Other" = "darkblue")) +  # Define colors for points
  scale_shape_manual(values = c("Significant" = 17, "Not Significant" = 16)) +  # Define shapes
  labs(shape = "Significance",  # Rename the shape legend
       fill = "Group") +  # Optionally rename the fill legend
  theme_minimal() +  # Minimal theme for clean plots
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels
    legend.position = "right"  # Position legend
  ) +
  xlab("Scaffold") +
  ylab("Log2(Fold Change)") +
  ggtitle("Log2(Fold Change) Per Scaffold for XY/XO Gene Expression")

```

![alt text](image-1.png)


## Gene coverage
```
hisat2-build ../mergedAssembly20241101.transcripts.fa DaffM_idx

hisat2 -x DaffM_idx -1 ../../../KUResearch/Datafiles/assemblies/shortreads/affinisM_frag_R1.fq -2 ../../../KUResearch/Datafiles/assemblies/shortreads/affinisM_frag_R2.fq -S aligned_reads.sam

hisat2 -x DaffM_idx -1 ../../../KUResearch/Datafiles/assemblies/shortreads/affinisF_frag_R1.fq.gz -2 ../../../KUResearch/Datafiles/assemblies/shortreads/affinisF_frag_R2.fq.gz -S fem_aligned_reads.sam

samtools view -bS aligned_reads.sam | samtools sort -o male_aligned_reads.sorted.bam
samtools view -bS fem_aligned_reads.sam | samtools sort -o fem_aligned_reads.sorted.bam

samtools index male_aligned_reads.sorted.bam
samtools index fem_aligned_reads.sorted.bam

gffread ../mergedAssembly20241101.gff -T -o mergedasm.gtf

samtools depth -a male_aligned_reads.sorted.bam > male_coverage.txt
samtools depth -a fem_aligned_reads.sorted.bam > fem_coverage.txt

samtools coverage male_aligned_reads.sorted.bam > male_aligned_reads.sorted.cov
samtools coverage fem_aligned_reads.sorted.bam > fem_aligned_reads.sorted.cov

awk -F'\t' '{gsub("ID=", "", $9); print}' /Users/conway/Desktop/CurrentWorkingDatasets/Daffinis/mergedAssembly20241101.gff > cleaned_gff.txt

awk -F'\t' 'NR==FNR {gsub("ID=", "", $9); $9 = $9 ".1"; map[$9] = $1; next} $1 in map {print $0, map[$1]}' OFS="\t" /Users/conway/Desktop/CurrentWorkingDatasets/Daffinis/mergedAssembly20241101.gff /Users/conway/Desktop/CurrentWorkingDatasets/Daffinis/male_coverage.txt > /Users/conway/Desktop/CurrentWorkingDatasets/Daffinis/male_depth_annotated.txt
```

![alt text](image-3.png)
![alt text](image-4.png)
![alt text](image-5.png)
![alt text](image-6.png)
![alt text](image-7.png)


## threshold checking 

I'm looking at what is on the red line. Making primers for those to see if Y linked. Primers below.

![alt text](image-2.png)

## Primers
### Scaffold 30
DaffY_ptg30_1_F - TCACCAAGTCAGAGCTTCATTGA

DaffY_ptg30_1_R - CGGTATGTTCTGATCTTTGAAGCA

### Scaffold 32
DaffY_ptg32_1_F - CGGATCCAGACTTCTTCGTCC

DaffY_ptg32_1_R - GACGGTTTCGCCAAATACGTC

## SexFindR

[SexFindR Link](https://sexfindr.readthedocs.io/en/latest/Step%200.%20Mapping%20and%20variant%20calling.html)

```
# create index
bowtie2-build ../../DaffM_hifi_nobact.fa fugu &> index_outerr.txt &

# align reads
bowtie2 -x fugu -1 ~/Desktop/KUResearch/Datafiles/assemblies/shortreads/affinisM_frag_R1.fq -2 ~/Desktop/KUResearch/Datafiles/assemblies/shortreads/affinisM_frag_R2.fq -S affinisM.sam &> bt2_affinisM_outerr.txt &

samtools faidx ../../DaffM_hifi_nobact.fa


```




## DNApipeTE

[DNApipeTE github](https://github.com/clemgoub/dnaPipeTE)

```
seqtk seq -a ../DaffM_hifi_nobact.fa | seqtk sample -s100 - 1000000 > simulated_reads.fastq


```
