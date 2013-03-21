# Alignment reference
AREF=../Genome

# Main directories
HUMAN_DIR=../Human_DNA
LIST_DIR=../Listeria_RNA
MIT_DIR=../Mitocondrial_DNA
CONTROL_DIR=../control.vs
REPLICATE_DIR=.

READS_DIR =/home/data/reads
BOWTIE_DIR=/home/data/joanacruz/bowtie2-2.0.0-beta7
SCRIPT_DIR=../Scripts
PLOTM_DIR=../Scripts/Graphics

## Auxiliary directories
# Results
TRIM_DIR=$(HUMAN_DIR)/trim
SORTED_DIR=$(REPLICATE_DIR)/alignsort
PAIRED_DIR=$(HUMAN_DIR)/paired_reads
COUNTS_DIR=$(HUMAN_DIR)/counts
COUNTS_CONTROL_DIR=$(CONTROL_DIR)/counts
COUNTS_REP_DIR=$(REPLICATE_DIR)/counts
DESEQ_DIR=$(REPLICATE_DIR)/DESeq
GO_DIR=$(REPLICATE_DIR)/GOStats
WEB_DIR=$(REPLICATE_DIR)/web.site


# Graphics
PLOT_DIR=$(REPLICATE_DIR)/Plots
PLOT_DESEQ=$(PLOT_DIR)/DESeq.analysis
PLOT_GO=$(PLOT_DIR)/Genes.expression

GTF_FILE := $(AREF)/hg19_refseq.gtf

all: plot-deseq

deseq-replicate : $(DESEQ_DIR)/DESeq_rep.vs.control_pvalue.csv.gz

counts: $(addprefix $(COUNTS_REP_DIR)/, $(notdir $(subst .fastq,.txt.gz, $(wildcard $(REPLICATE_DIR)/*.fastq)))))

bowtie: $(addprefix $(SORTED_DIR)/, $(notdir $(subst .fastq.gz,.sam.gz, $(wildcard $(REPLICATE_DIR)/*.fastq.gz)))))


# Make sequences alignment with the human genome (using Bowtie2 with --sensitive option and with all paired ends that aligned written at .sam)
$(SORTED_DIR)/rep%_Control_HeLa_L1.sam.gz : $(REPLICATE_DIR)/rep%_Control_HeLa_L1.fastq.gz | $(SORTED_DIR)/.d
	@ gunzip rep$*_Control_HeLa_L1.fastq.gz
	$(BOWTIE_DIR)/bowtie2 -x $(AREF)/hg19 --sensitive $(REPLICATE_DIR)/rep$*_Control_HeLa_L1.fastq | gzip > $(SORTED_DIR)/rep$*_Control_HeLa_L1.sam.gz
	@ echo "Alignment done."

sort: $(REPLICATE_DIR)/Control_HeLa_L1_replicate.bam
## Sort the BAM file by read name and convert to SAM (to use in HTSeq)
	@ samtools sort -n $(REPLICATE_DIR)/Control_HeLa_L1_replicate.bam $(SORTED_DIR)/Control_HeLa_L1_replicate_sorted
	@ samtools view $(SORTED_DIR)/Control_HeLa_L1_replicate_sorted.bam -h > $(SORTED_DIR)/Control_HeLa_L1_replicate_sorted.sam
	@ gzip $(SORTED_DIR)/Control_HeLa_L1_replicate_sorted.sam
	
	
# Uses HTSeq-count script to count how many reads map to each feature (being a feature a range of positions on a chromosome)
$(COUNTS_REP_DIR)/rep%_Control_HeLa_L1.txt.gz : $(SORTED_DIR)/rep%_Control_HeLa_L1.sam.gz $(COUNTS_REP_DIR)/.d
	@ gunzip $(SORTED_DIR)/rep$*_Control_HeLa_L1.sam.gz
	python -m HTSeq.scripts.count -q -m union -s no -i gene_id $(SORTED_DIR)/rep$*_Control_HeLa_L1.sam $(AREF)/hg19_refseq.gtf | gzip > $(COUNTS_REP_DIR)/rep$*_Control_HeLa_L1.txt.gz
	@ gzip $(SORTED_DIR)/rep$*_Control_HeLa_L1.sam
	

# concatenate public data with our control
$(COUNTS_REP_DIR)/rep%_Control.txt.gz : $(COUNTS_REP_DIR)/rep%_Control_HeLa_L1.txt.gz
	@ gunzip $(COUNTS_DIR)/Control_HeLa_L1.txt.gz $(COUNTS_REP_DIR)/rep$*_Control_HeLa_L1.txt.gz
	@ paste $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_REP_DIR)/rep$*_Control_HeLa_L1.txt | cut -f1,2,4- | sed '1 i\Control_HeLa\tControl_HeLa_rep'$*'' > $(COUNTS_REP_DIR)/rep$*_Control.txt
	@ gzip $(COUNTS_REP_DIR)/rep$*_Control.txt $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_REP_DIR)/rep$*_Control_HeLa_L1.txt


# concatenate all controls from public data with our control (one single file) & concatenate each public control with each other (to compare similarty between public controls)
$(COUNTS_REP_DIR)/ALL_rep_Control.txt.gz:
	@ gunzip $(COUNTS_DIR)/Control_HeLa_L1.txt.gz  $(COUNTS_REP_DIR)/rep*_Control_HeLa_L1.txt.gz
	@ paste $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_REP_DIR)/rep1_Control_HeLa_L1.txt $(COUNTS_REP_DIR)/rep2_Control_HeLa_L1.txt | cut -f1,2,4,6- | sed '1 i\Control_HeLa\tControl_HeLa_rep1\tControl_HeLa_rep2' > $(COUNTS_REP_DIR)/ALL_rep_Control.txt
	@ gzip $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_REP_DIR)/rep*_Control_HeLa_L1.txt $(COUNTS_REP_DIR)/ALL_rep_Control.txt	
	
	

# Do differential analysis between replicates and control
$(DESEQ_DIR)/DESeq_rep.vs.control_pvalue.csv.gz : $(COUNTS_REP_DIR)/ALL_rep_Control.txt.gz $(DESEQ_DIR)/.d
	@ gunzip $(COUNTS_REP_DIR)/ALL_rep_Control.txt.gz
	@ Rscript $(SCRIPT_DIR)/DEseq_replicate2.R $(COUNTS_REP_DIR)/ALL_rep_Control.txt $(DESEQ_DIR)/
	@ gzip $(DESEQ_DIR)/DESeq_* $(COUNTS_REP_DIR)/ALL_rep_Control.txt


# Plot the dispersion associated with DESeq analysis
plot-deseq : $(PLOT_DIR)/.d
	@ gunzip $(DESEQ_DIR)/DESeq_*_ALL.csv.gz
	@ Rscript $(SCRIPT_DIR)/analysis.deseq.rep.R $(DESEQ_DIR)/ $(PLOT_DIR)/
	@ gzip  $(DESEQ_DIR)/DESeq_*_ALL.csv

	
# Create a directory (use DIR/.d)
%/.d:
	@ mkdir -p $(@D)
	@ touch $@

