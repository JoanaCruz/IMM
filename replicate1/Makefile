# Alignment reference
AREF=../Genome

# Main directories
HUMAN_DIR=../Human_DNA
LIST_DIR=../Listeria_RNA
MIT_DIR=../Mitocondrial_DNA
CONTROL_DIR=../control.vs
REPLICANT_DIR=.

READS_DIR =/home/data/reads
BOWTIE_DIR=/home/data/joanacruz/bowtie2-2.0.0-beta7
SCRIPT_DIR=../Scripts
PLOTM_DIR=../Scripts/Graphics

## Auxiliary directories
# Results
TRIM_DIR=$(HUMAN_DIR)/trim
SORTED_DIR=$(REPLICANT_DIR)/alignsort
PAIRED_DIR=$(HUMAN_DIR)/paired_reads
COUNTS_DIR=$(HUMAN_DIR)/counts
COUNTS_CONTROL_DIR=$(CONTROL_DIR)/counts
COUNTS_REP_DIR=$(REPLICANT_DIR)/counts
DESEQ_DIR=$(REPLICANT_DIR)/DESeq
GO_DIR=$(REPLICANT_DIR)/GOStats
WEB_DIR=$(REPLICANT_DIR)/web.site


# Graphics
PLOT_DIR=$(REPLICANT_DIR)/Plots
PLOT_DESEQ=$(PLOT_DIR)/DESeq.analysis
PLOT_GO=$(PLOT_DIR)/Genes.expression

GTF_FILE := $(AREF)/hg19_refseq.gtf

all: gostats

gostats: $(addprefix $(GO_DIR)/GOsummary_, $(notdir $(subst _L1_1.fq.gz,.csv.gz, $(wildcard $(READS_DIR)/LM*_L1_1.fq.gz)))))

deseq: $(addprefix $(DESEQ_DIR)/DESeq_, $(notdir $(subst L1_1.fq.gz,pvalue.csv.gz, $(wildcard $(READS_DIR)/LM*_L1_1.fq.gz))))

concatenate: $(addprefix $(COUNTS_REP_DIR)/,$(notdir $(subst _L1_1.fq.gz,_control.txt.gz, $(wildcard $(READS_DIR)/LM*_L1_1.fq.gz))))


# Make sequences alignment with the human genome (using Bowtie2 with --sensitive option and with all paired ends that aligned written at .sam)
bowtie : $(REPLICANT_DIR)/Control_HeLa_L1_replicate_b.fastq | $(SORTED_DIR)/.d
	$(BOWTIE_DIR)/bowtie2 -x $(AREF)/hg19 --sensitive $(REPLICANT_DIR)/Control_HeLa_L1_replicate_b.fastq | gzip > $(SORTED_DIR)/Control_HeLa_L1_replicate_b.sam.gz
	@ echo "Alignment done."

sort: $(REPLICANT_DIR)/Control_HeLa_L1_replicate.bam
## Sort the BAM file by read name and convert to SAM (to use in HTSeq)
	@ samtools sort -n $(REPLICANT_DIR)/Control_HeLa_L1_replicate.bam $(SORTED_DIR)/Control_HeLa_L1_replicate_sorted
	@ samtools view $(SORTED_DIR)/Control_HeLa_L1_replicate_sorted.bam -h > $(SORTED_DIR)/Control_HeLa_L1_replicate_sorted.sam
	@ gzip $(SORTED_DIR)/Control_HeLa_L1_replicate_sorted.sam
	
	
# Uses HTSeq-count script to count how many reads map to each feature (being a feature a range of positions on a chromosome)
htseq : $(SORTED_DIR)/Control_HeLa_L1_replicate_b.sam
	#@ gunzip $(SORTED_DIR)/Control_HeLa_L1_replicate_b.sam.gz
	@ python -m HTSeq.scripts.count -q -m union -s no -i gene_id $(SORTED_DIR)/Control_HeLa_L1_replicate_b.sam $(AREF)/hg19_refseq.gtf | gzip > $(COUNTS_REP_DIR)/Control_HeLa_L1_replicate_b.txt.gz
	#@ gzip $(SORTED_DIR)/Control_HeLa_L1_replicate_sorted.sam $(COUNTS_DIR)/Control_HeLa_L1_replicate.txt

	
# Concatenate files from HTSeq-count (control vs (LM1 or LM2) )
$(COUNTS_REP_DIR)/%_control.txt.gz : $(COUNTS_DIR)/%_L1.txt.gz $(COUNTS_REP_DIR)/.d
	@ gunzip $(COUNTS_DIR)/$*_L1.txt.gz $(COUNTS_DIR)/Control_HeLa_L1.txt.gz $(COUNTS_DIR)/Control_HeLa_L1_replicate.txt.gz
	@ paste $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_DIR)/Control_HeLa_L1_replicate.txt $(COUNTS_DIR)/$*_L1.txt | cut -f1,2,4,6- | sed '1 i\Control_HeLa\tControl_HeLa_replicate\tLM1_HeLa_'$*''> $(COUNTS_REP_DIR)/$*_control.txt
	@ gzip $(COUNTS_REP_DIR)/$*_control.txt $(COUNTS_DIR)/$*_L1.txt $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_DIR)/Control_HeLa_L1_replicate.txt



concatenate-control : $(COUNTS_DIR)/Control_HeLa_L1.txt.gz $(COUNTS_DIR)/Control_HeLa_L1_replicate.txt.gz $(COUNTS_REP_DIR)/.d
	@ gunzip $(COUNTS_DIR)/Control_HeLa_L1.txt.gz $(COUNTS_DIR)/Control_HeLa_L1_replicate.txt.gz
	@ paste $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_DIR)/Control_HeLa_L1_replicate.txt | cut -f1,2,4- | sed '1 i\Control_HeLa\tControl_HeLa_replicate' > $(COUNTS_REP_DIR)/Control_control_replicate.txt
	@ gzip $(COUNTS_REP_DIR)/Control_control_replicate.txt $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_DIR)/Control_HeLa_L1_replicate.txt
	

# Do differencial expression analysis between control and wildtype & between control and mutant (using DESeq)
$(DESEQ_DIR)/DESeq_%_pvalue.csv.gz : $(COUNTS_REP_DIR)/%_control.txt.gz $(DESEQ_DIR)/.d
	@ gunzip $(COUNTS_REP_DIR)/$*_control.txt.gz
	@ Rscript $(SCRIPT_DIR)/DEseq_replicate1.R $(COUNTS_REP_DIR)/$*_control.txt $* $(DESEQ_DIR)/
	@ gzip $(DESEQ_DIR)/DESeq_$*_*.csv $(COUNTS_REP_DIR)/$*_control.txt

	
plot-deseq : $(PLOT_DIR)/.d
	@ gunzip $(DESEQ_DIR)/DESeq_*_ALL.csv.gz
	@ Rscript $(SCRIPT_DIR)/analysis.deseq.rep.R $(DESEQ_DIR)/ $(PLOT_DIR)/
	@ gzip  $(DESEQ_DIR)/DESeq_*_ALL.csv
	
	

deseq-control : $(COUNTS_REP_DIR)/Control_control_replicate.txt.gz $(DESEQ_DIR)/.d
	@ gunzip $(COUNTS_REP_DIR)/Control_control_replicate.txt.gz
	@ Rscript $(SCRIPT_DIR)/DEseq.R $(COUNTS_REP_DIR)/Control_control_replicate.txt Control $(DESEQ_DIR)/
	@ gzip $(DESEQ_DIR)/DESeq_Control_*.csv $(COUNTS_REP_DIR)/Control_control_replicate.txt

	
# Create a directory (use DIR/.d)
%/.d:
	@ mkdir -p $(@D)
	@ touch $@

