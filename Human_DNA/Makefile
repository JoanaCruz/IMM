# Alignment reference
AREF=../Genome

# Main directories
HUMAN_DIR=.
LIST_DIR=Listeria_RNA
MIT_DIR=Mitocondrial_DNA
READS_DIR =/data/reads
BOWTIE_DIR=../../bowtie2-2.0.0-beta7
SCRIPT_DIR=../Scripts
PLOTM_DIR=../Scripts/Graphics

## Auxiliary directories
# Results
TRIM_DIR=$(HUMAN_DIR)/trim
SORTED_DIR=$(HUMAN_DIR)/alignsort
PAIRED_DIR=$(HUMAN_DIR)/paired_reads
COUNTS_DIR=$(HUMAN_DIR)/counts
DESEQ_DIR=$(HUMAN_DIR)/DESeq
GO_DIR=$(HUMAN_DIR)/GOStats

# Graphics
PLOT_DIR=$(HUMAN_DIR)/Plots
PLOT_DESEQ_DIR=$(PLOT_DIR)/DESeq.analysis
PLOT_GO_DIR=$(PLOT_DIR)/Genes.per.goterm
PLOT_GO_LM1=$(PLOT_GO_DIR)/LM1
PLOT_GO_LM2=$(PLOT_GO_DIR)/LM2

GTF_FILE := $(AREF)/hg19_refseq.gtf

all: gostats

gostats: $(addprefix $(GO_DIR)/,$(subst LM1_HeLa, GOsummary, $(notdir $(subst _L1_1.fq.gz,.csv.gz, $(wildcard $(READS_DIR)/LM1_*_1.fq.gz)))))

deseq: $(addprefix $(DESEQ_DIR)/,$(subst LM1, DESeq_LM, $(notdir $(subst L1_1.fq.gz,pvalue.csv.gz, $(wildcard $(READS_DIR)/LM1_*_1.fq.gz)))))

concatenate: $(addprefix $(COUNTS_DIR)/,$(subst LM1, LM, $(notdir $(subst _L1_1.fq.gz,.txt.gz, $(wildcard $(READS_DIR)/LM1_*_1.fq.gz)))))

counts: $(addprefix $(COUNTS_DIR)/,$(addsuffix .txt.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

paired: $(addprefix $(PAIRED_DIR)/,$(addsuffix _paired.sam.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

sorted: $(addprefix $(SORTED_DIR)/,$(addsuffix _sorted.sam.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

trim: $(addprefix $(TRIM_DIR)/,$(addsuffix _trimmed.gz,$(notdir $(basename $(basename $(wildcard $(READS_DIR)/*.fq.gz))))))


# Make bases trim (delete first 15 bases from the beginning of the sequences)
$(TRIM_DIR)/%_trimmed.gz : $(READS_DIR)/%.fq.gz | $(TRIM_DIR)/.d
	@ python $(SCRIPT_DIR)/trim.py $(READS_DIR)/ $* 15 $(TRIM_DIR)/ 
	@ echo "Bases trimmed."


## Sort the SAM file by read name (to use in HTSeq)
# Convert .bam sorted to .sam
$(SORTED_DIR)/%_sorted.sam.gz : $(SORTED_DIR)/%_sorted.bam
	@ samtools view $(SORTED_DIR)/$*_sorted.bam -h > $(SORTED_DIR)/$*_sorted.sam
	@ gzip $(SORTED_DIR)/$*_sorted.sam
	@ echo "SAM file sorted."
# Convert .sam to .bam and sort
$(SORTED_DIR)/%_sorted.bam : $(SORTED_DIR)/%.sam.gz
	@ samtools view -bS $(SORTED_DIR)/$*.sam.gz | samtools sort -n - $(SORTED_DIR)/$*_sorted

# Make sequences alignment with the human genome (using Bowtie2 with --sensitive option and with all paired ends that aligned written at .sam)
$(SORTED_DIR)/%.sam.gz :  $(TRIM_DIR)/%_1_trimmed.gz $(TRIM_DIR)/%_2_trimmed.gz | $(SORTED_DIR)/.d
	$(BOWTIE_DIR)/bowtie2 -x $(AREF)/hg19 --sensitive -1 $(TRIM_DIR)/$*_1_trimmed.gz -2 $(TRIM_DIR)/$*_2_trimmed.gz | gzip > $(SORTED_DIR)/$*.sam.gz
	@ echo "Alignment done."


# Filter paired end reads which only one align
$(PAIRED_DIR)/%_paired.sam.gz : $(SORTED_DIR)/%_sorted.sam.gz | $(PAIRED_DIR)/.d
	@ python $(SCRIPT_DIR)/paired_end.py $(SORTED_DIR)/ $* $(PAIRED_DIR)/
	@ echo "Reads filtered."


# Uses HTSeq-count script to count how many reads map to each feature (being a feature a range of positions on a chromosome)
$(COUNTS_DIR)/%.txt.gz : $(PAIRED_DIR)/%_paired.sam.gz $(GTF_FILE) | $(COUNTS_DIR)/.d
	@ gunzip $(PAIRED_DIR)/*.sam.gz 
	@ python -m HTSeq.scripts.count -m union -s no -i gene_id $(PAIRED_DIR)/$*_paired.sam $(AREF)/hg19_refseq.gtf | gzip > $@
	@ gzip $(PAIRED_DIR)/*.sam


# Concatenate files from HTSeq-count (same timepoint diferent condition)
$(COUNTS_DIR)/LM_HeLa_%.txt.gz : $(COUNTS_DIR)/LM1_HeLa_%_L1.txt.gz $(COUNTS_DIR)/LM2_HeLa_%_L1.txt.gz 
	@ gunzip $(COUNTS_DIR)/LM*_HeLa_$*_L1.txt.gz
	@ join -t "`/bin/echo -e '\t'`" --head $(COUNTS_DIR)/LM1_HeLa_$*_L1.txt $(COUNTS_DIR)/LM2_HeLa_$*_L1.txt > $(COUNTS_DIR)/LM_HeLa_$*.txt
	@ sed -i '1s/^/LM1_HeLa_'$*'\tLM2_HeLa_'$*'\n/' $(COUNTS_DIR)/LM_HeLa_$*.txt
	@ gzip $(COUNTS_DIR)/*.txt


# Do differencial expression analysis (using DESeq)
$(DESEQ_DIR)/DESeq_LM_HeLa_%_pvalue.csv.gz: $(COUNTS_DIR)/LM_HeLa_%.txt.gz | $(DESEQ_DIR)/.d
	@ gunzip $(COUNTS_DIR)/LM_HeLa_$*.txt.gz
	@ Rscript $(SCRIPT_DIR)/DEseq.R $(COUNTS_DIR)/LM_HeLa_$*.txt $* $(DESEQ_DIR)/
	@ gzip $(DESEQ_DIR)/*.csv


# GO annotation (using GOstats)
$(GO_DIR)/GOsummary_%.csv.gz: $(DESEQ_DIR)/DESeq_LM_HeLa_%_pvalue.csv.gz | $(GO_DIR)/.d
	@ gunzip $(DESEQ_DIR)/DESeq_LM_HeLa_$*_pvalue.csv.gz
	@ Rscript $(SCRIPT_DIR)/GOStats_def.R $(DESEQ_DIR)/DESeq_LM_HeLa_$*_pvalue.csv $* $(GO_DIR)/
	@ gzip $(GO_DIR)/*.csv
	@ gzip $(DESEQ_DIR)/DESeq_LM_HeLa_$*_pvalue.csv


plot: $(addprefix $(COUNTS_DIR)/, $(notdir $(subst 20_L1_1.fq.gz,all_timep.txt.gz, $(wildcard $(READS_DIR)/LM*_HeLa_20_L1_1.fq.gz))))

GO_FILES := $(basename $(wildcard $(GO_DIR)/*.csv.gz))

# Concatenate files with the same condition (LM1 or LM2), do heatmap with the most differential expressed genes and do plots of gene expression variation for each one of the goterms.
$(COUNTS_DIR)/LM%_HeLa_all_timep.txt.gz : $(COUNTS_DIR)/LM%_HeLa_*_L1.txt.gz $(COUNTS_DIR)/Control_HeLa_L1.txt.gz | $(PLOT_GO_LM1)/.d $(PLOT_DESEQ_DIR)/.d $(PLOT_GO_LM2)/.d 
	@ gunzip $(COUNTS_DIR)/Control_HeLa_L1.txt.gz $(COUNTS_DIR)/LM$*_HeLa_*_L1.txt.gz
	@ paste $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_DIR)/LM$*_HeLa_20_L1.txt $(COUNTS_DIR)/LM$*_HeLa_60_L1.txt $(COUNTS_DIR)/LM$*_HeLa_120_L1.txt $(COUNTS_DIR)/LM$*_HeLa_240_L1.txt | cut -f1,2,4,6,8,10- | sed '1 i\Control_HeLa\tLM'$*'_HeLa_20\tLM'$*'_HeLa_60\tLM'$*'_HeLa_120\tLM'$*'_HeLa_240' > $(COUNTS_DIR)/LM$*_HeLa_all_timep.txt
	@ gzip $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_DIR)/LM$*_HeLa_*_L1.txt
# Gene expression graphics
	@ gunzip $(GO_DIR)/*.csv.gz
	@ Rscript $(PLOTM_DIR)/heatmap.R  $(COUNTS_DIR)/LM$*_HeLa_all_timep.txt $(PLOT_GO_DIR)/ LM$*
	@ Rscript $(PLOTM_DIR)/gene.expression.per.goterm.R $(COUNTS_DIR)/LM$*_HeLa_all_timep.txt $(PLOT_GO_DIR)/ LM$* $(GO_FILES)
	@ gzip $(COUNTS_DIR)/LM$*_HeLa_all_timep.txt $(GO_DIR)/*

# Create a directory (use DIR/.d)
%/.d:
	@ mkdir -p $(@D)
	@ touch $@


# Delete everything
clean:
	rm -r -f $(TRIM_DIR)
	rm -r -f $(SORTED_DIR)
	rm -r -f $(PAIRED_DIR)
	rm -r -f $(COUNTS_DIR)
	rm -r -f $(DESEQ_DIR)
	rm -r -f $(GO_DIR)
	rm -r -f $(PLOT_DIR)

