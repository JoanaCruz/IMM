# Alignment reference
AREF=../Genome

# Main directories
HUMAN_DIR=Human_DNA
LIST_DIR=Listeria_RNA
MIT_DIR=Mitocondrial_DNA
READS_DIR =/data/reads
BOWTIE_DIR=../../bowtie2-2.0.0-beta7

#Auxiliary directories
TRIM_DIR=$(HUMAN_DIR)/trim
SORTED_DIR=$(LIST_DIR)/align&sorted_list
PAIRED_DIR=$(LIST_DIR)/paired_reads_list
COUNTS_DIR=$(LIST_DIR)/counts_list
DESEQ_DIR=$(LIST_DIR)/DESeq_list
GO_DIR=$(LIST_DIR)/GOStats_list

GTF_FILE := $(AREF)/listeria_monocytogenes_j0161_8_transcripts.gtf

all: gostats

gostats: $(addprefix $(GO_DIR)/,$(subst LM1_HeLa, GOsummary, $(notdir $(subst _L1_1.fq.gz,.txt.gz, $(wildcard $(READS_DIR)/LM1_*_1.fq.gz)))))

deseq: $(addprefix $(DESEQ_DIR)/,$(subst LM1, DESeq_LM, $(notdir $(subst _L1_1.fq.gz,.csv.gz, $(wildcard $(READS_DIR)/LM1_*_1.fq.gz)))))

concatenate: $(addprefix $(COUNTS_DIR)/,$(subst LM1, LM, $(notdir $(subst L1_1.fq.gz,HeLa.txt.gz, $(wildcard $(READS_DIR)/LM1_*_1.fq.gz)))))

counts: $(addprefix $(COUNTS_DIR)/,$(addsuffix .txt.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

paired: $(addprefix $(PAIRED_DIR)/,$(addsuffix _paired.sam.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

sorted: $(addprefix $(SORTED_DIR)/,$(addsuffix _sorted.sam.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

trim: $(addprefix $(TRIM_DIR)/,$(addsuffix _trimmed.gz,$(notdir $(basename $(basename $(wildcard $(READS_DIR)/*.fq.gz))))))


# Make bases trim (delete first 15 bases from the beginning of the sequences)
$(TRIM_DIR)/%_trimmed.gz : $(READS_DIR)/%.fq.gz | $(TRIM_DIR)/.d
	@ python $(SCRIPT_DIR)/trim.py $(READS_DIR)/ $* 15 $(TRIM_DIR)/ 
	@ echo "Bases trimmed."


# Make sequences alignment with the human genome (using Bowtie2 with --sensitive option and with all paired ends that aligned written at .sam)
$(SORTED_DIR)/%.sam.gz :  $(TRIM_DIR)/%_1_trimmed.gz $(TRIM_DIR)/%_2_trimmed.gz | $(SORTED_DIR)/.d
	$(BOWTIE_DIR)/bowtie2 -x $(AREF)/listeria_transcripts_j0161 --sensitive -1 $(TRIM_DIR)/$*_1_trimmed.gz -2 $(TRIM_DIR)/$*_2_trimmed.gz | gzip > $(SORTED_DIR)/$*.sam.gz
	@ echo "Alignment done."

### Arrange SAM files.

## Sort the SAM file by read name (to use in HTSeq) unzip
# Convert .sam to .bam and sort
$(SORTED_DIR)/%_sorted.bam : $(SORTED_DIR)/%.sam.gz
	@ samtools view -bS $(SORTED_DIR)/$*.sam.gz | samtools sort -n - $(SORTED_DIR)/$*_sorted
# Convert .bam sorted to .sam
$(SORTED_DIR)/%_sorted.sam.gz : $(SORTED_DIR)/%_sorted.bam
	@ samtools view $(SORTED_DIR)/$*_sorted.bam -h > $(SORTED_DIR)/$*_sorted.sam
	@ gzip $(SORTED_DIR)/$*_sorted.sam
	@ echo "SAM file sorted."


# Filter paired end reads which only one align
$(PAIRED_DIR)/%_paired.sam.gz : $(SORTED_DIR)/%_sorted.sam.gz | $(PAIRED_DIR)/.d
	@ python $(SCRIPT_DIR)/paired_end.py $(SORTED_DIR)/ $*$(PAIRED_DIR)/
	@ echo "Reads filtered."


# Uses HTSeq-count script to count how many reads map to each feature (being a feature a range of positions on a chromosome)
$(COUNTS_DIR)/%.txt.gz : $(PAIRED_DIR)/%_paired.sam.gz $(GTF_FILE) | $(COUNTS_DIR)/.d
	@ gunzip $(PAIRED_DIR)/*.sam.gz 
	@ python -m HTSeq.scripts.count -m union -s no -i gene_id $(PAIRED_DIR)/$*_paired.sam $(AREF)/listeria_monocytogenes_j0161_8_transcripts.gtf | gzip > $@
	@ gzip $(PAIRED_DIR)/*.sam


# Concatenate files from HTSeq-count (same timepoint diferent condition)
$(COUNTS_DIR)/LM_HeLa_%.txt.gz : $(COUNTS_DIR)/LM1_HeLa_%_L1.txt.gz $(COUNTS_DIR)/LM2_HeLa_%_L1.txt.gz 
	@ gunzip $(COUNTS_DIR)/LM*_HeLa_$*_L1.txt.gz
	@ join -t "`/bin/echo -e '\t'`" --head $(COUNTS_DIR)/LM1_HeLa_$*_L1.txt $(COUNTS_DIR)/LM2_HeLa_$*_L1.txt > $(COUNTS_DIR)/LM_HeLa_$*.txt
	@ sed -i '1s/^/LM1_HeLa_'$*'\tLM2_HeLa_'$*'\n/' $(COUNTS_DIR)/LM_HeLa_$*.txt
	@ gzip $(COUNTS_DIR)/*.txt


# Do differencial expression analysis (using DESeq)
$(DESEQ_DIR)/DESeq_LM_HeLa_%.csv.gz: $(COUNTS_DIR)/LM_HeLa_%.txt.gz | $(DESEQ_DIR)/.d
	@ gunzip $(COUNTS_DIR)/LM_HeLa_$*.txt.gz
	@ Rscript $(SCRIPT_DIR)/DEseq.R $(COUNTS_DIR)/LM_HeLa_$*.txt $* $(DESEQ_DIR)/
	@ gzip $(DESEQ_DIR)/*.csv


# GO annotation (using GOstats)
$(GO_DIR)/GOsummary_%.csv.gz: $(DESEQ_DIR)/DESeq_LM_HeLa_%_pvalue.csv.gz | $(GO_DIR)/.d
	@ gunzip $(DESEQ_DIR)/DESeq_LM_HeLa_$*_pvalue.csv.gz
	@ Rscript $(SCRIPT_DIR)/GOStats_def.R $(DESEQ_DIR)/DESeq_LM_HeLa_$*_pvalue.csv $* $(GO_DIR)/
	@ gzip $(GO_DIR)/*.csv
	@ gzip $(DESEQ_DIR)/DESeq_LM_HeLa_$*_pvalue.csv


# Create a directory (use DIR/.d)
%/.d:
	@ mkdir -p $(@D)
	@ touch $@



# Delete everything
clean:
	rm -r -f $(SORTED_DIR)
	rm -r -f $(PAIRED_DIR)
	rm -r -f $(COUNTS_DIR)
	rm -r -f $(DESEQ_DIR)
	rm -r -f $(GO_DIR)

