# Alignment reference
AREF=../Genome

# Directories
MAKE_DIR=.
READS_DIR =/data/reads
TRIM_DIR=trim
SAM_DIR=sam_aligned
SORTED_DIR=sam_sorted
BOWTIE_DIR=../bowtie2-2.0.0-beta7
COUNTS_DIR=counted_reads
PAIRED_DIR=paired_reads
DESEQ_DIR=DESeq
GO_DIR=GOStats

GTF_FILE := $(AREF)/hg19_refseq.gtf

all: GO

GO: $(addprefix $(GO_DIR)/,$(subst LM1_HeLa, GOsummary, $(notdir $(subst _L1_1.fq.gz,.txt.gz, $(wildcard $(READS_DIR)/LM1_*_1.fq.gz)))))

R: $(addprefix $(DESEQ_DIR)/,$(subst LM1, DESeq_LM, $(notdir $(subst _L1_1.fq.gz,.txt.gz, $(wildcard $(READS_DIR)/LM1_*_1.fq.gz)))))

concatenate: $(addprefix $(COUNTS_DIR)/,$(subst LM1, LM, $(notdir $(subst L1_1.fq.gz,HeLa.txt.gz, $(wildcard $(READS_DIR)/LM1_*_1.fq.gz)))))

counts: $(addprefix $(COUNTS_DIR)/,$(addsuffix .txt.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

paired: $(addprefix $(PAIRED_DIR)/,$(addsuffix _paired.sam.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

sorted: $(addprefix $(SORTED_DIR)/,$(addsuffix _sorted.sam.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

trim: $(addprefix $(TRIM_DIR)/,$(addsuffix _trimmed.gz,$(notdir $(basename $(basename $(wildcard $(READS_DIR)/*.fq.gz))))))


# Make bases trim (delete first 15 bases from the beginning of the sequences)
$(TRIM_DIR)/%_trimmed.gz : $(READS_DIR)/%.fq.gz | $(TRIM_DIR)/.d
	@ python trim.py $(READS_DIR)/ $* 15 $(MAKE_DIR)/$(TRIM_DIR)/ 
	@ echo "Bases trimmed."


# Make sequences alignment with the human genome (using Bowtie2 with --sensitive option and with all paired ends that aligned written at .sam)
$(SORTED_DIR)/%.sam.gz :  $(TRIM_DIR)/%_1_trimmed.gz $(TRIM_DIR)/%_2_trimmed.gz | $(SORTED_DIR)/.d
	$(BOWTIE_DIR)/bowtie2 -x $(AREF)/hg19 --sensitive -1 $(TRIM_DIR)/$*_1_trimmed.gz -2 $(TRIM_DIR)/$*_2_trimmed.gz | gzip > $(MAKE_DIR)/$(SORTED_DIR)/$*.sam.gz
	@ echo "Alignment done."


### Arrange SAM files.

## Sort the SAM file by read name (to use in HTSeq) unzip
# Convert .sam to .bam and sort
$(SORTED_DIR)/%_sorted.bam : $(SAM_DIR)/%.sam.gz
	@ samtools view -bS $(SAM_DIR)/$*.sam.gz | samtools sort -n - $(SORTED_DIR)/$*_sorted
# Convert .bam sorted to .sam
$(SORTED_DIR)/%_sorted.sam.gz : $(SORTED_DIR)/%_sorted.bam
	@ samtools view $(SORTED_DIR)/$*_sorted.bam -h > $(SORTED_DIR)/$*_sorted.sam
	@ gzip $(SORTED_DIR)/$*_sorted.sam
	@ echo "SAM file sorted."


# Filter paired end reads which only one align
$(PAIRED_DIR)/%_paired.sam.gz : $(SORTED_DIR)/%_sorted.sam.gz | $(PAIRED_DIR)/.d
	@ python paired_end.py $(SORTED_DIR)/ $* $(MAKE_DIR)/$(PAIRED_DIR)/
	@ echo "Reads filtered."


# Uses HTSeq-count script to count how many reads map to each feature (being a feature a range of positions on a chromosome)
$(COUNTS_DIR)/%.txt.gz : $(PAIRED_DIR)/%_paired.sam.gz $(GTF_FILE) | $(COUNTS_DIR)/.d
	@ gunzip $(PAIRED_DIR)/*.sam.gz 
	@ python -m HTSeq.scripts.count -m union -s no -i gene_id $(PAIRED_DIR)/$*_paired.sam $(AREF)/hg19_refseq.gtf | gzip > $@
	@ gzip $(PAIRED_DIR)/*.sam


# Concatenate files from HTSeq-count (same timepoint diferetn condition)
$(COUNTS_DIR)/LM_HeLa_%.txt.gz : $(COUNTS_DIR)/LM1_HeLa_%_L1.txt.gz $(COUNTS_DIR)/LM2_HeLa_%_L1.txt.gz 
	@ gunzip $(COUNTS_DIR)/LM*_HeLa_$*_L1.txt.gz
	@ join -t "`/bin/echo -e '\t'`" --head $(COUNTS_DIR)/LM1_HeLa_$*_L1.txt $(COUNTS_DIR)/LM2_HeLa_$*_L1.txt > $(COUNTS_DIR)/LM_HeLa_$*.txt
	@ sed -i '1s/^/LM1_HeLa_'$*'\tLM2_HeLa_'$*'\n/' $(COUNTS_DIR)/LM_HeLa_$*.txt
	@ gzip $(COUNTS_DIR)/*.txt


# Do differencial expression analysis (using DESeq)
$(DESEQ_DIR)/DESeq_LM_HeLa_%.txt.gz: $(COUNTS_DIR)/LM_HeLa_%.txt.gz | $(DESEQ_DIR)/.d
	@ gunzip $(COUNTS_DIR)/LM_HeLa_$*.txt.gz
	@ Rscript $(MAKE_DIR)/DEseq.R $(COUNTS_DIR)/LM_HeLa_$*.txt $* $(DESEQ_DIR)/
	@ gzip $(DESEQ_DIR)/*.txt

# GO annotation (using GOstats)
$(GO_DIR)/GOsummary_%.txt.gz: $(DESEQ_DIR)/DESeq_LM_HeLa_%_pvalue.txt.gz | $(GO_DIR)/.d
	@ gunzip $(DESEQ_DIR)/DESeq_LM_HeLa_$*_pvalue.txt.gz
	@ Rscript $(MAKE_DIR)/GOStats_def.R $(DESEQ_DIR)/DESeq_LM_HeLa_$*_pvalue.txt $* $(GO_DIR)/
	@ gzip $(GO_DIR)/*.txt
	@ gzip $(DESEQ_DIR)/DESeq_LM_HeLa_$*_pvalue.txt

# Create a directory (use DIR/.d)
%/.d:
	@ mkdir -p $(@D)
	@ touch $@



# Delete everything
clean:
	rm -f *_trimmed.gz 
	rm -f *.sam
	rm -f *_sorted.sam
	rm -f *_sorted.bam
	rm -r -f $(SAM_DIR)
	rm -r -f $(TRIM_DIR)
	rm -r -f $(SORTED_DIR) 
	rm -r -f $(COUNTS_DIR)
	rm -r -f $(PAIRED_DIR)
	rm -r -f $(DESEQ_DIR)
	rm -r -f $(GO_DIR) 
