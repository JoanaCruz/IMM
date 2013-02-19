# Alignment reference
AREF=../Genome

# Main directories
HUMAN_DIR=../Human_DNA
LIST_DIR=../Listeria_RNA
MIT_DIR=../Mitocondrial_DNA
CONTROL_DIR=.
WILDTYPE_DIR=$(CONTROL_DIR)/control.vs.wildtype
MUTANT_DIR=$(CONTROL_DIR)/control.vs.mutant
READS_DIR =/home/data/reads
BOWTIE_DIR=/home/data/joanacruz/bowtie2-2.0.0-beta7
SCRIPT_DIR=../Scripts
PLOTM_DIR=../Scripts/Graphics

## Auxiliary directories
# Results
TRIM_DIR=$(HUMAN_DIR)/trim
SORTED_DIR=$(HUMAN_DIR)/alignsort
PAIRED_DIR=$(HUMAN_DIR)/paired_reads
COUNTS_DIR=$(HUMAN_DIR)/counts
COUNTS_CONTROL_DIR=$(CONTROL_DIR)/counts
DESEQ_DIR=$(CONTROL_DIR)/DESeq
GO_DIR=$(CONTROL_DIR)/GOStats
WEB_DIR=$(CONTROL_DIR)/web_site


# Graphics
PLOT_DIR=$(CONTROL_DIR)/Plots
PLOT_DESEQ=$(PLOT_DIR)/DESeq.analysis
PLOT_GO=$(PLOT_DIR)/Genes.expression

GTF_FILE := $(AREF)/hg19_refseq.gtf

all: gostats

gostats: $(addprefix $(GO_DIR)/GOsummary_,$(subst DESeq_,, $(notdir $(subst _ALL,, $(wildcard $(DESEQ_DIR)/*_ALL.csv.gz)))))

deseq: $(addprefix $(DESEQ_DIR)/DESeq_, $(notdir $(subst L1_1.fq.gz,pvalue.csv.gz, $(wildcard $(READS_DIR)/*_1.fq.gz))))

concatenate: $(addprefix $(COUNTS_CONTROL_DIR)/,$(notdir $(subst _L1_1.fq.gz,_control.txt.gz, $(wildcard $(READS_DIR)/*_1.fq.gz))))

counts: $(addprefix $(COUNTS_DIR)/, $(notdir $(subst _1.fq.gz,.txt.gz, $(wildcard $(READS_DIR)/*_1.fq.gz))))

paired: $(addprefix $(PAIRED_DIR)/,$(addsuffix _paired.sam.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

sorted: $(addprefix $(SORTED_DIR)/,$(addsuffix _sorted.sam.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

trim: $(addprefix $(TRIM_DIR)/,$(addsuffix _trimmed.gz,$(notdir $(basename $(basename $(wildcard $(READS_DIR)/*.fq.gz))))))


#Make bases trim (delete first 15 bases from the beginning of the sequences)
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


#Uses HTSeq-count script to count how many reads map to each feature (being a feature a range of positions on a chromosome)
$(COUNTS_DIR)/%.txt.gz : $(PAIRED_DIR)/%_paired.sam.gz $(GTF_FILE) | $(COUNTS_DIR)/.d
	echo "htseq"
	@ gunzip $(PAIRED_DIR)/*.sam.gz 
	@ python -m HTSeq.scripts.count -m union -s no -i gene_id $(PAIRED_DIR)/$*_paired.sam $(AREF)/hg19_refseq.gtf | gzip > $@
	@ gzip $(PAIRED_DIR)/*.sam


# Concatenate files from HTSeq-count (control vs (LM1 or LM2) )
$(COUNTS_CONTROL_DIR)/%_control.txt.gz : $(COUNTS_DIR)/%_L1.txt.gz $(COUNTS_CONTROL_DIR)/.d
	@ gunzip $(COUNTS_DIR)/$*_L1.txt.gz $(COUNTS_DIR)/Control_HeLa_L1.txt.gz
	@ paste $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_DIR)/$*_L1.txt | cut -f1,2,4- | sed '1 i\Control_HeLa\tLM1_HeLa_'$*''> $(COUNTS_CONTROL_DIR)/$*_control.txt
	@ gzip $(COUNTS_CONTROL_DIR)/$*_control.txt $(COUNTS_DIR)/$*_L1.txt $(COUNTS_DIR)/Control_HeLa_L1.txt


# Do differencial expression analysis between control and wildtype & between control and mutant (using DESeq)
$(DESEQ_DIR)/DESeq_%_pvalue.csv.gz : $(COUNTS_CONTROL_DIR)/%_control.txt.gz $(DESEQ_DIR)/.d
	@ gunzip $(COUNTS_CONTROL_DIR)/$*_control.txt.gz
	@ Rscript $(SCRIPT_DIR)/DEseq.R $(COUNTS_CONTROL_DIR)/$*_control.txt $* $(DESEQ_DIR)/
	@ gzip $(DESEQ_DIR)/DESeq_*.csv $(COUNTS_CONTROL_DIR)/$*_control.txt


# GO annotation (using GOstats)
$(GO_DIR)/GOsummary_%.csv.gz: $(DESEQ_DIR)/DESeq_%_pvalue.csv.gz $(GO_DIR)/.d
	@ gunzip $(DESEQ_DIR)/DESeq_$*_pvalue.csv.gz 
	@ Rscript $(SCRIPT_DIR)/GOStats_control.R $(DESEQ_DIR)/DESeq_$*_pvalue.csv $(GO_DIR)/
	@ gzip $(GO_DIR)/* $(DESEQ_DIR)/DESeq_$*_pvalue.csv

# Do web site to see gene expression
site : $(addprefix $(WEB_DIR)/page_,$(subst GOsummary_,, $(notdir $(subst .csv.gz,.html, $(wildcard $(GO_DIR)/GOsummary_*.csv.gz)))))
	@ echo  $(wildcard $(GO_DIR)/GOsummary_*.csv.gz)
	
$(WEB_DIR)/page_%.html: $(GO_DIR)/GOsummary_%.csv.gz $(WEB_DIR)/.d
	@ gunzip $(GO_DIR)/GOsummary_$*.csv.gz
	@ python $(SCRIPT_DIR)/html.py $(GO_DIR)/GOsummary_$*.csv $(WEB_DIR)/
	@ python $(SCRIPT_DIR)/web.py $(GO_DIR)/GOsummary_$*.csv $(WEB_DIR)/ $(PLOT_GO)/
	@ gzip $(GO_DIR)/GOsummary_$*.csv
	

plotdeseq: $(addprefix $(PLOT_DESEQ)/,$(subst DESeq_,, $(notdir $(subst 20_ALL.csv.gz,deseq.png, $(wildcard $(DESEQ_DIR)/DESeq_LM1_HeLa_20_ALL.csv.gz)))))


$(PLOT_DESEQ)/LM1_HeLa_deseq.png : $(DESEQ_WILDTYPE_DIR)/DESeq_LM1_HeLa_20_ALL.csv.gz $(DESEQ_MUTANT_DIR)/DESeq_LM2_HeLa_20_ALL.csv.gz $(PLOT_DESEQ)/.d
	@ gunzip $(DESEQ_DIR)/DESeq_LM1_HeLa_*_ALL.csv.gz $(DESEQ_DIR)/DESeq_LM2_HeLa_*_ALL.csv.gz
	@ Rscript $(PLOTM_DIR)/deseq.analysis.R $(DESEQ_WILDTYPE_DIR)/ LM1 $(PLOT_DESEQ)/
	@ Rscript $(PLOTM_DIR)/deseq.analysis.R $(DESEQ_MUTANT_DIR)/ LM2 $(PLOT_DESEQ)/
	@ gzip $(DESEQ_WILDTYPE_DIR)/DESeq_LM1_HeLa_*_ALL.csv $(DESEQ_MUTANT_DIR)/DESeq_LM2_HeLa_*_ALL.csv


plot: $(addprefix $(COUNTS_CONTROL_DIR)/, $(subst DESeq_,, $(notdir $(subst 20_ALL.csv.gz,all_timep.txt.gz, $(wildcard $(DESEQ_DIR)/DESeq_LM*_HeLa_20_ALL.csv.gz)))))
	@ echo  $(addprefix $(COUNTS_DIR)/, $(subst DESeq_,, $(notdir $(subst 20_ALL.csv.gz,all_timep.txt.gz, $(wildcard $(DESEQ_DIR)/DESeq_LM*_HeLa_20_ALL.csv.gz)))))

GO_FILES := $(basename $(wildcard $(GO_DIR)/*.csv.gz))

# Concatenate files with the same condition (LM1 or LM2), do heatmap with the most differential expressed genes and do plots of gene expression variation for each one of the goterms.
$(COUNTS_CONTROL_DIR)/LM%_HeLa_all_timep.txt.gz : $(COUNTS_DIR)/LM%_HeLa_*_L1.txt.gz $(COUNTS_DIR)/Control_HeLa_L1.txt.gz | $(PLOT_GO)/.d
	@ gunzip $(COUNTS_DIR)/Control_HeLa_L1.txt.gz $(COUNTS_DIR)/LM$*_HeLa_*_L1.txt.gz
	@ paste $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_DIR)/LM$*_HeLa_20_L1.txt $(COUNTS_DIR)/LM$*_HeLa_60_L1.txt $(COUNTS_DIR)/LM$*_HeLa_120_L1.txt $(COUNTS_DIR)/LM$*_HeLa_240_L1.txt | cut -f1,2,4,6,8,10- | sed '1 i\Control_HeLa\tLM'$*'_HeLa_20\tLM'$*'_HeLa_60\tLM'$*'_HeLa_120\tLM'$*'_HeLa_240' > $(COUNTS_CONTROL_DIR)/LM$*_HeLa_all_timep.txt
	@ gzip $(COUNTS_DIR)/Control_HeLa_L1.txt $(COUNTS_DIR)/LM$*_HeLa_*_L1.txt
# Gene expression graphics
	@ gunzip $(GO_DIR)/*.gz
	@ Rscript $(PLOTM_DIR)/heatmap.R  $(COUNTS_DIR)/LM$*_HeLa_all_timep.txt $(PLOT_GO_DIR)/ LM$*
	@ Rscript $(PLOTM_DIR)/gene.expression.per.goterm.R $(COUNTS_CONTROL_DIR)/LM$*_HeLa_all_timep.txt $(PLOT_GO)/ $(GO_FILES)
	@ gzip $(COUNTS_CONTROL_DIR)/LM$*_HeLa_all_timep.txt $(GO_DIR)/*


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

