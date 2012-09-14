# Alignment reference
AREF=Genome

# Directories
MAKE_DIR=.
READS_DIR =/home/luispedro/youtaro-data
TRIM_DIR=trim
SAM_DIR=sam_aligned
SORTED_DIR=sam_sorted
BOWTIE_DIR=bowtie2-2.0.0-beta7
COUTS_DIR=counted_reads

all: counts

counts: $(addprefix $(COUNTS_DIR)/,$(addsuffix .txt.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

sorted: $(addprefix $(SORTED_DIR)/,$(addsuffix _sorted.sam.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

sam: $(addprefix $(SAM_DIR)/,$(addsuffix .sam.gz, $(notdir $(subst _1.fq.gz,, $(wildcard $(READS_DIR)/*_1.fq.gz)))))

trim: $(addprefix $(TRIM_DIR)/,$(addsuffix _trimmed.gz,$(notdir $(basename $(basename $(wildcard $(READS_DIR)/*.fq.gz))))))


# Make bases trim (delete first 15 bases from the beginning of the sequences)
$(TRIM_DIR)/%_trimmed.gz : $(READS_DIR)/%.fq.gz | $(TRIM_DIR)/.d
	@ python trim.py $(READS_DIR)/ $* 15 $(MAKE_DIR)/$(TRIM_DIR)/ 
	@ echo "Bases trimmed."

# Make sequences alignment with the human genome (using Bowtie2 with --sensitive option and with all paired ends that aligned written at .sam)
$(SAM_DIR)/%.sam.gz :  $(TRIM_DIR)/%_1_trimmed.gz | $(SAM_DIR)/.d
	@ gunzip $(AREF)/*.bt2.gz
	@ gunzip $(AREF)/*.sh.gz
	$(BOWTIE_DIR)/bowtie2 -x $(AREF)/hg19 --sensitive -1 $(TRIM_DIR)/$*_1_trimmed.gz -2 $(TRIM_DIR)/$*_2_trimmed.gz -S $(MAKE_DIR)/$(SAM_DIR)/$*.sam
	@ gzip $(MAKE_DIR)/$(SAM_DIR)/$*.sam
	@ gzip $(AREF)/*.bt2
	@ gzip $(AREF)/*.sh
	@ echo "Alignment done."

## Sort the SAM file by read name (to use in HTSeq) unzip
# Convert .sam to .bam and sort
$(SORTED_DIR)/%_sorted.bam : $(SAM_DIR)/%.sam.gz | $(SORTED_DIR)/.d
	@ samtools view -bS $(SAM_DIR)/$*.sam.gz | samtools sort -n - $(SORTED_DIR)/$*_sorted

# Convert .bam sorted to .sam
$(SORTED_DIR)/%_sorted.sam.gz : $(SORTED_DIR)/%_sorted.bam
	@ samtools view $(SORTED_DIR)/$*_sorted.bam -h > $(SORTED_DIR)/$*_sorted.sam
	@ gzip $(SORTED_DIR)/$*_sorted.sam
	@ echo "SAM file sorted."

# Uses HTSeq-count script to count how many reads map to each feature (being a feature a range of positions on a chromosome)
$(COUNTS_DIR)/%.txt.gz : $(SORTED_DIR)/%_sorted.sam.gz | $(COUNTS_DIR)/.d
	@ gunzip $(SORTED_DIR)/*.sam.gz 
	@ python -m HTSeq.scripts.count -i ID $(SORTED_DIR)/$*_sorted.sam $(AREF)/*.gff
	@ gzip $(SORTED_DIR)/*.sam
	@ gzip $(COUNTS_DIR)/*.txt

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

