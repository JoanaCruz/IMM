# Alignment reference
AREF=Genome


# Directories
MAKE_DIR=.
READS_DIR =/home/luispedro/youtaro-data
TRIM_DIR=trim
SAM_DIR=sam_aligned
SORTED_DIR=sam_sorted
BOWTIE_DIR=bowtie2-2.0.0-beta7

all: sorted

sorted: $(addprefix $(SORTED_DIR)/,$(addsuffix _sorted.sam.gz, $(notdir $(basename $(basename $(wildcard $(READS_DIR)/*.fq.gz))))))
	
sam: $(addprefix $(SAM_DIR)/,$(addsuffix .sam.gz, $(notdir $(basename $(basename $(wildcard $(READS_DIR)/*.fq.gz))))))

trim: $(addprefix $(TRIM_DIR)/,$(addsuffix _trimmed.gz,$(notdir $(basename $(basename $(wildcard $(READS_DIR)/*.fq.gz))))))



# make bases trim (delete first 15 bases from the beginning of the sequences)
$(TRIM_DIR)/%_trimmed.gz : $(READS_DIR)/%.fq.gz | $(TRIM_DIR)/.d
	@ python trim.py $(READS_DIR)/ $* 15 $(MAKE_DIR)/$(TRIM_DIR)/ 
	@ echo "Bases trimmed."

# make sequences alignment with the human genome (using Bowtie2 with 0 number of mismatches allowed and with all paired ends that aligned written at %_all.sam)
$(SAM_DIR)/%.sam.gz :  $(TRIM_DIR)/%_trimmed.gz | $(SAM_DIR)/.d
	@ gunzip $(AREF)/*.bt2.gzip
	@ gunzip $(AREF)/*.sh.gzip
	$(BOWTIE_DIR)/bowtie2 -x $(AREF)/hg19 --sensitive -1 $(TRIM_DIR)/*1_trimmed.gz -2 $(TRIM_DIR)/*2_trimmed.gz -S $(MAKE_DIR)/$(SAM_DIR)/$*.sam
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

# create a directory (use DIR/.d)
%/.d:
	@ mkdir -p $(@D)
	@ touch $@

# delete everything
clean:
	rm -f *_trimmed.gz 
	rm -f *.sam
	rm -f *_sorted.sam
	rm -f *_sorted.bam
	rm -r -f $(SAM_DIR)
	rm -r -f $(TRIM_DIR)
	rm -r -f $(SORTED_DIR)


