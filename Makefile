# Alignment reference
AREF=indexes


# Directories
READS_DIR =.
TRIM_DIR=trim
SAM_DIR=sam
BOWTIE_DIR=bowtie2-2.0.0-beta7

all: sam

sam: $(addprefix $(SAM_DIR)/,$(addsuffix _T.sam, $(basename $(basename $(wildcard $(READS_DIR)/*.fq.gz)))))

trim: $(addprefix $(TRIM_DIR)/,$(addsuffix _trimmed.gz, $(basename $(basename $(wildcard $(READS_DIR)/*.fq.gz)))))

# make bases trim (delete first 15 bases from the beginning of the sequences)
$(TRIM_DIR)/%_trimmed.gz : $(READS_DIR)/%.fq.gz | $(TRIM_DIR)/.d
	python trim.py $(READS_DIR)/ $* 15 $(READS_DIR)/$(TRIM_DIR)/ 
	@echo "Bases trimmed."

# make sequences alignment with the human genome (using Bowtie2 with 0 number of mismatches allowed and with all paired ends that aligned written at %_all.sam)
$(SAM_DIR)/%_T.sam :  $(TRIM_DIR)/%_trimmed.gz | $(SAM_DIR)/.d
	$(BOWTIE_DIR)/bowtie2 -x $(AREF)/e_coli -N 0 -1 $(TRIM_DIR)/*1_trimmed.gz -2 $(TRIM_DIR)/*2_trimmed.gz -S $*_all.sam
	@echo "Alignment done."

# create a directory (use DIR/.d)
%/.d:
	mkdir -p $(@D)
	touch $@

# delete everything
clean:
	rm -f *_trimmed.gz 
	rm -f *.sam
	rm -r -f $(SAM_DIR)
	rm -r -f $(TRIM_DIR)


