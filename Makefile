# Alignment reference
AREF=indexes

# Directories
READS_DIR =
TRIM_DIR=trim
SAM_DIR=sam
BOWTIE_DIR=bowtie2-2.0.0-beta7

# all: e_coli_1000_1_trimmed.gz
all: sam

sam: $(addprefix $(SAM_DIR)/,$(addsuffix _T.sam, $(basename $(basename $(wildcard *.fq.gz)))))

#trim: $(addprefix $(TRIM_DIR)/,$(addsuffix _trimmed.gz, $(basename $(basename $(wildcard *.fq.gz)))))


# make bases trim (delete first 15 bases from the beginning of the sequences)
$(TRIM_DIR)/%_trimmed.gz : %.fq.gz
	python trim.py ./ $(notdir $*) 15 ./$(TRIM_DIR)/
	echo $(notdir $*)  
	echo "Bases trimmed."

# make sequences alignment with the human genome (using Bowtie2 with 0 number of mismatches allowed and with all paired ends that aligned written at %_all.sam)
$(SAM_DIR)/%_T.sam :  $(TRIM_DIR)/%_trimmed.gz
	$(BOWTIE_DIR)/bowtie2 -x $(AREF)/e_coli -N 0 -1 $(TRIM_DIR)/*1_trimmed.gz -2 $(TRIM_DIR)/*2_trimmed.gz -S $*_all.sam
	echo "Alignment done."

#bowtie2-2.0.0-beta7/./bowtie2 -x indexes/e_coli -N 0 --al-conc-gz ecoli_all.sam -1 e_coli_1000_1.fq.gz -2 e_coli_1000_2.fq.gz -S coli_all.sam


# delete everything
clean:
	rm -f *_trimmed.gz 
	rm -f *.sam
	rm -r -f $(SAMS_DIR)
	rm -r -f $(TRIM_DIR)


