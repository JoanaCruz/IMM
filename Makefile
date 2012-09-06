# Directories
READS_DIR =


# all: e_coli_1000_1_trimmed.gz
all: trim

trim: $(addsuffix _trimmed.gz, $(basename $(basename $(wildcard *.fq.gz))))

%_trimmed.gz : %.fq.gz
	python trim.py ./ $* 15 ./  
	echo "Bases trimmed."

# delete everything
clean:
	rm -f *_trimmed.gz 


