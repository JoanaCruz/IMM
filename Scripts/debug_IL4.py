### DEBUG
from sys import argv
import subprocess
import fileinput

#import HTSeq
# python Scripts/debug_IL4.py /mnt/hgfs/reads/debug/IL4.fasta /mnt/hgfs/Genome/Homo_sapiens.GRCh37.71 /mnt/hgfs/Genome/Homo_sapiens.GRCh37.71.gtf ~/Downloads/bowtie2-2.0.0-beta7/ ./

script, gene_file, genome_ref, anotation_gtf, bowtie_ref, out_dir =argv


BOWTIE_DIR = bowtie_ref
AREF1 = genome_ref
AREF2 = anotation_gtf

## Create ficticious fq file
#read .fq file
processed_fasta=[]
for line in fileinput.input([gene_file]):
	if fileinput.isfirstline():
		header=line.replace('\r\n','')
	else:
		processed_fasta+=str(line).replace('\r\n','')
		


with open('./IL4-cutted.fq','w+') as new_fasta_file:
	#write header
	new_fasta_file.write(header+'\n')
	# Cut fq file in portions with s size
	trim=10
	beggin=0
	for i in range(trim,len(processed_fasta)):
		new_fasta_file.write(''.join(processed_fasta[beggin:i])+'\n')
		beggin+=1

		
## Align with human genome
subprocess.call([BOWTIE_DIR+'/bowtie2', '-f', '-r', '--sensitive', '-x', AREF1, '-U', './IL4-cutted.fq', '-S', out_dir+'IL4-cutted.sam'])

## Map reads
subprocess.call(['python', '-m', 'HTSeq.scripts.count', '-m', 'union', '-s', 'no', '-i', 'gene_id', 'IL4-cutted.sam', AREF2, '>', out_dir+ 'IL4.txt'])
