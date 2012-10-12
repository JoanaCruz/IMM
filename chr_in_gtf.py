from sys import argv
import gzip

script, in_directory, input_file, out_directory =argv

inf = gzip.open(in_directory + input_file + ".gtf.gz", mode='rb')
out_file = gzip.open(out_directory + input_file + "_chr.gtf.gz", mode='wb',compresslevel=9) 

for (i,line) in enumerate(inf):
	string=line.split()
	ins = string[13] 
	string[13] = ins[0]+"chr"+ins[1:]
	
inf.close()
out_file.close()
