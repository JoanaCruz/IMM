from sys import argv
import gzip

script, in_directory, input_file, out_directory =argv

inf = gzip.open(in_directory + input_file + ".gtf.gz", mode='rb')
out_file = open(out_directory + input_file + "_chr.gtf", 'a') 

for (i,line) in enumerate(inf):
	string=line.split()
	ins = string[13] 
	string[13] = ins[0]+"chr"+ins[1:]
	string=" ".join(string)
	print string
	out_file.write(string +"\n")

inf.close()
out_file.close()
