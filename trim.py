from sys import argv
import gzip

script, in_directory, fast_file, delete, out_directory =argv

fastq_file=in_directory + fast_file + ".fq.gz"
delete=int(delete)

inf=gzip.open(fastq_file, mode='rb')
out_file = gzip.open(out_directory + fast_file + "_trimmed.gz", mode='wb',compresslevel=9)
for (i, line) in enumerate(inf):
	if i % 2 != 0:
		out_file.write(line[delete:])
		if i in range (0,10):
			print line[delete:]
			print line
	else:
		out_file.write(line)
		if i in range (0,10):
			print line

inf.close()
out_file.close()

