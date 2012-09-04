from sys import argv
import gzip

script, fast_file, delete =argv

fastq_file=fast_file + ".fq.gz"
delete=int(delete)

inf=gzip.open(fastq_file, mode='rb')
out_file = gzip.open(fast_file + "_trimmed.gz", mode='wb',compresslevel=9)
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

