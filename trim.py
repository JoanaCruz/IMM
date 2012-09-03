from sys import argv

script, fastq_file, delete , to_file=argv

delete=int(delete)

with open(fastq_file, "r") as inf:
	out_file = open(to_file, 'a')
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