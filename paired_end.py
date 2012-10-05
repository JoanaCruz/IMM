from sys import argv
import gzip

script, in_directory, sam_file, out_directory, =argv

inf=gzip.open(in_directory + sam_file + "_sorted.sam.gz", mode='rb')
out_file = gzip.open(out_directory + sam_file + "_paired.sam.gz", mode='wb',compresslevel=9)

j=0
y=0
o=0
for (i,line) in enumerate(inf):
	if line[0] == "@":
		out_file.write(line)
		j=j+1
	else:
		o=o+1
		if (i+j)%2 == 0:
			l = line.split()[5]
			name = line
		else:
			print l, line.split()[5]
			if l == line.split()[5]:
				out_file.write(name)
				out_file.write(line)
				y=y+2
			
print o
print y
inf.close()
out_file.close()
