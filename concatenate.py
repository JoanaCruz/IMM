from sys import argv

script, f1, f2 = argv

file_1 = open(f1, "r")
file_2 = open(f2, "r")
    
for line in file_2:
	print file_1.readline() + " " + line.split()[1]


