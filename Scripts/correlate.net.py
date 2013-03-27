from sys import argv
import linecache
import glob
import os


# extract mrna data related with gene from NCBI file
def extract_mrna(input_file_data):
	mrna_data=[]
	source=[]
	gene_name=[]
	with open(input_file_data) as read_data_file:
		for (i, line) in enumerate(read_data_file):
			if i==0:
				source.append(line)
			if line[0:10]=="<pre><pre>":
				gene_name.append(line[10:-1])
			if line[0:4]=="\tNM_":
				mrna_data.append(line.split()[0].split(".")[0])
	return mrna_data, source, gene_name


# write html page with gene expression profile (taking into account mrna data extracted from NCBI file)
def write_plot_html(gene_name, mrna_data, source, plot_dir, out_directory, rownum):
	#create folder for plot html
	if not os.path.exists(out_directory + "html_plots/"):
		os.makedirs(out_directory + "html_plots/")
	with open(out_directory + "html_plots/column" + str(rownum) + ".html", "w+") as html_plot:
		html_plot.write('<html><head><title>NETWORK ANALYSIS</title></head><body><font size="5">Network of immune response - 27th March 2013</font>\n')
		html_plot.write('<br><br><font size="4">Gene expression profile</font>')
		
		html_plot.write('<p align="right"> <a href="../network_index.html"> Previous page </a> </p>')	
		
		html_plot.write('<table rules="all" CELLSPACING="0" COLS="5" BORDER="1"><tr>' + "\n" )
		
		header=["Gene name", "RefSeq", "NCBI acession"]
		for head in header:
			html_plot.write('<th>' + head + '</th>' + "\n" ) # write header columns
		html_plot.write('</tr>'+ "\n") # write </tr> tag
		
		html_plot.write('<tr class="color1">' + "\n")
		column=[gene_name, mrna_data, source]
		colnum=1
		for i in column:
			if i==source:
				html_plot.write('<td class="column_' + str(colnum) + '">' + '<p><a href=' + str(source[0]) +'>GO to NCBI site</a> </p></td>\n')
			elif len(mrna_data) > 1:
				html_plot.write('<td class="column_' + str(colnum) + '">')
				for j in range(len(mrna_data)):
					html_plot.write(str(i[j]) + "  ")	
				html_plot.write('</td>\n')
			else:
				html_plot.write('<td class="column_' + str(colnum) + '">' + str(i[0]) + '</td>'+ "\n")
			
		html_plot.write('</tr></table>')

		i=0
		for i in range(len(mrna_data)):
			html_plot.write('<br><table><tr><td><img src='+plot_dir+mrna_data[i]+".png"+ '></a></td></tr></table></body></html>')


def write_html(mrna_data, source, gene_name, rownum, out_directory, plot_dir):
	with open(out_directory+"network_index.html","w+") as index:
		print(mrna_data)
		if rownum == 0:
			index.write('<html><head><title>NETWORK ANALYSIS</title></head><body><font size="5">Network of immune response - 27th March 2013</font></body></html><br>\n')
			index.write('<br><table rules="all" CELLSPACING="0" COLS="5" BORDER="1">' + "\n" )
			index.write('<tr>'+ "\n") # write <tr> tag
			header=["Gene name", "RefSeq", "NCBI acession", "Gene expression profile"]
			for head in header:
				index.write('<th>' + head + '</th>' + "\n" ) # write header columns
			index.write('</tr>'+ "\n") # write </tr> tag
		
		if rownum % 2 == 0:
			index.write('<tr class="color1">' + "\n")
		else:
			index.write('<tr class="color2">' + "\n")
		column=[ gene_name, mrna_data, source, "See plot" ]
		colnum=0
		for i in column:
			if i=="See plot":
				index.write('<td class="column_' + str(colnum) + '">' + '<p><a href="html_plots/column'+ str(rownum) +'.html">See plot</a> </p>' + '</td>\n')
				write_plot_html(gene_name, mrna_data, source, plot_dir, out_directory, rownum)
			elif i==source:
				index.write('<td class="column_' + str(colnum) + '">' + '<p><a href=' + str(source[0]) +'>GO to NCBI site</a> </p></td>\n')
			elif len(mrna_data) > 1:
				index.write('<td class="column_' + str(colnum) + '">')
				for j in range(len(mrna_data)):
					index.write(str(i[j]) + "  ")	
				index.write('</td>\n')
			else:
				index.write('<td class="column_' + str(colnum) + '">' + str(i[0]) + '</td>'+ "\n")
			colnum += 1
		index.write('</tr></table>' + "\n")


def do_analysis(data_directory, plot_directory, out_directory):
	os.chdir(data_directory)
	rownum=0
	for input_file_data in glob.glob("*.txt"):
		mrna_data, source, gene_name = extract_mrna(input_file_data)
		write_html(mrna_data, source, gene_name, rownum, out_directory, plot_directory)
		rownum += 1
	

script, data_directory, plot_directory, out_directory = argv

#call function
do_analysis(argv[1], argv[2], argv[3])
 

	
	






