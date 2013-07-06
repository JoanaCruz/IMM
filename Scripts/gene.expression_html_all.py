import subprocess
import linecache
import csv

R_dir='/home/joana/results/Scripts/'
count_tab='/home/joana/results/control.vs_9/counts/LM1_HeLa_all_timep.txt'
#subprocess.call(["Rscript", R_dir+"normalize_count.tab.R", count_tab])
count_table = count_tab + '_normalized'

## Write index html
# Get gene name, description, acession from Emsembl site
def get_emsembl_data(ensembl_id):
	print ensembl_id
	source = "http://www.ensembl.org/id/"+ensembl_id
	subprocess.call(["wget", source])
	gene_name = linecache.getline(ensembl_id, 86).split(" ")[12] #line 86
	gene_description = linecache.getline(ensembl_id, 89).split("<p>")[1].split("</p>")[0] #89
	subprocess.call(["rm","-f",ensembl_id])
	return source, gene_name, gene_description
	
def write_plot_html(file_name,colnum, gene_name, ensembl_id, gene_description, web_page):
	f1=open(file_name, 'w+')
	write.f1('<html><head><title>Gene Expression profile</title></head><body>')
	write.f1('<font size="5">Gene Expression profile</font> <br><br>')
	header=["Gene name", "Ensembl gene ID", "Description", "Emsembl acession", "Gene expression profile"]
	
	## Start table (write header)
	f1.write('<table rules="all" CELLSPACING="0" COLS="5" BORDER="1">\n' )
	for head in header:
		f1.write('<th>' + head + '</th>\n' ) # write header columns
	f1.write('</tr>\n') # write </tr> tag
	column = [gene_name, ensembl_id, gene_description, web_page]
	for i in column:
		if i=="See plot":
			f1.write('<td class="column_' + str(colnum) + '">' + '<div align="center"><p><a href="html_plots/column'+ str(rownum) +'.html">See plot</a> </p></div></td>\n')
		elif i==web_page:
			f1.write('<td class="column_' + str(colnum) + '"><div align="center"><p><a href=' + str(web_page) +'>GO to Ensembl site</a></div></p></td>\n')
		elif i==gene_description:
			f1.write('<td class="column_' + str(colnum) + '">' + str(i) + '</td>\n')		  
		else:
			f1.write('<td class="column_' + str(colnum) + '"><div align="center">' + str(i) + '</div></td>\n')

	
	
	
def write_index(R_dir, count_table):
	f1=open('index.html', 'w+')
	f1.write('<html><head><title>Gene Expression profile</title></head><body>\n')
	f1.write('<font size="5">Gene Expression profile</font> <br><br>')
	header=["Gene name", "Ensembl gene ID", "Description", "Emsembl acession", "Gene expression profile"]
	
	## Start table (write header)
	f1.write('<table rules="all" CELLSPACING="0" COLS="5" BORDER="1">\n' )
	for head in header:
		f1.write('<th>' + head + '</th>\n' ) # write header columns
	f1.write('</tr>\n') # write </tr> tag
	
	## Write gene entries 
	# Read csv file with normalized counts
	reader = open(count_table,'r')
	# initialize rownum variable
	rownum = 0
	# write <table> tag
	firstline = True
	for row in reader: # Read a single row count table
		if firstline:    #skip first line
			firstline = False
		else:
			ensembl_id=row[1:16]
			print ensembl_id
			web_page, gene_name, gene_description = get_emsembl_data(ensembl_id)  
			column = [gene_name, ensembl_id, gene_description, web_page, 'See plot']
			colnum = 0
			for i in column:
				if i=="See plot":
					f1.write('<td class="column_' + str(colnum) + '">' + '<div align="center"><p><a href="html_plots/column'+ str(rownum) +'.html">See plot</a> </p></div></td>\n')
					#write_plot_html('column_'+str(colnum)+'.html', colnum, gene_name, ensembl_id, gene_description, web_page)
				elif i==web_page:
					f1.write('<td class="column_' + str(colnum) + '"><div align="center"><p><a href=' + str(web_page) +'>GO to Ensembl site</a></div></p></td>\n')
				elif i==gene_description:
					f1.write('<td class="column_' + str(colnum) + '">' + str(i) + '</td>\n')		  
				else:
					f1.write('<td class="column_' + str(colnum) + '"><div align="center">' + str(i) + '</div></td>\n')
				colnum += 1
			f1.write('</tr>\n')
			rownum += 1
		# write </table> tag
	f1.write('</table>')
	#subprocess.call(["rm","-f",count_table])
	
write_index(R_dir, count_table)
## Write html with plots for each gene
# Write description data
# Extract data from DESeq tables and write
# Do plots and write table with plots
 