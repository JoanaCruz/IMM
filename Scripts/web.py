from sys import argv
from bs4 import BeautifulSoup
import re
import linecache
script, file_read, out_dir,tp, ont, plot_dir =argv

#file_read: Table from excell with the information of one time point, one ontology (in html format)
#out_dir: Directory where html with plots will be saved.
#plot_dir: Directory where are the gene expression plots.

#web.py ./ pageBP_20_semplot.html ./ 20 BP ./
#python web.py pageMF_60.html 60_timep/ 60 MF ../../genes/

f = open( file_read , 'r')
new_file= open ( "New_"+file_read, 'w')

# 0. get goterm
# 1. get genes
# 2. hiperlink to the plot page
# 3. write plots in plot html page

length_f=0
for length_f,l in enumerate(f):
	length_f+=1
f.close()

f = open( file_read , 'r')

i=0
j=0
var="FALSE"

for i, line in enumerate(f):
	if i==(length_f-6):
		break
	if i==21:
		if ont=="BP":
			ontology="Biological processes"
		elif ont=="CC":
			ontology="Celular components"
		elif ont=="MF":
			ontology="Molecular function"
		#go to previous page
		head="""<font size="5">Gene Expression profile - 7th February 2013</font> <br><br><p><font size="4">GO ontology: """+ontology+" -> Timepoint:"+tp+""" </font></p>
<p align="right"> <a href="page"""+ont+"""_timepoints.html"><font size="3"> Previous page</font></a></p><br><br>"""
		new_file.write(head)
	if i==27:
		#define width of gene expression column
		width="""<COLGROUP WIDTH="100"></COLGROUP>"""
		new_file.write(width)
	if i==33:
		# add column to see plots and respective link
		column="""<TD STYLE="border-top: 1px solid #000000; border-bottom: 1px solid #000000; border-left: 1px solid #000000; border-right: 1px solid #000000" ALIGN="CENTER" VALIGN=MIDDLE><B>Genes expression profile</B></TD>"""
		new_file.write(column)
	if i == 37+j:
		#get GO term
		goterm_line=linecache.getline(file_read, 37+j)
		r = re.compile('>(.*?)<')
		m = r.search(goterm_line)
		goterm_name=m.group(1)
		goterm_file=goterm_name.replace(":","-")
	elif i == 40+j:
		#get genes
		genes_line=linecache.getline(file_read, 40+j)
		r = re.compile('>(.*?)<')
		m = r.search(genes_line)
		genes_group = m.group(1)
		# Write "see plot" and respective hyperlink
		html='''<TD STYLE="border-top: 1px solid #000000; border-bottom: 1px solid #000000; border-left: 1px solid #000000; border-right: 1px solid #000000" ALIGN="CENTER" VALIGN=MIDDLE><A HREF="''' +out_dir+ "page" +ont+"_"+tp+"_"+goterm_file+""".html">See plot</A></TD>"""
		new_file.write(html)
		j+=7
		var="TRUE"
	new_file.write(line)
		
	if var == "TRUE":
		genes=genes_group.split(", ")
		out_file=open(out_dir+"page"+ont+"_"+tp+"_"+goterm_file+".html", 'w+')
		out_file.write("""<html><head><title>Gene Expression profile</title></head><body><p><font size="5">Gene Expression profile - 7th February 2013</font> <br></p>
	<p><font size="4">GO ontology: Biological Processes -> Timepoint:20 -> """+goterm_name+ """</font></p><p align="right"> <a href="../New_page""" +ont+"_"+tp+ """.html"> Previous page </a> </p>
	<table><tr>""")
		for i in range(len(genes)):
			out_file.write("""<td><img src="""+plot_dir+genes[i]+".png"+ """></a></td></tr></table></body></html>""")
		out_file.close()
	

new_file.close()
f.close()

