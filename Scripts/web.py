from sys import argv
from bs4 import BeautifulSoup
import re
import linecache
import os
script, file_read, web_dir, plot_dir =argv

#file_read: Table from excell with the information of one time point, one ontology (in html format)
#out_dir: Directory where html with plots will be saved.
#plot_dir: Directory where are the gene expression plots.

#web.py ./ pageBP_20_semplot.html ./ 20 BP ./
#python web.py pageMF_60.html 60_timep/ 60 MF ../../genes/

name_in_file=file_read.split("_")

tp=name_in_file[2]
ont=name_in_file[1]
out_dir=web_dir+ont+"/"+tp+"_timep"

f = open( web_dir +ont+"/"+ file_read + ".html", 'r')
new_file= open ( web_dir +ont+ "/New_" + file_read +".html", 'w+')

# 0. get goterm
# 1. get genes
# 2. hiperlink to the plot page
# 3. write plots in plot html page

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

length_f=0
for length_f,l in enumerate(f):
	length_f+=1
f.close()

f = open( web_dir +ont+"/"+ file_read + ".html", 'r')

i=0
j=0
var="FALSE"

for i, line in enumerate(f):
	if i==(length_f-2):
		break
	if i==0:
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
	if i==7:
		# add column to see plots
		column="""<TD STYLE="border-top: 1px solid #000000; border-bottom: 1px solid #000000; border-left: 1px solid #000000; border-right: 1px solid #000000" ALIGN="CENTER" VALIGN=MIDDLE><B>Genes expression profile</B></TD>"""
		new_file.write(column)
	if i == 10+j:
		#get GO term
		goterm_line=linecache.getline(web_dir +ont+"/"+ file_read + ".html", 11+j)
		r = re.compile('>(.*?)<')
		m = r.search(goterm_line)
		goterm_name=m.group(1)
		goterm_file=goterm_name.replace(":","-")
	elif i == 14+j:
		#get genes
		genes_line=linecache.getline(web_dir +ont+"/"+ file_read + ".html", 14+j)
		r = re.compile('>(.*?)<')
		m = r.search(genes_line)
		genes_group = m.group(1)
		# Write "see plot" and respective hyperlink
		html='''<TD STYLE="border-top: 1px solid #000000; border-bottom: 1px solid #000000; border-left: 1px solid #000000; border-right: 1px solid #000000" ALIGN="CENTER" VALIGN=MIDDLE><A HREF="''' + web_dir + tp + "_timep/page_"+ont+"_"+tp+"_"+goterm_file+""".html">See plot</A></TD>"""
		print html
		new_file.write(html)
		j+=7
		var="TRUE"
	new_file.write(line)
		
	if var == "TRUE":
		genes=genes_group.split(", ")
		#write html with plots
		out_file=open(web_dir+ont+"/"+tp+"_timep/page_"+ont+"_"+tp+"_"+goterm_file+".html", 'w+')
		out_file.write("""<html><head><title>Gene Expression profile</title></head><body><p><font size="5">Gene Expression profile - 7th February 2013</font> <br></p>
	<p><font size="4">GO ontology: """ +ontology+" -> Timepoint: " +tp+ " -> " +goterm_name+ """</font></p><p align="right"> <a href="../New_""" +ont+"_"+tp+ """.html"> Previous page </a> </p>
	<table><tr>""")
		for i in range(len(genes)):
			out_file.write("""<td><img src="""+plot_dir+genes[i]+".png"+ """></a></td></tr></table></body></html>""")
		out_file.close()
	

new_file.close()
f.close()
#os.remove(web_dir +ont+"/"+ file_read + ".html")
