import sys
import csv
import os

if len(sys.argv) < 3:
	print "Usage: csvToTable.py csv_file html_file"
	exit(1)

# Open the CSV file for reading
reader = csv.reader(open(sys.argv[1]))

web_dir=sys.argv[2]
#create first page1
f1=open(web_dir+"page1.html", 'w+')

f1.write('<html><head><title>Gene Expression profile</title></head><body><font size="5">Gene Expression profile - 7th February 2013</font> <br><br><font size="4">GO ontologies:</font><p> <a href="'+web_dir+'BP/pageBP_timepoints.html">Biological Processes</a></p><p> <a href="'+web_dir+'CC/pageCC_timepoints.html">Celular components</a></p><p><a href="'+web_dir+'MF/pageMF_timepoints.html">Molecular function</a></p></body></html>')

f1.close()

# Create the HTML file for output
file_dir_name=sys.argv[1].split("/")
file_out_dir=sys.argv[2]
file_name=file_dir_name[1]
ont=file_name.split("_")[1]
tp=file_name.split("_")[2].split(".")[0]
file_out=file_out_dir+ont

if not os.path.exists(file_out):
    os.makedirs(file_out)

tp_html=open(file_out+"/page"+ont+"_"+"timepoints.html","w+")

if ont=="BP":
	ontology="Biological processes"
elif ont=="CC":
	ontology="Celular components"
elif ont=="MF":
	ontology="Molecular function"


tp_html.write("""<html>
<head>
<title>Gene Expression profile</title>
</head>
<body>
<font size="5">Gene Expression profile - 7th February 2013</font> <br>
<br>
<font size="4">GO ontologies: """+ontology+""" -> Timepoint: </font>
<p> <a href="New_page_"""+ont+"""_20.html">20 minutes</a> </p>
<p> <a href="New_page_"""+ont+"""_60.html">60 minutes</a> </p>
<p> <a href="New_page_"""+ont+"""_120.html">120 minutes</a> </p>
<p> <a href="New_page_"""+ont+"""_240.html">240 minutes</a> </p>
<br>
<p> <a href="../page1.html"> Previous page </a> </p>
</body>
</html>""")

tp_html.close()


htmlfile = open(file_out+"/page"+"_"+ont+"_"+tp+".html","w")

# initialize rownum variable
rownum = 0

# write <table> tag
htmlfile.write('<table rules="all" CELLSPACING="0" COLS="5" BORDER="1">' + "\n" )

for row in reader: # Read a single row from the CSV file
	if rownum == 0:
		htmlfile.write('<tr>'+ "\n") # write <tr> tag
  		for (i,column) in enumerate(row):
			if i==0 or i==1 or i==2 or i==7 or i==8:
				htmlfile.write('<th>' + column + '</th>' + "\n" ) # write header columns
		htmlfile.write('</tr>'+ "\n") # write </tr> tag
  	else: # write all other rows
  		colnum = 1
		if rownum % 2 == 0:
			htmlfile.write('<tr class="color1">' + "\n")
		else:
			htmlfile.write('<tr class="color2">' + "\n")
		for (i,column) in enumerate(row):
			if i==0 or i==1 or i==2 or i==7 or i==8:
				htmlfile.write('<td class="column_' + str(colnum) + '">' + column + '</td>'+ "\n")
				colnum += 1
		htmlfile.write('</tr>' + "\n")

	rownum += 1

# write </table> tag
htmlfile.write('</table>')
exit(0)
