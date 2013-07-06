# Create html pages to see DESEQ's tables
import sys
import os
import glob
import csv

# Page with conditions
cond_pag=open('cond_pag.html', 'w+')

cond_pag.write('<html><head><title>Differential expressed genes</title></head><body>')
cond_pag.write('<font size="5">Differential expressed genes - DESEQ output</font>')
cond_pag.write('<br><br><font size="4">Conditions:</font>')
cond_pag.write('<p><a href="LM1/timepoints.html">Control vs LM1 (Unifected vs Infected-wildtype)</a></p>')
cond_pag.write('<p><a href="LM2/timepoints.html">Control vs LM2 (Unifected vs Infected-mutant) </a></p>')
cond_pag.write('<p><a href="LM1.vs.LM2/timepoints.html">LM1 vs LM2 (Infected-wildtype vs Infected-mutant)</a></p>')
cond_pag.write('<p><a href="LM_adjacent/timepoints.html">Adjacent time points</a></p>')
cond_pag.write('</body></html>')
cond_pag.close()

# Page with timepoints
i=['LM1', 'LM2', 'LM1.vs.LM2', 'LM_adjacent']

for cond in i:
	if not os.path.exists(cond):
		os.makedirs(cond)
	if cond == 'LM1': condition = 'Unifected vs Infected-wildtype (Control vs LM1)'
	elif cond == 'LM2': condition = 'Unifected vs Infected-mutant (Control vs LM2)'
	elif cond == 'LM1.vs.LM2': condition = 'Infected-wildtype vs Infected-mutant (LM1 vs LM2)'
	elif cond == 'LM_adjacent': condition ='Adjacent time points'
	tp_pag=open(cond + 'timepoints.html', 'w+')
	tp_pag.write('<html><head><title>Differential expressed genes</title></head><body>')
	tp_pag.write('<font size="5">Differential expressed genes - DESEQ output</font>')
	tp_pag.write('<br><br><font size="4">Condition: '+ condition+ '</font>')
	if i == 'LM_adjacent':
		samples=['LM1', 'LM2']
		for var in samples:
			tp_pag.write('<p><a href="LM_adjacent/'+var+'_20.vs.60.html">'+var+': 20 vs 60 TP</a> </p>')
			tp_pag.write('<p><a href="LM_adjacent/'+var+'_60.vs.120.html">'+var+': 60 vs 120 TP</a> </p>')
			tp_pag.write('<p><a href="LM_adjacent/'+var+'_120.vs.240.html">'+var+': 120 vs 240 TP</a> </p>')
	else:
		tp_pag.write('<p><a href="'+cond+'/'+cond+'_20.html">'+cond+': Timepoint 20</a> </p>')
		tp_pag.write('<p><a href="'+cond+'/'+cond+'_60.html">'+cond+': Timepoint 60</a> </p>')
		tp_pag.write('<p><a href="'+cond+'/'+cond+'_120.html">'+cond+': Timepoint 120</a> </p>')
		tp_pag.write('<p><a href="'+cond+'/'+cond+'_240.html">'+cond+': Timepoint 240</a> </p>')
tp_pag.close()



# Convert csv tables to html
def csv2html(directory):
	files=glob.glob('./*.csv')
	for f in files:
		reader = csv.reader(open(f))
		htmlfile=open(cond+'/.html', 'w+')
		# initialize rownum variable
		rownum = 0
		# write <table> tag
		htmlfile.write('<table rules="all" CELLSPACING="0" COLS="5" BORDER="1">\n' )
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
		
		
csv2html('./')