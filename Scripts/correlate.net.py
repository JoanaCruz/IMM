from sys import argv
import linecache


#extract go terms related with gene from NCBI file
def extract_go(input_file_data, nb_line):
	i=nb_line+3
	bol=True
	go_data=[]
	while bol:
		go_line=linecache.getline(input_file_data, i)
		go_line_find=go_line.find("GO:")
		if go_line_find!=-1:
			go_data.append(go_line[go_line_find:go_line_find+10])
		else:
			bol=False
		i=i+1
	return go_data,i

# extract mrna data related with gene from NCBI file
def extract_mrna(input_file_data, nb_line):
	i=nb_line+1
	bol=True
	mrna_data=[]
	while bol:
		mrna_line=linecache.getline(input_file_data, i)
		if mrna_line[0:4]=="\tNM_":
			line_with_mrna=mrna_line.split()[0].split(".")[0]
			mrna_data.append(line_with_mrna)
		else:
			bol=False
		i=i+1
	return mrna_data,i

# mRNA and goterms associated with the gene being analyzed (summary downloaded from NCBI site)
def get_data(input_file_data, read_data_file):
	#initialize variables
	parsing_go=False
	parsing_mrna=False
	go_data_MF, go_data_BP, go_data_CC, mrna_data = [[]]*4

	for (i, line) in enumerate(read_data_file):
		if line=="Gene Ontology Provided by GOA (www.ebi.ac.uk/GOA)\n":
			parsing_go=True
		if line=="GENERAL PROTEIN INFORMATION\n":
			parsing_go=False	
		if line==" mRNA and Protein(s)\n":
			parsing_mrna=True
		if line=="RELATED SEQUENCES\n":
			parsing_mrna=False
		if parsing_go:
			if len(line.split())==4 and line.split()[0]=="Function":
				go_data_MF, new_i = extract_go(input_file_data, i)
				i=new_i
			if len(line.split())==4 and line.split()[0]=="Processs":
				go_data_BP, new_i = extract_go(input_file_data, i)
				i=new_i
			if len(line.split())==4 and line.split()[0]=="Component":
				go_data_CC, new_i = extract_go(input_file_data, i)
				i=new_i
		if parsing_mrna:
			mrna_data, new_i = extract_mrna(input_file_data, i)
			i=new_i
	return go_data_MF,go_data_BP, go_data_CC, mrna_data
			
# search data in GOStats files
def find_terms(input_file_data, read_data_file, tp):
	go_data_MF,go_data_BP, go_data_CC, mrna_data=get_data(input_file_data, read_data_file)
	for ont in ["MF","BP","CC"]:
		print(ont)
		read_go_file=open(go_directory+"/GOsummary_"+ont+"_"+tp+".csv")
		search_result=open(out_directory+"network_seach_"+tp+"_"+ont+".txt", "w")
		for (j,line) in enumerate(read_go_file):
			go_line=line.split(",")
			go_term=go_line[1].strip('"')
			mrna=go_line[7]
			go_data_ont=eval("go_data_"+ont)
			found=0
			for i in range(0,len(go_data_ont)):
				if j==1 and i==0:
					print(go_term==go_data_ont[i])
				if go_term==go_data_ont[i]:
					search_result.write(line)
					found=found+1
	read_go_file.close()
	search_result.close()
	return found

script, input_file_data, input_go_file, out_directory =argv

read_data_file=open(input_file_data)
go_directory=input_go_file.split("/")[0]
tp=input_go_file.split("_")[2].split(".")[0]

# Do table with terms found
find_terms(input_file_data, read_data_file, tp)

read_data_file.close()
