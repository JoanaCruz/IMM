## Load count data
args <- commandArgs(TRUE)
library( "DESeq" )

args = c("LM1_HeLa_all_timep.txt", "LM2_HeLa_all_timep.txt")
count_table_LM1 = read.table(file = args[1], header=TRUE, row.names=1)
count_table_LM2 = read.table(file = args[2], header=TRUE, row.names=1)


## Define data design
# LM1
dataDesign_LM1 = data.frame(
        row.names = colnames( count_table_LM1 ),
        condition = colnames( count_table_LM1 ),
        libType = rep("paired-end",ncol(count_table_LM1)))
conditions_LM1=dataDesign_LM1$condition
libraryType_LM1=dataDesign_LM1$libType

data_LM1 = newCountDataSet(count_table_LM1, conditions_LM1)
data_LM1 = estimateSizeFactors(data_LM1)
normalized_LM1 = counts(data_LM1, normalized=TRUE)


# LM2
dataDesign_LM2 = data.frame(
        row.names = colnames( count_table_LM2 ),
        condition = colnames( count_table_LM2 ),
        libType = rep("paired-end",ncol(count_table_LM2)))
conditions_LM2=dataDesign_LM2$condition
libraryType_LM2=dataDesign_LM2$libType

data_LM2 = newCountDataSet(count_table_LM2, conditions_LM2)
data_LM2 = estimateSizeFactors(data_LM2)
normalized_LM2 = counts(data_LM2, normalized=TRUE)

## Get genes of interest

get_ensembl_genes = function(gene_array){
	ensembl_array=c()
	for( i in 1:length(gene_array)){
		print(gene_array[i])
		#system(sprintf("wget http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s -O %s", gene_array[i], gene_array[i]))
		found=FALSE
		f=file(gene_array[i], open = "r")
		while (length(oneLine <- readLines(f, n = 1, warn = FALSE)) > 0 && !found) {
			if (grepl('>ENSG', oneLine)){
				ensembl_id = strsplit(strsplit(oneLine, split=">")[[1]][2], "<")[[1]][1]	
				ensembl_array=c(ensembl_array, ensembl_id)
				found=TRUE
			}
		}
		print(i)
		print(ensembl_array)
		close(f)
		#system(sprintf("rm %s", gene_array[i]))
		
	}
	return(ensembl_array)
}

get_ensembl_genes(gene_array)

# x1 NDGR1
# 
# ATF3
# IL8
# IL4
# CCL2
# IFNG
# IL5
# CHOP "DDIT3"CHOP
# TNF-Alpha TNF
# IL12
# 
# IKKB IKBKB
# NFKB NFKB1
# IL8

get_goi = function(count_table){
	gene_array=c("NDRG1", "ATF3", "IL8", "IL4", "CCL2", "IFNG", "IL5", "DDIT3", "TNF", "IL13", "IL12A", "IKBKB", "NFKB1")
	#goi=get_ensembl_genes(gene_array)
	goi=c("ENSG00000104419", "ENSG00000162772", "ENSG00000169429", "ENSG00000113520", "ENSG00000108691", "ENSG00000111537", "ENSG00000113525", "ENSG00000175197", "ENSG00000232810", "ENSG00000169194", "ENSG00000168811", "ENSG00000104365", "ENSG00000109320")

	tab_goi=c()
	genes_col=rownames(count_table)
	for( j in 1:length(goi) ){
		pos_gene =  goi[j] == row.names(count_table)
		line_count_table = as.matrix(count_table[pos_gene,])
		tab_goi = rbind(tab_goi,line_count_table)
	}
	return(t(tab_goi))
}

get_X_dot = function(count_table){
	X_dot = c()
	for (i in 1:nrow(count_table)){
		row = t(as.matrix(count_table[i,]))
		colnames(row) = NULL
		tp=c(0, 20, 60, 120, 240)
		linear_interpol_vec = c()
		for (j in 1:(length(row)-1)){
			adj_tp = c(tp[j], tp[j+1])
			adj_exp = c(row[1,j],row[1,j+1])
			reg = lm( adj_exp ~ adj_tp )$coefficients[2]
			linear_interpol_vec=c(linear_interpol_vec, reg)
		}
		linear_interpol_vec=as.matrix(linear_interpol_vec)
		row.names(linear_interpol_vec) = NULL
		X_dot = rbind(X_dot, t(linear_interpol_vec))
	}
	return(t(X_dot))
}


## Get interaction matrix
get_w = function(goi_120, X_dot_mat){
	w = c()
	for( i in 1:13){
		if( i==1 ){
			# select only column with non zero w
			X_select = cbind(goi_120[,i], goi_120[,8])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			w=cbind(w, c(w_selected[1], rep(0,6), w_selected[2]))
		} else if ( i==2 || i==6 ) {
			# select only column with non zero w
			X_select = cbind(goi_120[,1], goi_120[,i], goi_120[,8])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			if(i==2) { w=cbind(w, c(w_selected[1], w_selected[2], rep(0,5), w_selected[3])) }
			if (i==6) { w=cbind(w, c(w_selected[1], rep(0,4), w_selected[2], 0, w_selected[3])) }
		} else if ( i==3 ) {
			# select only column with non zero w
			X_select = cbind(goi_120[,2], goi_120[,i], goi_120[,7], goi_120[,8])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			w = cbind(w, c(0, w_selected[1], w_selected[2], rep(0,3), w_selected[3], w_selected[4]))
		} else if ( i==7 ) {
			# select only column with non zero w
			X_select = cbind(goi_120[,6], goi_120[,i], goi_120[,8])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			w = cbind(w, c(rep(0,5), w_selected[1], w_selected[2], w_selected[3]))
		} else if ( i==4 || i==5) {
			# select only column with non zero w
			X_select = cbind(goi_120[,2], goi_120[,i], goi_120[,8])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			if(i==4) { w=cbind(w, c(0, w_selected[1], 0, w_selected[2], rep(0,3), w_selected[3])) }
			if (i==5) { w=cbind(w, c(0, w_selected[1], 0, 0, w_selected[2], rep(0,2), w_selected[3])) }
			#w=cbind(w, c(0, w_selected[1], rep(0,i-3), w_selected[2], rep(0,13-i), w_selected[3]))
		}
	}
	return(w)
}

## Simulate kinetics
# Set functions

func1 = function(x1){
	W[1,1]*x1+W[8,1]
}

func2 = function(x1, x2){
	W[1,2]*x1 + W[2,2]*x2 + W[8,2]
}

func3  = function(x2, x3, x7){
	W[2,3]*x2 + W[3,3]*x3 + W[7,3]*x7 + W[8,3]
}

func4  = function(x2, x4){
	W[2,4]*x2 + W[4,4]*x4 + W[8,4]
}

func5  = function(x2, x5){
	W[2,5]*x2 + W[5,5]*x5 + W[8,5]
}

func6 = function(x1, x6){
	W[1,6]*x1 + W[6,6]*x6 + W[8,6]
}

func7  = function(x6, x7){
	W[6,7]*x6 + W[7,7]*x7 + W[8,7]
}

# Numerical integration. Each integral is calculated from div to div.
find_kinetics = function(div, goi_mat){
	xvec=seq(0, 240, div)
	func_array = c(func1, func2, func4, func5, func6, func7,  func3)
	expression_list=list()
	order=c(1, 2, 4, 5, 6, 7, 3)
	for (j in 1:length(func_array)){
		func = func_array[[j]]
		x_exp = goi_mat[1,order[j]]
		if( order[j] == 1 ){
			# define x(0+div) or x[i]
			for(i in 2:length(xvec)){
				x_exp=c(x_exp, x_exp[i-1] + func(x_exp[i-1])*div)
			}
			expression_list[[order[j]]] = x_exp
		} else if ( order[j] == 2 || order[j] == 5) {
			for(i in 2:length(xvec)){
				x_exp=c(x_exp, x_exp[i-1] + func(expression_list[[1]][i-1], x_exp[i-1])*div)
			}
			expression_list[[order[j]]] = x_exp
		} else if ( order[j] == 3 ) {
			for(i in 2:length(xvec)){
				x_exp=c(x_exp, x_exp[i-1] + func(expression_list[[2]][i-1],  x_exp[i-1], expression_list[[7]][i-1])*div)
			}
			expression_list[[order[j]]] = x_exp
		} else if ( order[j] == 7 ) {
			for(i in 2:length(xvec)){
				x_exp=c(x_exp, x_exp[i-1] + func(expression_list[[6]][i-1],  x_exp[i-1])*div)
			}
			expression_list[[order[j]]] = x_exp
		} else {
			for(i in 2:length(xvec)){
				x_exp=c(x_exp, x_exp[i-1] + func(expression_list[[2]][i-1], x_exp[i-1])*div)
			}
			expression_list[[order[j]]] = x_exp
		}
	}
	expression_list[[j+1]]=xvec
	#expression_list_ordered=
	return(expression_list)
}

expression_list = find_kinetics(10, goi_LM1_mat)

## Plot simulated expression kinetics
plot_kinetics = function(expression_list, goi_mat, div, cond){
	genes=c("NDRG1", "ATF3", "IL8", "CCL2", "DDIT3", "IKBKB", "NFKB1")
	tp = c(0, 20, 60, 120, 240)
	png(filename=sprintf('%s_kinetics_NDRG1_%s.png', cond, div))
	par(mfrow=c(3,3))
	for( i in 1:7 ){
	      exp=goi_mat[,i]
	      #png(filename= sprintf("%skinetics_plot_sub.ctrl_log_gene%s_%s.png",cond, i, div))
	      plot( tp, exp, main=genes[i], ylim=c(min(min(expression_list[[i]]), min(exp)), max(max(expression_list[[i]]), max(exp))) )
	      lines(expression_list[[8]], expression_list[[i]])
	      #dev.off()
	}
	dev.off()
}


##Permutation test
mse = function(sim, obs) {mean( (sim - obs)^2, na.rm = FALSE)}


find_goi = function(count_mat, used_index, cutoff){
# investigate if each column have solution and if variation of genes set is > 10
	genes_sol = c()
	found=TRUE
	while ( found==TRUE ){
	# get genes which entries are all non-zero
		for( i in 1:ncol(count_mat) ){
			have_zero = length(which(count_mat[2:5,i]==0))>0
			if (!have_zero){
				genes_sol = cbind(genes_sol,count_mat[,i])
			}
		}
	#get genes which variance is greater than 10
		acept_var=cutoff
		not_found=TRUE
		while(not_found==TRUE){
			genes_index = sample(1:ncol(genes_sol), 7)
			goi = genes_sol[,genes_index]
			var_array = c()
			for(i in 1:ncol(goi)){
				var_array = c(var_array, var(goi[,i]))
			}
			if(sum(var_array)>acept_var){not_found=FALSE}
		}
	# Investigate if this gene set already was tested
	# se o primeiro for diferente do que esta a ser testado vai dar false e caba
	# o 2o pode ser igual (dá true)
		if(length(used_index)==0){
			found=FALSE
		} else {
			found_array=c()
			for(i in 1:ncol(used_index)){
				found_array = c(found_array, identical(genes_index, used_index[,i]))
			}
			if(length(unique(found_array))==1){
				found=FALSE
				used_index=cbind(used_index, genes_index)
			}
		}
	}
	return(list(goi, used_index))
}

# mat1=c(16702, 11122, 18739, 8642, 17141, 6591, 22223)
# mat2=c(37214, 39293, 4617, 25342, 23984, 12197, 21151)
# used_index=cbind(mat1, mat2)
# genes_index=mat1

permutation_test = function(count_mat, nb_resample, div, cutoff){
	tp=c(0, 20, 60, 120, 240)
	mse_sample=c()
	var_sample=c()
	mse1_sample=c()
	used_index=c()
	for( j in 1:nb_resample ){
		print(j)
		#c(genes_index, goi, used_index)
		list = find_goi(count_mat, used_index, cutoff)
		# get goi
		goi=cbind(list[[1]], rep(1,5))
		used_index=list[[2]]
		#get X_dot
		X_dot = get_X_dot(t(goi))
		X_dot_mat = as.matrix(X_dot)
		# get W
		goi_120 = goi[1:4,]
		W = get_w(goi_120, X_dot_mat)
		# Integrate
		expression_list = find_kinetics(div, goi)
		# measeure mse
		known_expression_id = c()
		for( i in 1:length(tp)){
			known_expression_id = c(known_expression_id, match(tp[i], expression_list[[8]]))
		}
		mse_matrix = c()
		mse_array = c()
		var_array = c()
		for(i in 1:7){
			mse_matrix = rbind(mse_matrix, expression_list[[i]][known_expression_id])
			mse_array=c(mse_array, mse(goi[,i],mse_matrix[i,]))
			var_array = c(var_array, var(goi[,i]))
		}
	mse_sample=c(mse_sample, sum(mse_array/var_array))
	var_sample=c(var_sample, sum(var_array))
	mse1_sample=c(mse1_sample, sum(mse_array))
	}
	return(list(mse_sample, var_sample, mse1_sample))
}

measure_goi_mse = function(goi_mat, div){
	tp=c(0, 20, 60, 120, 240)
	mse_sample=c()
	var_sample=c()
	mse1_sample=c()
	# get X_dot
	X_dot = get_X_dot(t(goi_mat))
	X_dot_mat = as.matrix(X_dot)
	# get W
	goi_120 = goi_mat[1:4,]
	W = get_w(goi_120, X_dot_mat)
	# Integrate
	expression_list = find_kinetics(div, goi_mat)
	# measeure mse
	known_expression_id = c()
	for( i in 1:length(tp)){
		known_expression_id = c(known_expression_id, match(tp[i], expression_list[[8]]))
	}
	mse_matrix = c()
	mse_array = c()
	var_array = c()
	for(i in 1:7){
		mse_matrix = rbind(mse_matrix, expression_list[[i]][known_expression_id])
		mse_array=c(mse_array, mse(goi_mat[,i],mse_matrix[i,]))
		var_array = c(var_array, var(goi_mat[,i]))
	}
	mse_sample= sum(mse_array/var_array)
	var_sample= sum(var_array)
	mse1_sample= sum(mse_array)
	return(list(mse_sample, var_sample, mse1_sample))
}



### Resolve gene-gene interaction matrix (Matrix W)
## LM1

# get genes of interest from LM1 condition (Matrix X)
goi_LM1 = get_goi(as.data.frame(normalized_LM1))
goi_LM1_mat = as.matrix(goi_LM1)
goi_LM1_mat = log2(as.matrix(goi_LM1_mat+1))
#colnames(goi_LM1_mat) = row.names(goi_LM1_mat) = NULL
#subtract control
zero=c()
for(i in 1:13){
	goi_LM1_mat[,i] = goi_LM1_mat[,i] - goi_LM1_mat[1,i]
	if(goi_LM1_mat[2,i]==0){
		zero=c(zero, i)
	}
}
goi_LM1_mat=goi_LM1_mat[,-zero]
goi_LM1_mat = cbind(goi_LM1_mat, rep(1,5))
goi_120=goi_LM1_mat[1:4,]

# get X dot from linear interpolation (Matrix X.)
X_dot = get_X_dot(t(goi_LM1_mat))
X_dot_mat = as.matrix(X_dot)


# Define known interaction in W matrix
W = get_w(goi_120, X_dot_mat)


## Simulate kinetics

expression_list = find_kinetics(1, goi_LM1_mat)
plot_kinetics(expression_list, goi_LM1_mat, 0.1, "LM1")

tp = c(0, 20, 60, 120, 240)
     plot( tp, exp, ylim=c(min(min(expression_list[[i]]), min(exp)), max(max(expression_list[[i]]), max(exp))) )
	      lines(expression_list[[7]], expression_list[[i]])
	      dev.off()


# # gene 1
exp = goi_LM1_mat[,1]
tp = c(0, 20, 60, 120, 240)

plot(tp, exp,ylab="log2(expression)", xlab="Time [min]")
lines(expression_list[[7]], expression_list[[1]])



tp = c(0, 20, 60, 120, 240)
goi_mat=goi_LM1_mat
gene=c('IL1A', 'HLA-DMA', 'NFKBIE', 'CD59', 'STAT1', 'STAT5A')
par(mfrow=c(2,3))
	for( i in 1:6 ){
	      print(i)
	      exp=goi_mat[,i]
	      plot( tp, exp, col="blue", pch=19,main=gene[i], xlab=NA, ylab=NA, ylim=c(min(min(expression_list[[i]]), min(exp)), max(max(expression_list[[i]]), max(exp))) )
	      lines(expression_list[[7]], expression_list[[i]])
	      mtext("Time [min]", side=1, outer=TRUE, line=-2)
	mtext("log2(expression)", side=2, outer=TRUE, line=-1.5)
	  }
	      dev.off()

## Permutation test

# Set all genes as log and subtract control value
count_mat = as.matrix(normalized_LM1)
count_mat = t(log2(count_mat+1))
#colnames(count_mat) = row.names(count_mat) = NULL
#subtract control
for(i in 1:ncol(count_mat)){
	count_mat[,i] = count_mat[,i] - count_mat[1,i] 
}

cutoff=0
nb_resample=100
div=1
mse_sample = permutation_test(count_mat, nb_resample, div,10)
mse_goi = measure_goi_mse(goi_LM1_mat, div)

ratio1 = length(mse_sample[[1]][mse_sample[[1]] < mse_goi[[1]]])
p_value = ratio1/nb_resample

png(sprintf("hist.unexplained_%ssamples_%sdev_LM1_MSEvar_cutoff=%s.png", nb_resample, div, cutoff))
hist(c(mse_goi[[1]],mse_sample[[1]]), breaks=50, main=sprintf("Permutation test\nLM1 sample with p-value=%s", p_value),xlab="MSE/Var")
abline(v=mse_goi[[1]], col="red")
dev.off()

# Investigate variance distribution
png(sprintf("hist.unexplained_%ssamples_%sdev_LM1_Var_cutoff=%s.png", nb_resample, div, cutoff))
hist(c(mse_goi[[2]],mse_sample[[2]]), breaks=50, main=sprintf("Permutation test\nLM1 sample", nb_resample, div),xlab="Var")
abline(v=mse_goi[[2]], col="red")
dev.off()

#
ratio2 = length(mse_sample[[3]][mse_sample[[3]] < mse_goi[[3]]])
p_value2 = ratio2/nb_resample

png(sprintf("hist.unexplained_%ssamples_%sdev_LM1_MSE_cutoff=%s.png", nb_resample, div, cutoff))
hist(c(mse_goi[[3]],mse_sample[[3]]), breaks=50, main=sprintf("Permutation test\nLM1 sample with p-value=%s", p_value2),xlab="MSE") 
abline(v=mse_goi[[3]], col="red")
dev.off()



### LM2
# get genes of interest from LM1 condition (Matrix X)
goi_LM2 = get_goi(as.data.frame(normalized_LM2))
goi_LM2_mat = as.matrix(goi_LM2)
goi_LM2_mat = log2(as.matrix(goi_LM2_mat+1))
#colnames(goi_LM2_mat) = row.names(goi_LM2_mat) = NULL
#subtract control
zero=c()
for(i in 1:13){
	goi_LM2_mat[,i] = goi_LM2_mat[,i] - goi_LM2_mat[1,i]
	if(goi_LM2_mat[2,i]==0){
		zero=c(zero, i)
	}
}
goi_LM2_mat=goi_LM2_mat[,-zero]
goi_LM2_mat = cbind(goi_LM2_mat, rep(1,5))
goi_120=goi_LM2_mat[1:4,]

 
# get X dot from linear interpolation (Matrix X.)
X_dot = get_X_dot(t(goi_LM2_mat))
X_dot_mat = as.matrix(X_dot)

# Define known interaction in W matrix
W = get_w(goi_120, X_dot_mat)

## Simulate kinetics

expression_list = find_kinetics(0.01, goi_LM2_mat)
plot_kinetics(expression_list, goi_LM2_mat, 0.1, "LM2")

# # gene 1
# exp = goi_LM1_mat[,1]
# tp = c(0, 20, 60, 120, 240)
# 
# plot(tp, exp)
# lines(expression_list[[7]], expression_list[[1]])

tp = c(0, 20, 60, 120, 240)
goi_mat=goi_LM2_mat
gene=c('IL1A', 'HLA-DMA', 'NFKBIE', 'CD59', 'STAT1', 'STAT5A')
par(mfrow=c(2,3))
	for( i in 1:6 ){
	      print(i)
	      exp=goi_mat[,i]
	      plot( tp, exp, col="blue", pch=19,main=gene[i], xlab=NA, ylab=NA, ylim=c(min(min(expression_list[[i]]), min(exp)), max(max(expression_list[[i]]), max(exp))) )
	      lines(expression_list[[7]], expression_list[[i]])
	      mtext("Time [min]", side=1, outer=TRUE, line=-2)
	mtext("log2(expression)", side=2, outer=TRUE, line=-1.5)
	  }

## Permutation test

# Set all genes as log and subtract control value
count_mat = as.matrix(normalized_LM2)
count_mat = t(log2(count_mat+1))
#colnames(count_mat) = row.names(count_mat) = NULL
#subtract control
for(i in 1:ncol(count_mat)){
	count_mat[,i] = count_mat[,i] - count_mat[1,i] 
}

nb_resample=100
div=1
cutoff=0
mse_sample = permutation_test(count_mat, nb_resample, div, cutoff)
mse_goi = measure_goi_mse(goi_LM2_mat, div)

ratio1 = length(mse_sample[[1]][mse_sample[[1]] < mse_goi[[1]]])
p_value = ratio1/nb_resample

png(sprintf("hist.unexplained_%ssamples_%sdev_LM2_MSEvar_cutoff=%s.png", nb_resample, div, cutoff))
hist(c(mse_goi[[1]],mse_sample[[1]]), breaks=50, main=sprintf("Permutation test\nLM1 sample with p-value=%s", p_value),xlab="MSE/Var")
abline(v=mse_goi[[1]], col="red")
dev.off()

# Investigate variance distribution
png(sprintf("hist.unexplained_%ssamples_%sdev_LM2_Var_cutoff=%s.png", nb_resample, div, cutoff))
hist(c(mse_goi[[2]],mse_sample[[2]]), breaks=50, main=sprintf("Permutation test\nLM1 sample", nb_resample, div),xlab="Var")
abline(v=mse_goi[[2]], col="red")
dev.off()

#
ratio2 = length(mse_sample[[3]][mse_sample[[3]] < mse_goi[[3]]])
p_value2 = ratio2/nb_resample

png(sprintf("hist.unexplained_%ssamples_%sdev_LM2_MSE_cutoff=%s.png", nb_resample, div, cutoff))
hist(c(mse_goi[[3]],mse_sample[[3]]), breaks=50, main=sprintf("Permutation test\nLM1 sample with p-value=%s", p_value2),xlab="MSE") 
abline(v=mse_goi[[3]], col="red")
dev.off()
