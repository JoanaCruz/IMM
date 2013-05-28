## Load count data
args <- commandArgs(TRUE)
library( "DESeq" )

args = "LM1_HeLa_all_timep.txt"
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
        row.names = colnames( count_table_read_LM2 ),
        condition = colnames( count_table_read_LM2 ),
        libType = rep("paired-end",ncol(count_table_read_LM2)))
conditions_LM2=dataDesign$condition
libraryType_LM2=dataDesign$libType

data_LM2 = newCountDataSet(count_table_LM2, conditions_LM2)
normalized_LM2 = counts(data_LM2, normalized=TRUE)

## Get genes of interest

# 1 IL1A ENSG00000115008
# 2 HLA-DMA ENSG00000204257
# 3 NFKBIE ENSG00000146232
# 4 CD59 ENSG00000085063
# 5 STAT1 ENSG00000115415
# 6 STAT5A ENSG00000126561

get_goi = function(count_table){
	goi=c("ENSG00000115008", "ENSG00000204257", "ENSG00000146232", "ENSG00000085063", "ENSG00000115415", "ENSG00000126561")
	tab_goi=c()
	genes_col=rownames(count_table)
	for( j in 1:length(goi) ){
		pos_gene =  goi[j] == row.names(count_table)
		line_count_table = as.matrix(count_table[pos_gene,])
		tab_goi = rbind(tab_goi,line_count_table)
	}
	tab_goi = rbind(tab_goi,rep(1,6))
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


### Resolve gene-gene interaction matrix (Matrix W)
## LM1

# get genes of interest from LM1 condition (Matrix X)
goi_LM1 = get_goi(as.data.frame(normalized_LM1))
goi_LM1_mat = as.matrix(goi_LM1)
# goi_LM1_mat = log2(as.matrix(goi_LM1))
# goi_LM1_mat[goi_LM1_mat=="-Inf"] = -10
colnames(goi_LM1_mat) = row.names(goi_LM1_mat) = NULL
goi_120=goi_LM1_mat[1:4,]
 
# get X dot from linear interpolation (Matrix X.)
X_dot = get_X_dot(t(goi_LM1_mat))
X_dot_mat = as.matrix(X_dot)


# Define known interaction in W matrix
# 1 -> 2,3,4
# 2 -> NULL
# 3 -> 5
# Infect,1 -> 4 -> NULL
# 5 -> 6
# 6 -> NULL

get_w = function(goi_120, X_dot_mat){
	#because the first gene is influenced by the infection pseudo-gene
	w = matrix(0 , nrow=7)
	for( i in 1:6){
		if(i==2 || i==3 || i==4 ){
			# select only column with non zero w
			X_select = cbind(goi_120[,1], goi_120[,7])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			w=cbind(w, c(w_selected[1], rep(0,5), w_selected[2]))
		} else if ( i==5 ) {
			# select only column with non zero w
			X_select = cbind(goi_120[,3], goi_120[,7])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			w = cbind(w, c(0, 0, w_selected[1], 0, 0, 0, w_selected[2]))
		} else if ( i==6 ) {
			# select only column with non zero w
			X_select = cbind(goi_120[,5], goi_120[,7])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			w = cbind(w, c(rep(0,4), w_selected[1], 0, w_selected[2]))
		}
	}
	return(w)
	
}


get_w = function(goi_120, X_dot_mat){
	#because the first gene is influenced by the infection pseudo-gene
	w = c()
	for( i in 1:6){
		if( i==1 ){
			# select only column with non zero w
			X_select = cbind(goi_120[,i], goi_120[,7])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			w=cbind(w, c(w_selected[1], rep(0,5), w_selected[2]))
		} else if(i==2 || i==3 || i==4 ){
			# select only column with non zero w
			X_select = cbind(goi_120[,1], goi_120[,i], goi_120[,7])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			if(i==2) { w=cbind(w, c(w_selected[1], w_selected[2], rep(0,4), w_selected[3])) }
			else if (i==3) { w=cbind(w, c(w_selected[1], 0, w_selected[2], rep(0,3), w_selected[3])) }
			else if (i==4) { w=cbind(w, c(w_selected[1], 0, 0, w_selected[2], 0, 0, w_selected[3])) }
		} else if ( i==5 ) {
			# select only column with non zero w
			X_select = cbind(goi_120[,3], goi_120[,i], goi_120[,7])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			w = cbind(w, c(0, 0, w_selected[1], 0, w_selected[2], 0, w_selected[3]))
		} else if ( i==6 ) {
			# select only column with non zero w
			X_select = cbind(goi_120[,5], goi_120[,7])
			X_dot_select = X_dot_mat[,i]
			# resolve
			w_selected = qr.coef(qr(X_select), X_dot_select)
			w = cbind(w, c(rep(0,4), w_selected[1], 0, w_selected[2]))
		}
	}
	return(w)
	
}

W = get_w(goi_120, X_dot_mat)

# #ALL equal
# pinv = function(A){
# 	s <- svd(A)
# 	D <- diag(s$d)
# 	Dinv <- diag(1/s$d)
# 	U <- s$u; V <- s$v
# 	X = V%*%Dinv%*%t(U)
# 	return(X)
# }
#  pinv(X_select)%*%X_dot_select
#            [,1]
# [1,] -0.1080768
# [2,]  2.6945242
# > lm.fit(X_select, X_dot_select)$coefficients
#         x1         x2 
# -0.1080768  2.6945242 
# > qr.coef(qr(X_select), X_dot_select)
# [1] -0.1080768  2.6945242


## Simulate kinetics

# Set functions

func1 = function(x1){
	W[1,1]*x1+W[7,1]
}

func2 = function(x1, x2){
	W[1,2]*x1 + W[2,2]*x2 + W[7,2]
}

func3 = function(x1, x3){
	W[1,3]*x1 + W[3,3]*x3 + W[7,3]
}

func4  = function(x1, x3){
	W[1,4]*x1 + W[4,4]*x3 + W[7,4]
}

func5  = function(x1, x5){
	W[3,5]*x1 + W[5,5]*x5 + W[7,5]
}

func6  = function(x5){
	W[5,6]*x5 + W[7,6]
}

# Numerical integration. Each integral is calculated from div to div.
find_kinetics = function(div, goi_mat){
	xvec=seq(0, 240, div)
	func_array = c(func1, func2, func3, func4, func5, func6)
	expression_list=list()
	for (j in 1:length(func_array)){
		func = func_array[[j]]
		x_exp = goi_mat[1,j]
		if( j == 1 ){
			# define x(0+div) or x[i]
			for(i in 2:length(xvec)){
				x_exp=c(x_exp, x_exp[i-1] + func(x_exp[i-1])*div)
			}
			expression_list[[j]] = x_exp
		} else if ( j == 5 ) {
			for(i in 2:length(xvec)){
				x_exp=c(x_exp, x_exp[i-1] + func(expression_list[[3]][i-1], x_exp[i-1])*div)
			}
			expression_list[[j]] = x_exp
		} else if ( j == 6 ) {
			for(i in 2:length(xvec)){
				x_exp=c(x_exp, x_exp[i-1] + func(expression_list[[5]][i-1])*div)
			}
			expression_list[[j]] = x_exp
		} else {
			for(i in 2:length(xvec)){
				x_exp=c(x_exp, x_exp[i-1] + func(expression_list[[1]][i-1], x_exp[i-1])*div)
			}
			expression_list[[j]] = x_exp
		}
	}
	expression_list[[j+1]]=xvec
	return(expression_list)
}

expression_list = find_kinetics(0.01, goi_LM1_mat)


## Plot simulated expression kinetics
plot_kinetics = function(expression_list, goi_mat, div){
	tp = c(0, 20, 60, 120, 240)
	for( i in 1:6 ){
	      exp=goi_mat[,i]
	      png(filename= sprintf("kinetics_plot_gene%s_%s.png", i, div))
	      plot(tp, exp)
	      lines(expression_list[[7]], expression_list[[i]])
	      dev.off()
	}
}

plot_kinetics(expression_list, goi_LM1_mat, 0.01)


# gene 1
exp = goi_LM1_mat[,1]
tp = c(0, 20, 60, 120, 240)

plot(tp, exp)
lines(expression_list[[7]], expression_list[[1]])

