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
	tab_goi = rbind(tab_goi,rep(1,5))
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

func4  = function(x1, x4){
	W[1,4]*x1 + W[4,4]*x4 + W[7,4]
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

## Plot simulated expression kinetics
plot_kinetics = function(expression_list, goi_mat, div, cond){
	tp = c(0, 20, 60, 120, 240)
	for( i in 1:6 ){
	      exp=goi_mat[,i]
	      png(filename= sprintf("%skinetics_plot_sub.ctrl_log_gene%s_%s.png",cond, i, div))
	      plot( tp, exp, ylim=c(min(min(expression_list[[i]]), min(exp)), max(max(expression_list[[i]]), max(exp))) )
	      lines(expression_list[[7]], expression_list[[i]])
	      dev.off()
	}
}


##Permutation test
mse = function(sim, obs) {mean( (sim - obs)^2, na.rm = FALSE)}


permutation_test = function(count_mat, nb_resample, div){
	tp=c(0, 20, 60, 120, 240)
	mse_sample=c()
	var_sample=c()
	for( j in 1:nb_resample ){
		# random goi selection
		genes_id = sample(1:ncol(count_mat), 6)
		goi = count_mat[,genes_id]
		goi = cbind(goi, rep(1,5))
		# investigate if each column have solution
		used_id = genes_id
		for( i in 1:ncol(goi) ){
			goi_col=goi[,i]
			while(length(which(goi_col[2:5]==0))>0){
				genes_id_new = sample(1:ncol(count_mat), 1)
				if(!(genes_id_new %in% used_id)){
				goi_col = count_mat[,genes_id_new]
				}
			}
			goi[, i] = goi_col
		}
		# get X_dot
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
			known_expression_id = c(known_expression_id, match(tp[i], expression_list[[7]]))
		}
		mse_matrix = c()
		mse_array = c()
		var_array = c()
		for(i in 1:6){
			mse_matrix = rbind(mse_matrix, expression_list[[i]][known_expression_id])
			mse_array=c(mse_array, mse(goi[,i],mse_matrix[i,]))
			var_array = c(var_array, var(goi[,i]))
		}
	mse_sample=c(mse_sample, sum(mse_array/var_array))
	var_sample=c(var_sample, sum(var_array))
	}
	
	return(list(mse_sample, var_sample))
}

measure_goi_mse = function(goi_mat, div){
	tp=c(0, 20, 60, 120, 240)
	mse_sample=c()
	var_sample=c()
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
		known_expression_id = c(known_expression_id, match(tp[i], expression_list[[7]]))
	}
	mse_matrix = c()
	mse_array = c()
	var_array = c()
	for(i in 1:6){
		mse_matrix = rbind(mse_matrix, expression_list[[i]][known_expression_id])
		mse_array=c(mse_array, mse(goi_mat[,i],mse_matrix[i,]))
		var_array = c(var_array, var(goi_mat[,i]))
	}
	mse_sample=c(mse_sample, sum(mse_array/var_array))
	var_sample=c(var_sample, sum(var_array))
	return(list(mse_sample, var_sample))
}


### Resolve gene-gene interaction matrix (Matrix W)
## LM1

# get genes of interest from LM1 condition (Matrix X)
goi_LM1 = get_goi(as.data.frame(normalized_LM1))
goi_LM1_mat = as.matrix(goi_LM1)
goi_LM1_mat = log2(as.matrix(goi_LM1_mat+1))
#colnames(goi_LM1_mat) = row.names(goi_LM1_mat) = NULL
#subtract control
for(i in 1:6){
	goi_LM1_mat[,i] = goi_LM1_mat[,i] - goi_LM1_mat[1,i] 
}
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

expression_list = find_kinetics(0.1, goi_LM1_mat)
plot_kinetics(expression_list, goi_LM1_mat, 0.01, "LM1")

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

nb_resample=1000
div=1
mse_sample = permutation_test(count_mat, nb_resample, div)
mse_goi = measure_goi_mse(goi_LM1_mat, div)
ratio = length(mse_sample[[1]][mse_sample[[1]] < mse_goi[[1]]])
p_value = ratio/nb_resample

png(sprintf("hist.unexplained_%ssamples_%sdev_LM1.png", nb_resample, div))
hist(c(mse_goi[[1]],mse_sample[[1]]), breaks=100, main=sprintf("Permutation test\nLM1 sample with p-value=%s",  p_value),xlab="MSE/Var")
abline(v=mse_goi[[1]], col="red")
dev.off()


hist(c(mse_goi[[1]],mse_sample[[1]]), breaks=400, xlab="MSE/Var", xlim=c(0,250))
abline(v=mse_goi[[1]], col="red", lwd=2)


# Investigate variance distribution
hist(c(mse_goi[[2]],mse_sample[[2]]), breaks=100, main=sprintf("Unexplained variation\nLM1; %s samples; %s dev; p-value=%s", nb_resample, div, p_value),xlab="unexplained variation")
abline(v=mse_goi[[2]], col="red")



### LM2
# get genes of interest from LM1 condition (Matrix X)
goi_LM2 = get_goi(as.data.frame(normalized_LM2))
goi_LM2_mat = as.matrix(goi_LM2)
goi_LM2_mat = log2(as.matrix(goi_LM2_mat+1))
colnames(goi_LM2_mat) = row.names(goi_LM2_mat) = NULL
#subtract control
for(i in 1:6){
	goi_LM2_mat[,i] = goi_LM2_mat[,i] - goi_LM2_mat[1,i] 
}
goi_120=goi_LM2_mat[1:4,]
 
# get X dot from linear interpolation (Matrix X.)
X_dot = get_X_dot(t(goi_LM2_mat))
X_dot_mat = as.matrix(X_dot)

# Define known interaction in W matrix
# 1 -> 2,3,4
# 2 -> NULL
# 3 -> 5
# Infect,1 -> 4 -> NULL
# 5 -> 6
# 6 -> NULL

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

expression_list = find_kinetics(0.01, goi_LM2_mat)
plot_kinetics(expression_list, goi_LM2_mat, 0.01, "LM2")

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

nb_resample=1000
div=1
mse_sample = permutation_test(count_mat, nb_resample, div)
mse_goi = measure_goi_mse(goi_LM2_mat, div)

ratio = length(mse_sample[[1]][mse_sample[[1]] < mse_goi[[1]]])
p_value = ratio/nb_resample

png(sprintf("hist_%ssamples_%sdev_LM2.png", nb_resample, div))
hist(c(mse_goi[[1]],mse_sample[[1]]), breaks=100, main=sprintf("LM2; %s samples; %s dev; p-value=%s", nb_resample, div, p_value))
abline(v=mse_goi, col="red")
#dev.off()

png(sprintf("hist.unexplained_%ssamples_%sdev_LM2.png", nb_resample, div))
hist(c(mse_goi[[1]],mse_sample[[1]]), breaks=100, main=sprintf("Unexplained variation\nLM2; %s samples; %s dev; p-value=%s", nb_resample, div, p_value), xlab="unexplained variation")
abline(v=mse_goi, col="red")
dev.off()

# Investigate variance distribution
hist(c(mse_goi[[2]],mse_sample[[2]]), breaks=100, main=sprintf("Unexplained variation\nLM1; %s samples; %s dev; p-value=%s", nb_resample, div, p_value),xlab="unexplained variation")
abline(v=mse_goi[[2]], col="red")



hist(c(mse_goi[[1]],mse_sample[[1]]), breaks=500, xlab="MSE/Var", xlim=c(0,400))
abline(v=mse_goi[[1]], col="red", lwd=2)
