x=c(0, 0.5, 1, 2, 4)
y=c(0, 3.1, 2.6, 2.9, 1.5)
y2=c(0,4.1,4.5,5.3,4.7)
d=-2.41*y-2.20*y2+14.9

splinefunH <-function (x, y, m)                                                                                                   
{                                                                                                                    
    n <- length(x)                                                                                                   
    stopifnot(is.numeric(x), is.numeric(y), is.numeric(m), length(y) ==                                              
        n, length(m) == n)                                                                                           
    if (is.unsorted(x)) {                                                                                            
        i <- sort.list(x)                                                                                            
        x <- x[i]                                                                                                    
        y <- y[i]                                                                                                    
        m <- m[i]                                                                                                    
    }                                                                                                                
    dx <- x[-1L] - x[-n]                                                                                             
    if (any(is.na(dx)) || any(dx == 0))                                                                              
        stop("'x' must be *strictly* increasing (non - NA)")                                                         
    splinefunH0(x, y, m)                                                                                    
}                                                     

splinefunH0 <- function(x, y, m, dx = x[-1L] - x[-length(x)])
{
	u=seq(min(x), max(x), 0.01)
	deriv=0
	extrapol = c("cubic")
	u = seq(min(x), max(x), 0.01)
	i = findInterval(u, x, all.inside = (extrapol == "cubic"))

	interp <- function(u, i) {
		h <- dx[i]
		t <- (u - x[i]) / h
		## Compute the 4 Hermite (cubic) polynomials h00, h01,h10, h11
		t1 <- t-1
		h01 <- t*t*(3 - 2*t)
		h00 <- 1 - h01
		tt1 <- t*t1
		h10 <- tt1 * t1
		h11 <- tt1 * t
		y[i]  * h00 + h*m[i]  * h10 +
		y[i+1]* h01 + h*m[i+1]* h11
	}

	plot(x, y, xlim=c(0,120))
	lines(u, interp(u, i))
}


splinefunH(x,y,d)

func1 = function(x1){
	-2.99*x1+14.8
}

func2 = function(x1, x2){
	-2.41*x1-2.20*x2+14.9
}


find_expression = function(div){
	xvec=seq(0, 4, div)
	func_array = c(func1, func2, func3, func4, fun5, fun6)
	for (j in 1:length(func_array)){
		func = func_array[[j]]
		x_exp=0
		if( j == 1 || j == 5 || j == 6){
			for(i in 2:length(xvec)){
				x_exp=c(x_exp, x_exp[i-1] + func(x_exp[i-1])*div)
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

expression_list = find_expression(0.1)

plot(x, y2, xlim=c(0,4), ylim=c(0,6))
lines(expression_list[[3]], expression_list[[1]])

plot(x, y, xlim=c(0,4), ylim=c(0,6))
lines(expression_list[[3]], expression_list[[2]])




