library(survcomp)
library(genefu)

Get_DI_continous = function( surv , time , time_censor , variable ){
	data = data.frame( surv=surv , time=time , variable=variable )
	data$time = as.numeric(as.character(data$time))
	data$variable = as.numeric( as.character(data$variable) )

  	for(i in 1:nrow(data)){
	    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time_censor ){
	      data[ i , "time" ] = time_censor
	      data[ i , "surv" ] = 0
    	}
  	}

	d_i = D.index( x=data$variable , surv.time=data$time , surv.event=data$surv , na.rm=TRUE )

	coef = d_i$coef
	se = d_i$se 

	c( round( coef , 2 ) , round( se , 2 ) , round( coef - (1.96 * se ) , 2 ) , round( coef + (1.96 * se ) , 2 ) , d_i$p.value )
}