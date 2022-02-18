library(survcomp)
library(genefu)


###########################################################################################
###########################################################################################

Get_HR_dicho = function( surv , time , time_censor , variable , cutoff , title , ylab , dir ){
	data = data.frame( surv=surv , time=time , variable=variable )
	data$time = as.numeric(as.character(data$time))
	data$variable = ifelse( as.numeric(as.character(data$variable)) >= cutoff , 1 , 0 )

  	for(i in 1:nrow(data)){
	    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time_censor ){
	      data[ i , "time" ] = time_censor
	      data[ i , "surv" ] = 0
    	}
  	}

	if( length( data$variable[ data$variable == 1 ] )>=5 & length( data$variable[ data$variable == 0 ] ) >= 5 ){

		cox = coxph( formula= Surv( time , surv ) ~ variable , data=data )
		mean <- round( coef( cox )[1] , 2 )
		se = round( summary(cox)$coefficients[3] , 2 )
		low <- round( confint( cox , level=.95 )[ 1 , 1 ] , 2 )
		up <- round( confint( cox , level=.95 )[ 1 , 2 ] , 2 )
		pval <- summary(cox)$coef[1,5]

		pdf( paste( dir , "/KMplot_", title , ".pdf",sep=""),height=5,width=5.5,bg="transparent", onefile=FALSE)
			km.coxph.plot( formula.s=Surv( time, surv ) ~ variable ,data.s=data, x.label="Time (Months)", y.label=ylab, main.title=paste( title , "\n(cutoff=" , round( cutoff , 2 ) , ")" , sep="" ) ,sub.title="", 
				leg.text=c( "Low" , "High"), 
				leg.pos="topright", .col=c( "black","#e53935"),  show.n.risk=TRUE, n.risk.step=6, n.risk.cex=0.85, ylim=c(0,1), leg.inset=0,.lwd=3 , verbose=FALSE )
		dev.off()

	} else{
		mean <- NA
		se = NA
		low <- NA
		up <- NA
		pval <- NA	
	}	

	c( mean , se , low , up , pval )
}

###########################################################################################
###########################################################################################

Get_HR_tertile = function( surv , time , time_censor , variable , title , ylab , dir ){

	data = data.frame( surv=surv , time=time , variable=variable )
	qt = quantile( as.numeric(as.character(data$variable)) , na.rm=TRUE , probs= c( .33 , .66 ) )
	data$time = as.numeric(as.character(data$time))
	data$variable = ifelse( as.numeric(as.character(data$variable)) <= qt[1] , 0 , 
						ifelse( as.numeric(as.character(data$variable)) <= qt[2] , 1 , 
						ifelse( as.numeric(as.character(data$variable)) > qt[2] , 2 , NA ) ) ) 

  	for(i in 1:nrow(data)){
	    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time_censor ){
	      data[ i , "time" ] = time_censor
	      data[ i , "surv" ] = 0
    	}
  	}

	if( length( data$variable[ data$variable == 1 ] )>=5 & length( data$variable[ data$variable == 0 ] ) >= 5 ){

		cox = coxph( formula= Surv( time , surv ) ~ variable , data=data )
		mean <- round( coef( cox )[1] , 2 )
		se = round( summary(cox)$coefficients[3] , 2 )
		low <- round( confint( cox , level=.95 )[ 1 , 1 ] , 2 )
		up <- round( confint( cox , level=.95 )[ 1 , 2 ] , 2 )
		pval <- summary(cox)$coef[1,5]

		pdf( paste( dir , "/KMplot_", title , ".pdf",sep=""),height=5,width=5.5,bg="transparent", onefile=FALSE)
			km.coxph.plot( formula.s=Surv( time, surv ) ~ variable ,data.s=data, x.label="Time (Months)", y.label=ylab, main.title=paste( title , "\n(Tertile)" , sep="" ) ,sub.title="", 
				leg.text=c( "Low" , "Int" , "High"), 
				leg.pos="topright", .col=c( "black" , "#42a5f5" , "#e53935" ),  show.n.risk=TRUE, n.risk.step=6, n.risk.cex=0.85, ylim=c(0,1), leg.inset=0,.lwd=3 , verbose=FALSE )
		dev.off()

	} else{
		mean <- NA
		se = NA
		low <- NA
		up <- NA
		pval <- NA	
	}	

	c( mean , se , low , up , pval )
}

###########################################################################################
###########################################################################################

Get_HR_quartile = function( surv , time , time_censor , variable , title , ylab , dir ){

	data = data.frame( surv=surv , time=time , variable=variable )
	qt = quantile( as.numeric(as.character(data$variable)) , na.rm=TRUE , probs= c( .25 , .5 , .75 ) )
	data$time = as.numeric(as.character(data$time))
	data$variable = ifelse( as.numeric(as.character(data$variable)) <= qt[1] , 0 , 
						ifelse( as.numeric(as.character(data$variable)) <= qt[2] , 1 , 
						ifelse( as.numeric(as.character(data$variable)) <= qt[3] , 2 , 
						ifelse( as.numeric(as.character(data$variable)) > qt[3] , 3 , NA ) ) ) )

  	for(i in 1:nrow(data)){
	    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time_censor ){
	      data[ i , "time" ] = time_censor
	      data[ i , "surv" ] = 0
    	}
  	}

	if( length( data$variable[ data$variable == 1 ] )>=5 & length( data$variable[ data$variable == 0 ] ) >= 5 ){

		cox = coxph( formula= Surv( time , surv ) ~ variable , data=data )
		mean <- round( coef( cox )[1] , 2 )
		se = round( summary(cox)$coefficients[3] , 2 )
		low <- round( confint( cox , level=.95 )[ 1 , 1 ] , 2 )
		up <- round( confint( cox , level=.95 )[ 1 , 2 ] , 2 )
		pval <- summary(cox)$coef[1,5]

		pdf( paste( dir , "/KMplot_", title , ".pdf",sep=""),height=6,width=6.5,bg="transparent", onefile=FALSE)
			km.coxph.plot( formula.s=Surv( time, surv ) ~ variable ,data.s=data, x.label="Time (Months)", y.label=ylab, main.title=paste( title , "\n(Quartile)" , sep="" ) ,sub.title="", 
				leg.text=c( "1st (Lowest)" , "2nd" , "3rd" , "4th(Highest)"), 
				leg.pos="topright", .col=c( "black" , "#42a5f5" , "#66bb6a" , "#e53935" ),  show.n.risk=TRUE, n.risk.step=6, n.risk.cex=0.85, ylim=c(0,1), leg.inset=0,.lwd=3 , verbose=FALSE )
		dev.off()

	} else{
		mean <- NA
		se = NA
		low <- NA
		up <- NA
		pval <- NA	
	}	

	c( mean , se , low , up , pval )
}


###########################################################################################
###########################################################################################

Get_HR_continous = function( surv , time , time_censor , variable ){
	data = data.frame( surv=surv , time=time , variable=variable )
	data$time = as.numeric(as.character(data$time))
	data$variable = as.numeric( as.character(data$variable) )

  	for(i in 1:nrow(data)){
	    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time_censor ){
	      data[ i , "time" ] = time_censor
	      data[ i , "surv" ] = 0
    	}
  	}

	cox = coxph( formula= Surv( time , surv ) ~ variable , data=data )
	mean <- round( coef( cox )[1] , 2 )
	low <- round( confint( cox , level=.95 )[ 1 , 1 ] , 2 )
	up <- round( confint( cox , level=.95 )[ 1 , 2 ] , 2 )
	pval <- summary(cox)$coef[1,5]

	c( mean , round( summary(cox)$coefficients[3] , 2 ) , low , up , pval )
}
