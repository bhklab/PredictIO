
source("meta/Meta-Analysis.R")

library(meta)
library(metafor)
library(genefu)
library(forestplot)
library(dmetar)

########################################################################################################################
########################################################################################################################

get_Meta_OR = function( data , signature , prefix , height , width ){
	
	data = as.data.frame( data )

	data$study <- as.character( data$study )
	data$study[ grep( "_" , data$study ) ] <- sapply( data$study[ grep( "_" , data$study ) ] , function(x){ unlist( strsplit( x , "_" , fixed=TRUE ) )[2] } )
	data$Sequencing <- as.character( data$Sequencing )
	data$Primary <- as.character( data$Primary )
	data$N <- as.numeric(as.character( data$N ))
	data$coef <- as.numeric(as.character( data$coef ))
	data$SE <- as.numeric(as.character( data$SE ))
	data$Pval <- as.numeric(as.character( data$Pval )) 

	data = data[ !is.na( data$study ) , ]
	data = data[ order( data$coef ) , ]

	if( nrow(data) >= 3 ){

		data$study <- paste( data$study , ", " , data$Primary , ", n = " , data$N , sep= "" ) 

		meta <- metagen( TE = coef,
	                  seTE = SE,
	                  data = data,
	                  studlab = study ,
	                  fixed = FALSE ,
	                  random = TRUE ,
	                  control = list( maxiter = 10000 , stepadj=0.5 ) )

		######################################################################
		######################################################################

		pdf( prefix , height= height, width= width , bg="transparent" )

			m <- c( min( c( data$coef , -1 ) , na.rm=TRUE) - .5 , ( max( c( 1 , abs(data$coef) ) , na.rm=TRUE) ) + .5 )

			forest( meta , 
	            leftcols = c("studlab", "effect.ci" , "Pval"),
				leftlabs= c( "Study" , "logOR [95%CI]" , "P-value" ) , 
	   			xlab = "Estimated logOR",
				digits.se = 2,
	   			colgap.forest=unit(10, "mm") ,
		      	plotwidth = unit( 30 , "mm") ,
		       	pooled.totals = TRUE,
		       	smlab = " ",
		       	comb.random =TRUE,
		       	comb.fixed = FALSE,
		       	text.fixed.w = FALSE,
			    layout = "JAMA",
			    print.I2.ci = TRUE,
			    print.Q = FALSE,
			    print.pval.Q = TRUE,
			    print.I2 = TRUE,
			    print.tau2 = FALSE,
			    resid.hetstat = FALSE,
		       	test.overall.random = TRUE,
		       	test.overall.fixed = FALSE,
		       	xlim = m , 
		       	col.square= "black" ,  
		       	col.study= "black" ,  
		       	col.square.lines = "black" ,
		       	col.diamond.random  = "#1565c0"  ,
		       	col.diamond.lines.random  ="#1565c0" ,
		       	col.by = "#1565c0" ,
			    addrow.subgroups=TRUE
			)
		dev.off()
	}
}

##########################################
##########################################
## Response
get_LogReg_Response = function( cohort , geneID ){
	out = NULL
	for(i in 1:length( cohort ) ){
		log_RData= paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , cohort[i] , "/" , geneID , "_LogReg_result.RData"  , sep="")
		load( log_RData )

		out = rbind( out , log_response)
	}
	out
}
get_LogReg_Dicho_Response = function( cohort , geneID ){
	out = NULL
	for(i in 1:length( cohort ) ){
		log_RData= paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , cohort[i] , "/" , geneID , "_LogReg_result.RData"  , sep="")
		load( log_RData )

		out = rbind( out , log_dicho_response)
	}
	out
}

########################################################################################################################
########################################################################################################################

Get_Meta_Analysis <- function( cohort , geneID ){

	########################################################################################################################
	########################################################################################################################
	## Meta-analysis of the Log Regression models (Response) Continous
	log_response = get_LogReg_Response( cohort=cohort , geneID= geneID )
	get_Meta_OR( data=log_response , signature= geneID , 
				prefix=paste( "../results/Resistance_Gene_Validation/" , geneID , "/meta_analysis/Meta_analysis_" , geneID , "_Response_Continuous.pdf" , sep="") ,
				height= 6 ,
				width= 7 )

	## Meta-analysis of the Log Regression models (Response) High vs Low
	log_dicho_response = get_LogReg_Dicho_Response( cohort=cohort , geneID= geneID )
	get_Meta_OR( data=log_dicho_response , signature= IO_MetaScoreID[k] , 
				prefix= paste( "../results/Resistance_Gene_Validation/" , geneID , "/meta_analysis/Meta_analysis_" , geneID , "_Response_Dicho.pdf" , sep="") ,
				height= 6 ,
				width= 7 )

}

