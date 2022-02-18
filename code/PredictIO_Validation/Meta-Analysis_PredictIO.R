#################################################################
#################################################################
source("meta/Meta-Analysis.R")

library(meta)
library(metafor)
library(genefu)
library(forestplot)
library(dmetar)

#################################################################
#################################################################

get_Meta_HR = function( data , signature , prefix , height , width ){

	data = as.data.frame( data )

	data$study <- as.character( data$study )
	data$study[ grep( "_" , data$study ) ] <- sapply( data$study[ grep( "_" , data$study ) ] , function(x){ unlist( strsplit( x , "_" , fixed=TRUE ) )[2] } )
	data$Sequencing <- as.character( data$Sequencing )
	data$Primary <- as.character( data$Primary )
	data$N <- as.numeric(as.character( data$N ))
	data$HR <- as.numeric(as.character( data$HR ))
	data$SE <- as.numeric(as.character( data$SE ))
	data$Pval <- as.numeric(as.character( data$Pval )) 

	data = data[ !is.na( data$study ) , ]
	data = data[ order( data$HR ) , ]

	if( nrow(data) >= 3 ){

		data$study <- paste( data$study , ", " , data$Primary , ", n = " , data$N , sep= "" ) 

		meta <- metagen( TE = HR,
	                  seTE = SE,
	                  data = data,
	                  studlab = study,
	                  comb.fixed = FALSE,
	                  comb.random = TRUE)

		######################################################################
		######################################################################

		pdf( prefix , height= height , width= width , bg="transparent" )

			m <- c( min( c( data$HR , -1 ) , na.rm=TRUE) - .5 , ( max( c( 1 , abs(data$HR) ) , na.rm=TRUE) ) + .5 )

			forest( meta , 
	            leftcols = c("studlab", "effect.ci" , "Pval" ),
				leftlabs= c( "Study" , "logHR [95%CI]" , "P-value" ) , 
	   			xlab = "logHR estimate",
				digits.se = 2 ,
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
		       	col.by = "#1565c0",
			    addrow.subgroups=TRUE 
		    )
		dev.off()
	}
}

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
	                  studlab = study,
	                  comb.fixed = FALSE,
	                  comb.random = TRUE)

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

########################################################################################################################
########################################################################################################################


##########################################
##########################################
## OS

get_Cox_OS = function( cohort , IO_MetaScoreID ){
	out = NULL
	for(i in 1:length( cohort ) ){
		cox_RData= paste( "../results/PredictIO_Validation/" , cohort[i] , "/" , IO_MetaScoreID , "/" , IO_MetaScoreID , "_COX_result.RData"  , sep="")
		load( cox_RData )

		out = rbind( out , cox_os)
	}
	out
}

get_Cox_Dicho_OS = function( cohort , IO_MetaScoreID ){
	out = NULL
	for(i in 1:length( cohort ) ){
		cox_RData= paste( "../results/PredictIO_Validation/" , cohort[i] , "/" , IO_MetaScoreID , "/" , IO_MetaScoreID , "_COX_result.RData"  , sep="")
		load( cox_RData )

		out = rbind( out , cox_dicho_os)
	}
	out
}

##########################################
##########################################
## PFS

get_Cox_PFS = function( cohort , IO_MetaScoreID ){
	out = NULL
	for(i in 1:length( cohort ) ){
		cox_RData= paste( "../results/PredictIO_Validation/" , cohort[i] , "/" , IO_MetaScoreID , "/" , IO_MetaScoreID , "_COX_result.RData"  , sep="")
		load( cox_RData )

		out = rbind( out , cox_pfs)
	}
	out
}

get_Cox_Dicho_PFS = function( cohort , IO_MetaScoreID ){
	out = NULL
	for(i in 1:length( cohort ) ){
		cox_RData= paste( "../results/PredictIO_Validation/" , cohort[i] , "/" , IO_MetaScoreID , "/" , IO_MetaScoreID , "_COX_result.RData"  , sep="")
		load( cox_RData )

		out = rbind( out , cox_dicho_pfs)
	}
	out
}

##########################################
##########################################
## Response

get_LogReg_Response = function( cohort , IO_MetaScoreID ){
	out = NULL
	for(i in 1:length( cohort ) ){
		log_RData= paste( "../results/PredictIO_Validation/" , cohort[i] , "/" , IO_MetaScoreID , "/" , IO_MetaScoreID , "_LogReg_result.RData"  , sep="")
		load( log_RData )

		out = rbind( out , log_response)
	}
	out
}

get_LogReg_Dicho_Response = function( cohort , IO_MetaScoreID ){
	out = NULL
	for(i in 1:length( cohort ) ){
		log_RData= paste( "../results/PredictIO_Validation/" , cohort[i] , "/" , IO_MetaScoreID , "/" , IO_MetaScoreID , "_LogReg_result.RData"  , sep="")
		load( log_RData )

		out = rbind( out , log_dicho_response)
	}
	out
}
	
########################################################################################################################
########################################################################################################################

Get_Meta_Analysis_PredictIO <- function( cohort , threshold ){

	IO_MetaScoreID = c( "IO_Resistance" , "IO_Sensitive" , "PredictIO" )

	for(k in 1:length(IO_MetaScoreID) ){

		###########################################################
		###########################################################
		## Meta-analysis of the COX models (OS) Continous
		cox_os = get_Cox_OS( cohort=cohort , IO_MetaScoreID= IO_MetaScoreID[k] )
		get_Meta_HR( data=cox_os , signature= IO_MetaScoreID[k] ,
					prefix=paste( "../results/PredictIO_Validation/meta_analysis/OS/Meta_analysis_" , IO_MetaScoreID[k] , "_OS_Continuous.pdf" , sep="") ,
					height= 6 ,
					width= 7 )

		## Meta-analysis of the COX models (OS) High vs Low
		cox_dicho_os = get_Cox_Dicho_OS( cohort=cohort , IO_MetaScoreID= IO_MetaScoreID[k] )
		get_Meta_HR( data=cox_dicho_os , signature= IO_MetaScoreID[k] , 
					prefix=paste( "../results/PredictIO_Validation/meta_analysis/OS/Meta_analysis_" , IO_MetaScoreID[k] , "_OS_Dicho.pdf" , sep="") ,
					height= 6 ,
					width= 7 )

		#############################################################
		#############################################################
		## Meta-analysis of the COX models (PFS) Continous
		cox_pfs = get_Cox_PFS( cohort=cohort , IO_MetaScoreID= IO_MetaScoreID[k] )
		get_Meta_HR( data=cox_pfs , signature= IO_MetaScoreID[k] , 
					prefix= paste( "../results/PredictIO_Validation/meta_analysis/PFS/Meta_analysis_" , IO_MetaScoreID[k] , "_PFS_Continuous.pdf" , sep="") ,
					height= 6 ,
					width= 7 )

		## Meta-analysis of the COX models (PFS) High vs Low
		cox_dicho_pfs = get_Cox_Dicho_PFS( cohort=cohort , IO_MetaScoreID= IO_MetaScoreID[k] )
		get_Meta_HR( data= cox_dicho_pfs , signature= IO_MetaScoreID[k] , 
					prefix= paste( "../results/PredictIO_Validation/meta_analysis/PFS/Meta_analysis_" , IO_MetaScoreID[k] , "_PFS_Dicho.pdf" , sep="") ,
					height= 6 ,
					width= 7 )

		#####################################################
		#####################################################
		## Meta-analysis of the Log Regression models (Response) Continous
		log_response = get_LogReg_Response( cohort=cohort , IO_MetaScoreID= IO_MetaScoreID[k] )
		get_Meta_OR( data=log_response , signature= IO_MetaScoreID[k] , 
					prefix=paste( "../results/PredictIO_Validation/meta_analysis/Response/Meta_analysis_" , IO_MetaScoreID[k] , "_Response_Continuous.pdf" , sep="") ,
					height= 6 ,
					width= 7 )

		## Meta-analysis of the Log Regression models (Response) High vs Low
		log_dicho_response = get_LogReg_Dicho_Response( cohort=cohort , IO_MetaScoreID= IO_MetaScoreID[k] )
		get_Meta_OR( data=log_dicho_response , signature= IO_MetaScoreID[k] , 
					prefix= paste( "../results/PredictIO_Validation/meta_analysis/Response/Meta_analysis_" , IO_MetaScoreID[k] , "_Response_Dicho.pdf" , sep="") ,
					height= 6 ,
					width= 7 )



	}
}

