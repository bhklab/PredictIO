##########################################################################################################################################################
##########################################################################################################################################################

library(survcomp)
library(genefu)
library(data.table)
library(Biobase)

library(beeswarm)
library(plotROC)
#Remove rows if count is < zero in 50% of sample
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0) ))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}

get_Directory <- function( dir , cohort ){

	file = paste( "../results/signature_Validation/" , cohort , "/" , dir , sep="" )
	
	if( file.exists( file ) ) {
		unlink( file )
	}

	dir.create( paste( file ) )

	dir.create( paste( file , "/KMPlot" , sep="" ) )
	dir.create( paste( file , "/KMPlot/OS" , sep="" ) )
	dir.create( paste( file , "/KMPlot/PFS" , sep="" ) )

	dir.create( paste( file , "" , sep="" ) )

	dir.create( paste( file , "/Funnel" , sep="" ) )
	dir.create( paste( file , "/Funnel/OS" , sep="" ) )
	dir.create( paste( file , "/Funnel/PFS" , sep="" ) )
	dir.create( paste( file , "/Funnel/Response" , sep="" ) )
	
	dir.create( paste( file , "/Overall" , sep="" ) )
	dir.create( paste( file , "/Overall/OS" , sep="" ) )
	dir.create( paste( file , "/Overall/PFS" , sep="" ) )
	dir.create( paste( file , "/Overall/Response" , sep="" ) )
	
	dir.create( paste( file , "/PerCancer" , sep="" ) )
	dir.create( paste( file , "/PerCancer/OS" , sep="" ) )
	dir.create( paste( file , "/PerCancer/PFS" , sep="" ) )
	dir.create( paste( file , "/PerCancer/Response" , sep="" ) )
	
	dir.create( paste( file , "/PerSequencing" , sep="" ) )
	dir.create( paste( file , "/PerSequencing/OS" , sep="" ) )
	dir.create( paste( file , "/PerSequencing/PFS" , sep="" ) )
	dir.create( paste( file , "/PerSequencing/Response" , sep="" ) )

}

if_NULL = function( x , colnames , study ){
	if( is.null( x ) ){
		x = t( as.data.frame( rep( NA , length( colnames ) ) ) )	
		rownames(x) = study	
	}

	x
}

#################################################
#################################################

Get_TIDE = function( expr , primary ){
	options("encoding" = "UTF-8")
	
	primary = ifelse( primary %in% "Melanoma" , "Melanoma" , 
				ifelse( primary %in% "Lung" , "NSCLC" , "Other" ) )

	expr_output = paste( tempfile() , "_EXPR" , sep="" )
	tide_output = paste( tempfile() , "_TIDE" , sep="" )

	mean_expr = apply( expr , 1 , mean , na.rm=TRUE )

	output =  expr - mean_expr
	write.table( output , file= expr_output , quote=FALSE , sep="\t" , 
				col.names=TRUE , row.names=TRUE , fileEncoding= "UTF-8" )
	
	print( paste( "tidepy " , expr_output , " -o " , tide_output , " -c " , primary , " --ignore_norm" , sep="" ) )
	system( paste( "tidepy " , expr_output , " -o " , tide_output , " -c " , primary , " --ignore_norm" , sep="" ) )

	tide = NULL
	if( file.exists( tide_output ) ) {
		tide = read.table( tide_output , header=TRUE, sep="\t" , stringsAsFactors=FALSE ) 
		names = tide[,1]

		tide = tide$TIDE
		tide = as.numeric( scale( tide ) )

		names( tide ) = names
	}
	tide
}

#################################################
#################################################


Get_TIDE_signature = function( cohort , sigType ){

	source('meta/Get_HR.R')
	source('meta/Get_DI.R')
	get_Directory( dir= "TIDE" , cohort= cohort )


	load( paste( "../data/validation_cohort/", cohort , ".RData" , sep="" ) )
	remove <- rem(expr)
	if( length(remove) ){
		expr <- expr[-remove,]
	}
	data = expr

	primary = unique( clin$primary )
	primary = ifelse( length( primary ) > 1 , "Pan-Cancer" , primary )

	cox_os = cox_pfs = cox_dicho_os = cox_dicho_pfs = NULL

	di_os = di_pfs = NULL
	log_response = log_dicho_response = NULL
	auc = NULL

	TIDE=NULL
	if( ncol(data) ){ 
		TIDE = Get_TIDE( expr=data , primary= primary )
	}
	
	if( ! is.null( length( TIDE ) ) & length( TIDE ) > 20 ){

		if( length( TIDE[ !is.na( clin$os ) ] ) >= 10 ){
			## Compute the association ADO signature with OS (36 month cutoff) using Cox Regression Model
			hr = Get_HR_continous( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= TIDE )
			cox_os = rbind( cox_os , 
					c( cohort , primary , "TPM" , length( TIDE[ !is.na( clin$os ) ] ) ,  hr ) )

			## Compute the association ADO signature with OS (36 month cutoff) using Concordence Index
			cci = Get_DI_continous( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= TIDE )
			di_os = rbind( di_os , 
					c( cohort , primary , "TPM" , length( TIDE[ !is.na( clin$os ) ] ) ,  cci ) )

		
			## Compute the association ADO signature (HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
			hr = Get_HR_dicho( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= TIDE , cutoff= median( TIDE , na.rm=TRUE ) ,
									 title = paste( cohort , primary , 36 , "OS_Overall" , sep="_" ) , ylab="Overall Survival" , dir= paste( "../results/signature_Validation/" , cohort , "/TIDE/KMPlot/OS" , sep="" ) )
			cox_dicho_os = rbind( cox_dicho_os , 
					c( cohort , primary , "TPM" , length( TIDE[ !is.na( clin$os ) ] ) ,  hr ) )

		

		}

		if( length( TIDE[ !is.na( clin$pfs ) ]  ) >= 10 ){
			## Compute the association of ADO signature with PFS using Cox Regression Model
			hr = Get_HR_continous( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= TIDE )
			cox_pfs = rbind( cox_pfs , 
					c( cohort , primary , "TPM" , length( TIDE[ !is.na( clin$pfs ) ] ) ,  hr ) )

			## Compute the association of ADO signature with PFS using Concordence Index
			cci = Get_DI_continous( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= TIDE )
			di_pfs = rbind( di_pfs , 
					c( cohort , primary , "TPM" , length( TIDE[ !is.na( clin$pfs ) ] ) ,  cci ) )
		
			
			## Compute the association of ADO signature (HIGH vs LOW) with PFS using Cox Regression Model
			hr = Get_HR_dicho( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= TIDE , cutoff= median( TIDE , na.rm=TRUE ) ,
									 title = paste( cohort , primary , 36 , "PFS_Overall" , sep="_" ) , ylab="Progression-Free Survival" , dir= paste( "../results/signature_Validation/" , cohort , "/TIDE/KMPlot/PFS" , sep="" ) )
			cox_dicho_pfs = rbind( cox_dicho_pfs , 
					c( cohort , primary , "TPM" , length( TIDE[ !is.na( clin$pfs ) ] ) ,  hr ) )
		}

		
		if( length( TIDE[ !is.na( clin$response ) ]  ) >= 10 ){

			## Association with Response (Continous)
			x = ifelse( clin$response %in% "R" , 0 , ifelse( clin$response %in% "NR" , 1 , NA ) )

			fit = glm( x ~ TIDE , family=binomial( link="logit" ) )

			log_response = rbind( log_response , c( cohort , primary , 
						"TPM" , 
						length( TIDE[ !is.na( clin$response ) ] ) , 
						round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
						round( confint(fit)[ 2 , ] , 2 ) , 
						summary(fit)$coefficients[ 2 , 4 ] ) )



			################################################################
			################################################################

	  		d = as.data.frame( cbind( x , TIDE ) )
	  		colnames(d) = c( "response" , "sig" )
	  		d = d[ !is.na( d$response ) , ]
	  		d$response = as.numeric( as.character( d$response ) )
	  		d$sig = as.numeric( as.character( d$sig ) )

			pdf( paste( "../results/signature_Validation/" , cohort , "/TIDE/ICB_Response.pdf" , sep="" ) , height=3.5,width=3,bg="transparent")

		  		xLabels <- paste( c( "R" , "NR" ) , "\n(N=" , table(d$response) , ")" , sep = "" )
			    yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=.25 ) 
			    boxplot( sig ~ response , data=d , ylab= "TIDE" , xlab="ICB Response" , main="" , 
			    	col=adjustcolor( c( "#00bcd4" , "#f44336" ), alpha.f = .4), 	    	 
			    	boxlty = 1 ,outline=FALSE, axes=FALSE, ylim=c(min( d$sig ) ,  max( d$sig )))

				beeswarm( sig ~ response , data=d , 
				        pch = 19, corral="wrap", cex=.3, 
				        col=adjustcolor( c( "#00bcd4" , "#f44336" ), alpha.f = .6), 
				        add = TRUE)

			    axis(side = 2, at=yLabels, labels= as.character( yLabels ), las= 2,
			             cex.axis=1,tick=1,col="black")

		    	axis(side = 1, at=seq( 1 , length(xLabels) , 1 ) , padj=.5, labels=xLabels, las= 1,
			      cex.axis=1,tick=1,col="black")

				mtext( paste( "LogOR=" , round( summary(fit)$coefficients[ 2 , 1 ] , 2 ) , 
						" (95CI= [" , round( confint(fit)[ 2 , 1 ] , 2 ) , "; " , 
						round( confint(fit)[ 2 , 2 ] , 2 ) , "])\nP" , 
						ifelse( summary(fit)$coefficients[ 2 , 4 ] <= .001 , "≤ 0.001" , paste( "=" , format.pval( summary(fit)$coefficients[ 2 , 4 ] , 1 ) ) ) , sep="" ) , 
						col="#6D6D6D" )
			dev.off()



			################################################################
			################################################################

			if( sigType %in% "sensitive" ){
				x = ifelse( clin$response %in% "NR" , 0 , ifelse( clin$response %in% "R" , 1 , NA ) )
			} else{			
				x = ifelse( clin$response %in% "R" , 0 , ifelse( clin$response %in% "NR" , 1 , NA ) )
			}

	  		d = as.data.frame( cbind( x , TIDE ) )
	  		colnames(d) = c( "response" , "sig" )
	  		d = d[ !is.na( d$response ) , ]
	  		d$response = as.numeric( as.character( d$response ) )
	  		d$sig = as.numeric( as.character( d$sig ) )

			basicplot <- ggplot(d, aes(d = response, m = sig )) + geom_roc(n.cuts = 100, labels = FALSE)
			ROCPlot = basicplot + 
			  style_roc() +
			  theme(axis.text = element_text(colour = "black")) +
			  ggtitle( "TIDE" ) + 
			  annotate("text", x = .75, y = .25, 
			           label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2))) +
			  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))

			auc = calc_auc(basicplot)$AUC

			ggsave(filename= paste( "../results/signature_Validation/" , cohort , "/TIDE/ICB_Response_ROC.pdf" , sep="" ), plot=ROCPlot, device="pdf", height=4 , width=4 )


			################################################################
			################################################################
			## Association with Response (HIGH vs Low)
			m = 1
			y = ifelse( TIDE >= m , 1 ,0 )
			x = ifelse( clin$response %in% "R" , 0 , ifelse( clin$response %in% "NR" , 1 , NA ) )

			if( length(TIDE[ !is.na( clin$response ) ] )>= 40 ){
				
				fit = glm( x ~ y , family=binomial( link="logit" ) )

				log_dicho_response = rbind( log_dicho_response , c( cohort , primary ,
							"TPM" , 
							length( TIDE[ !is.na( clin$response ) ] ) , 
							round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
							round( confint(fit)[ 2 , ] , 2 ) , 
							summary(fit)$coefficients[ 2 , 4 ] ) )
			}


		}
	}	

	colnames = c( "study" , "Primary" , "Sequencing" , "N" , "HR" , "SE" , "95di_low" , "95di_high"  , "Pval" )
	cox_os = if_NULL( x= cox_os , colnames= colnames , study= cohort )
	cox_pfs = if_NULL( x= cox_pfs , colnames= colnames , study= cohort )
	cox_dicho_os = if_NULL( x= cox_dicho_os , colnames= colnames , study= cohort )
	cox_dicho_pfs = if_NULL( x= cox_dicho_pfs , colnames= colnames , study= cohort )

	di_os = if_NULL( x= di_os , colnames= colnames , study= cohort )
	di_pfs = if_NULL( x= di_pfs , colnames= colnames , study= cohort )

	log_response = if_NULL( x= log_response , colnames= colnames , study= cohort )
	log_dicho_response = if_NULL( x= log_dicho_response , colnames= colnames , study= cohort )


	colnames(cox_os) = colnames(cox_pfs) = colnames(cox_dicho_os) = colnames(cox_dicho_pfs) = c( "study" , "Primary" , "Sequencing" , "N" , "HR" , "SE" , "95di_low" , "95di_high"  , "Pval" )
	colnames(di_os) = colnames(di_pfs) = c( "study" , "Primary" , "Sequencing" , "N" , "DI" , "SE" , "95di_low" , "95di_high"  , "Pval" )
	colnames(log_response) = colnames(log_dicho_response) = c( "study" , "Primary" , "Sequencing" , "N", "coef" , "SE" , "95di_low" , "95di_high" , "Pval" )


	cox_os = as.data.frame( cox_os )
	di_os = as.data.frame( di_os )

	cox_pfs = as.data.frame( cox_pfs )
	di_pfs = as.data.frame( di_pfs )

	cox_dicho_os = as.data.frame( cox_dicho_os )


	cox_dicho_pfs = as.data.frame( cox_dicho_pfs )

	log_response = as.data.frame( log_response )
	log_dicho_response = as.data.frame( log_dicho_response )

	save( log_response , log_dicho_response , file= paste( "../results/signature_Validation/" , cohort , "/TIDE/TIDE_LogReg_result.RData" , sep="" ) )
	save( cox_os , cox_pfs , cox_dicho_os , cox_dicho_pfs , file= paste( "../results/signature_Validation/" , cohort , "/TIDE/TIDE_COX_result.RData" , sep="" ) ) 
	save(  di_os , di_pfs , file= paste( "../results/signature_Validation/" , cohort , "/TIDE/TIDE_DI_result.RData" , sep="" ) ) 
	save(  auc , file=paste( "../results/signature_Validation/" , cohort , "/TIDE/TIDE_AUC.RData" , sep="" ) ) 

}