##########################################################################################################################################################
##########################################################################################################################################################

require(survcomp)
require(genefu)
require(GSVA)

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
		unlist( file )
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

get_Scale = function( x ){
	rid = rownames(x)
	cid = colnames(x)
	out = t( apply( x , 1 , scale ) )
	rownames(out) = rid
	colnames(out) = cid
	out
}

if_NULL = function( x , colnames , study ){
	if( is.null( x ) ){
		x = t( as.data.frame( rep( NA , length( colnames ) ) ) )	
		rownames(x) = study	
	}

	x
}

Mean_Signature = function( signature_name , cohort , sigType ){

	source( 'meta/Get_HR.R' )
	source( 'meta/Get_DI.R' )

	get_Directory( dir=signature_name , cohort= cohort )

	sig = read.csv( file= paste( "../data/signatures/" , signature_name , ".csv" , sep="" ) , sep=';' , header=TRUE , stringsAsFactor=FALSE)
	rownames(sig) = sig$gene
	sig$coef <- as.numeric( as.character( sig$coef ) )
	sig$gene <- as.character( sig$gene )

	load( paste( "../data/validation_cohort/" , cohort , ".RData" , sep="" ) )

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
				
	if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene , ] ) ) / nrow( sig ) > 0.8 & ncol(data) >= 10 ){

		gene = intersect( rownames(data) , sig$gene) 
		s = sig[ gene , ]
		geneSig <- apply( get_Scale( x=data[ gene , ] ) , 2 , function(x) ( sum( ( x * s$coef ) ) /  nrow( s ) ) )
		names( geneSig ) = colnames(data)
		
		if( length( geneSig[ !is.na( clin$os ) ] ) >= 10 ){
			## Compute the association Weighted-mean signature with OS (36 month cutoff) using Cox Regression Model
			hr = Get_HR_continous( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= geneSig )
			cox_os = rbind( cox_os , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$os ) ] ) ,  hr ) )

			## Compute the association Weighted-mean signature with OS (36 month cutoff) using Concordence Index
			cci = Get_DI_continous( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= geneSig )
			di_os = rbind( di_os , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$os ) ] ) ,  cci ) )

		
			## Compute the association Weighted-mean signature (HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
			hr = Get_HR_dicho( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= geneSig , cutoff= median( geneSig , na.rm=TRUE ) ,
									 title = paste( cohort , primary , 36 , "OS_Overall" , sep="_" ) , ylab="Overall Survival" , dir= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/KMPlot/OS" , sep="" ) )
			cox_dicho_os = rbind( cox_dicho_os , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$os ) ] ) ,  hr ) )



		}

		if( length( geneSig[ !is.na( clin$pfs ) ]  ) >= 10 ){
			## Compute the association of Weighted-mean signature with PFS using Cox Regression Model
			hr = Get_HR_continous( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= geneSig )
			cox_pfs = rbind( cox_pfs , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$pfs ) ] ) ,  hr ) )

			## Compute the association of Weighted-mean signature with PFS using Concordence Index
			cci = Get_DI_continous( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= geneSig )
			di_pfs = rbind( di_pfs , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$pfs ) ] ) ,  cci ) )
		
			
			## Compute the association of Weighted-mean signature (HIGH vs LOW) with PFS using Cox Regression Model
			hr = Get_HR_dicho( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= geneSig , cutoff= median( geneSig , na.rm=TRUE ) ,
									 title = paste( cohort , primary , 36 , "PFS_Overall" , sep="_" ) , ylab="Progression-Free Survival" , dir= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/KMPlot/PFS" , sep="" ) )
			cox_dicho_pfs = rbind( cox_dicho_pfs , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$pfs ) ] ) ,  hr ) )


		}

		
		if( length( geneSig[ !is.na( clin$response ) ]  ) >= 10 ){

			## Association with Response (Continous)
			x = ifelse( clin$response %in% "R" , 0 , ifelse( clin$response %in% "NR" , 1 , NA ) )

			fit = glm( x ~ geneSig , family=binomial( link="logit" ) )

			log_response = rbind( log_response , c( cohort , primary , 
						"TPM" , 
						length( geneSig[ !is.na( clin$response ) ]  ) , 
						round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
						round( confint(fit)[ 2 , ] , 2 ) , 
						summary(fit)$coefficients[ 2 , 4 ] ) )


			################################################################
			################################################################

	  		d = as.data.frame( cbind( x , geneSig ) )
	  		colnames(d) = c( "response" , "sig" )
	  		d = d[ !is.na( d$response ) , ]
	  		d$response = as.numeric( as.character( d$response ) )
	  		d$sig = as.numeric( as.character( d$sig ) )

			pdf( paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/ICB_Response.pdf" , sep="" ) , height=3.5,width=3,bg="transparent")

		  		xLabels <- paste( c( "R" , "NR" ) , "\n(N=" , table(d$response) , ")" , sep = "" )
			    yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=.25 ) 
			    boxplot( sig ~ response , data=d , ylab= signature_name , xlab="ICB Response" , main="" , 
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

	  		d = as.data.frame( cbind( x , geneSig ) )
	  		colnames(d) = c( "response" , "sig" )
	  		d = d[ !is.na( d$response ) , ]
	  		d$response = as.numeric( as.character( d$response ) )
	  		d$sig = as.numeric( as.character( d$sig ) )

			basicplot <- ggplot(d, aes(d = response, m = sig )) + geom_roc(n.cuts = 100, labels = FALSE)
			ROCPlot = basicplot + 
			  style_roc() +
			  theme(axis.text = element_text(colour = "black")) +
			  ggtitle( signature_name ) + 
			  annotate("text", x = .75, y = .25, 
			           label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2))) +
			  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))

			auc = calc_auc(basicplot)$AUC

			ggsave(filename= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/ICB_Response_ROC.pdf" , sep="" ), plot=ROCPlot, device="pdf", height=4 , width=4 )


			################################################################
			################################################################
			## Association with Response (HIGH vs Low)
			m = median( geneSig , na.rm=TRUE )
			y = ifelse( geneSig >= m , 1 ,0 )
			x = ifelse( clin$response %in% "R" , 0 , ifelse( clin$response %in% "NR" , 1 , NA ) )

			if( length(geneSig[ !is.na( clin$response ) ] ) >= 40 ){
				
				fit = glm( x ~ y , family=binomial( link="logit" ) )

				log_dicho_response = rbind( log_dicho_response , c( cohort , primary ,
							"TPM" , 
							length( geneSig[ !is.na( clin$response ) ]  ) , 
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

	save( log_response , log_dicho_response , file= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/" , signature_name , "_LogReg_result.RData" , sep="" ) )
	save( cox_os , cox_pfs , cox_dicho_os , cox_dicho_pfs , file= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/" , signature_name , "_COX_result.RData" , sep="" ) ) 
	save(  di_os , di_pfs , file= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/" , signature_name , "_DI_result.RData" , sep="" ) ) 
	save(  auc , file=paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/" , signature_name , "_AUC.RData" , sep="" ) ) 

}


GSVA_Signature = function( signature_name , cohort , sigType ){

	source('meta/Get_HR.R')
	source('meta/Get_DI.R')

	get_Directory( dir= signature_name , cohort= cohort )

	sig = as.character( read.csv( file= paste( "../data/signatures/" , signature_name , ".csv" , sep="" ) , sep=';' , header=TRUE , stringsAsFactor=FALSE)[,1] )


	load( paste( "../data/validation_cohort/" , cohort , ".RData" , sep="" ) )
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
				
	if( ifelse( is.null( nrow( data[ rownames(data) %in% sig , ]) ) , 1 , nrow( data[ rownames(data) %in% sig , ] ) ) / length( sig ) > 0.8 & ncol(data) >= 10 ){

		geneSig <- as.numeric(gsva( get_Scale( x=data ) , list(sig) , verbose=FALSE ) )
		names( geneSig ) = colnames(data)

		if( length( geneSig[ !is.na( clin$os ) ] ) >= 10 ){
			## Compute the association GSVA signature with OS (36 month cutoff) using Cox Regression Model
			hr = Get_HR_continous( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= geneSig )
			cox_os = rbind( cox_os , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$os ) ] ) ,  hr ) )

			## Compute the association GSVA signature with OS (36 month cutoff) using Concordence Index
			cci = Get_DI_continous( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= geneSig )
			di_os = rbind( di_os , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$os ) ] ) ,  cci ) )

		
			## Compute the association GSVA signature (HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
			hr = Get_HR_dicho( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= geneSig , cutoff= median( geneSig , na.rm=TRUE ) ,
									 title = paste( cohort , primary , 36 , "OS_Overall" , sep="_" ) , ylab="Overall Survival" , dir= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/KMPlot/OS" , sep="" ) )
			cox_dicho_os = rbind( cox_dicho_os , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$os ) ] ) ,  hr ) )



		}

		if( length( geneSig[ !is.na( clin$pfs ) ]  ) >= 10 ){
			## Compute the association of GSVA signature with PFS using Cox Regression Model
			hr = Get_HR_continous( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= geneSig )
			cox_pfs = rbind( cox_pfs , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$pfs ) ] ) ,  hr ) )

			## Compute the association of GSVA signature with PFS using Concordence Index
			cci = Get_DI_continous( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= geneSig )
			di_pfs = rbind( di_pfs , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$pfs ) ] ) ,  cci ) )
		
			
			## Compute the association of GSVA signature (HIGH vs LOW) with PFS using Cox Regression Model
			hr = Get_HR_dicho( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= geneSig , cutoff= median( geneSig , na.rm=TRUE ) ,
									 title = paste( cohort , primary , 36 , "PFS_Overall" , sep="_" ) , ylab="Progression-Free Survival" , dir= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/KMPlot/PFS" , sep="" ) )
			cox_dicho_pfs = rbind( cox_dicho_pfs , 
					c( cohort , primary , "TPM" , length( geneSig[ !is.na( clin$pfs ) ] ) ,  hr ) )


		}

		
		if( length( geneSig[ !is.na( clin$response ) ]  ) >= 10 ){

			## Association with Response (Continous)
			x = ifelse( clin$response %in% "R" , 0 , ifelse( clin$response %in% "NR" , 1 , NA ) )

			fit = glm( x ~ geneSig , family=binomial( link="logit" ) )

			log_response = rbind( log_response , c( cohort , primary , 
						"TPM" , 
						length( geneSig[ !is.na( clin$response ) ]  ) , 
						round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
						round( confint(fit)[ 2 , ] , 2 ) , 
						summary(fit)$coefficients[ 2 , 4 ] ) )


			################################################################
			################################################################

	  		d = as.data.frame( cbind( x , geneSig ) )
	  		colnames(d) = c( "response" , "sig" )
	  		d = d[ !is.na( d$response ) , ]
	  		d$response = as.numeric( as.character( d$response ) )
	  		d$sig = as.numeric( as.character( d$sig ) )

			pdf( paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/ICB_Response.pdf" , sep="" ) , height=3.5,width=3,bg="transparent")

		  		xLabels <- paste( c( "R" , "NR" ) , "\n(N=" , table(d$response) , ")" , sep = "" )
			    yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=.25 ) 
			    boxplot( sig ~ response , data=d , ylab= signature_name , xlab="ICB Response" , main="" , 
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

	  		d = as.data.frame( cbind( x , geneSig ) )
	  		colnames(d) = c( "response" , "sig" )
	  		d = d[ !is.na( d$response ) , ]
	  		d$response = as.numeric( as.character( d$response ) )
	  		d$sig = as.numeric( as.character( d$sig ) )

			basicplot <- ggplot(d, aes(d = response, m = sig )) + geom_roc(n.cuts = 100, labels = FALSE)
			ROCPlot = basicplot + 
			  style_roc() +
			  theme(axis.text = element_text(colour = "black")) +
			  ggtitle( signature_name ) + 
			  annotate("text", x = .75, y = .25, 
			           label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2))) +
			  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))

			auc = calc_auc(basicplot)$AUC

			ggsave(filename= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/ICB_Response_ROC.pdf" , sep="" ), plot=ROCPlot, device="pdf", height=4 , width=4 )

			################################################################
			################################################################
			## Association with Response (HIGH vs Low)
			m = median( geneSig , na.rm=TRUE )
			y = ifelse( geneSig >= m , 1 ,0 )
			x = ifelse( clin$response %in% "R" , 0 , ifelse( clin$response %in% "NR" , 1 , NA ) )

			if( length(geneSig[ !is.na( clin$response ) ] )>= 40 ){
				
				fit = glm( x ~ y , family=binomial( link="logit" ) )

				log_dicho_response = rbind( log_dicho_response , c( cohort , primary ,
							"TPM" , 
							length( geneSig[ !is.na( clin$response ) ] ) , 
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

	save( log_response , log_dicho_response , file= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/" , signature_name , "_LogReg_result.RData" , sep="" ) )
	save( cox_os , cox_pfs , cox_dicho_os , cox_dicho_pfs , file= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/" , signature_name , "_COX_result.RData" , sep="" ) ) 
	save(  di_os , di_pfs , file= paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/" , signature_name , "_DI_result.RData" , sep="" ) ) 
	save(  auc , file=paste( "../results/signature_Validation/" , cohort , "/" , signature_name , "/" , signature_name , "_AUC.RData" , sep="" ) ) 

}