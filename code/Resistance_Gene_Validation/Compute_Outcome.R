##########################################################################################################################################################
##########################################################################################################################################################

library(survcomp)
library(genefu)
library(GSVA)
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


Get_Gene <- function( exprID , geneID ){

	source('meta/Get_HR.R')
	source('meta/Get_DI.R')

	load( paste( "../data/validation_cohort/" , exprID , ".RData" , sep="" ) )

	cox_os = cox_pfs = cox_dicho_os = cox_dicho_pfs = cox_tertile_os = cox_tertile_pfs = NULL
	di_os = di_pfs = NULL
	log_response = log_dicho_response = log_tertile_response = NULL
	auc = NULL

	remove <- rem(expr)
	if( length(remove) ){
		expr <- expr[-remove,]
	}
	data = expr

	primary = unique( clin$primary )
	primary = ifelse( length( primary ) > 1 , "Pan-Cancer" , primary )

	if( geneID %in% rownames(data) & ncol(data) > 20 ){
		geneSig <- as.numeric( scale( as.numeric( as.character( data[ geneID , ] ) ) ) )
		names( geneSig) = colnames( data )

		if( length( geneSig[ !is.na( clin$os ) ] ) >= 20 ){
			## Compute the association ADO signature with OS (36 month cutoff) using Cox Regression Model
			hr = Get_HR_continous( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= geneSig )
			cox_os = rbind( cox_os , 
					c( exprID , primary , "TPM" , length( geneSig[ !is.na( clin$os ) ] ) ,  hr ) )

			## Compute the association ADO signature with OS (36 month cutoff) using Concordence Index
			cci = Get_DI_continous( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= geneSig )
			di_os = rbind( di_os , 
					c( exprID , primary , "TPM" , length( geneSig[ !is.na( clin$os ) ] ) ,  cci ) )

			## Compute the association ADO signature (HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
			hr = Get_HR_dicho( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= geneSig , cutoff= median( geneSig , na.rm=TRUE ) ,
									 title = paste( exprID , 36 , "OS_Median" , sep="_" ) , ylab="Overall Survival" , dir= paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , exprID , sep="" ) )
			cox_dicho_os = rbind( cox_dicho_os , 
					c( exprID , primary , "TPM" , length( geneSig[ !is.na( clin$os ) ] ) ,  hr ) )

			## Compute the association ADO signature (HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
			hr = Get_HR_dicho( surv=clin$os , time=clin$t.os ,
									 time_censor=36 , variable= geneSig , cutoff= quantile(geneSig, probs=.66) ,
									 title = paste( exprID , 36 , "OS_Tertile" , sep="_" ) , ylab="Overall Survival" , dir= paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , exprID , sep="" ) )
			cox_tertile_os = rbind( cox_tertile_os , 
					c( exprID , primary , "TPM" , length( geneSig[ !is.na( clin$os ) ] ) ,  hr ) )
		}

		if( length( geneSig[ !is.na( clin$pfs ) ]  ) >= 20 ){
			## Compute the association of ADO signature with PFS using Cox Regression Model
			hr = Get_HR_continous( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= geneSig )
			cox_pfs = rbind( cox_pfs , 
					c( exprID , primary , "TPM" , length( geneSig[ !is.na( clin$pfs ) ] ) ,  hr ) )

			## Compute the association of ADO signature with PFS using Concordence Index
			cci = Get_DI_continous( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= geneSig )
			di_pfs = rbind( di_pfs , 
					c( exprID , primary , "TPM" , length( geneSig[ !is.na( clin$pfs ) ] ) ,  cci ) )
		
			
			## Compute the association of ADO signature (HIGH vs LOW) with PFS using Cox Regression Model
			hr = Get_HR_dicho( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= geneSig , cutoff= median( geneSig , na.rm=TRUE ) ,
									 title = paste( exprID , 36 , "PFS_Median" , sep="_" ) , ylab="Progression-Free Survival" , dir= paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , exprID , sep="" ) )
			cox_dicho_pfs = rbind( cox_dicho_pfs , 
					c( exprID , primary , "TPM" , length( geneSig[ !is.na( clin$pfs ) ] ) ,  hr ) )
			
			## Compute the association of ADO signature (HIGH vs LOW) with PFS using Cox Regression Model
			hr = Get_HR_dicho( surv=clin$pfs , time=clin$t.pfs ,
									 time_censor=24 , variable= geneSig , cutoff= quantile(geneSig, probs=.66) ,
									 title = paste( exprID , 36 , "PFS_Tertile" , sep="_" ) , ylab="Progression-Free Survival" , dir= paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , exprID , sep="" ) )
			cox_tertile_pfs = rbind( cox_tertile_pfs , 
					c( exprID , primary , "TPM" , length( geneSig[ !is.na( clin$pfs ) ] ) ,  hr ) )
		}

		
		if( length( geneSig[ !is.na( clin$response ) ]  ) >= 20 ){

			## Association with Response (Continous)
			x = ifelse( clin$response %in% "R" , 0 , ifelse( clin$response %in% "NR" , 1 , NA ) )

			fit = glm( x ~ geneSig , family=binomial( link="logit" ) )

			log_response = rbind( log_response , c( exprID , primary , 
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

			pdf( paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , exprID , "/" , geneID , "_ICB_Response.pdf" , sep="" ) , height=3.5,width=3,bg="transparent")

		  		xLabels <- paste( c( "R" , "NR" ) , "\n(N=" , table(d$response) , ")" , sep = "" )
			    yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=.25 ) 
			    boxplot( sig ~ response , data=d , ylab= paste( geneID , "expression" ) , xlab="ICB Response" , main="" , 
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

			basicplot <- ggplot(d, aes(d = response, m = sig )) + geom_roc(n.cuts = 100, labels = FALSE)
			ROCPlot = basicplot + 
			  style_roc() +
			  theme(axis.text = element_text(colour = "black")) +
			  annotate("text", x = .75, y = .25, 
			           label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2))) +
			  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))

			auc = calc_auc(basicplot)$AUC

			ggsave(filename= paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , exprID , "/" , geneID , "_ICB_Response_ROC.pdf" , sep="" ) , plot=ROCPlot, device="pdf", height=4 , width=4 )


			################################################################
			################################################################
			## Association with Response (HIGH vs Low)
			m = median( geneSig , na.rm=TRUE )
			y = ifelse( geneSig >= m , 1 ,0 )
			x = ifelse( clin$response %in% "R" , 0 , ifelse( clin$response %in% "NR" , 1 , NA ) )

			if( length(geneSig[ !is.na( clin$response ) ] )>= 20 ){
				
				fit = glm( x ~ y , family=binomial( link="logit" ) )

				log_dicho_response = rbind( log_dicho_response , c( exprID , primary , 
							"TPM" , 
							length( geneSig[ !is.na( clin$response ) ]  ) , 
							round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
							round( confint(fit)[ 2 , ] , 2 ) , 
							summary(fit)$coefficients[ 2 , 4 ] ) )
			}

			## Association with Response (HIGH vs Low)
			m = quantile(geneSig, probs=.66)
			y = ifelse( geneSig >= m , 1 ,0 )
			x = ifelse( clin$response %in% "R" , 0 , ifelse( clin$response %in% "NR" , 1 , NA ) )

			if( length(geneSig[ !is.na( clin$response ) ] )>= 20 ){
				
				fit = glm( x ~ y , family=binomial( link="logit" ) )

				log_tertile_response = rbind( log_tertile_response , c( exprID , primary , 
							"TPM" , 
							length( geneSig[ !is.na( clin$response ) ]  ) , 
							round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
							round( confint(fit)[ 2 , ] , 2 ) , 
							summary(fit)$coefficients[ 2 , 4 ] ) )
			}
		}
	}

	colnames = c( "study" , "Primary" , "Sequencing" , "N" , "HR" , "SE" , "95di_low" , "95di_high"  , "Pval" )
	cox_os = if_NULL( x= cox_os , colnames= colnames , study= exprID )
	cox_pfs = if_NULL( x= cox_pfs , colnames= colnames , study= exprID )
	cox_dicho_os = if_NULL( x= cox_dicho_os , colnames= colnames , study= exprID )
	cox_dicho_pfs = if_NULL( x= cox_dicho_pfs , colnames= colnames , study= exprID )
	cox_tertile_os = if_NULL( x= cox_tertile_os , colnames= colnames , study= exprID )
	cox_tertile_pfs = if_NULL( x= cox_tertile_pfs , colnames= colnames , study= exprID )

	di_os = if_NULL( x= di_os , colnames= colnames , study= exprID )
	di_pfs = if_NULL( x= di_pfs , colnames= colnames , study= exprID )

	log_response = if_NULL( x= log_response , colnames= colnames , study= exprID )
	log_dicho_response = if_NULL( x= log_dicho_response , colnames= colnames , study= exprID )
	log_tertile_response = if_NULL( x= log_tertile_response , colnames= colnames , study= exprID )

	colnames(cox_os) = colnames(cox_pfs) = colnames(cox_dicho_os) = colnames(cox_dicho_pfs) = colnames(cox_tertile_os) = colnames(cox_tertile_pfs) = c( "study" , "Primary" , "Sequencing" , "N" , "HR" , "SE" , "95di_low" , "95di_high"  , "Pval" )
	colnames(di_os) = colnames(di_pfs) = c( "study" , "Primary" , "Sequencing" , "N" , "DI" , "SE" , "95di_low" , "95di_high"  , "Pval" )
	colnames(log_response) = colnames(log_dicho_response) = colnames(log_tertile_response) = c( "study" , "Primary" , "Sequencing" , "N", "coef" , "SE" , "95di_low" , "95di_high" , "Pval" )


	cox_os = as.data.frame( cox_os )
	di_os = as.data.frame( di_os )

	cox_pfs = as.data.frame( cox_pfs )
	di_pfs = as.data.frame( di_pfs )

	cox_dicho_os = as.data.frame( cox_dicho_os )
	cox_dicho_pfs = as.data.frame( cox_dicho_pfs )
	cox_tertile_os = as.data.frame( cox_tertile_os )
	cox_tertile_pfs = as.data.frame( cox_tertile_pfs )

	log_response = as.data.frame( log_response )
	log_dicho_response = as.data.frame( log_dicho_response )
	log_tertile_response = as.data.frame( log_tertile_response )

	save( log_response , log_dicho_response , log_tertile_response , file = paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , exprID , "/" , geneID , "_LogReg_result.RData" , sep="" ) )
	save( cox_os , cox_pfs , cox_dicho_os , cox_dicho_pfs , cox_tertile_os , cox_tertile_pfs , file = paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , exprID , "/" , geneID , "_COX_result.RData" , sep="" ) ) 
	save(  di_os , di_pfs , file = paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , exprID , "/" , geneID , "_DI_result.RData" , sep="" ) ) 
	save(  auc , file = paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , exprID , "/" , geneID , "_AUC.RData" , sep="" ) ) 
}