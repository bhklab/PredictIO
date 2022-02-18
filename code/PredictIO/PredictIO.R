##########################################################################################################################################################
##########################################################################################################################################################

library(survcomp)
library(genefu)
library(GSVA)


#################################################
#################################################


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

getIOscore <- function( data , meta_res ){
	
	remove <- rem(data)
	if( length(remove) ){
		data <- data[-remove,]
	}		
	
	sensitive = meta_res[ meta_res$coef < 0 & meta_res$include %in% 1 , ]$gene
	resistance = meta_res[ meta_res$coef > 0 & meta_res$include %in% 1 , ]$gene
	
	IO_resistance = NULL
	if( ifelse( is.null( nrow( data[ rownames(data) %in% resistance , ]) ) , 1 , nrow( data[ rownames(data) %in% resistance , ] ) ) / length( resistance ) > 0.8 ){
		IO_resistance = as.numeric( gsva( get_Scale( x= data ) , list(resistance) , verbose=FALSE ) )
	}
	
	IO_sensitive = NULL
	if( ifelse( is.null( nrow( data[ rownames(data) %in% sensitive , ]) ) , 1 , nrow( data[ rownames(data) %in% sensitive , ] ) ) / length( sensitive ) > 0.8 ){
		IO_sensitive = as.numeric( gsva( get_Scale( x= data ) , list(sensitive) , verbose=FALSE ) )
	}

	####################################################################################
	####################################################################################
	## Compte PredictIO
	signature = NULL
	if( !is.null( IO_resistance ) & !is.null( IO_sensitive ) ){

		signature = IO_sensitive - IO_resistance
		names(signature) = colnames(data)
	}
	signature
}

Get_PredictIO <- function( meta_res ){

	source('meta/Get_HR.R')
	source('meta/Get_DI.R')


	load("../data/process_data/ICB_exp_filtered.RData")

	expr
	study = names(expr)

	cox_os = cox_pfs = cox_dicho_os = cox_dicho_pfs = cox_tertile_os = cox_tertile_pfs = NULL
	di_os = di_pfs = NULL
	log_response = log_dicho_response = log_tertile_response = NULL

	for( i in 1:length(study)){

		tumor = names( table( phenoData( expr[[i]] )$primary )[ table( phenoData( expr[[i]] )$primary ) >= 20 ] )

		if( length(tumor) > 0 ){
			for( j in 1:length(tumor)){

				data = exprs(expr[[i]])[ ,  phenoData(expr[[i]])$primary %in% tumor[j] & phenoData(expr[[i]])$rna %in% c( "fpkm" , "tpm" ) ]
				
				IOscore = NULL
				if( ncol(data)){ IOscore <- getIOscore( data = data , meta_res ) }

				if( !is.null( IOscore ) ){
					if( length( IOscore[ !is.na( phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ] ) >= 20 ){
						## Compute the association PredictIO signature with OS (36 month cutoff) using Cox Regression Model
						hr = Get_HR_continous( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= IOscore )
						cox_os = rbind( cox_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IOscore[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

						## Compute the association PredictIO signature with OS (36 month cutoff) using Concordence Index
						cci = Get_DI_continous( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= IOscore )
						di_os = rbind( di_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IOscore[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  cci ) )

					
						## Compute the association PredictIO signature (HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= IOscore , cutoff= median( IOscore , na.rm=TRUE ) ,
												 title = paste( study[i] , tumor[j] , 36 , "OS_Overall" , sep="_" ) , ylab="Overall Survival" , dir= "../results/PredictIO/PredictIO/KMPlot/OS" )
						cox_dicho_os = rbind( cox_dicho_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IOscore[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

						## Compute the association PredictIO signature (HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= IOscore , cutoff= quantile( IOscore , probs=.66 ) ,
												 title = paste( study[i] , tumor[j] , 36 , "OS_Overall" , sep="_" ) , ylab="Overall Survival" , dir= "../results/PredictIO/PredictIO/KMPlot/OS" )
						cox_tertile_os = rbind( cox_tertile_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IOscore[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )



					}

					if( length( IOscore[ !is.na( phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) >= 20 ){
						## Compute the association of PredictIO signature with PFS using Cox Regression Model
						hr = Get_HR_continous( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= IOscore )
						cox_pfs = rbind( cox_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IOscore[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

						## Compute the association of PredictIO signature with PFS using Concordence Index
						cci = Get_DI_continous( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= IOscore )
						di_pfs = rbind( di_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IOscore[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  cci ) )
					
						
						## Compute the association of PredictIO signature (HIGH vs LOW) with PFS using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= IOscore , cutoff= median( IOscore , na.rm=TRUE ) ,
												 title = paste( study[i] , tumor[j] , 36 , "PFS_Overall" , sep="_" ) , ylab="Progression-Free Survival" , dir= "../results/PredictIO/PredictIO/KMPlot/PFS" )
						cox_dicho_pfs = rbind( cox_dicho_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IOscore[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

						
						## Compute the association of PredictIO signature (HIGH vs LOW) with PFS using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= IOscore , cutoff= quantile( IOscore , probs=.66 ) ,
												 title = paste( study[i] , tumor[j] , 36 , "PFS_Overall" , sep="_" ) , ylab="Progression-Free Survival" , dir= "../results/PredictIO/PredictIO/KMPlot/PFS" )
						cox_tertile_pfs = rbind( cox_tertile_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IOscore[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )


					}

					
					if( length( IOscore[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) >= 20 ){

						## Association of PredictIO signature with Response (Continous)
						x = ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

						fit = glm( x ~ IOscore , family=binomial( link="logit" ) )

						log_response = rbind( log_response , c( study[i] , 
									tumor[j] , 
									toupper( phenoData(expr[[i]])$rna[1] ) , 
									length( IOscore[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) , 
									round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
									round( confint(fit)[ 2 , ] , 2 ) , 
									summary(fit)$coefficients[ 2 , 4 ] ) )

						## Association of PredictIO signature with Response (HIGH vs Low)
						m = median( IOscore , na.rm=TRUE )
						y = ifelse( IOscore >= m , 1 ,0 )
						x = ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

						if( length(IOscore[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ] )>= 20 ){
							
							fit = glm( x ~ y , family=binomial( link="logit" ) )

							log_dicho_response = rbind( log_dicho_response , c( study[i] , 
										tumor[j] , 
										toupper( phenoData(expr[[i]])$rna[1] ) , 
										length( IOscore[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) , 
										round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
										round( confint(fit)[ 2 , ] , 2 ) , 
										summary(fit)$coefficients[ 2 , 4 ] ) )
						}
						## Association of PredictIO signature with Response (HIGH vs Low)
						m = quantile( IOscore , probs=.66 )
						y = ifelse( IOscore >= m , 1 ,0 )
						x = ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

						if( length(IOscore[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ] )>= 20 ){
							
							fit = glm( x ~ y , family=binomial( link="logit" ) )

							log_tertile_response = rbind( log_tertile_response , c( study[i] , 
										tumor[j] , 
										toupper( phenoData(expr[[i]])$rna[1] ) , 
										length( IOscore[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) , 
										round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
										round( confint(fit)[ 2 , ] , 2 ) , 
										summary(fit)$coefficients[ 2 , 4 ] ) )
						}
					}
				}
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

	save( log_response , log_dicho_response , log_tertile_response , file= "../results/PredictIO/PredictIO/PredictIO_LogReg_result.RData" )
	save( cox_os , cox_pfs , cox_dicho_os , cox_dicho_pfs , cox_tertile_os , cox_tertile_pfs , file= "../results/PredictIO/PredictIO/PredictIO_COX_result.RData" ) 
	save( di_os , di_pfs , file= "../results/PredictIO/PredictIO/PredictIO_DI_result.RData" ) 
}