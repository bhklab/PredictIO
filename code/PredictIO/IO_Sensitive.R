###############################################################################
###############################################################################

library(survcomp)
library(genefu)
library(GSVA)

#Remove rows if count is < zero in 50% of sample
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0) ))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}


Get_IO_sensitive <- function( meta_res ){

	sig = meta_res[ meta_res$coef < 0 & meta_res$include %in% 1 , ]$gene


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
				remove <- rem(data)
				if( length(remove) ){
					data <- data[-remove,]
				}

				if( ifelse( is.null( nrow(data[ rownames(data) %in% sig,]) ) , 1 , nrow(data[ rownames(data) %in% sig,]) ) / length( sig ) > 0.7 & ncol(data) > 20 ){
					geneSig <- (genefu::rescale( as.numeric( gsva( data, list(sig) , verbose=FALSE ) ) , na.rm=TRUE , q=0.05 ) - 0.5) * 2

					if( length( geneSig[ !is.na( phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ] ) >= 20 ){
						## Compute the association IO_Sensitive signature with OS (36 month cutoff) using Cox Regression Model
						hr = Get_HR_continous( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= geneSig )
						cox_os = rbind( cox_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( geneSig[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

						## Compute the association IO_Sensitive signature with OS (36 month cutoff) using Concordence Index
						cci = Get_DI_continous( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= geneSig )
						di_os = rbind( di_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( geneSig[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  cci ) )

					
						## Compute the association IO_Sensitive signature (HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= geneSig , cutoff= median( geneSig , na.rm=TRUE ) ,
												 title = paste( study[i] , tumor[j] , 36 , "OS_Overall" , sep="_" ) , ylab="Overall Survival" , dir= "../results/PredictIO/IO_Sensitive/KMPlot/OS" )
						cox_dicho_os = rbind( cox_dicho_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( geneSig[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )


						## Compute the association IO_Sensitive signature (HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= geneSig , cutoff= quantile(geneSig, probs=.66) ,
												 title = paste( study[i] , tumor[j] , 36 , "OS_Overall" , sep="_" ) , ylab="Overall Survival" , dir= "../results/PredictIO/IO_Sensitive/KMPlot/OS" )
						cox_tertile_os = rbind( cox_tertile_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( geneSig[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )



					}

					if( length( geneSig[ !is.na( phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) >= 20 ){
						## Compute the association of IO_Sensitive signature with PFS using Cox Regression Model
						hr = Get_HR_continous( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= geneSig )
						cox_pfs = rbind( cox_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( geneSig[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

						## Compute the association of IO_Sensitive signature with PFS using Concordence Index
						cci = Get_DI_continous( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= geneSig )
						di_pfs = rbind( di_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( geneSig[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  cci ) )
					
						
						## Compute the association of IO_Sensitive signature (HIGH vs LOW) with PFS using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= geneSig , cutoff= median( geneSig , na.rm=TRUE ) ,
												 title = paste( study[i] , tumor[j] , 36 , "PFS_Overall" , sep="_" ) , ylab="Progression-Free Survival" , dir= "../results/PredictIO/IO_Sensitive/KMPlot/PFS" )
						cox_dicho_pfs = rbind( cox_dicho_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( geneSig[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

						## Compute the association of IO_Sensitive signature (HIGH vs LOW) with PFS using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= geneSig , cutoff= quantile(geneSig, probs=.66) ,
												 title = paste( study[i] , tumor[j] , 36 , "PFS_Overall" , sep="_" ) , ylab="Progression-Free Survival" , dir= "../results/PredictIO/IO_Sensitive/KMPlot/PFS" )
						cox_tertile_pfs = rbind( cox_tertile_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( geneSig[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )


					}

					
					if( length( geneSig[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) >= 20 ){

						## Association of IO_Sensitive signature with Response (Continous)
						x = ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

						fit = glm( x ~ geneSig , family=binomial( link="logit" ) )

						log_response = rbind( log_response , c( study[i] , 
									tumor[j] , 
									toupper( phenoData(expr[[i]])$rna[1] ) , 
									length( geneSig[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) , 
									round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
									round( confint(fit)[ 2 , ] , 2 ) , 
									summary(fit)$coefficients[ 2 , 4 ] ) )

						## Association of IO_Sensitive signature with Response (HIGH vs Low)
						m = median( geneSig , na.rm=TRUE )
						y = ifelse( geneSig >= m , 1 ,0 )
						x = ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

						if( length(geneSig[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ] )>= 20 ){
							
							fit = glm( x ~ y , family=binomial( link="logit" ) )

							log_dicho_response = rbind( log_dicho_response , c( study[i] , 
										tumor[j] , 
										toupper( phenoData(expr[[i]])$rna[1] ) , 
										length( geneSig[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) , 
										round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
										round( confint(fit)[ 2 , ] , 2 ) , 
										summary(fit)$coefficients[ 2 , 4 ] ) )
						}

						## Association of IO_Sensitive signature with Response (HIGH vs Low)
						m = quantile(geneSig, probs=.66)
						y = ifelse( geneSig >= m , 1 ,0 )
						x = ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

						if( length(geneSig[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ] )>= 20 ){
							
							fit = glm( x ~ y , family=binomial( link="logit" ) )

							log_tertile_response = rbind( log_tertile_response , c( study[i] , 
										tumor[j] , 
										toupper( phenoData(expr[[i]])$rna[1] ) , 
										length( geneSig[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) , 
										round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
										round( confint(fit)[ 2 , ] , 2 ) , 
										summary(fit)$coefficients[ 2 , 4 ] ) )
						}


					}
				}	

			}
		}
	}

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

	save( log_response , log_dicho_response , log_tertile_response , file= "../results/PredictIO/IO_Sensitive/IO_Sensitive_LogReg_result.RData" )
	save( cox_os , cox_pfs , cox_dicho_os , cox_dicho_pfs , cox_tertile_os , cox_tertile_pfs , file= "../results/PredictIO/IO_Sensitive/IO_Sensitive_COX_result.RData" ) 
	save(  di_os , di_pfs , file= "../results/PredictIO/IO_Sensitive/IO_Sensitive_DI_result.RData" ) 
}