##########################################################################################################################################################
##########################################################################################################################################################

library(survcomp)
library(genefu)
library(data.table)
library(Biobase)

#Remove rows if count is < zero in 50% of sample
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0) ))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}

if_NULL = function( x , colnames , study ){
	if( is.null( x ) ){
		x = as.data.frame( cbind( study , matrix( NA , nrow= length(study) , ncol=length( colnames )-1 ) ) )
	}

	x
}

get_Directory <- function( dir ){

	file = paste( "../results/Signature/" , dir , sep="" )
	
	if( file.exists( file ) ) {
		unlink( file , recursive = TRUE )
	}

	dir.create( file , recursive = TRUE )

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


Get_TIDE_signature = function( expr ){

	source('meta/Get_HR.R')
	source('meta/Get_DI.R')
	get_Directory( dir= "TIDE" )

	study = names(expr)

	cox_os = cox_pfs = cox_dicho_os = cox_dicho_pfs = NULL

	di_os = di_pfs = NULL
	log_response = log_dicho_response = NULL

	for( i in 1:length(study)){

		tumor = names( table( phenoData( expr[[i]] )$primary )[ table( phenoData( expr[[i]] )$primary ) >= 20 ] )

		if( length(tumor) > 0 ){
			for( j in 1:length(tumor)){

				data = exprs(expr[[i]])[ ,  phenoData(expr[[i]])$primary %in% tumor[j] & phenoData(expr[[i]])$rna %in% c( "fpkm" , "tpm" ) ]
				remove <- rem(data)
				if( length(remove) ){
					data <- data[-remove,]
				}
				
				TIDE=NULL
				if( ncol(data) ){ 
					TIDE = Get_TIDE( expr=data , primary= tumor )
				}
				
				if( ! is.null( length( TIDE ) ) & length( TIDE ) > 20 ){

					if( length( TIDE[ !is.na( phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ] ) >= 20 ){
						## Compute the association ADO signature with OS (36 month cutoff) using Cox Regression Model
						hr = Get_HR_continous( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= TIDE )
						cox_os = rbind( cox_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( TIDE[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

						## Compute the association ADO signature with OS (36 month cutoff) using Concordence Index
						cci = Get_DI_continous( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= TIDE )
						di_os = rbind( di_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( TIDE[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  cci ) )

					
						## Compute the association ADO signature (HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= TIDE , cutoff= median( TIDE , na.rm=TRUE ) ,
												 title = paste( study[i] , tumor[j] , 36 , "OS_Overall" , sep="_" ) , ylab="Overall Survival" , dir="../results/Signature/TIDE/KMPlot/OS" )
						cox_dicho_os = rbind( cox_dicho_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( TIDE[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

					

					}

					if( length( TIDE[ !is.na( phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) >= 20 ){
						## Compute the association of ADO signature with PFS using Cox Regression Model
						hr = Get_HR_continous( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= TIDE )
						cox_pfs = rbind( cox_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( TIDE[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

						## Compute the association of ADO signature with PFS using Concordence Index
						cci = Get_DI_continous( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= TIDE )
						di_pfs = rbind( di_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( TIDE[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  cci ) )
					
						
						## Compute the association of ADO signature (HIGH vs LOW) with PFS using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= TIDE , cutoff= median( TIDE , na.rm=TRUE ) ,
												 title = paste( study[i] , tumor[j] , 36 , "PFS_Overall" , sep="_" ) , ylab="Progression-Free Survival" , dir="../results/Signature/TIDE/KMPlot/PFS" )
						cox_dicho_pfs = rbind( cox_dicho_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( TIDE[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )
					}

					
					if( length( TIDE[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) >= 20 ){

						## Association with Response (Continous)
						x = ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

						fit = glm( x ~ TIDE , family=binomial( link="logit" ) )

						log_response = rbind( log_response , c( study[i] , 
									tumor[j] , 
									toupper( phenoData(expr[[i]])$rna[1] ) , 
									length( TIDE[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) , 
									round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
									round( confint(fit)[ 2 , ] , 2 ) , 
									summary(fit)$coefficients[ 2 , 4 ] ) )

						## Association with Response (HIGH vs Low)
						m = 1
						y = ifelse( TIDE >= m , 1 ,0 )
						x = ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

						if( length(TIDE[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ] )>= 40 ){
							
							fit = glm( x ~ y , family=binomial( link="logit" ) )

							log_dicho_response = rbind( log_dicho_response , c( study[i] , 
										tumor[j] , 
										toupper( phenoData(expr[[i]])$rna[1] ) , 
										length( TIDE[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) , 
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
	cox_os = if_NULL( x= cox_os , colnames= colnames , study= study )
	cox_pfs = if_NULL( x= cox_pfs , colnames= colnames , study= study )
	cox_dicho_os = if_NULL( x= cox_dicho_os , colnames= colnames , study= study )
	cox_dicho_pfs = if_NULL( x= cox_dicho_pfs , colnames= colnames , study= study )

	di_os = if_NULL( x= di_os , colnames= colnames , study= study )
	di_pfs = if_NULL( x= di_pfs , colnames= colnames , study= study )

	log_response = if_NULL( x= log_response , colnames= colnames , study= study )
	log_dicho_response = if_NULL( x= log_dicho_response , colnames= colnames , study= study )

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

	save( log_response , log_dicho_response , file="../results/Signature/TIDE/TIDE_LogReg_result.RData" )
	save( cox_os , cox_pfs , cox_dicho_os , cox_dicho_pfs , file="../results/Signature/TIDE/TIDE_COX_result.RData" ) 
	save(  di_os , di_pfs , file="../results/Signature/TIDE/TIDE_DI_result.RData" ) 
}