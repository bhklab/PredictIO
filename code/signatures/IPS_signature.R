##########################################################################################################################################################
##########################################################################################################################################################

library(survcomp)
library(genefu)

#################################################################
#################################################################
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

#################################################################
#################################################################

Get_IPS = function( expr ){
	####################################################
	##
	##   This R-script can be used to calculate Immunophenoscore (IPS) and generate Immunophenogram from "EXPR.txt" and "IPS_genes.txt"
	##   (C) ICBI, Medical University of Innsbruck, Biocenter, Division of Bioinformatics
	##   Version 1.0 08.07.2016
	##   Needs packages ggplot2,grid,gridExtra
	##
	####################################################

	## calculate Immunophenoscore
	ipsmap<- function (x) {
		if (x<=0) {
			ips<-0
		} else {
			if (x>=3) {
			 ips<-10
			} else {
				ips<-round(x*10/3, digits=0)
			}
		}
		return(ips)
	}

	expr = as.data.frame(expr)
	sample_names <- names(expr)


	## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
	# For different 
	IPSG<-read.table("../data/signatures/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
	unique_ips_genes<-as.vector(unique(IPSG$NAME))

	IPS<-NULL
	MHC<-NULL
	CP<-NULL
	EC<-NULL
	SC<-NULL
	AZ<-NULL

	# Gene names in expression file
	GVEC <- row.names( expr )
	# Genes names in IPS genes file
	VEC <- as.vector( IPSG$GENE )
	# Match IPS genes with genes in expression file
	ind <- which( is.na( match( VEC , GVEC ) ) )
	if( length( ind ) ){
		IPSG <- IPSG[-ind,]
	}

	for (i in 1:length(sample_names)) {	
		GE<-expr[[i]]
		mGE<-mean(GE, na.rm=T)
		sGE<-sd(GE, na.rm=T)
		Z1<-(expr[as.vector(IPSG$GENE),i]-mGE)/sGE
		W1<-IPSG$WEIGHT
		WEIGHT<-NULL
		MIG<-NULL
		k<-1
		for (gen in unique_ips_genes) {
			MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
			WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)], na.rm=T)
			k<-k+1
		}
		WG<-MIG*WEIGHT
		MHC[i]<-mean(WG[1:10],na.rm=T)
		CP[i]<-mean(WG[11:20],na.rm=T)
		EC[i]<-mean(WG[21:24],na.rm=T)
		SC[i]<-mean(WG[25:26],na.rm=T)
		AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
		IPS[i]<-ipsmap(AZ[i])
	}
	AZ
}

#################################################
#################################################

Get_IPS_signature = function( expr ){

	source('meta/Get_HR.R')
	source('meta/Get_DI.R')
	get_Directory( dir= "IPS" )


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
				IPS=NULL
				if( ncol(data) & nrow(data)>10000 ){ 
					IPS = as.numeric( scale( Get_IPS( expr=data ) ) )
					names( IPS ) = colnames(data)
				}
				
				if( ! is.null( length( IPS ) ) & length( IPS ) > 20 ){

					if( length( IPS[ !is.na( phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ] ) >= 20 ){
						## Compute the association ADO signature with OS (36 month cutoff) using Cox Regression Model
						hr = Get_HR_continous( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= IPS )
						cox_os = rbind( cox_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IPS[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

						## Compute the association ADO signature with OS (36 month cutoff) using Concordence Index
						cci = Get_DI_continous( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= IPS )
						di_os = rbind( di_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IPS[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  cci ) )

					
						## Compute the association ADO signature (HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=36 , variable= IPS , cutoff= median( IPS , na.rm=TRUE ) ,
												 title = paste( study[i] , tumor[j] , 36 , "OS_Overall" , sep="_" ) , ylab="Overall Survival" , dir="../results/Signature/IPS/KMPlot/OS" )
						cox_dicho_os = rbind( cox_dicho_os , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IPS[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )


					}

					if( length( IPS[ !is.na( phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) >= 20 ){
						## Compute the association of ADO signature with PFS using Cox Regression Model
						hr = Get_HR_continous( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= IPS )
						cox_pfs = rbind( cox_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IPS[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

						## Compute the association of ADO signature with PFS using Concordence Index
						cci = Get_DI_continous( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= IPS )
						di_pfs = rbind( di_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IPS[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  cci ) )
					
						
						## Compute the association of ADO signature (HIGH vs LOW) with PFS using Cox Regression Model
						hr = Get_HR_dicho( surv=phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] , time=phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% tumor[j] ] ,
												 time_censor=24 , variable= IPS , cutoff= median( IPS , na.rm=TRUE ) ,
												 title = paste( study[i] , tumor[j] , 36 , "PFS_Overall" , sep="_" ) , ylab="Progression-Free Survival" , dir="../results/Signature/IPS/KMPlot/PFS" )
						cox_dicho_pfs = rbind( cox_dicho_pfs , 
								c( study[i] , tumor[j] , toupper( phenoData(expr[[i]])$rna[1] ) , length( IPS[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

					}

					
					if( length( IPS[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) >= 20 ){

						## Association with Response (Continous)
						x = ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

						fit = glm( x ~ IPS , family=binomial( link="logit" ) )

						log_response = rbind( log_response , c( study[i] , 
									tumor[j] , 
									toupper( phenoData(expr[[i]])$rna[1] ) , 
									length( IPS[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) , 
									round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
									round( confint(fit)[ 2 , ] , 2 ) , 
									summary(fit)$coefficients[ 2 , 4 ] ) )

						## Association with Response (HIGH vs Low)
						m = median( IPS , na.rm=TRUE )
						y = ifelse( IPS >= m , 1 ,0 )
						x = ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

						if( length(IPS[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ] )>= 40 ){
							
							fit = glm( x ~ y , family=binomial( link="logit" ) )

							log_dicho_response = rbind( log_dicho_response , c( study[i] , 
										tumor[j] , 
										toupper( phenoData(expr[[i]])$rna[1] ) , 
										length( IPS[ !is.na( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% tumor[j] ] ) ]  ) , 
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

	save( log_response , log_dicho_response , file="../results/Signature/IPS/IPS_LogReg_result.RData" )
	save( cox_os , cox_pfs , cox_dicho_os , cox_dicho_pfs , file="../results/Signature/IPS/IPS_COX_result.RData" ) 
	save(  di_os , di_pfs , file="../results/Signature/IPS/IPS_DI_result.RData" ) 
}