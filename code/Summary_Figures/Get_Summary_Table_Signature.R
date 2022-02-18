########################################################################################################################
########################################################################################################################


library(meta)
library(metafor)
library(genefu)
library(dmetar)


get_Meta_HR = function( data , signature ){

	data$study <- as.character( data$study )
	data$Sequencing <- as.character( data$Sequencing )
	data$Primary <- as.character( data$Primary )
	data$HR <- as.numeric(as.character( data$HR ))
	data$SE <- as.numeric(as.character( data$SE ))
	data$Pval <- as.numeric(as.character( data$Pval )) 

	meta <- metagen( TE = HR ,
					seTE = SE ,
					data = data ,
					studlab = study ,
					fixed = FALSE ,
					random = TRUE ,
					control = list( maxiter = 10000 , stepadj=0.5 ) )

	######################################################################
	######################################################################
	## Save the merged HR and Pvalue
	meta_res <- c( signature , nrow(data) ,
					round( meta$TE.random , 2 ) ,  
					round( meta$seTE.random , 2 ) ,   
					round( meta$lower.random , 2 ) ,   
					round( meta$upper.random , 2 ) ,
					meta$pval.random , 
					round( meta$I2 , 2 ) ,
					meta$pval.Q )
	names(meta_res) <- c( "study" , "N" , "logHR" , "se_logHR" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" )
	meta_res
}

get_Meta_OR = function( data , signature ){

	data$study <- as.character( data$study )
	data$Sequencing <- as.character( data$Sequencing )
	data$Primary <- as.character( data$Primary )
	data$coef <- as.numeric(as.character( data$coef ))
	data$SE <- as.numeric(as.character( data$SE ))
	data$Pval <- as.numeric(as.character( data$Pval )) 

	meta <- metagen( TE = coef,
					seTE = SE,
					data = data,
					studlab = study ,
					fixed = FALSE ,
					random = TRUE ,
					control = list( maxiter = 10000 , stepadj=0.5 ) )

	######################################################################
	######################################################################
	## Save the merged coef and Pvalue
	meta_res <- c( signature , nrow(data) ,
					round( meta$TE.random , 2 ) ,  
					round( meta$seTE.random , 2 ) ,   
					round( meta$lower.random , 2 ) ,   
					round( meta$upper.random , 2 ) ,
					meta$pval.random , 
					round( meta$I2 , 2 ) ,
					meta$pval.Q )
	names(meta_res) <- c( "study" , "N" , "logOR" , "se_logOR" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" )
	meta_res
}

get_Meta_DI = function( data , signature ){

	data$study <- as.character( data$study )
	data$Sequencing <- as.character( data$Sequencing )
	data$Primary <- as.character( data$Primary )
	data$DI <- as.numeric(as.character( data$DI ))
	data$SE <- as.numeric(as.character( data$SE ))
	data$Pval <- as.numeric(as.character( data$Pval )) 

	meta <- metagen( TE = DI,
					seTE = SE,
					data = data,
					studlab = study ,
					fixed = FALSE ,
					random = TRUE ,
					control = list( maxiter = 10000 , stepadj=0.5 ) )

	######################################################################
	######################################################################
	## Save the merged DI and Pvalue
	meta_res <- c( signature , nrow(data) ,
					round( meta$TE.random , 2 ) ,  
					round( meta$seTE.random , 2 ) ,   
					round( meta$lower.random , 2 ) ,   
					round( meta$upper.random , 2 ) ,
					meta$pval.random , 
					round( meta$I2 , 2 ) ,
					meta$pval.Q )
	names(meta_res) <- c( "study" , "N" , "DI" , "se_DI" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" )
	meta_res
}

########################################################################################################################
########################################################################################################################

tumorID = c( "Melanoma" , "Lung" , "Kidney" )
sequencingID = c( "TPM" , "FPKM" )

signature = read.table( file= "../data/signatures/signature_INFO.txt" , sep="\t" , header=TRUE , stringsAsFactor=FALSE)
signature$Signature <- as.character( signature$Signature )
signature$method <- as.character( signature$method )
signature$association <- as.character( signature$association )

SignatureID = signature$Signature

res_indi = res_meta = NULL

for(k in 1:length(SignatureID)){

	cox_RData = paste( "../results/Signature/" , SignatureID[k] , "/" , SignatureID[k] , "_COX_result.RData"  , sep="")
	di_RData = paste( "../results/Signature/" , SignatureID[k] , "/" , SignatureID[k] , "_DI_result.RData" , sep="")
	log_RData = paste( "../results/Signature/" , SignatureID[k] , "/" , SignatureID[k] , "_LogReg_result.RData"  , sep="")

	sigID = SignatureID[k]


	########################################################################################################################
	########################################################################################################################

	load( cox_RData )

	##################################################
	## Meta-analysis of the COX models (OS) Continous
	cox_os = cox_os[ cox_os$SE <= 10 & !is.na( cox_os$Pval ) , ]

	if( nrow( cox_os ) ){

		#######################
		## Global meta-analysis
		meta_res = c( "ALL" , get_Meta_HR( data=cox_os , signature = "ALL" ) )


		#######################################
		## Subgroup meta-analysis :  Tumor Type
		for( j in 1:length(tumorID) ){
			data = cox_os[ cox_os$Primary %in% tumorID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_HR( data=data , signature = tumorID[j] ) 
				meta_res = rbind( meta_res , c( "Tumor" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Tumor" , tumorID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		## Subgroup meta-analysis : Sequencing Type
		for( j in 1:length(sequencingID) ){
			data = cox_os[ cox_os$Sequencing %in% sequencingID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_HR( data=data , signature = sequencingID[j] ) 
				meta_res = rbind( meta_res , c( "Sequencing" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Sequencing" , sequencingID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		###########################################
		## Result merging

		out = cbind( sigID , "Signature" , "OS" , "COX" , cox_os )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
		res_indi = rbind( res_indi , out )

		out = cbind( sigID , "Signature" , "OS" , "COX" , meta_res )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
		res_meta = rbind( res_meta , out )

	}
	##################################################
	## Meta-analysis of the COX models (PFS) Continous
	cox_pfs = cox_pfs[ cox_pfs$SE <= 10 & !is.na( cox_pfs$Pval ) , ]

	if( nrow( cox_pfs ) ){
		
		#######################
		## Global meta-analysis
		meta_res = c( "ALL" , get_Meta_HR( data=cox_pfs , signature = "ALL" ) )


		#######################################
		## Subgroup meta-analysis :  Tumor Type
		for( j in 1:length(tumorID) ){
			data = cox_pfs[ cox_pfs$Primary %in% tumorID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_HR( data=data , signature = tumorID[j] ) 
				meta_res = rbind( meta_res , c( "Tumor" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Tumor" , tumorID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		## Subgroup meta-analysis : Sequencing Type
		for( j in 1:length(sequencingID) ){
			data = cox_pfs[ cox_pfs$Sequencing %in% sequencingID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_HR( data=data , signature = sequencingID[j] ) 
				meta_res = rbind( meta_res , c( "Sequencing" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Sequencing" , sequencingID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		###########################################
		## Result merging

		out = cbind( sigID , "Signature" , "PFS" , "COX" , cox_pfs )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
		res_indi = rbind( res_indi , out )

		out = cbind( sigID , "Signature" , "PFS" , "COX" , meta_res )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
		res_meta = rbind( res_meta , out )
	}
	########################################################################################################################
	########################################################################################################################

	load( di_RData )

	##################################################
	## Meta-analysis of the DI models (OS) Continous
	di_os = di_os[ di_os$SE <= 10 & !is.na( di_os$Pval ) , ]

	if( nrow( di_os ) ){
	
		#######################
		## Global meta-analysis
		meta_res = c( "ALL" , get_Meta_DI( data=di_os , signature = "ALL" ) )


		#######################################
		## Subgroup meta-analysis :  Tumor Type
		for( j in 1:length(tumorID) ){
			data = di_os[ di_os$Primary %in% tumorID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_DI( data=data , signature = tumorID[j] ) 
				meta_res = rbind( meta_res , c( "Tumor" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Tumor" , tumorID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		## Subgroup meta-analysis : Sequencing Type
		for( j in 1:length(sequencingID) ){
			data = di_os[ di_os$Sequencing %in% sequencingID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_DI( data=data , signature = sequencingID[j] ) 
				meta_res = rbind( meta_res , c( "Sequencing" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Sequencing" , sequencingID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		###########################################
		## Result merging

		out = cbind( sigID , "Signature" , "OS" , "DI" , di_os )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
		res_indi = rbind( res_indi , out )

		out = cbind( sigID , "Signature" , "OS" , "DI" , meta_res )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
		res_meta = rbind( res_meta , out )
	}
	##################################################
	## Meta-analysis of the DI models (PFS) Continous
	di_pfs = di_pfs[ di_pfs$SE <= 10 & !is.na( di_pfs$Pval ) , ]

	if( nrow( di_pfs ) ){
	
		#######################
		## Global meta-analysis
		meta_res = c( "ALL" , get_Meta_DI( data=di_pfs , signature = "ALL" ) )


		#######################################
		## Subgroup meta-analysis :  Tumor Type
		for( j in 1:length(tumorID) ){
			data = di_pfs[ di_pfs$Primary %in% tumorID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_DI( data=data , signature = tumorID[j] ) 
				meta_res = rbind( meta_res , c( "Tumor" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Tumor" , tumorID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		## Subgroup meta-analysis : Sequencing Type
		for( j in 1:length(sequencingID) ){
			data = di_pfs[ di_pfs$Sequencing %in% sequencingID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_DI( data=data , signature = sequencingID[j] ) 
				meta_res = rbind( meta_res , c( "Sequencing" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Sequencing" , sequencingID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		###########################################
		## Result merging

		out = cbind( sigID , "Signature" , "PFS" , "DI" , di_pfs )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
		res_indi = rbind( res_indi , out )

		out = cbind( sigID , "Signature" , "PFS" , "DI" , meta_res )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
		res_meta = rbind( res_meta , out )
	}
	########################################################################################################################
	########################################################################################################################

	load( log_RData )

	##################################################
	## Meta-analysis of the Log regression models (OS) Continous
	log_response = log_response[ log_response$SE <= 10 & !is.na( log_response$Pval ) , ]

	if( nrow( log_response ) ){
		
		#######################
		## Global meta-analysis
		meta_res = c( "ALL" , get_Meta_OR( data=log_response , signature = "ALL" ) )


		#######################################
		## Subgroup meta-analysis :  Tumor Type
		for( j in 1:length(tumorID) ){
			data = log_response[ log_response$Primary %in% tumorID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_OR( data=data , signature = tumorID[j] ) 
				meta_res = rbind( meta_res , c( "Tumor" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Tumor" , tumorID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		## Subgroup meta-analysis : Sequencing Type
		for( j in 1:length(sequencingID) ){
			data = log_response[ log_response$Sequencing %in% sequencingID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_OR( data=data , signature = sequencingID[j] ) 
				meta_res = rbind( meta_res , c( "Sequencing" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Sequencing" , sequencingID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		###########################################
		## Result merging

		out = cbind( sigID , "Signature" , "Response" , "Log_regression" , log_response )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
		res_indi = rbind( res_indi , out )

		out = cbind( sigID , "Signature" , "Response" , "Log_regression" , meta_res )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
		res_meta = rbind( res_meta , out )
	}
}


########################################################################################################################
########################################################################################################################


tumorID = c( "Melanoma" , "Lung" , "Kidney" )
sequencingID = c( "TPM" , "FPKM" )

SignatureID = "PredictIO"

for(k in 1:length(SignatureID)){

	cox_RData = paste( "../results/PredictIO/" , SignatureID[k] , "/" , SignatureID[k] , "_COX_result.RData"  , sep="")
	di_RData = paste( "../results/PredictIO/" , SignatureID[k] , "/" , SignatureID[k] , "_DI_result.RData" , sep="")
	log_RData = paste( "../results/PredictIO/" , SignatureID[k] , "/" , SignatureID[k] , "_LogReg_result.RData"  , sep="")

	sigID = SignatureID[k]


	########################################################################################################################
	########################################################################################################################

	load( cox_RData )

	##################################################
	## Meta-analysis of the COX models (OS) Continous
	cox_os = cox_os[ cox_os$SE <= 10 & !is.na( cox_os$Pval ) , ]

	if( nrow( cox_os ) ){
	
		#######################
		## Global meta-analysis
		meta_res = c( "ALL" , get_Meta_HR( data=cox_os , signature = "ALL" ) )


		#######################################
		## Subgroup meta-analysis :  Tumor Type
		for( j in 1:length(tumorID) ){
			data = cox_os[ cox_os$Primary %in% tumorID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_HR( data=data , signature = tumorID[j] ) 
				meta_res = rbind( meta_res , c( "Tumor" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Tumor" , tumorID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		## Subgroup meta-analysis : Sequencing Type
		for( j in 1:length(sequencingID) ){
			data = cox_os[ cox_os$Sequencing %in% sequencingID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_HR( data=data , signature = sequencingID[j] ) 
				meta_res = rbind( meta_res , c( "Sequencing" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Sequencing" , sequencingID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		###########################################
		## Result merging

		out = cbind( sigID , "MetaScore" , "OS" , "COX" , cox_os )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
		res_indi = rbind( res_indi , out )

		out = cbind( sigID , "MetaScore" , "OS" , "COX" , meta_res )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
		res_meta = rbind( res_meta , out )
	}
	##################################################
	## Meta-analysis of the COX models (PFS) Continous
	cox_pfs = cox_pfs[ cox_pfs$SE <= 10 & !is.na( cox_pfs$Pval ) , ]

	if( nrow( cox_pfs ) ){
	
		#######################
		## Global meta-analysis
		meta_res = c( "ALL" , get_Meta_HR( data=cox_pfs , signature = "ALL" ) )


		#######################################
		## Subgroup meta-analysis :  Tumor Type
		for( j in 1:length(tumorID) ){
			data = cox_pfs[ cox_pfs$Primary %in% tumorID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_HR( data=data , signature = tumorID[j] ) 
				meta_res = rbind( meta_res , c( "Tumor" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Tumor" , tumorID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		## Subgroup meta-analysis : Sequencing Type
		for( j in 1:length(sequencingID) ){
			data = cox_pfs[ cox_pfs$Sequencing %in% sequencingID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_HR( data=data , signature = sequencingID[j] ) 
				meta_res = rbind( meta_res , c( "Sequencing" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Sequencing" , sequencingID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		###########################################
		## Result merging

		out = cbind( sigID , "MetaScore" , "PFS" , "COX" , cox_pfs )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
		res_indi = rbind( res_indi , out )

		out = cbind( sigID , "MetaScore" , "PFS" , "COX" , meta_res )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
		res_meta = rbind( res_meta , out )
	}
	########################################################################################################################
	########################################################################################################################

	load( di_RData )

	##################################################
	## Meta-analysis of the DI models (OS) Continous
	di_os = di_os[ di_os$SE <= 10 & !is.na( di_os$Pval ) , ]

	if( nrow( di_os ) ){
	
		#######################
		## Global meta-analysis
		meta_res = c( "ALL" , get_Meta_DI( data=di_os , signature = "ALL" ) )


		#######################################
		## Subgroup meta-analysis :  Tumor Type
		for( j in 1:length(tumorID) ){
			data = di_os[ di_os$Primary %in% tumorID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_DI( data=data , signature = tumorID[j] ) 
				meta_res = rbind( meta_res , c( "Tumor" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Tumor" , tumorID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		## Subgroup meta-analysis : Sequencing Type
		for( j in 1:length(sequencingID) ){
			data = di_os[ di_os$Sequencing %in% sequencingID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_DI( data=data , signature = sequencingID[j] ) 
				meta_res = rbind( meta_res , c( "Sequencing" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Sequencing" , sequencingID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		###########################################
		## Result merging

		out = cbind( sigID , "MetaScore" , "OS" , "DI" , di_os )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
		res_indi = rbind( res_indi , out )

		out = cbind( sigID , "MetaScore" , "OS" , "DI" , meta_res )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
		res_meta = rbind( res_meta , out )
	}
	##################################################
	## Meta-analysis of the DI models (PFS) Continous
	di_pfs = di_pfs[ di_pfs$SE <= 10 & !is.na( di_pfs$Pval ) , ]

	if( nrow( di_pfs ) ){
	
		#######################
		## Global meta-analysis
		meta_res = c( "ALL" , get_Meta_DI( data=di_pfs , signature = "ALL" ) )


		#######################################
		## Subgroup meta-analysis :  Tumor Type
		for( j in 1:length(tumorID) ){
			data = di_pfs[ di_pfs$Primary %in% tumorID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_DI( data=data , signature = tumorID[j] ) 
				meta_res = rbind( meta_res , c( "Tumor" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Tumor" , tumorID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		## Subgroup meta-analysis : Sequencing Type
		for( j in 1:length(sequencingID) ){
			data = di_pfs[ di_pfs$Sequencing %in% sequencingID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_DI( data=data , signature = sequencingID[j] ) 
				meta_res = rbind( meta_res , c( "Sequencing" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Sequencing" , sequencingID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		###########################################
		## Result merging

		out = cbind( sigID , "MetaScore" , "PFS" , "DI" , di_pfs )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
		res_indi = rbind( res_indi , out )

		out = cbind( sigID , "MetaScore" , "PFS" , "DI" , meta_res )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
		res_meta = rbind( res_meta , out )
	}
		########################################################################################################################
	########################################################################################################################

	load( log_RData )

	##################################################
	## Meta-analysis of the Log regression models (OS) Continous
	log_response = log_response[ log_response$SE <= 10 & !is.na( log_response$Pval ) , ]

	if( nrow( log_response ) ){
	
		#######################
		## Global meta-analysis
		meta_res = c( "ALL" , get_Meta_OR( data=log_response , signature = "ALL" ) )


		#######################################
		## Subgroup meta-analysis :  Tumor Type
		for( j in 1:length(tumorID) ){
			data = log_response[ log_response$Primary %in% tumorID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_OR( data=data , signature = tumorID[j] ) 
				meta_res = rbind( meta_res , c( "Tumor" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Tumor" , tumorID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		## Subgroup meta-analysis : Sequencing Type
		for( j in 1:length(sequencingID) ){
			data = log_response[ log_response$Sequencing %in% sequencingID[j] , ]
			if( nrow( data ) >= 3 ){
				res = get_Meta_OR( data=data , signature = sequencingID[j] ) 
				meta_res = rbind( meta_res , c( "Sequencing" , res ) )
			} else{
				meta_res = rbind( meta_res , c( "Sequencing" , sequencingID[j] , nrow(data) , rep( NA , 7 ) ) )
			}
		}

		###########################################
		###########################################
		## Result merging

		out = cbind( sigID , "MetaScore" , "Response" , "Log_regression" , log_response )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
		res_indi = rbind( res_indi , out )

		out = cbind( sigID , "MetaScore" , "Response" , "Log_regression" , meta_res )
		colnames(out) = c( "signature" , "signatureType" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
		res_meta = rbind( res_meta , out )
	}
}

########################################################################################################################
########################################################################################################################

write.table( res_indi , file= "../results/Summary_Figure/Table/Signature/Signature_Individual.txt" , sep="\t" , row.names=FALSE , col.names=TRUE , quote=FALSE )
write.table( res_meta , file= "../results/Summary_Figure/Table/Signature/Signature_Meta_analysis.txt" , sep="\t" , row.names=FALSE , col.names=TRUE , quote=FALSE )
