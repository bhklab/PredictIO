########################################################################################################################
########################################################################################################################

library(meta)
library(metafor)
library(genefu)
library(forestplot)
library(dmetar)


get_Meta_HR = function( data , signature ){

	data$study <- as.character( data$study )
	data$Sequencing <- as.character( data$Sequencing )
	data$Primary <- as.character( data$Primary )
	data$HR <- as.numeric(as.character( data$HR ))
	data$SE <- as.numeric(as.character( data$SE ))
	data$Pval <- as.numeric(as.character( data$Pval )) 

	meta <- metagen( TE = HR,
					seTE = SE,
					data = data,
					studlab = study ,
					fixed = FALSE ,
					random = TRUE ,
					control = list( maxiter = 10000 , stepadj=0.5 ) )

	######################################################################
	######################################################################
	## Save the merged HR and Pvalue
	meta_res <- c( signature , nrow(data) ,
					meta$TE.random ,  
					meta$seTE.random ,   
					meta$lower.random ,   
					meta$upper.random ,
					meta$pval.random , 
					meta$I2 ,
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
					meta$TE.random ,  
					meta$seTE.random ,   
					meta$lower.random ,   
					meta$upper.random ,
					meta$pval.random , 
					meta$I2 ,
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
					meta$TE.random ,  
					meta$seTE.random ,   
					meta$lower.random ,   
					meta$upper.random ,
					meta$pval.random , 
					meta$I2 ,
					meta$pval.Q )
	names(meta_res) <- c( "study" , "N" , "DI" , "se_DI" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" )
	meta_res
}

########################################################################################################################
########################################################################################################################

tumorID = c( "Melanoma" , "Lung" , "Kidney" )
sequencingID = c( "WES" , "TGS" )

res_indi = res_meta = NULL

cox_RData = paste( "../results/nsTMB/nsTMB_COX_result.RData"  , sep="")
di_RData = paste( "../results/nsTMB/nsTMB_DI_result.RData" , sep="")
log_RData = paste( "../results/nsTMB/nsTMB_LogReg_result.RData"  , sep="")

sigID = "nsTMB"


########################################################################################################################
########################################################################################################################

load( cox_RData )

##################################################
##################################################
## Meta-analysis of the COX models (OS) Continuous
cox_os = cox_os[ as.numeric( as.character( cox_os$SE ) ) <= 10 & !is.na( cox_os$Pval ) , ]


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

out = cbind( sigID , "continuous" , "OS" , "COX" , cox_os )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
res_indi = rbind( res_indi , out )

out = cbind( sigID , "continuous" , "OS" , "COX" , meta_res )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
res_meta = rbind( res_meta , out )

##################################################
##################################################
## Meta-analysis of the COX models (OS) Continuous
cox_dicho_median_os = cox_dicho_median_os[ as.numeric( as.character( cox_dicho_median_os$SE ) ) <= 10 & !is.na( cox_dicho_median_os$Pval ) , ]


#######################
## Global meta-analysis
meta_res = c( "ALL" , get_Meta_HR( data=cox_dicho_median_os , signature = "ALL" ) )


#######################################
## Subgroup meta-analysis :  Tumor Type
for( j in 1:length(tumorID) ){
	data = cox_dicho_median_os[ cox_dicho_median_os$Primary %in% tumorID[j] , ]
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
	data = cox_dicho_median_os[ cox_dicho_median_os$Sequencing %in% sequencingID[j] , ]
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

out = cbind( sigID , "median" , "OS" , "COX" , cox_dicho_median_os )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
res_indi = rbind( res_indi , out )

out = cbind( sigID , "median" , "OS" , "COX" , meta_res )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
res_meta = rbind( res_meta , out )


##################################################
##################################################
## Meta-analysis of the COX models (OS) Continuous
cox_dicho_10TMB_os = cox_dicho_10TMB_os[ as.numeric( as.character( cox_dicho_10TMB_os$SE ) ) <= 10 & !is.na( cox_dicho_10TMB_os$Pval ) , ]


#######################
## Global meta-analysis
meta_res = c( "ALL" , get_Meta_HR( data=cox_dicho_10TMB_os , signature = "ALL" ) )


#######################################
## Subgroup meta-analysis :  Tumor Type
for( j in 1:length(tumorID) ){
	data = cox_dicho_10TMB_os[ cox_dicho_10TMB_os$Primary %in% tumorID[j] , ]
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
	data = cox_dicho_10TMB_os[ cox_dicho_10TMB_os$Sequencing %in% sequencingID[j] , ]
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

out = cbind( sigID , "10TMB" , "OS" , "COX" , cox_dicho_10TMB_os )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
res_indi = rbind( res_indi , out )

out = cbind( sigID , "10TMB" , "OS" , "COX" , meta_res )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
res_meta = rbind( res_meta , out )


######################################################################################################################################################
######################################################################################################################################################

##################################################
## Meta-analysis of the COX models (PFS) Continuous
cox_pfs = cox_pfs[ as.numeric( as.character( cox_pfs$SE ) ) <= 10 & !is.na( cox_pfs$Pval ) , ]


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

out = cbind( sigID , "continuous" , "PFS" , "COX" , cox_pfs )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
res_indi = rbind( res_indi , out )

out = cbind( sigID , "continuous" , "PFS" , "COX" , meta_res )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
res_meta = rbind( res_meta , out )


##################################################
## Meta-analysis of the COX models (PFS) Median
cox_dicho_median_pfs = cox_dicho_median_pfs[ as.numeric( as.character( cox_dicho_median_pfs$SE ) ) <= 10 & !is.na( cox_dicho_median_pfs$Pval ) , ]


#######################
## Global meta-analysis
meta_res = c( "ALL" , get_Meta_HR( data=cox_dicho_median_pfs , signature = "ALL" ) )


#######################################
## Subgroup meta-analysis :  Tumor Type
for( j in 1:length(tumorID) ){
	data = cox_dicho_median_pfs[ cox_dicho_median_pfs$Primary %in% tumorID[j] , ]
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
	data = cox_dicho_median_pfs[ cox_dicho_median_pfs$Sequencing %in% sequencingID[j] , ]
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

out = cbind( sigID , "median" , "PFS" , "COX" , cox_dicho_median_pfs )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
res_indi = rbind( res_indi , out )

out = cbind( sigID , "median" , "PFS" , "COX" , meta_res )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
res_meta = rbind( res_meta , out )

##################################################
## Meta-analysis of the COX models (PFS) 10TMB
cox_dicho_10TMB_pfs = cox_dicho_10TMB_pfs[ as.numeric( as.character( cox_dicho_10TMB_pfs$SE ) ) <= 10 & !is.na( cox_dicho_10TMB_pfs$Pval ) , ]


#######################
## Global meta-analysis
meta_res = c( "ALL" , get_Meta_HR( data=cox_dicho_10TMB_pfs , signature = "ALL" ) )


#######################################
## Subgroup meta-analysis :  Tumor Type
for( j in 1:length(tumorID) ){
	data = cox_dicho_10TMB_pfs[ cox_dicho_10TMB_pfs$Primary %in% tumorID[j] , ]
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
	data = cox_dicho_10TMB_pfs[ cox_dicho_10TMB_pfs$Sequencing %in% sequencingID[j] , ]
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

out = cbind( sigID , "10TMB" , "PFS" , "COX" , cox_dicho_10TMB_pfs )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
res_indi = rbind( res_indi , out )

out = cbind( sigID , "10TMB" , "PFS" , "COX" , meta_res )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
res_meta = rbind( res_meta , out )

########################################################################################################################
########################################################################################################################

load( di_RData )

##################################################
## Meta-analysis of the DI models (OS) Continuous
di_os = di_os[ as.numeric( as.character( di_os$SE ) ) <= 10 & !is.na( di_os$Pval ) , ]


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

out = cbind( sigID , "continuous" , "OS" , "DI" , di_os )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
res_indi = rbind( res_indi , out )

out = cbind( sigID , "continuous" , "OS" , "DI" , meta_res )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
res_meta = rbind( res_meta , out )

##################################################
## Meta-analysis of the DI models (PFS) Continuous
di_pfs = di_pfs[ as.numeric( as.character( di_pfs$SE ) ) <= 10 & !is.na( di_pfs$Pval ) , ]


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

out = cbind( sigID , "continuous" , "PFS" , "DI" , di_pfs )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
res_indi = rbind( res_indi , out )

out = cbind( sigID , "continuous" , "PFS" , "DI" , meta_res )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
res_meta = rbind( res_meta , out )

######################################################################################################################################################
######################################################################################################################################################

load( log_RData )

##################################################
## Meta-analysis of the Log regression models (Response) Continuous
log_response = log_response[ as.numeric( as.character( log_response$SE ) ) <= 10 & !is.na( log_response$Pval ) , ]


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

out = cbind( sigID , "continuous" , "Response" , "Log_regression" , log_response )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
res_indi = rbind( res_indi , out )

out = cbind( sigID , "continuous" , "Response" , "Log_regression" , meta_res )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
res_meta = rbind( res_meta , out )


##################################################
##################################################
## Meta-analysis of the Log regression models (Response) Median
log_dicho_median_response = log_dicho_median_response[ as.numeric( as.character( log_dicho_median_response$SE ) ) <= 10 & !is.na( log_dicho_median_response$Pval ) , ]


#######################
## Global meta-analysis
meta_res = c( "ALL" , get_Meta_OR( data=log_dicho_median_response , signature = "ALL" ) )


#######################################
## Subgroup meta-analysis :  Tumor Type
for( j in 1:length(tumorID) ){
	data = log_dicho_median_response[ log_dicho_median_response$Primary %in% tumorID[j] , ]
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
	data = log_dicho_median_response[ log_dicho_median_response$Sequencing %in% sequencingID[j] , ]
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

out = cbind( sigID , "median" , "Response" , "Log_regression" , log_dicho_median_response )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
res_indi = rbind( res_indi , out )

out = cbind( sigID , "median" , "Response" , "Log_regression" , meta_res )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
res_meta = rbind( res_meta , out )

##################################################
##################################################
## Meta-analysis of the Log regression models (Response) 10TMB
log_dicho_10TMB_response = log_dicho_10TMB_response[ as.numeric( as.character( log_dicho_10TMB_response$SE ) ) <= 10 & !is.na( log_dicho_10TMB_response$Pval ) , ]


#######################
## Global meta-analysis
meta_res = c( "ALL" , get_Meta_OR( data=log_dicho_10TMB_response , signature = "ALL" ) )


#######################################
## Subgroup meta-analysis :  Tumor Type
for( j in 1:length(tumorID) ){
	data = log_dicho_10TMB_response[ log_dicho_10TMB_response$Primary %in% tumorID[j] , ]
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
	data = log_dicho_10TMB_response[ log_dicho_10TMB_response$Sequencing %in% sequencingID[j] , ]
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

out = cbind( sigID , "10TMB" , "Response" , "Log_regression" , log_dicho_10TMB_response )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "study" , "Primary" , "Sequencing" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" )
res_indi = rbind( res_indi , out )

out = cbind( sigID , "10TMB" , "Response" , "Log_regression" , meta_res )
colnames(out) = c( "signature" , "variable" , "outcome" , "model" , "Subgroup" , "Type" , "N" , "Effect_size" , "SE" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
res_meta = rbind( res_meta , out )


######################################################################################################################################################
######################################################################################################################################################

write.table( res_indi , file= "../results/Summary_Figure/Table/nsTMB/nsTMB_Individual.txt" , sep="\t" , row.names=FALSE , col.names=TRUE , quote=FALSE )
write.table( res_meta , file= "../results/Summary_Figure/Table/nsTMB/nsTMB_Meta_analysis.txt" , sep="\t" , row.names=FALSE , col.names=TRUE , quote=FALSE )
