
########################################################################################################################
########################################################################################################################
source("meta/Meta-Analysis.R")


Get_Meta_analysis = function( SignatureID ){

	cox_RData= paste( "../results/Signature/" , SignatureID , "/" , SignatureID , "_COX_result.RData"  , sep="")
	di_RData= paste( "../results/Signature/" , SignatureID , "/" , SignatureID , "_DI_result.RData" , sep="")
	log_RData= paste( "../results/Signature/" , SignatureID , "/" , SignatureID , "_LogReg_result.RData"  , sep="")
	sigID= SignatureID

	height_cox_1 = 6
	width_cox_1 = 10
	height_cox_2 = 9
	height_cox_3 = 9
	width_cox_3 = 11
	width_cox_2 = 10

	height_DI_1 = 7
	width_DI_1 = 10
	height_DI_2 = 9
	width_DI_2 = 10
	height_DI_3 = 9
	width_DI_3 = 11

	height_log_1 = 6
	width_log_1 = 10
	height_log_2 = 9
	width_log_2 = 10
	height_log_3 = 9
	width_log_3 = 11

	height_funnel_4 = 5
	width_funnel_4 = 8

	########################################################################################################################
	########################################################################################################################
	Get_Cancer = function( cancer ){

		if( length( grep( "HR" , colnames( cancer ) ) ) ){
				cancer$HR <- as.numeric( as.character( cancer$HR ) )
				
		} else{
			if( length( grep( "coef" , colnames( cancer ) ) ) ){
					cancer$coef <- as.numeric( as.character( cancer$coef ) )
			} else{
			cancer$DI <- as.numeric( as.character( cancer$DI ) )
			}
		}

		cancer$study <- as.character( cancer$study )
		cancer$Sequencing <- as.character( cancer$Sequencing )
		cancer$Primary <- as.character( cancer$Primary )
		cancer$Pval <- as.numeric(as.character( cancer$Pval )) 
		cancer$SE <- as.numeric(as.character( cancer$SE )) 

		tab <- table( cancer$Primary)[ table(cancer$Primary) %in% c(1,2) ]

		cancer$study[ cancer$Primary %in% names(tab) ] <- paste( cancer$study[ cancer$Primary %in% names(tab) ] , ", " , 
															cancer$Primary[ cancer$Primary %in% names(tab) ] , sep= "" ) 
		cancer$study <- paste( cancer$study , ", n = " , cancer$N , sep= "" ) 
		cancer$Primary[ cancer$Primary %in% names(tab) ] <- "Other"

		cancer
	}

	Get_Seq = function( seq ){

		if( length( grep( "HR" , colnames( seq ) ) ) ){
				seq$HR <- as.numeric( as.character( seq$HR ) )
		} else{
			if( length( grep( "coef" , colnames( seq ) ) ) ){
					seq$coef <- as.numeric( as.character( seq$coef ) )
			} else{
				seq$DI <- as.numeric( as.character( seq$DI ) )
			}
		}
		
		seq$study <- as.character( seq$study )
		seq$Sequencing <- as.character( seq$Sequencing )
		seq$Primary <- as.character( seq$Primary )
		seq$SE <- as.numeric( as.character( seq$SE ) )
		seq$Pval <- as.numeric( as.character( seq$Pval ) ) 

		seq$study <- paste( seq$study , ", " , seq$Primary , ", n = " , seq$N , sep= "" ) 

		seq
	}

	########################################################################################################################
	########################################################################################################################

	load( cox_RData )

	########################################
	## Meta-analysis of the COX models (OS) Continous
	cox_os = cox_os[ cox_os$SE <= 10 & !is.na( cox_os$Pval ) , ]
	if( nrow( cox_os[ !is.na( cox_os[ , 1 ] ) , ] ) ){
		cancer <- Get_Cancer( cancer=cox_os )
		seq <- Get_Seq( seq=cox_os )
		Get_Cox_Forestplot( data = cox_os , cancer = cancer , seq = seq , prefix= paste( "../results/Signature/" , sigID , sep="") , 
							label= paste( "OS_" , sigID , "_Continous_Cox" , sep="" ) , dir="OS" , 
							height_1 = height_cox_1 , width_1 = width_cox_1 , height_2 = height_cox_2 , width_2 = width_cox_2 , 
							height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
	}

	## Meta-analysis of the COX models (OS) High vs Low
	cox_dicho_os = cox_dicho_os[ cox_dicho_os$SE <= 10 & !is.na( cox_dicho_os$Pval ) , ]
	if( nrow( cox_dicho_os[ !is.na( cox_dicho_os[ , 1 ] ) , ] ) ){
		cancer <- Get_Cancer( cancer=cox_dicho_os )
		seq <- Get_Seq( seq=cox_dicho_os )
		Get_Cox_Forestplot( data = cox_dicho_os , cancer = cancer , seq = seq , prefix=paste( "../results/Signature/" , sigID , sep="") , 
							label= paste( "OS_" , sigID , "_Dicho_Cox" , sep="" ) , dir="OS" , 
							height_1 = height_cox_1 , width_1 = width_cox_1 , height_2 = height_cox_2 , width_2 = width_cox_2 , 
							height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
	}
	## Meta-analysis of the COX models (PFS) Continous
	cox_pfs = cox_pfs[ cox_pfs$SE <= 10 & !is.na( cox_pfs$Pval ) , ]
	if( nrow( cox_pfs[ !is.na( cox_pfs[ , 1 ] ) , ] ) ){
		cancer <- Get_Cancer( cancer=cox_pfs )
		seq <- Get_Seq( seq=cox_pfs )
		Get_Cox_Forestplot( data = cox_pfs , cancer = cancer , seq = seq , prefix= paste( "../results/Signature/" , sigID , sep="") , 
							label= paste( "PFS_" , sigID , "_Continous_Cox" , sep="" ) , dir="PFS" , 
							height_1 = height_cox_1 , width_1 = width_cox_1 , height_2 = height_cox_2 , width_2 = width_cox_2 , 
							height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
	}
	## Meta-analysis of the COX models (PFS) High vs Low
	cox_dicho_pfs = cox_dicho_pfs[ cox_dicho_pfs$SE <= 10 & !is.na( cox_dicho_pfs$Pval ) , ]
	if( nrow( cox_dicho_pfs[ !is.na( cox_dicho_pfs[ , 1 ] ) , ] ) ){
		cancer <- Get_Cancer( cancer=cox_dicho_pfs )
		seq <- Get_Seq( seq=cox_dicho_pfs )
		Get_Cox_Forestplot( data = cox_dicho_pfs , cancer = cancer , seq = seq , prefix= paste( "../results/Signature/" , sigID , sep="") , 
							label= paste( "PFS_" , sigID , "_Dicho_Cox" , sep="" ) , dir="PFS" , 
							height_1 = height_cox_1 , width_1 = width_cox_1 , height_2 = height_cox_2 , width_2 = width_cox_2 , 
							height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
	}

	load( di_RData )

	########################################
	## Meta-analysis of the DI models (OS) Continous
	di_os = di_os[ di_os$SE <= 10 & !is.na( di_os$Pval ) , ]
	if( nrow( di_os[ !is.na( di_os[ , 1 ] ) , ] ) ){
		cancer <- Get_Cancer( cancer=di_os )
		seq <- Get_Seq( seq=di_os )
		Get_DI_Forestplot( data = di_os , cancer = cancer , seq = seq , prefix= paste( "../results/Signature/" , sigID , sep="") , 
							label= paste( "OS_" , sigID , "_Continous_DI" , sep="" ) , dir="OS" , 
							height_1 = height_DI_1 , width_1 = width_DI_1 , height_2 = height_DI_2 , width_2 = width_DI_2 , 
							height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
	}

	## Meta-analysis of the DI models (PFS) Continous
	di_pfs = di_pfs[ di_pfs$SE <= 10 & !is.na( di_pfs$Pval ) , ]
	if( nrow( di_pfs[ !is.na( di_pfs[ , 1 ] ) , ] ) ){
		cancer <- Get_Cancer( cancer=di_pfs )
		seq <- Get_Seq( seq=di_pfs )
		Get_DI_Forestplot( data = di_pfs , cancer = cancer , seq = seq , prefix= paste( "../results/Signature/" , sigID , sep="") , 
							label= paste( "PFS_" , sigID , "_Continous_DI" , sep="" ) , dir="PFS" , 
							height_1 = height_DI_1 , width_1 = width_DI_1 , height_2 = height_DI_2 , width_2 = width_DI_2 , 
							height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
	}


	load( log_RData )
	## Meta-analysis of the Log Regression models (Response) Continous
	log_response = log_response[ log_response$SE <= 10 & !is.na( log_response$Pval ) , ]
	if( nrow( log_response[ !is.na( log_response[ , 1 ] ) , ] ) ){
		cancer <- Get_Cancer( cancer=log_response )
		seq <- Get_Seq( seq=log_response )
		Get_LogReg_Forestplot( data = log_response , cancer = cancer , seq = seq , prefix= paste( "../results/Signature/" , sigID , sep="") , 
							label= paste( "Response_" , sigID , "_Continous_LogReg" , sep="" ) , dir="Response" , 
							height_1 = height_log_1 , width_1 = width_log_1 , height_2 = height_log_2 , width_2 = width_log_2 , 
							height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
	}

	## Meta-analysis of the Log Regression models (Response) High vs Low
	log_dicho_response = log_dicho_response[ log_dicho_response$SE <= 10 & !is.na( log_dicho_response$Pval ) , ]
	if( nrow( log_dicho_response[ !is.na( log_dicho_response[ , 1 ] ) , ] ) ){
		cancer <- Get_Cancer( cancer=log_dicho_response )
		seq <- Get_Seq( seq=log_dicho_response )
		Get_LogReg_Forestplot( data = log_dicho_response , cancer = cancer , seq = seq , prefix= paste( "../results/Signature/" , sigID , sep="") , 
							label= paste( "Response_" , sigID , "_Dicho_LogReg" , sep="" ) , dir="Response" , 
							height_1 = height_log_1 , width_1 = width_log_1 , height_2 = height_log_2 , width_2 = width_log_2 , 
							height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
	}

}
