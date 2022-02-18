########################################################################################################################
########################################################################################################################

source("meta/meta-Analysis.R")

########################################################################################################################
########################################################################################################################

cox_RData= paste( "../results/nsTMB/nsTMB_COX_result.RData"  , sep="")
di_RData= paste( "../results/nsTMB/nsTMB_DI_result.RData" , sep="")
log_RData= paste( "../results/nsTMB/nsTMB_LogReg_result.RData"  , sep="")
sigID= "nsTMB"

height_cox_1 = 7
width_cox_1 = 7
height_cox_2 = 11
width_cox_2 = 7
height_cox_3 = 9
width_cox_3 = 7

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


################################################################################
################################################################################
## Meta-analysis of the COX models (OS) Continous
cox_os = cox_os[ as.numeric( as.character( cox_os$SE)) < 10 & !is.na(cox_os$SE) , ]
dim(cox_os)
sum(as.numeric(as.character(cox_os$N)))

cancer <- Get_Cancer( cancer=cox_os )
seq <- Get_Seq( seq=cox_os )
Get_Cox_Forestplot( data = cox_os , cancer = cancer , seq = seq , prefix= paste( "../results/" , sigID , sep="") , 
					label= paste( "OS_" , sigID , "_Continous_Cox" , sep="" ) , dir="OS" , 
					height_1 = 7 , width_1 = 7 , height_2 = 11 , width_2 = 7 , 
					height_3 = 9 , width_3 = 7 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 


################################################################################
################################################################################
## Meta-analysis of the COX models (OS) High vs Low (Median)
cox_dicho_median_os = cox_dicho_median_os[ as.numeric( as.character( cox_dicho_median_os$SE)) < 10 & !is.na(cox_dicho_median_os$SE) , ]
dim(cox_dicho_median_os)
sum(as.numeric(as.character(cox_dicho_median_os$N)))

cancer <- Get_Cancer( cancer=cox_dicho_median_os )
seq <- Get_Seq( seq=cox_dicho_median_os )
Get_Cox_Forestplot( data = cox_dicho_median_os , cancer = cancer , seq = seq , prefix=paste( "../results/" , sigID , sep="") , 
					label= paste( "OS_" , sigID , "_Dicho_Median_Cox" , sep="" ) , dir="OS" , 
					height_1 = 7 , width_1 = 7 , height_2 = 11 , width_2 = 7 , 
					height_3 = 9 , width_3 = 7 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 

## Meta-analysis of the COX models (OS) High vs Low (10TMB)
cox_dicho_10TMB_os = cox_dicho_10TMB_os[ as.numeric( as.character( cox_dicho_10TMB_os$SE)) < 10 & !is.na(cox_dicho_10TMB_os$SE) , ]
dim(cox_dicho_10TMB_os)
sum(as.numeric(as.character(cox_dicho_10TMB_os$N)))

cancer <- Get_Cancer( cancer=cox_dicho_10TMB_os )
seq <- Get_Seq( seq=cox_dicho_10TMB_os )
Get_Cox_Forestplot( data = cox_dicho_10TMB_os , cancer = cancer , seq = seq , prefix=paste( "../results/" , sigID , sep="") , 
					label= paste( "OS_" , sigID , "_Dicho_10TMB_Cox" , sep="" ) , dir="OS" , 
					height_1 = 6 , width_1 = 7 , height_2 = 8 , width_2 = 7 , 
					height_3 = 7 , width_3 = 7 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 

################################################################################
################################################################################
## Meta-analysis of the COX models (PFS) Continous
cox_pfs = cox_pfs[ as.numeric( as.character( cox_pfs$SE)) < 10 & !is.na(cox_pfs$SE) , ]
dim(cox_pfs)
sum(as.numeric(as.character(cox_pfs$N)))
cancer <- Get_Cancer( cancer=cox_pfs )
seq <- Get_Seq( seq=cox_pfs )
Get_Cox_Forestplot( data = cox_pfs , cancer = cancer , seq = seq , prefix= paste( "../results/" , sigID , sep="") , 
					label= paste( "PFS_" , sigID , "_Continous_Cox" , sep="" ) , dir="PFS" , 
					height_1 = 5 , width_1 = 7 , height_2 = 7 , width_2 = 7 , 
					height_3 = 6 , width_3 = 7 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 


################################################################################
################################################################################
## Meta-analysis of the COX models (PFS) High vs Low (Median)
cox_dicho_median_pfs = cox_dicho_median_pfs[ as.numeric( as.character( cox_dicho_median_pfs$SE)) < 10 & !is.na(cox_dicho_median_pfs$SE) , ]
dim(cox_dicho_median_pfs)
sum(as.numeric(as.character(cox_dicho_median_pfs$N)))
cancer <- Get_Cancer( cancer=cox_dicho_median_pfs )
seq <- Get_Seq( seq=cox_dicho_median_pfs )
Get_Cox_Forestplot( data = cox_dicho_median_pfs , cancer = cancer , seq = seq , prefix= paste( "../results/" , sigID , sep="") , 
					label= paste( "PFS_" , sigID , "_Dicho_Median_Cox" , sep="" ) , dir="PFS" , 
					height_1 = 4 , width_1 = 7 , height_2 = 6 , width_2 = 7 , 
					height_3 = 5.5 , width_3 = 7 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 


## Meta-analysis of the COX models (PFS) High vs Low (10TMB)
cox_dicho_10TMB_pfs = cox_dicho_10TMB_pfs[ as.numeric( as.character( cox_dicho_10TMB_pfs$SE)) < 10 & !is.na(cox_dicho_10TMB_pfs$SE) , ]
dim(cox_dicho_10TMB_pfs)
sum(as.numeric(as.character(cox_dicho_10TMB_pfs$N)))
cancer <- Get_Cancer( cancer=cox_dicho_10TMB_pfs )
seq <- Get_Seq( seq=cox_dicho_10TMB_pfs )
Get_Cox_Forestplot( data = cox_dicho_10TMB_pfs , cancer = cancer , seq = seq , prefix= paste( "../results/" , sigID , sep="") , 
					label= paste( "PFS_" , sigID , "_Dicho_10TMB_Cox" , sep="" ) , dir="PFS" , 
					height_1 = 4 , width_1 = 7 , height_2 = 5 , width_2 = 7 , 
					height_3 = 5 , width_3 = 7 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 


################################################################################
################################################################################

load( di_RData )

################################################################################
################################################################################
## Meta-analysis of the CI models (OS) Continous
di_os = di_os[ as.numeric( as.character( di_os$SE)) < 10 & !is.na(di_os$SE) , ]
dim(di_os)
sum(as.numeric(as.character(di_os$N)))
cancer <- Get_Cancer( cancer=di_os )
seq <- Get_Seq( seq=di_os )
Get_DI_Forestplot( data = di_os , cancer = cancer , seq = seq , prefix= paste( "../results/" , sigID , sep="") , 
					label= paste( "OS_" , sigID , "_Continous_DI" , sep="" ) , dir="OS" , 
					height_1 = 7 , width_1 = 7 , height_2 = 11 , width_2 = 7 , 
					height_3 = 8.5 , width_3 = 7 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 


################################################################################
################################################################################
## Meta-analysis of the CI models (PFS) Continous
di_pfs = di_pfs[ as.numeric( as.character( di_pfs$SE)) < 10 & !is.na(di_pfs$SE) , ]
dim(di_pfs)
sum(as.numeric(as.character(di_pfs$N)))
cancer <- Get_Cancer( cancer=di_pfs )
seq <- Get_Seq( seq=di_pfs )
Get_DI_Forestplot( data = di_pfs , cancer = cancer , seq = seq , prefix= paste( "../results/" , sigID , sep="") , 
					label= paste( "PFS_" , sigID , "_Continous_DI" , sep="" ) , dir="PFS" , 
					height_1 = 5 , width_1 = 7 , height_2 = 7 , width_2 = 7 , 
					height_3 = 6 , width_3 = 7 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 



################################################################################
################################################################################

load( log_RData )

################################################################################
################################################################################
## Meta-analysis of the Log Regression models (Response) Continous
log_response = log_response[ as.numeric( as.character( log_response$SE)) < 10 & !is.na(log_response$SE) , ]
dim(log_response)
sum(as.numeric(as.character(log_response$N)))
cancer <- Get_Cancer( cancer=log_response )
seq <- Get_Seq( seq=log_response )
Get_LogReg_Forestplot( data = log_response , cancer = cancer , seq = seq , prefix= paste( "../results/" , sigID , sep="") , 
					label= paste( "Response_" , sigID , "_Continous_LogReg" , sep="" ) , dir="Response" , 
					height_1 = 6 , width_1 = 7 , height_2 = 8 , width_2 = 7 , 
					height_3 = 7 , width_3 = 7 , height_4 = height_funnel_4 , width_4 = width_funnel_4 )  


################################################################################
################################################################################
## Meta-analysis of the Log Regression models (Response) High vs Low (Median)
log_dicho_median_response = log_dicho_median_response[ as.numeric( as.character( log_dicho_median_response$SE)) < 10 & !is.na(log_dicho_median_response$SE) , ]
dim(log_dicho_median_response)
sum(as.numeric(as.character(log_dicho_median_response$N)))
cancer <- Get_Cancer( cancer=log_dicho_median_response )
seq <- Get_Seq( seq=log_dicho_median_response )
Get_LogReg_Forestplot( data = log_dicho_median_response , cancer = cancer , seq = seq , prefix= paste( "../results/" , sigID , sep="") , 
					label= paste( "Response_" , sigID , "_Dicho_Median_LogReg" , sep="" ) , dir="Response" , 
					height_1 = 5 , width_1 = 7 , height_2 = 6.5 , width_2 = 7 , 
					height_3 = 6.5 , width_3 = 7 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 

################################################################################
################################################################################
## Meta-analysis of the Log Regression models (Response) High vs Low (10TMB)
log_dicho_10TMB_response = log_dicho_10TMB_response[ as.numeric( as.character( log_dicho_10TMB_response$SE)) < 10 & !is.na(log_dicho_10TMB_response$SE) , ]
dim(log_dicho_10TMB_response)
sum(as.numeric(as.character(log_dicho_10TMB_response$N)))
cancer <- Get_Cancer( cancer=log_dicho_10TMB_response )
seq <- Get_Seq( seq=log_dicho_10TMB_response )
Get_LogReg_Forestplot( data = log_dicho_10TMB_response , cancer = cancer , seq = seq , prefix= paste( "../results/" , sigID , sep="") , 
					label= paste( "Response_" , sigID , "_Dicho_10TMB_LogReg" , sep="" ) , dir="Response" , 
					height_1 = 4 , width_1 = 7 , height_2 = 5.5 , width_2 = 7 , 
					height_3 = 5.5 , width_3 = 7 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 



