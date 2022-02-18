########################################################################################################################
########################################################################################################################

source("meta/Meta-Analysis.R")

########################################################################################################################
########################################################################################################################

get_Directory <- function( dir ){

	tumorID = c( "Melanoma" , "Lung" , "Kidney" )
	file = paste( "../results/Per_TumorType/Signature/" , dir , sep="" )
	
	if( file.exists( file ) ) {
		unlink( file , recursive = TRUE )
	}

	dir.create( file , recursive = TRUE )

	for( i in 1:length(tumorID) ){
		dir.create( paste( file , "/" , tumorID[i] , sep="" ) )

		dir.create( paste( file , "/" , tumorID[i] , "/KMPlot" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/KMPlot/OS" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/KMPlot/PFS" , sep="" ) )

		dir.create( paste( file , "/" , tumorID[i] , "/Funnel" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/Funnel/OS" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/Funnel/PFS" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/Funnel/Response" , sep="" ) )
		
		dir.create( paste( file , "/" , tumorID[i] , "/Overall" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/Overall/OS" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/Overall/PFS" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/Overall/Response" , sep="" ) )
		
		dir.create( paste( file , "/" , tumorID[i] , "/PerCancer" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/PerCancer/OS" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/PerCancer/PFS" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/PerCancer/Response" , sep="" ) )
		
		dir.create( paste( file , "/" , tumorID[i] , "/PerSequencing" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/PerSequencing/OS" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/PerSequencing/PFS" , sep="" ) )
		dir.create( paste( file , "/" , tumorID[i] , "/PerSequencing/Response" , sep="" ) )
	}

}

########################################################################################################################
########################################################################################################################

signature = read.table( file= "../data/signatures/signature_INFO.txt" , sep="\t" , header=TRUE , stringsAsFactor=FALSE)
signature$Signature <- as.character( signature$Signature )
signature$method <- as.character( signature$method )
signature$association <- as.character( signature$association )

########################################################################################################################
########################################################################################################################

SignatureID = signature$Signature

tumorID = c( "Melanoma" , "Lung" , "Kidney" )

for(k in 1:length(SignatureID)){

	get_Directory( dir= SignatureID[k] )

	for(l in 1:length(tumorID)){

		cox_RData= paste( "../results/Signature/" , SignatureID[k] , "/" , SignatureID[k] , "_COX_result.RData"  , sep="")
		di_RData= paste( "../results/Signature/" , SignatureID[k] , "/" , SignatureID[k] , "_DI_result.RData" , sep="")
		log_RData= paste( "../results/Signature/" , SignatureID[k] , "/" , SignatureID[k] , "_LogReg_result.RData"  , sep="")
		sigID= SignatureID[k]


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
		cox_os = cox_os[ cox_os$Primary %in% tumorID[l] , ]


		if( nrow( cox_os ) >= 3 ){
			cancer <- Get_Cancer( cancer=cox_os )
			seq <- Get_Seq( seq=cox_os )
			Get_Cox_Forestplot( data = cox_os , cancer = cancer , seq = seq , prefix= paste( "../results/Per_TumorType/Signature/" , sigID , "/" , tumorID[l]  , sep="") , 
								label= paste( "OS_" , sigID , "_Continous_Cox" , sep="" ) , dir="OS" , 
								height_1 = height_cox_1 , width_1 = width_cox_1 , height_2 = height_cox_2 , width_2 = width_cox_2 , 
								height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
		}
		## Meta-analysis of the COX models (OS) High vs Low
		cox_dicho_os = cox_dicho_os[ cox_dicho_os$Primary %in% tumorID[l] , ]

		if( nrow( cox_dicho_os ) >= 3 ){
			cancer <- Get_Cancer( cancer=cox_dicho_os )
			seq <- Get_Seq( seq=cox_dicho_os )
			Get_Cox_Forestplot( data = cox_dicho_os , cancer = cancer , seq = seq , prefix=paste( "../results/Per_TumorType/Signature/" , sigID , "/" , tumorID[l]  , sep="") , 
								label= paste( "OS_" , sigID , "_Dicho_Cox" , sep="" ) , dir="OS" , 
								height_1 = height_cox_1 , width_1 = width_cox_1 , height_2 = height_cox_2 , width_2 = width_cox_2 , 
								height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
		}
		## Meta-analysis of the COX models (PFS) Continous
		cox_pfs = cox_pfs[ cox_pfs$Primary %in% tumorID[l] , ]

		if( nrow( cox_pfs ) >= 3 ){
			cancer <- Get_Cancer( cancer=cox_pfs )
			seq <- Get_Seq( seq=cox_pfs )
			Get_Cox_Forestplot( data = cox_pfs , cancer = cancer , seq = seq , prefix= paste( "../results/Per_TumorType/Signature/" , sigID , "/" , tumorID[l]  , sep="") , 
								label= paste( "PFS_" , sigID , "_Continous_Cox" , sep="" ) , dir="PFS" , 
								height_1 = height_cox_1 , width_1 = width_cox_1 , height_2 = height_cox_2 , width_2 = width_cox_2 , 
								height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
		}
		
		## Meta-analysis of the COX models (PFS) High vs Low
		cox_dicho_pfs = cox_dicho_pfs[ cox_dicho_pfs$Primary %in% tumorID[l] , ]
		if( nrow( cox_dicho_pfs ) >= 3 ){
			cancer <- Get_Cancer( cancer=cox_dicho_pfs )
			seq <- Get_Seq( seq=cox_dicho_pfs )
			Get_Cox_Forestplot( data = cox_dicho_pfs , cancer = cancer , seq = seq , prefix= paste( "../results/Per_TumorType/Signature/" , sigID , "/" , tumorID[l]  , sep="") , 
								label= paste( "PFS_" , sigID , "_Dicho_Cox" , sep="" ) , dir="PFS" , 
								height_1 = height_cox_1 , width_1 = width_cox_1 , height_2 = height_cox_2 , width_2 = width_cox_2 , 
								height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
		}
		
		load( di_RData )

		########################################
		## Meta-analysis of the DI models (OS) Continous
		di_os <- di_os[ di_os$Primary %in% tumorID[l] , ]
		if( nrow( di_os ) >= 3 ){
			cancer <- Get_Cancer( cancer=di_os )
			seq <- Get_Seq( seq=di_os )
			Get_DI_Forestplot( data = di_os , cancer = cancer , seq = seq , prefix= paste( "../results/Per_TumorType/Signature/" , sigID , "/" , tumorID[l]  , sep="") , 
								label= paste( "OS_" , sigID , "_Continous_DI" , sep="" ) , dir="OS" , 
								height_1 = height_DI_1 , width_1 = width_DI_1 , height_2 = height_DI_2 , width_2 = width_DI_2 , 
								height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
		}
		

		## Meta-analysis of the DI models (PFS) Continous
		di_pfs = di_pfs[ di_pfs$Primary %in% tumorID[l] , ]
		if( nrow( di_pfs ) >= 3 ){
			cancer <- Get_Cancer( cancer=di_pfs )
			seq <- Get_Seq( seq=di_pfs )
			Get_DI_Forestplot( data = di_pfs , cancer = cancer , seq = seq , prefix= paste( "../results/Per_TumorType/Signature/" , sigID , "/" , tumorID[l]  , sep="") , 
								label= paste( "PFS_" , sigID , "_Continous_DI" , sep="" ) , dir="PFS" , 
								height_1 = height_DI_1 , width_1 = width_DI_1 , height_2 = height_DI_2 , width_2 = width_DI_2 , 
								height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
		}


		load( log_RData )
		## Meta-analysis of the Log Regression models (Response) Continous
		log_response = log_response[ log_response$Primary %in% tumorID[l] , ]
		if( nrow( log_response ) >= 3 ){
			cancer <- Get_Cancer( cancer=log_response )
			seq <- Get_Seq( seq=log_response )
			Get_LogReg_Forestplot( data = log_response , cancer = cancer , seq = seq , prefix= paste( "../results/Per_TumorType/Signature/" , sigID , "/" , tumorID[l]  , sep="") , 
								label= paste( "Response_" , sigID , "_Continous_LogReg" , sep="" ) , dir="Response" , 
								height_1 = height_log_1 , width_1 = width_log_1 , height_2 = height_log_2 , width_2 = width_log_2 , 
								height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
		}

		## Meta-analysis of the Log Regression models (Response) High vs Low
		log_dicho_response = log_dicho_response[ log_dicho_response$Primary %in% tumorID[l] , ]
		if( nrow( log_dicho_response ) >= 3 ){
			cancer <- Get_Cancer( cancer=log_dicho_response )
			seq <- Get_Seq( seq=log_dicho_response )
			Get_LogReg_Forestplot( data = log_dicho_response , cancer = cancer , seq = seq , prefix= paste( "../results/Per_TumorType/Signature/" , sigID , "/" , tumorID[l]  , sep="") , 
								label= paste( "Response_" , sigID , "_Dicho_LogReg" , sep="" ) , dir="Response" , 
								height_1 = height_log_1 , width_1 = width_log_1 , height_2 = height_log_2 , width_2 = width_log_2 , 
								height_3 = height_cox_3 , width_3 = width_cox_3 , height_4 = height_funnel_4 , width_4 = width_funnel_4 ) 
		}
	}
}