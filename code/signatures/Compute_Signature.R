##########################################################################################################################################################
##########################################################################################################################################################

library(survcomp)
library(genefu)
library(data.table)
library(Biobase)


load("../data/process_data/ICB_exp_filtered.RData")
expr

source( "signatures/Get_Meta-analysis.R")
source( "signatures/Get_Signature.R")
source( "signatures/TIDE_signature.R")
source( "signatures/IPS_signature.R")
source( "signatures/IPRES_signature.R")
source( "signatures/COX_IS_signature.R")

##########################################################################################################################################################
##########################################################################################################################################################
signature = read.table( file= "../data/signatures/signature_INFO.txt" , sep="\t" , header=TRUE , stringsAsFactor=FALSE)
signature$Signature <- as.character( signature$Signature )
signature$method <- as.character( signature$method )

for(i in 1:nrow( signature ) ){

	if( signature$method[i] %in% "GSVA" ){

		print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )

		GSVA_Signature( signature_name= signature$Signature[i] , expr= expr )
		Get_Meta_analysis( SignatureID= signature$Signature[i] )

	}
	if( signature$method[i] %in% "weighted_mean" ){

		print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
		
		Mean_Signature( signature_name= signature$Signature[i] , expr= expr )
		Get_Meta_analysis( SignatureID= signature$Signature[i] )

	}
	if( signature$method[i] %in% "specific" ){

		if( signature$Signature[i] %in% "IPS" ){
		
			print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
		
			Get_IPS_signature( expr= expr )
			Get_Meta_analysis( SignatureID= "IPS" )

		}		
		if( signature$Signature[i] %in% "COX_IS" ){
		
			print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
		
			Get_COX_IS_signature( expr= expr )
			Get_Meta_analysis( SignatureID= "COX_IS" )

		}
		if( signature$Signature[i] %in% "IPRES" ){

			print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )

			Get_IPRES_signature( expr= expr )
			Get_Meta_analysis( SignatureID= "IPRES" )

		}
		if( signature$Signature[i] %in% "TIDE" ){

			print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )

			Get_TIDE_signature( expr= expr )
			Get_Meta_analysis( SignatureID= "TIDE" )

		}	
	}
}
##########################################################################################################################################################
##########################################################################################################################################################

source( "signatures/PerTumorType_Signature_Meta-analysis.R" )


##########################################################################################################################################################
##########################################################################################################################################################

