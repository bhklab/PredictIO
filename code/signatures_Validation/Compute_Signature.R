##########################################################################################################################################################
##########################################################################################################################################################

library(survcomp)
library(genefu)
library(data.table)
library(Biobase)

source( "signatures_Validation/Get_Signature.R")
source( "signatures_Validation/TIDE_signature.R")
source( "signatures_Validation/IPS_signature.R")
source( "signatures_Validation/IPRES_signature.R")
source( "signatures_Validation/COX_IS_signature.R")

dir.create( "../results/signature_Validation/" )

get_Directory_cohort <- function( cohort ){
	file = paste( "../results/signature_Validation/" , cohort , sep="" )
	if( file.exists( file ) ) {
		unlink( file )
	}
	dir.create( file )
}

signature = read.table( file= "../data/signatures/signature_INFO.txt" , sep="\t" , header=TRUE , stringsAsFactor=FALSE)
signature$Signature <- as.character( signature$Signature )
signature$method <- as.character( signature$method )


cohort = c( "Kim" , "INSPIRE" , "Gide" , "Puch" , "VanDenEnde" , "Shiuan" )

for(z in 1:length( cohort ) ){

	get_Directory_cohort( cohort= cohort[z] )

	for(i in 1:nrow( signature ) ){

		if( signature$method[i] %in% "GSVA" ){

			print( paste( cohort[z] , "|" , signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )

			GSVA_Signature( signature_name= signature$Signature[i] , cohort= cohort[z] , sigType= signature$association[i] )
		}
		if( signature$method[i] %in% "weighted_mean" ){

			print( paste( cohort[z] , "|" , signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
			
			Mean_Signature( signature_name= signature$Signature[i] , cohort= cohort[z] , sigType= signature$association[i] )
		}
		if( signature$method[i] %in% "specific" ){

			if( signature$Signature[i] %in% "COX_IS" ){
			
				print( paste( cohort[z] , "|" , signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
			
				Get_COX_IS_signature( cohort= cohort[z] , sigType= signature$association[i] )
			}
			if( signature$Signature[i] %in% "IPS" ){
			
				print( paste( cohort[z] , "|" , signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
			
				Get_IPS_signature( cohort= cohort[z] , sigType= signature$association[i] )
			}
			if( signature$Signature[i] %in% "IPRES" ){

				print( paste( cohort[z] , "|" , signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )

				Get_IPRES_signature( cohort= cohort[z] , sigType= signature$association[i] )
			}
			if( signature$Signature[i] %in% "TIDE" ){

				print( paste( cohort[z] , "|" , signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )

				Get_TIDE_signature( cohort= cohort[z] , sigType= signature$association[i] )
			}
			
		}

	}
}