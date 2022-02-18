##########################################################################################################################################################
##########################################################################################################################################################

library(survcomp)
library(genefu)
library(GSVA)
library(reshape2)

#################################################
#################################################


#Remove rows if count is < zero in 50% of sample
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0) ))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}

#################################################
#################################################

Get_GSVA = function( signature_name , expr ){

	geneSig = rep( NA , ncol(expr) )
	if( ncol(expr) > 20 ){
		sig = as.character( read.csv( file= paste( "../data/signatures/" , signature_name , ".csv" , sep="" ) , sep=';' , header=TRUE , stringsAsFactor=FALSE)[,1] )

		if( ifelse( is.null( nrow( expr[ rownames(expr) %in% sig , ]) ) , 1 , nrow( expr[ rownames(expr) %in% sig , ] ) ) / length( sig ) > 0.8 & ncol(expr) > 20 ){

			geneSig <- as.numeric( gsva( expr , list(sig) , verbose=FALSE ) )
			names( geneSig ) = colnames(expr)
		}
	}
	geneSig
}

Get_Mean = function( signature_name , expr ){

	geneSig = rep( NA , ncol(expr) )
	if( ncol(expr) > 20 ){	
		sig = read.csv( file= paste( "../data/signatures/" , signature_name , ".csv" , sep="" ) , sep=';' , header=TRUE , stringsAsFactor=FALSE)
		rownames(sig) = sig$gene
		sig$coef <- as.numeric( as.character( sig$coef ) )
		sig$gene <- as.character( sig$gene )

		
		if( ifelse( is.null( nrow( expr[ rownames(expr) %in% sig$gene , ]) ) , 1 , nrow( expr[ rownames(expr) %in% sig$gene , ] ) ) / nrow( sig ) > 0.8 & ncol(expr) > 20 ){

			gene = intersect( rownames(expr) , sig$gene) 
			s = sig[ gene , ]
			geneSig <- apply( expr[ gene , ] , 2 , function(x) ( sum( ( x * s$coef ) ) /  nrow( s ) ) )
			names( geneSig ) = colnames(expr)
		}
	}
	geneSig
}

Get_COX_IS = function( expr ){

	geneSig = rep( NA , ncol(expr) )
	if( ncol(expr) > 20 ){
		sig = read.csv( file= paste( "../data/signatures/COX_IS.csv" , sep="" ) , sep=';' , header=TRUE , stringsAsFactor=FALSE)
		sig$gene <- as.character( sig$gene )
		sig$coef <- as.character( sig$coef )

		if( ifelse( is.null( nrow( expr[ rownames(expr) %in% sig$gene , ]) ) , 1 , nrow( expr[ rownames(expr) %in% sig$gene , ] ) ) / nrow( sig ) > 0.8 & ncol(expr) > 20 ){
			pos = apply( expr[ rownames( expr ) %in% sig[ sig$coef %in% 1 , ]$gene , ] , 2 , function(x){ ( sum( x ) /  length( x ) ) } )
			neg = apply( expr[ rownames( expr ) %in% sig[ sig$coef %in% -1 , ]$gene , ] , 2 , function(x){ ( sum( x ) /  length( x ) ) } )

			geneSig <- as.numeric( pos / neg )
			names( geneSig ) = colnames( expr )
			
			if( length( geneSig[ !is.nan( geneSig ) ] ) < 20){
				geneSig = NULL
			}
		}
	}
	geneSig

}

Get_TIDE = function( expr , primary ){

	tide = rep( NA , ncol(expr) )
	if( ncol(expr) > 20 ){
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

		if( file.exists( tide_output ) ) {
			tide = read.table( tide_output , header=TRUE, sep="\t" , stringsAsFactors=FALSE ) 
			names = tide[,1]

			tide = tide$TIDE
			tide = as.numeric( scale( tide ) )

			names( tide ) = names
		}
	}

	tide
}

Get_IPRES = function( expr ){

	geneSig = rep( NA , ncol(expr) )
	if( ncol(expr) > 20 ){
		sig = read.csv( file= paste( "../data/signatures/IPRES.csv" , sep="" ) , sep=';' , header=TRUE , stringsAsFactor=FALSE)
		gene.list <- as.character( sig$gene.list )
		Geneset <- as.character( sig$Geneset )

		sig = list()
		for( i in 1:length(gene.list)){
				sig[ Geneset[i] ] =  strsplit( gene.list[i] , "," , fixed=TRUE )
		}

		geneSig <- apply( gsva( expr , sig , verbose=FALSE ) , 2 , mean , na.rm=TRUE )
		names(geneSig) = colnames( expr )
	}
	geneSig

}

Get_IPS = function( expr ){
	####################################################
	##
	##   This R-script can be used to calculate Immunophenoscore (IPS) and generate Immunophenogram from "EXPR.txt" and "IPS_genes.txt"
	##   (C) ICBI, Medical University of Innsbruck, Biocenter, Division of Bioinformatics
	##   Version 1.0 08.07.2016
	##   Needs packages ggplot2,grid,gridExtra
	##
	####################################################

	AZ = rep( NA , ncol(expr) )
	if( ncol(expr) & nrow(expr)>10000 ){ 

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
		}
	}

	AZ
}

##########################################################################################################
##########################################################################################################

getCOR <- function( data , dataset , primary ){
	
	remove <- rem(data)
	if( length(remove) ){
		data <- data[-remove,]
	}		
		
	IO_signature = NULL

	#################################################
	#################################################
	signature = read.table( file= "../data/signatures/signature_INFO.txt" , sep="\t" , header=TRUE , stringsAsFactor=FALSE)
	signature$Signature <- as.character( signature$Signature )
	signature$method <- as.character( signature$method )

	for(i in 1:nrow( signature ) ){

		if( signature$method[i] %in% "GSVA" ){

			print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
			IO_signature = cbind( IO_signature , Get_GSVA( signature_name= signature$Signature[i] , expr= data ) )

		}
		if( signature$method[i] %in% "weighted_mean" ){	

			print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
			IO_signature = cbind( IO_signature , Get_Mean( signature_name= signature$Signature[i] , expr= data ) )

		}
		if( signature$method[i] %in% "specific" ){

			if( signature$Signature[i] %in% "IPS" ){

				print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
				IO_signature = cbind( IO_signature , Get_IPS( expr= data ) )

			}		
			if( signature$Signature[i] %in% "COX_IS" ){

				print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
				IO_signature = cbind( IO_signature , Get_COX_IS( expr= data ) )

			}
			if( signature$Signature[i] %in% "IPRES" ){

				print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
				IO_signature = cbind( IO_signature , Get_IPRES( expr= data ) )

			}
			if( signature$Signature[i] %in% "TIDE" ){

				print( paste( signature$Signature[i] , "|" , signature$method[i] , sep=" " ) )
				IO_signature = cbind( IO_signature , Get_TIDE( expr= data , primary= primary ) )

			}	
		}
	}
	colnames(IO_signature) = signature$Signature
	####################################################################################
	####################################################################################
	## Correlation matrix
	cor = cor( IO_signature , method="s" )
	cor = cbind( dataset , primary , nrow( IO_signature ) , melt( cor ) )
	cor = cor[ cor$Var1 != cor$Var2 , ]
	cor
}

##########################################################################################################
##########################################################################################################

load("../data/process_data/ICB_exp_filtered.RData")

expr
study = names(expr)

cor = NULL

for( i in 1:length(study)){

	tumor = names( table( phenoData( expr[[i]] )$primary )[ table( phenoData( expr[[i]] )$primary ) >= 20 ] )

	if( length(tumor) > 0 ){
		for( j in 1:length(tumor)){

			print( paste( "############## " , study[i] , ": " , tumor[j] , " ##############" , sep="" ) )

			data = exprs(expr[[i]])[ ,  phenoData(expr[[i]])$primary %in% tumor[j] & phenoData(expr[[i]])$rna %in% c( "fpkm" , "tpm" ) ]

			corM = NULL
			if( ncol( data ) ){ corM <- getCOR( data = data , dataset = study[i] , primary= tumor[j] ) }
			if( !is.null( corM ) ){ cor = rbind( cor , corM ) }
		}
	}
}	
colnames( cor ) = c( "dataset" , "primary" , "N" , "var1" , "var2" , "value" )

cor = cbind( cor , apply( cor , 1 , function(x){ paste( x[ "var1" ] , x[ "var2" ] , sep="_" ) } ) )
colnames( cor ) = c( "dataset" , "primary" , "N" , "var1" , "var2" , "value" , "id" )
cor = cor[ , c( "dataset" , "primary" , "N" , "var1" , "var2" , "id" , "value" ) ]


save( cor , file= "../results/Summary_Figure/CorrPlot/CorrPlot.RData" )
################################################
################################################

library(corrplot)
load( "../results/Summary_Figure/CorrPlot/CorrPlot.RData" )

signature = c( "B_cell_Helmink" , "APM_Thompson" , "M1" , "Response_Ivy" , "STAT1" , "Inflammatory" , "TIS" , "IRG_Ayers" , 
	"TLS" , "CD8_SF" , "T_cell_inflamed" , "PDL1" , "IFNG" , "peri_Tcell" , "ADO" , "CD8_Jiang" , "APM_Wang" , 
	"B_cell_Budczies" , "CYT" , "IPRES" , "Blood_Response" , "Myeloid_DC" , "IRG_Yang" , "IPS" , 
	"KDM5A" , "TIDE" , "TGFB_Mariathasan" , "NonResponse_Ivy" , "EMT_Thompson" , "ImmuneCells" , "LRRC15_CAF" , 
	"COX_IS" , "CRMA" , "IMPRES" , "MPS" , "C_ECM" , "T_cell_exclusion" , "PTEN_MITF" )


mat= matrix( NA , nrow= length(signature) , ncol= length(signature) )
colnames(mat) = rownames(mat) = signature

for( i in 1:length( signature ) ){
	for( j in 1:length( signature ) ){
		if( i != j ){
			c = cor[ cor$var1 %in% signature[i] & cor$var2 %in% signature[j] , ]

			mat[ signature[i] , signature[j] ] = round( median( c$value[ !is.na( c$value ) ] , na.rm=TRUE ) , 2 )

		} else{
			mat[ signature[i] , signature[j] ] = 1
		}
	}
}

pdf( "../results/Summary_Figure/CorrPlot/CorrPlot_Overall.pdf" , height=10 , width=10 , bg="transparent" )
	corrplot(mat , method = "color" , type = "lower" , lower.col = "black", number.cex = .4 ,  tl.col = "black" , tl.srt = 45 , diag = FALSE )
dev.off()

################################################
################################################

tumor = c( "Melanoma" , "Lung" , "Kidney")

for( k in 1:length(tumor) ){

	mat= matrix( NA , nrow= length(signature) , ncol= length(signature) )
	colnames(mat) = rownames(mat) = signature

	for( i in 1:length( signature ) ){
		for( j in 1:length( signature ) ){
			if( i != j ){
				c = cor[ cor$var1 %in% signature[i] & cor$var2 %in% signature[j] & cor$primary %in% tumor[k] , ]

				mat[ signature[i] , signature[j] ] = round( median( c$value[ !is.na( c$value ) ] , na.rm=TRUE ) , 2 )

			} else{
				mat[ signature[i] , signature[j] ] = 1
			}
		}
	}

	pdf( paste( "../results/Summary_Figure/CorrPlot/CorrPlot_", tumor[k] , ".pdf" , sep="" ) , height=10 , width=10 , bg="transparent" )
		corrplot(mat , method = "color" , type = "lower" , lower.col = "black", number.cex = .4 , tl.col = "black", tl.srt = 45 , diag = FALSE )
	dev.off()
}

