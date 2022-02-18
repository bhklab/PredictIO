
file = "../results/denovo_Single_Gene"

if( dir.exists( file ) ) {
	unlink( file , recursive = TRUE )
}

dir.create( file )

#Remove rows if count is < zero in 50% of sample
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0) ))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}

##################################################################################
##################################################################################

load("../data/process_data/ICB_exp_filtered.RData")

library(genefu)

expr
study = names(expr)

genes = tumor = NULL
for( i in 1:length(study)){

	primary = names( table( phenoData( expr[[i]] )$primary )[ table( phenoData( expr[[i]] )$primary ) >= 20 ] )
	tumor = c( tumor , sapply( primary , function(x){ paste( study[i] , x , sep = "__" ) } ) )

	data = exprs( expr[[ study[i] ]] )
	
	remove <- rem(data)
	if( length(remove) ){
		data <- data[-remove,]
	}
	genes = unique( sort( c( genes , rownames( data ) ) ) )
}

numberG = matrix( 0 , nrow=length(genes) , ncol=length(study) )
colnames(numberG) = study
rownames(numberG) = genes

for( i in 1:length(study)){

	data = exprs( expr[[ study[i] ]] )
	
	remove <- rem(data)
	if( length(remove) ){
		data <- data[-remove,]
	}

	numberG[ rownames(data) , study[i] ] = 1
}

genes = names( rowSums( numberG )[ rowSums( numberG ) %in% c( 12 , 13 , 14 ) ] )

##################################################################################
##################################################################################

pval = matrix( NA , nrow=length(genes) , ncol=length(tumor) )
colnames(pval) = tumor
rownames(pval) = genes

coef = matrix( NA , nrow=length(genes) , ncol=length(tumor) )
colnames(coef) = tumor
rownames(coef) = genes

se = matrix( NA , nrow=length(genes) , ncol=length(tumor) )
colnames(se) = tumor
rownames(se) = genes


for( i in 1:length(study)){

	primary = names( table( phenoData( expr[[i]] )$primary )[ table( phenoData( expr[[i]] )$primary ) >= 20 ] )

	if( length(primary) > 0 ){
		for( j in 1:length(primary)){

			print( paste( "study:" , study[i] , " | Tumor:" , primary[j] ) )

			data = exprs(expr[[i]])[ ,  phenoData(expr[[i]])$primary %in% primary[j] ]
			remove <- rem(data)
			if( length(remove) ){
				data <- data[-remove,]
			}

			data = as.matrix( data[ rownames(data) %in% genes , ] )
			if( nrow(data) & !( primary %in% "Lymph_node" ) ){
				for( k in 1:nrow(data) ){
						
					x = ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% primary[j] ] %in% "R" , 0 , 
						ifelse( phenoData(expr[[i]])$response[ phenoData(expr[[i]])$primary %in% primary[j] ] %in% "NR" , 1 , NA ) )

					g = as.numeric( scale( data[k , ] ) )

					fit = glm( x ~ g , family=binomial( link="logit" ) )

					pval[ rownames(data)[k] , paste( study[i] , primary[j] , sep="__" ) ] = summary(fit)$coefficients[ 2 , 4 ]
					coef[ rownames(data)[k] , paste( study[i] , primary[j] , sep="__" ) ] = round( summary(fit)$coefficients[ 2 , 1  ] , 2 )
					se[ rownames(data)[k] , paste( study[i] , primary[j] , sep="__" ) ] = round( summary(fit)$coefficients[ 2 , 2 ] , 2 )
				}
			}
		}
	}
}

fdr <- matrix( p.adjust( pval , method= "fdr" ) , ncol= ncol( pval ) , nrow= nrow( pval ) , dimnames= dimnames( pval ) )

save( pval , coef , se , fdr , file= "../results/denovo_Single_Gene/Single_Gene_LogReg_Response.RData" ) 


####################################################################################
####################################################################################

source('meta/Get_HR.R')

pval = matrix( NA , nrow=length(genes) , ncol=length(tumor) )
colnames(pval) = tumor
rownames(pval) = genes

hr = matrix( NA , nrow=length(genes) , ncol=length(tumor) )
colnames(hr) = tumor
rownames(hr) = genes

se = matrix( NA , nrow=length(genes) , ncol=length(tumor) )
colnames(se) = tumor
rownames(se) = genes


for( i in 1:length(study)){

	primary = names( table( phenoData( expr[[i]] )$primary )[ table( phenoData( expr[[i]] )$primary ) >= 20 ] )

	if( length(primary) > 0 ){
		for( j in 1:length(primary)){

			print( paste( "study:" , study[i] , " | Tumor:" , primary[j] ) )

			data = exprs(expr[[i]])[ , phenoData(expr[[i]])$primary %in% primary[j] ]
			remove <- rem(data)
			if( length(remove) ){
				data <- data[-remove,]
			}

			data = as.matrix( data[ rownames(data) %in% genes , ] )
			if( nrow(data) & !( primary %in% "Lymph_node" ) & length( phenoData(expr[[i]])$os[ !is.na( phenoData(expr[[i]])$os ) & phenoData(expr[[i]])$primary %in% primary[j] ] ) >= 20 ){
				for( k in 1:nrow(data) ){
						
					print( paste( "study:" , study[i] , " | Tumor:" , primary[j] , " | " , rownames(data)[k] , " | " , k , "/" , nrow( data ) ) )

					g = as.numeric( scale( data[k , ] ) )
					names( g ) = colnames( data )

					cox = Get_HR_continous( surv = phenoData(expr[[i]])$os[ phenoData(expr[[i]])$primary %in% primary[j] ] , 
											time = phenoData(expr[[i]])$t.os[ phenoData(expr[[i]])$primary %in% primary[j] ] ,
										time_censor=36 , variable= g )

					pval[ rownames(data)[k] , paste( study[i] , primary[j] , sep="__" ) ] = cox[5]
					hr[ rownames(data)[k] , paste( study[i] , primary[j] , sep="__" ) ] = cox[1]
					se[ rownames(data)[k] , paste( study[i] , primary[j] , sep="__" ) ] = cox[2]
				}
			}
		}
	}
}

fdr <- matrix( p.adjust( pval , method= "fdr" ) , ncol= ncol( pval ) , nrow= nrow( pval ) , dimnames= dimnames( pval ) )

save( pval , hr , se , fdr , file= "../results/denovo_Single_Gene/Single_Gene_COX_OS.RData" ) 

############################################################################################################
############################################################################################################

source('meta/Get_HR.R')

pval = matrix( NA , nrow=length(genes) , ncol=length(tumor) )
colnames(pval) = tumor
rownames(pval) = genes

hr = matrix( NA , nrow=length(genes) , ncol=length(tumor) )
colnames(hr) = tumor
rownames(hr) = genes

se = matrix( NA , nrow=length(genes) , ncol=length(tumor) )
colnames(se) = tumor
rownames(se) = genes


for( i in 1:length(study)){

	primary = names( table( phenoData( expr[[i]] )$primary )[ table( phenoData( expr[[i]] )$primary ) >= 20 ] )

	if( length(primary) > 0 ){
		for( j in 1:length(primary)){

			print( paste( "study:" , study[i] , " | Tumor:" , primary[j] ) )

			data = exprs(expr[[i]])[ , phenoData(expr[[i]])$primary %in% primary[j] ]
			remove <- rem(data)
			if( length(remove) ){
				data <- data[-remove,]
			}

			data = as.matrix( data[ rownames(data) %in% genes , ] )
			if( nrow(data) & !( primary %in% "Lymph_node" ) & length( phenoData(expr[[i]])$pfs[ !is.na( phenoData(expr[[i]])$pfs ) & phenoData(expr[[i]])$primary %in% primary[j] ] ) >= 20 ){
				for( k in 1:nrow(data) ){
						
					print( paste( "study:" , study[i] , " | Tumor:" , primary[j] , " | " , rownames(data)[k] , " | " , k , "/" , nrow( data ) ) )
					
					g = as.numeric( scale( data[k , ] ) )
					names( g ) = colnames( data )

					cox = Get_HR_continous( surv = phenoData(expr[[i]])$pfs[ phenoData(expr[[i]])$primary %in% primary[j] ] , 
											time = phenoData(expr[[i]])$t.pfs[ phenoData(expr[[i]])$primary %in% primary[j] ] ,
										time_censor=24 , variable= g )

					pval[ rownames(data)[k] , paste( study[i] , primary[j] , sep="__" ) ] = cox[5]
					hr[ rownames(data)[k] , paste( study[i] , primary[j] , sep="__" ) ] = cox[1]
					se[ rownames(data)[k] , paste( study[i] , primary[j] , sep="__" ) ] = cox[2]
				}
			}
		}
	}
}

fdr <- matrix( p.adjust( pval , method= "fdr" ) , ncol= ncol( pval ) , nrow= nrow( pval ) , dimnames= dimnames( pval ) )

save( pval , hr , se , fdr , file= "../results/denovo_Single_Gene/Single_Gene_COX_PFS.RData" ) 


############################################################################################################
############################################################################################################

