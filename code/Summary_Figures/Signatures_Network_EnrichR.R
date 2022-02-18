######################################################
######################################################
load( "../results/Summary_Figure/Signature_Network/Cluster.RData" )
library(reshape2)

######################################################
######################################################
signature = as.data.frame( t( read.csv( file= "../data/signatures/ALL_sig.csv" , sep=";" , header=FALSE , stringsAsFactor=FALSE) ) )

data = melt( signature, id.vars= "V1" )
data = data[ !( data$value %in% "" ) , ]
colnames( data ) = c( "signature" , "id" , "value" )

data$signature = as.character( data$signature )
data$value = as.character( data$value )

data = data[ order( data$signature ) , c( "signature" , "value" ) ]

clust = sort( names( table( cluster$cluster )[ table( cluster$cluster ) > 1 ] ) )

######################################################
######################################################
library(enrichR)
setEnrichrSite("Enrichr") # Human genes
dbs <- c( "KEGG_2016" )

kegg = NULL
for( i in 1:length( clust ) ){

	sig = cluster[ cluster$cluster %in% clust[ i ] , ]$signature

	gene = sort( unique( data[ data$signature %in% sig , ]$value ) )
	enriched <- enrichr( gene , dbs )

	kegg = rbind( kegg , cbind( clust[ i ] , enriched[[1]][ 1:5 , 1 ] ) )
}

save( kegg , file = "../results/Summary_Figure/Signature_Network/KEGG_Cluster.RData" )



