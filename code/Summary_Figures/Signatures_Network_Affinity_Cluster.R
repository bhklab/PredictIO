
library(reshape2)
library(apcluster)
library(RColorBrewer)
library(ggfortify)

##################################################################
##################################################################

signature = read.table( file= "../data/signatures/signature_INFO.txt" , sep="\t" , header=TRUE , stringsAsFactor=FALSE)
signature$Signature <- as.character( signature$Signature )
signature$method <- as.character( signature$method )
signature$association <- as.character( signature$association )


association = signature$association
names( association ) = signature$Signature
method = signature$method
names( method ) = signature$Signature

##################################################################
##################################################################

signature = as.data.frame( t( read.csv( file= "../data/signatures/ALL_sig.csv" , sep=";" , header=FALSE , stringsAsFactor=FALSE) ) )

data = melt( signature, id.vars= "V1" )
data = data[ !( data$value %in% "" ) , ]
colnames( data ) = c( "signature" , "id" , "value" )

data$signature = as.character( data$signature )
data$value = as.character( data$value )

data = data[ order( data$signature ) , c( "signature" , "value" ) ]

##################################################################
##################################################################

signature = sort( unique( data$signature ) )

association = association[ signature ]
method = method[ signature ]

overlap = matrix( nrow = length( signature ) , ncol = length( signature ) , 0 )
colnames( overlap ) = rownames( overlap ) = signature

for( i in 1:length( signature ) ){
	for( j in 1:length( signature ) ){

		s1 = data[ data$signature %in% signature[i] , ]$value
		s2 = data[ data$signature %in% signature[j] , ]$value

		int = intersect( s1 , s2 )

		overlap[ i , j ] = length( int )

	}
}

pca=prcomp( overlap , scale = TRUE)
x1 = pca$x[ , 1:2]
apres = apcluster( negDistMat( r = 2 ) , x1 )


col = brewer.pal( n = 8 , name ="Dark2" )

pdf( "../results/Summary_Figure/Signature_Network/APcluster_Signature.pdf" , height = 6 , width = 9 , bg = "transparent" )
	plot( apres , x1 , xlab = "PC1 (18.62%)" , ylab = "PC2 (5.97%)" )
	text( x1, row.names(x1) , cex = 0.6 , pos = 4 , col = "black" ) 
dev.off()

cl = apres@clusters
cluster = NULL
for( i in 1:length( cl ) ){
	c = unlist( cl[ i ] )
	cluster = rbind( cluster , cbind( names( c ) , i ) )
}
cluster = as.data.frame( cluster )
colnames(cluster ) = c( "signature" , "cluster" )

cluster$signature = as.character( cluster$signature )
cluster$cluster = as.character( cluster$cluster )

save( cluster , file = "../results/Summary_Figure/Signature_Network/Cluster.RData" )


