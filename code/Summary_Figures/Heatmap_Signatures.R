
signature = read.table( file= "../data/signatures/signature_INFO.txt" , sep="\t" , header=TRUE , stringsAsFactor=FALSE)
signature$Signature <- as.character( signature$Signature )
signature$method <- as.character( signature$method )
signature$association <- as.character( signature$association )

association = signature$association
names( association ) = signature$Signature

#####################################################################
#####################################################################
## Continous Response

res = NULL

#################
# Gene Signature
sig = signature$Signature

for( i in 1:length(sig)){
	file = paste( "../results/Signature/" , sig[ i ] , "/Overall/Response/Response_" , sig[ i ] , "_Continous_LogReg.RData" , sep="" )
	if( file.exists( file ) ){
		load( file )
		res = rbind( res , c( meta_res , "sig" ) )
	}
}

##################################
##################################
res = as.data.frame( cbind( res , NA ) )
colnames(res) = c( "feature" , "coef" , "se_coef" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" , "type" , "Padj_BH" )

res$feature = as.character( res$feature )
res$type = as.character( res$type )
res$coef = round( as.numeric( as.character( res$coef ) ) , 2 )
res$se_coef = round( ifelse( as.numeric( as.character( res$se_coef ) ) >= 0.3 , 0.3 , as.numeric( as.character( res$se_coef ) ) ) , 2 )
res$CI95_low = round( as.numeric( as.character( res$CI95_low ) ) , 2 )
res$CI95_high = round( as.numeric( as.character( res$CI95_high ) ) , 2 )
res$Pval = as.numeric( as.character( res$Pval ) )

res$Padj_BH = round( p.adjust( res$Pval , method="BH" ) , 2)

response = res

#####################################################################
#####################################################################
## Continous PFS

res = NULL

#################
# Gene Signature
sig = signature$Signature

for( i in 1:length(sig)){
	file = paste( "../results/Signature/" , sig[ i ] , "/Overall/PFS/PFS_" , sig[ i ] , "_Continous_Cox.RData" , sep="" )
	if( file.exists( file ) ){
		load( file )
		res = rbind( res , c( meta_res , "sig" ) )
	}
}

##################################
##################################
res = as.data.frame( cbind( res , NA ) )
colnames(res) = c( "feature" , "coef" , "se_coef" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" , "type" , "Padj_BH" )

res$feature = as.character( res$feature )
res$type = as.character( res$type )
res$coef = round( as.numeric( as.character( res$coef ) ) , 2 )
res$se_coef = round( ifelse( as.numeric( as.character( res$se_coef ) ) >= 0.3 , 0.3 , as.numeric( as.character( res$se_coef ) ) ) , 2 )
res$CI95_low = round( as.numeric( as.character( res$CI95_low ) ) , 2 )
res$CI95_high = round( as.numeric( as.character( res$CI95_high ) ) , 2 )
res$Pval = as.numeric( as.character( res$Pval ) )

res$Padj_BH = round( p.adjust( res$Pval , method="BH" ) , 2)


pfs = res

#####################################################################
#####################################################################
## Continous OS

res = NULL
#################
# Gene Signature
sig = signature$Signature

for( i in 1:length(sig)){
	file = paste( "../results/Signature/" , sig[ i ] , "/Overall/OS/OS_" , sig[ i ] , "_Continous_Cox.RData" , sep="" )
	if( file.exists( file ) ){
		load( file )
		res = rbind( res , c( meta_res , "sig" ) )
	}
}

##################################
##################################
res = as.data.frame( cbind( res , NA ) )
colnames(res) = c( "feature" , "coef" , "se_coef" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" , "type" , "Padj_BH" )

res$feature = as.character( res$feature )
res$type = as.character( res$type )
res$coef = round( as.numeric( as.character( res$coef ) ) , 2 )
res$se_coef = round( ifelse( as.numeric( as.character( res$se_coef ) ) >= 0.3 , 0.3 , as.numeric( as.character( res$se_coef ) ) ) , 2 )
res$CI95_low = round( as.numeric( as.character( res$CI95_low ) ) , 2 )
res$CI95_high = round( as.numeric( as.character( res$CI95_high ) ) , 2 )
res$Pval = as.numeric( as.character( res$Pval ) )

res$Padj_BH = round( p.adjust( res$Pval , method="BH" ) , 2)

os = res


#####################################################################
#####################################################################
## Generate Heatmap

library(pheatmap)
library(RColorBrewer)

data = cbind( response$coef , pfs$coef , os$coef )
rownames(data) = response$feature
colnames(data) = c( "Response" , "PFS" , "OS" )

pval = cbind( response$Pval , pfs$Pval , os$Pval )
rownames(pval) = response$feature
colnames(pval) = c( "Response" , "PFS" , "OS" )
pval[ is.na(pval) ] = 1

padj = matrix( p.adjust( as.vector( as.matrix(pval) ) , method='fdr') , ncol=3 )
rownames(padj) = response$feature
colnames(padj) = c( "Response" , "PFS" , "OS" )
padj[ is.na(padj) ] = 1

annotation_col = data.frame(  Signature_Type= factor( association[ rownames(data) ] ) , 
					OS_Sig = factor( ifelse( round( padj[ , 'OS' ] , 2 ) <= .05 , 'FDR' , ifelse( round( pval[ , "OS" ] , 2 ) <= 0.05 , "Pvalue.Only" ,  'NS' ) ) ) ,  
					PFS_Sig = factor( ifelse( round( padj[ , 'PFS' ] , 2 ) <= .05 , 'FDR' , ifelse( round( pval[ , "PFS" ] , 2 ) <= 0.05 , "Pvalue.Only" , 'NS' ) ) ) ,
					Response_Sig = factor( ifelse( round( padj[ , 'Response' ] , 2 ) <= .05 , 'FDR' , ifelse( round( pval[ , "Response" ] , 2 ) <= 0.05 , "Pvalue.Only" , 'NS' ) ) )
					)
		
	
rownames(annotation_col) = rownames(data)

ann_colors = list(  Signature_Type = c( sensitive = "#4caf50" , resistance = "#8e44ad" ) ,
					Response_Sig = c( FDR = "#222021", Pvalue.Only = "#828282", NS = "white" ) ,
					PFS_Sig = c( FDR = "#222021", Pvalue.Only = "#828282", NS = "white" ),
					OS_Sig = c( FDR = "#222021", Pvalue.Only = "#828282", NS = "white" ) )

neg = seq( round( min( data , na.rm=TRUE ) , 1 ) , 0 , by=.05 )
neg = neg[ -length(neg)]
pos = seq( 0 , round( max( data , na.rm=TRUE ) , 1 ) , by=.05 )

col = c( colorRampPalette( rev( brewer.pal(4, "Blues") ) )( length(neg) ) ,
		colorRampPalette( brewer.pal(4, "OrRd") )( length(pos) ) 
	)


pheatmap( t( data[ order( rowSums(data)) , ] ) , cluster_rows=FALSE , cluster_cols=FALSE , scale="none" , annotation_col = annotation_col, annotation_colors = ann_colors, 
    col = col , breaks = c( neg , pos ) , na_col="white" , border_color="#424242", 
    #display_numbers = data , number_color="black" , fontsize_number= 10 , 
    filename="../results/Summary_Figure/Heatmap/Heatmap_Summary.pdf", height=2.8, width=9.5 )  

pheatmap( t( data[ order( rowSums(data)) , ] ) , cluster_rows=FALSE , cluster_cols=FALSE , annotation_col = annotation_col, annotation_colors = ann_colors, 
    col = col , breaks = c( neg , pos ) , na_col="white" , border_color="#424242",
    #display_numbers = data , number_color="black" , fontsize_number= 10 , 
    filename="../results/Summary_Figure/Heatmap/Heatmap_Summary_legend.pdf", height=4, width=9.5 )  



