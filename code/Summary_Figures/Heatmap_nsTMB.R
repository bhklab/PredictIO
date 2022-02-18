
nsTMB_indi = read.table( file= "../results/Summary_Figure/Table/nsTMB/nsTMB_Individual.txt" , sep="\t" , header=TRUE , stringsAsFactor=FALSE)

sequencing = unique( cbind( nsTMB_indi$study , nsTMB_indi$Sequencing ) )
colnames(sequencing) = c( "study" , "sequencing" )
rownames(sequencing) =  sequencing[ , "study" ]

######################################################################
######################################################################

feature = unique( nsTMB_indi$variable )
study = sort( unique( paste( nsTMB_indi$study , nsTMB_indi$Primary , sep = "__" ) ) )

mat_response = mat_pfs = mat_os = matrix( NA , nrow = length( feature ) , ncol = length( study ) )
colnames( mat_response ) = colnames( mat_pfs ) = colnames( mat_os ) = study
rownames( mat_response ) = rownames( mat_pfs ) = rownames( mat_os ) = feature

pmat_response = pmat_pfs = pmat_os = matrix( NA , nrow = length( feature ) , ncol = length( study ) )
colnames( pmat_response ) = colnames( pmat_pfs ) = colnames( pmat_os ) = study
rownames( pmat_response ) = rownames( pmat_pfs ) = rownames( pmat_os ) = feature


for( i in 1:ncol( mat_response ) ){

	################################################
	################################################
	## Response
	
	cohort = unlist( strsplit( colnames(mat_response)[i] , "__" , fixed = TRUE ) )[1]
	primary = unlist( strsplit( colnames(mat_response)[i] , "__" , fixed = TRUE ) )[2]

	dat = nsTMB_indi[ nsTMB_indi$model %in% c( "COX" , "Log_regression" ) & nsTMB_indi$study %in% cohort & nsTMB_indi$Primary %in% primary & nsTMB_indi$outcome %in% "Response" , ]

	mat_response[ dat$variable , i ] = dat$Effect_size
	pmat_response[ dat$variable , i ] = dat$Pval

	################################################
	################################################
	## PFS
	
	cohort = unlist( strsplit( colnames(mat_response)[i] , "__" , fixed = TRUE ) )[1]
	primary = unlist( strsplit( colnames(mat_response)[i] , "__" , fixed = TRUE ) )[2]

	dat = nsTMB_indi[ nsTMB_indi$model %in% c( "COX" , "Log_regression" ) & nsTMB_indi$study %in% cohort & nsTMB_indi$Primary %in% primary & nsTMB_indi$outcome %in% "PFS" , ]

	mat_pfs[ dat$variable , i ] = dat$Effect_size
	pmat_pfs[ dat$variable , i ] = dat$Pval

	################################################
	################################################
	## OS
	
	cohort = unlist( strsplit( colnames(mat_response)[i] , "__" , fixed = TRUE ) )[1]
	primary = unlist( strsplit( colnames(mat_response)[i] , "__" , fixed = TRUE ) )[2]

	dat = nsTMB_indi[ nsTMB_indi$model %in% c( "COX" , "Log_regression" ) & nsTMB_indi$study %in% cohort & nsTMB_indi$Primary %in% primary & nsTMB_indi$outcome %in% "OS" , ]

	mat_os[ dat$variable , i ] = dat$Effect_size
	pmat_os[ dat$variable , i ] = dat$Pval

}

######################################################################
######################################################################


library(pheatmap)
library(RColorBrewer)


################################################
################################################

annotation_col = data.frame(  
					Primary = factor( sapply( colnames( mat_response ) , function(x){ unlist( strsplit( x , "__" , fixed = TRUE ))[2] } ) ) ,
					Sequencing = factor( sequencing[ sapply( colnames( mat_response ) , function(x){ unlist( strsplit( x , "__" , fixed = TRUE ))[1] } ) , "sequencing" ] ) 
				)
rownames(annotation_col) = colnames( mat_response ) 
annotation_col = annotation_col[ order( annotation_col$Primary , annotation_col$Sequencing ) , ]

primary = table( annotation_col$Primary )
primary = names( primary[ rev( order( primary ) ) ] )
ann = NULL
for( i in 1:length(primary) ){
	ann = rbind( ann , annotation_col[ annotation_col$Primary %in% primary[i] , ] )
}
annotation_col = ann

col = sort( c( brewer.pal( 9 , "Set1" ) , brewer.pal( 4 , "Accent" ) ) )

ann_colors = list(  Primary = c( 
						Bladder = "#377EB8" , 
						Brain = "#4DAF4A" , 
						Breast = "#F781BF" , 
						Colon = "#984EA3" , 
						Esophagus = "#999999" , 
						Eye = "#A65628" , 
						HNC = "#BEAED4" , 
						Kidney = "#E41A1C" , 
						Lung = "#7FC97F" , 
						Melanoma = "#FDC086" , 
						Stomach = "#FF7F00" , 
						Unknown = "#3c3c3c" , 
						Ureteral = "#FFFF99"
						) ,
					Sequencing = c( TGS = "#283747", WES = "#bfc9ca" )
			)



################################################
################################################

neg = seq( -8 , 0 , by=.1 )
neg = neg[ -length(neg) ]
pos = seq( 0 , 6 , by=.1 )

col = c( colorRampPalette( rev( brewer.pal(4, "Blues") ) )( length(neg) ) ,
		colorRampPalette( brewer.pal(4, "OrRd") )( length(pos) ) 
	)

pmat_response[ pmat_response > 0.05 ] = "" 
pmat_response[ pmat_response <= 0.05 & !( pmat_response %in% c( NA , "" ) ) ] = "*" 
pmat_response[ is.na( pmat_response) ] = "" 

mat_response = mat_response[ , rownames( annotation_col ) ]
pmat_response = pmat_response[ , rownames( annotation_col ) ]

studies = sapply( colnames( mat_response ) , function( x ){ unlist( strsplit( x , "__" , fixed = TRUE ) )[ 1 ] } )
pheatmap( mat_response , cluster_rows=FALSE , cluster_cols=FALSE , scale="none", labels_col = studies ,
    col = col , breaks = c( neg , pos ) , na_col="white" , border_color="#424242" ,
    display_numbers =  pmat_response , annotation_col = annotation_col, annotation_colors = ann_colors ,
    filename="../results/Summary_Figure/Genomic/Heatmap_Genomic_nsTMB_Response_legend.pdf", height = 5 , width = 8 )  

studies = sapply( colnames( mat_response ) , function( x ){ unlist( strsplit( x , "__" , fixed = TRUE ) )[ 1 ] } )
pheatmap( mat_response , cluster_rows=FALSE , cluster_cols=FALSE , scale="none", labels_col = studies ,
    col = col , breaks = c( neg , pos ) , na_col="white" , border_color="#424242" ,
    display_numbers =  pmat_response , annotation_col = annotation_col, annotation_colors = ann_colors ,
    filename="../results/Summary_Figure/Genomic/Heatmap_Genomic_nsTMB_Response.pdf", height = 1.8 , width = 8 )  


################################################
################################################

neg = seq( -8 , 0 , by=.1 )
neg = neg[ -length(neg) ]
pos = seq( 0 , 6 , by=.1 )

col = c( colorRampPalette( rev( brewer.pal(4, "Blues") ) )( length(neg) ) ,
		colorRampPalette( brewer.pal(4, "OrRd") )( length(pos) ) 
	)

pmat_pfs[ pmat_pfs > 0.05 ] = "" 
pmat_pfs[ pmat_pfs <= 0.05 & !( pmat_pfs %in% c( NA , "" ) ) ] = "*" 
pmat_pfs[ is.na( pmat_pfs) ] = "" 

mat_pfs = mat_pfs[ , rownames( annotation_col ) ]
pmat_pfs = pmat_pfs[ , rownames( annotation_col ) ]


studies = sapply( colnames( mat_pfs ) , function( x ){ unlist( strsplit( x , "__" , fixed = TRUE ) )[ 1 ] } )
pheatmap( mat_pfs , cluster_rows=FALSE , cluster_cols=FALSE , scale="none", labels_col = studies ,
    col = col , breaks = c( neg , pos ) , na_col="white" , border_color="#424242" ,
    display_numbers =  pmat_pfs , annotation_col = annotation_col, annotation_colors = ann_colors ,
    filename="../results/Summary_Figure/Genomic/Heatmap_Genomic_nsTMB_PFS.pdf", height = 1.8 , width = 8 )  


################################################
################################################

neg = seq( -8 , 0 , by=.1 )
neg = neg[ -length(neg) ]
pos = seq( 0 , 6 , by=.1 )

col = c( colorRampPalette( rev( brewer.pal(4, "Blues") ) )( length(neg) ) ,
		colorRampPalette( brewer.pal(4, "OrRd") )( length(pos) ) 
	)

pmat_os[ pmat_os > 0.05 ] = "" 
pmat_os[ pmat_os <= 0.05 & !( pmat_os %in% c( NA , "" ) ) ] = "*" 
pmat_os[ is.na( pmat_os) ] = "" 

mat_os = mat_os[ , rownames( annotation_col ) ]
pmat_os = pmat_os[ , rownames( annotation_col ) ]

studies = sapply( colnames( mat_os ) , function( x ){ unlist( strsplit( x , "__" , fixed = TRUE ) )[ 1 ] } )
pheatmap( mat_os , cluster_rows=FALSE , cluster_cols=FALSE , scale="none" , labels_col = studies ,
    col = col , breaks = c( neg , pos ) , na_col="white" , border_color="#424242" ,
    display_numbers =  pmat_os , annotation_col = annotation_col, annotation_colors = ann_colors ,
    filename="../results/Summary_Figure/Genomic/Heatmap_Genomic_nsTMB_OS.pdf", height = 1.8 , width = 8 )  

