
nsTMB_meta = read.table( file= "../results/Summary_Figure/Table/nsTMB/nsTMB_Meta_analysis.txt" , sep="\t" , header=TRUE , stringsAsFactor=FALSE)

######################################################################
######################################################################

feature = unique( nsTMB_meta$variable )

nsTMB_response = nsTMB_pfs = nsTMB_os = NULL

for( i in 1:length( feature ) ){

	################################################
	################################################
	## Response
	
	dat = nsTMB_meta[ nsTMB_meta$model %in% c( "COX" , "Log_regression" ) & nsTMB_meta$Subgroup %in% "ALL" & nsTMB_meta$variable %in% feature[i] & nsTMB_meta$outcome %in% "Response" , ]

	nsTMB_response = rbind( nsTMB_response , dat[ , c( "Effect_size" , "X95CI_low" , "X95CI_high" , "Pval" , "I2" , "Pval_I2" ) ] )

	################################################
	################################################
	## PFS

	dat = nsTMB_meta[ nsTMB_meta$model %in% c( "COX" , "Log_regression" ) & nsTMB_meta$Subgroup %in% "ALL" & nsTMB_meta$variable %in% feature[i] & nsTMB_meta$outcome %in% "PFS" , ]

	nsTMB_pfs = rbind( nsTMB_pfs , dat[ , c( "Effect_size" , "X95CI_low" , "X95CI_high" , "Pval" , "I2" , "Pval_I2" ) ] )


	################################################
	################################################
	## OS

	dat = nsTMB_meta[ nsTMB_meta$model %in% c( "COX" , "Log_regression" ) & nsTMB_meta$Subgroup %in% "ALL" & nsTMB_meta$variable %in% feature[i] & nsTMB_meta$outcome %in% "OS" , ]

	nsTMB_os = rbind( nsTMB_os , dat[ , c( "Effect_size" , "X95CI_low" , "X95CI_high" , "Pval" , "I2" , "Pval_I2" ) ] )


}

rownames( nsTMB_response ) = rownames( nsTMB_pfs ) = rownames( nsTMB_os ) = feature

######################################################################
######################################################################

library(forestplot)

pdf( "../results/Summary_Figure/Genomic/Forestplot_Meta-analysis_nsTMB.pdf" , height= 3, width= 7.5 , bg="transparent" , onefile=FALSE)
	res = rbind( nsTMB_response , nsTMB_pfs , nsTMB_os )
	res = res[ c( "continuous" , "median" , "10TMB" , "continuous1" , "median1" , "10TMB1" , "continuous2" , "median2" , "10TMB2" ) , ]
	rownames(res) = c( "continuous.res" , "median.res" , "10TMB.res" , "continuous.pfs" , "median.pfs" , "10TMB.pfs" , "continuous.os" , "median.os" , "10TMB.os" )

	res$Effect_size = round( as.numeric( as.character( res$Effect_size )) , 2 )
	res$X95CI_low = round( as.numeric( as.character( res$X95CI_low )) , 2 )
	res$X95CI_high = round( as.numeric( as.character( res$X95CI_high )) , 2 )
	res$Pval = round( as.numeric( as.character( res$Pval )) , 3 )
	res$I2 = round( as.numeric( as.character( res$I2 )) , 2 )
	res$Pval_I2 = round( as.numeric( as.character( res$Pval_I2 )) , 3 )

	result <- cbind( rownames( res ) ,
			res$Effect_size ,
			paste( "[" , res$X95CI_low , "; " , res$X95CI_high , "]" , sep="" ) ,
			res[ , c( "Pval" , "I2" , "Pval_I2" ) ] )
		
	result <- as.matrix(rbind(
				c("Variable" , "Effect_size" , "95% CI" , "P-value" , "Heterogeneity" , "Het. P-value"),
				result ))

	m.ns <- NA
	m.s <- NA
	l.ns <- NA
	l.s <- NA
	u.ns <- NA
	u.s <- NA
	for(i in 1:nrow(res) ){
		if( res$Pval[i] <= .05 ){
			m.ns <- c( m.ns , NA )
			m.s <- c( m.s , res$Effect_size[i] )
			l.ns <- c( l.ns , NA )
			l.s <- c( l.s , res$X95CI_low[i] )
			u.ns <- c( u.ns , NA )
			u.s <- c( u.s , res$X95CI_high[i] )
		} else {
			m.s <- c( m.s , NA )
			m.ns <- c( m.ns , res$Effect_size[i] )
			l.s <- c( l.s , NA )
			l.ns <- c( l.ns , res$X95CI_low[i] )
			u.s <- c( u.s , NA)
			u.ns <- c( u.ns , res$X95CI_high[i] )
		}
	}

	xlim = round( c( min( res$Effect_size , na.rm=TRUE)  - .2 , 0 , ( max( abs(res$Effect_size) , na.rm=TRUE) ) + .2 ) , 1 )

	forestplot(
		result,
		new_page = TRUE,
		colgap = unit( 2 , "mm") ,
		mean=cbind( m.ns , m.s ) ,
		lower=cbind( l.ns , l.s ) ,
		upper=cbind( u.ns , u.s ) ,
		is.summary= c( TRUE , rep( FALSE , nrow(res) ) ),
		hrzl_lines = list("2" = gpar(col = "#273746") ),
    	zero=0,
    	xlog=FALSE,
    	xlab="Combined LogOR (95% CI)",
    	boxsize=0.21,
    	xticks= xlim ,
		col=fpColors(
    		box=c("black","#D20000"),
    		lines=c("black","#D20000"), 
    		summary=c("black","#D20000")
    	)
	)	
dev.off()