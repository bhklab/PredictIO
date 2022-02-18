
source("Resistance_Gene_Validation/Compute_Outcome.R")
source("Resistance_Gene_Validation/Get_Meta-Analysis.R")

########################################################################################################################
########################################################################################################################

dir.create( "../results/Resistance_Gene_Validation" )

get_Directory <- function( dir , geneID ){
	file = paste( "../results/Resistance_Gene_Validation/" , geneID , sep="" )
	
	dir.create( file )
	dir.create( paste( file , "/meta_analysis" , sep = "" ) )

	file = paste( "../results/Resistance_Gene_Validation/" , geneID , "/" , dir , sep="" )

	if( file.exists( file ) ) {
		unlink( file )
	}

	dir.create( file )

}
########################################################################################################################
########################################################################################################################

geneID = c( "F2RL1" , "RBFOX2" )
cohort = c( "Shiuan" , "VanDenEnde" , "INSPIRE" , "Kim" , "Gide" , "Puch" )

for( i in 1:length( geneID ) ){
	for( z in 1:length( cohort ) ){
		get_Directory( dir = cohort[z] , geneID = geneID[i] )
		Get_Gene( exprID = cohort[z] , geneID = geneID[i] )
	}
	Get_Meta_Analysis( cohort = cohort , geneID = geneID[i] )
}

########################################################################################################################
########################################################################################################################

load( "../results/denovo_Single_Gene/Single_Gene_LogReg_Response.RData" ) 

response = NULL

for( i in 1:length( geneID ) ){

	data = as.data.frame( cbind( colnames( coef ) , colnames( coef ) , coef[ geneID[ i ] , ] , se[ geneID[ i ] , ] , pval[ geneID[ i ] , ] ) )
	colnames(data) = c( "study" , "primary" , "coef" , "se" , "pval" )
	data$study = sapply( as.character( data$study ) , function( x ){ unlist( strsplit( x , "__" , fixed = TRUE ) )[1] } )
	data$primary = sapply( as.character( data$primary ) , function( x ){ unlist( strsplit( x , "__" , fixed = TRUE ) )[2] } )
	data$coef = as.numeric(as.character( data$coef ))
	data$se = as.numeric(as.character( data$se ))
	data$pval = as.numeric(as.character( data$pval )) 
	data = data[ !is.na(data$coef) , ]
	data = data[ order( data$coef ) , ]
	
	data$study <- paste( data$study , ", " , data$primary , sep= "" ) 

	meta <- metagen( TE = coef,
					seTE = se,
					data = data,
					studlab = study ,
					fixed = FALSE ,
					random = TRUE ,
					control = list( maxiter = 10000 , stepadj=0.5 ) )

	######################################################################
	######################################################################

	response = rbind( response , c( meta$TE.fixed , meta$lower.fixed , meta$upper.fixed , meta$pval.random , meta$I2 , meta$pval.Q ) )

	######################################################################
	######################################################################

	pdf( paste( "../results/Resistance_Gene_Validation/" , geneID[i] , "/Forestplot_Response_" , geneID[i] , ".pdf" , sep = "" ) , height= 6 , width= 7  , bg="transparent" )

		m <- c( min( data$coef , na.rm=TRUE ) - .1 , max( abs(data$coef) , na.rm=TRUE) + .1 )

		forest( meta , 
            leftcols = c("studlab", "effect.ci" , "pval"),
			leftlabs= c( "Study" , "logOR [95%CI]" , "P-value" ) , 
   			xlab = "Estimated logOR",
			digits.se = 2,
   			colgap.forest=unit(10, "mm") ,
	      	plotwidth = unit( 30 , "mm") ,
	       	pooled.totals = TRUE,
	       	smlab = " ",
	       	comb.random =TRUE,
	       	comb.fixed = FALSE,
	       	text.fixed.w = FALSE,
		    layout = "JAMA",
		    print.I2.ci = TRUE,
		    print.Q = FALSE,
		    print.pval.Q = TRUE,
		    print.I2 = TRUE,
		    print.tau2 = FALSE,
		    resid.hetstat = FALSE,
	       	test.overall.random = TRUE,
	       	test.overall.fixed = FALSE,
	       	xlim = m , 
	       	col.square= "black" ,  
	       	col.study= "black" ,  
	       	col.square.lines = "black" ,
	       	col.diamond.random  = "#1565c0"  ,
	       	col.diamond.lines.random  ="#1565c0" ,
	       	col.by = "#1565c0" ,
		    addrow.subgroups=TRUE
		)
	dev.off()
}

########################################################################################################################
########################################################################################################################

geneID = c( "F2RL1" , "RBFOX2" )

load( "../results/denovo_Single_Gene/Single_Gene_COX_PFS.RData" ) 

pfs = NULL

for( i in 1:length( geneID ) ){

	data = as.data.frame( cbind( colnames( hr ) , colnames( hr ) , hr[ geneID[ i ] , ] , se[ geneID[ i ] , ] , pval[ geneID[ i ] , ] ) )
	colnames(data) = c( "study" , "primary" , "coef" , "se" , "pval" )
	data$study = sapply( as.character( data$study ) , function( x ){ unlist( strsplit( x , "__" , fixed = TRUE ) )[1] } )
	data$primary = sapply( as.character( data$primary ) , function( x ){ unlist( strsplit( x , "__" , fixed = TRUE ) )[2] } )
	data$coef = as.numeric(as.character( data$coef ))
	data$se = as.numeric(as.character( data$se ))
	data$pval = as.numeric(as.character( data$pval )) 
	data = data[ !is.na(data$coef) , ]
	data = data[ order( data$coef ) , ]
	
	data$study <- paste( data$study , ", " , data$primary , sep= "" ) 

	meta <- metagen( TE = coef,
					seTE = se,
					data = data,
					studlab = study ,
					fixed = FALSE ,
					random = TRUE ,
					control = list( maxiter = 10000 , stepadj=0.5 ) )

	######################################################################
	######################################################################

	pfs = rbind( pfs , c( meta$TE.fixed , meta$lower.fixed , meta$upper.fixed , meta$pval.random , meta$I2 , meta$pval.Q ) )

	######################################################################
	######################################################################

	pdf( paste( "../results/Resistance_Gene_Validation/" , geneID[i] , "/Forestplot_PFS_" , geneID[i] , ".pdf" , sep = "" ) , height= 6 , width= 7  , bg="transparent" )

		m <- c( min( data$coef , na.rm=TRUE ) - .1 , max( abs(data$coef) , na.rm=TRUE) + .1 )

		forest( meta , 
            leftcols = c("studlab", "effect.ci" , "pval"),
			leftlabs= c( "Study" , "logHR [95%CI]" , "P-value" ) , 
   			xlab = "Estimated logHR",
			digits.se = 2,
   			colgap.forest=unit(10, "mm") ,
	      	plotwidth = unit( 30 , "mm") ,
	       	pooled.totals = TRUE,
	       	smlab = " ",
	       	comb.random =TRUE,
	       	comb.fixed = FALSE,
	       	text.fixed.w = FALSE,
		    layout = "JAMA",
		    print.I2.ci = TRUE,
		    print.Q = FALSE,
		    print.pval.Q = TRUE,
		    print.I2 = TRUE,
		    print.tau2 = FALSE,
		    resid.hetstat = FALSE,
	       	test.overall.random = TRUE,
	       	test.overall.fixed = FALSE,
	       	xlim = m , 
	       	col.square= "black" ,  
	       	col.study= "black" ,  
	       	col.square.lines = "black" ,
	       	col.diamond.random  = "#1565c0"  ,
	       	col.diamond.lines.random  ="#1565c0" ,
	       	col.by = "#1565c0" ,
		    addrow.subgroups=TRUE
		)
	dev.off()
}
########################################################################################################################
########################################################################################################################

geneID = c( "F2RL1" , "RBFOX2" )

load( "../results/denovo_Single_Gene/Single_Gene_COX_OS.RData" ) 

os = NULL

for( i in 1:length( geneID ) ){

	data = as.data.frame( cbind( colnames( hr ) , colnames( hr ) , hr[ geneID[ i ] , ] , se[ geneID[ i ] , ] , pval[ geneID[ i ] , ] ) )
	colnames(data) = c( "study" , "primary" , "coef" , "se" , "pval" )
	data$study = sapply( as.character( data$study ) , function( x ){ unlist( strsplit( x , "__" , fixed = TRUE ) )[1] } )
	data$primary = sapply( as.character( data$primary ) , function( x ){ unlist( strsplit( x , "__" , fixed = TRUE ) )[2] } )
	data$coef = as.numeric(as.character( data$coef ))
	data$se = as.numeric(as.character( data$se ))
	data$pval = as.numeric(as.character( data$pval )) 
	data = data[ !is.na(data$coef) , ]
	data = data[ order( data$coef ) , ]
	
	data$study <- paste( data$study , ", " , data$primary , sep= "" ) 

	meta <- metagen( TE = coef,
					seTE = se,
					data = data,
					studlab = study ,
					fixed = FALSE ,
					random = TRUE ,
					control = list( maxiter = 10000 , stepadj=0.5 ) )

	######################################################################
	######################################################################

	os = rbind( os , c( meta$TE.fixed , meta$lower.fixed , meta$upper.fixed , meta$pval.random , meta$I2 , meta$pval.Q ) ) 

	######################################################################
	######################################################################

	pdf( paste( "../results/Resistance_Gene_Validation/" , geneID[i] , "/Forestplot_OS_" , geneID[i] , ".pdf" , sep = "" ) , height= 6 , width= 7  , bg="transparent" )

		m <- c( min( data$coef , na.rm=TRUE ) - .1 , max( abs(data$coef) , na.rm=TRUE) + .1 )

		forest( meta , 
            leftcols = c("studlab", "effect.ci" , "pval"),
			leftlabs= c( "Study" , "logHR [95%CI]" , "P-value" ) , 
   			xlab = "Estimated logHR",
			digits.se = 2,
   			colgap.forest=unit(10, "mm") ,
	      	plotwidth = unit( 30 , "mm") ,
	       	pooled.totals = TRUE,
	       	smlab = " ",
	       	comb.random =TRUE,
	       	comb.fixed = FALSE,
	       	text.fixed.w = FALSE,
		    layout = "JAMA",
		    print.I2.ci = TRUE,
		    print.Q = FALSE,
		    print.pval.Q = TRUE,
		    print.I2 = TRUE,
		    print.tau2 = FALSE,
		    resid.hetstat = FALSE,
	       	test.overall.random = TRUE,
	       	test.overall.fixed = FALSE,
	       	xlim = m , 
	       	col.square= "black" ,  
	       	col.study= "black" ,  
	       	col.square.lines = "black" ,
	       	col.diamond.random  = "#1565c0"  ,
	       	col.diamond.lines.random  ="#1565c0" ,
	       	col.by = "#1565c0" ,
		    addrow.subgroups=TRUE
		)
	dev.off()
}


########################################################################################################################
########################################################################################################################

rownames(response) = rownames(os) = rownames(pfs) = geneID
colnames(response) = colnames(os) = colnames(pfs) = c( "coef" , "CI95_low" , "CI95_high" , "pval" , "I2" , "I2_pval" ) 

geneID = c( "F2RL1" , "RBFOX2" )

library(forestplot)

for( i in 1:length( geneID ) ){

	pdf( paste( "../results/Resistance_Gene_Validation/" , geneID[i] , "/Forestplot_" , geneID[i] , ".pdf" , sep = "" ) , height= 2, width= 6 , bg="transparent" , onefile=FALSE)
		
		res = as.data.frame( rbind( response[ geneID[i] , ] ,  pfs[ geneID[i] , ] ,  os[ geneID[i] , ] ) )
		res$coef = round( as.numeric( as.character( res$coef ) ) , 2 )
		res$CI95_low = round( as.numeric( as.character( res$CI95_low ) ) , 2 )
		res$CI95_high = round( as.numeric( as.character( res$CI95_high ) ) , 2 )
		res$pval = round( as.numeric( as.character( res$pval ) ) , 3 )
		res$I2 = round( as.numeric( as.character( res$I2 ) ) , 2 )
		res$I2_pval = round( as.numeric( as.character( res$I2_pval ) ) , 3 )

		result <- cbind(
				c( "Response" , "PFS" , "OS" ) , res$coef ,
				paste( "[" , res$CI95_low , "; " , res$CI95_high , "]" , sep="" ) ,
				res[ , 4:6] )
			
		result <- as.matrix(rbind(
					c("Variable" , "coef" , "95% CI" , "P-value" , "Het." , "Het. P-value"),
					result ))

		m <- NA
		l <- NA
		u <- NA


		for(i in 1:nrow(res)){
			m <- c( m , res$coef[i])
			l <- c( l , res$CI95_low[i])
			u <- c( u , res$CI95_high[i])
		}

		xlim = round( c( min( res$coef , na.rm=TRUE)  - .2 , 0 , ( max( abs(res$coef) , na.rm=TRUE) ) + .05 ) , 1 )

		forestplot(
			result,
			new_page = TRUE,
			colgap = unit( 2 , "mm") ,
			mean= m ,
			lower= l ,
			upper= u ,
			is.summary= c( TRUE , rep( FALSE , nrow(res) ) ),
			hrzl_lines = list("2" = gpar(col = "#273746") ),
	    	zero=0,
	    	xlog=FALSE,
	    	xlab="Combined Effect-size (95% CI)",
	    	boxsize=0.21,
	    	xticks= xlim ,
			col=fpColors(
	    		box= "black" ,
	    		lines= "black" , 
	    		summary= "black" 
	    	)
		)	
	dev.off()

}
