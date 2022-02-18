
load("../data/mouseModel_Zemek/Zemek.RData" )

library(beeswarm)
library(plotROC)

######################################################################################
######################################################################################
load( "../results/denovo_Single_Gene/Meta-analysis_Single_Gene_Response.RData" )
m = meta_res[ meta_res$pval <= 0.05 & meta_res$I2_pval > 0.05 , ]
m = m[ rev( order( abs( m$coef) ) ) , ]

for(i in 1:nrow(meta_res)){

	if( meta_res$pval[i] <= 0.05 ){
		meta_res$include[i] = ifelse( abs( meta_res$coef[i] ) >= abs( m$coef[ 100 ] ) , 1 , 0 )
	} else{
		meta_res$include[i] = 0
	}
}
meta_res$cutoff = abs( m$coef[ 100 ] )

meta_res = meta_res[ meta_res$include %in% 1 , ]
######################################################################################
######################################################################################

gene = meta_res$gene
rownames( tpm ) = toupper( rownames( tpm ) )

gene = gene[ gene %in% rownames( tpm ) ]

ab1 = renca = NULL

for( i in 1:length( gene ) ){

	file = paste( "../results/MouseModel/Gene/" , gene[ i ] , sep = "" )
	if( file.exists( file ) ) {
		unlink( file )
	}
	dir.create( file )

	d = as.data.frame( cbind( annot[ colnames( tpm ) , ] ,  tpm[ gene[ i ] , ] ) )
	colnames(d) = c( colnames( annot ) , "gene" )
	d$gene = as.numeric( as.character( d$gene ) )
	d$patient = as.character( d$patient )
	d$tumorCells = as.character( d$tumorCells )
	d$tumorType = as.character( d$tumorType )
	d$response = as.character( d$response )

	pdf( paste( "../results/MouseModel/Gene/", gene[ i ] ,"/Boxplot_" , gene[ i ] , "_ALL.pdf" , sep = "" ) , height=3.5,width=3,bg="transparent")
		
		xLabels <- paste( names(table(d$response)) , "\n(N=" , table(d$response) , ")" , sep = "" )
		# yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=1 ) 
		yLabels <- seq( round( min( d$gene , na.rm=TRUE ) , 1 ) ,  round( max( d$gene , na.rm=TRUE ) , 1 ) , by=.25 ) 
		boxplot( gene ~ response , data=d , ylab= paste( gene[ i ] , "mRNA expression" ) , xlab="ICB Response" , main="" , 
			col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .4), 	    	 
			boxlty = 1 ,outline=FALSE, axes=FALSE, ylim=c(min( d$gene ) ,  max( d$gene )))

		beeswarm( gene ~ response , data=d , 
		        pch = 19, corral="wrap", cex=.3, 
		        col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .6), 
		        add = TRUE)

		axis(side = 2, at=yLabels, labels=yLabels, las= 2,
		         cex.axis=1,tick=1,col="black")

		axis(side = 1, at=seq(1,length(xLabels),1) , padj=.5, labels=xLabels, las= 1,
		  cex.axis=1,tick=1,col="black")

		wil = wilcox.test( gene ~ response , data=d )
		mtext( paste( "Mann Whitney U P" , 
			ifelse( wil$p.value <= .001 , "≤ 0.001" , paste( "=" , format.pval( wil$p.value , 1 ) ) ) , sep="" ) , 
			col="#6D6D6D" )
	dev.off()

	################################
	################################
	d$response = as.numeric( ifelse( as.character( d$response ) %in% "R" , 0 , 1 ) )

	basicplot <- ggplot( d, aes( d = response , m = gene )) + geom_roc( n.cuts = 100, labels = FALSE )
	ROCPlot = basicplot + 
	style_roc() +
	theme( axis.text = element_text( colour = "black" ) ) +
	ggtitle( paste( gene[ i ] , "mRNA expression" ) ) + 
	annotate( "text", x = .75, y = .25 , 
	       label = paste( "AUC =" , round( calc_auc( basicplot )$AUC , 2 ) ) ) +
	scale_x_continuous( "1 - Specificity" , breaks = seq( 0 , 1 , by = .1 ) )

	auc = calc_auc( basicplot )$AUC

	ggsave( filename = paste( "../results/MouseModel/Gene/", gene[ i ] ,"/AUC_" , gene[ i ] , "_ALL.pdf" , sep = "" ) , plot = ROCPlot , device = "pdf" , height = 4 , width = 4 )


	######################################################################################
	######################################################################################
	d = as.data.frame( cbind( annot[ colnames( tpm ) , ] ,  tpm[ gene[ i ] , ] ) )
	colnames(d) = c( colnames( annot ) , "gene" )
	d$gene = as.numeric( as.character( d$gene ) )
	d$patient = as.character( d$patient )
	d$tumorCells = as.character( d$tumorCells )
	d$tumorType = as.character( d$tumorType )
	d$response = as.character( d$response )

	d = d[ d$tumorCells %in% "AB1" , ]

	pdf( paste( "../results/MouseModel/Gene/", gene[ i ] ,"/Boxplot_" , gene[ i ] , "_AB1.pdf" , sep = "" ) , height=3.5,width=3,bg="transparent")
		
		xLabels <- paste( names(table(d$response)) , "\n(N=" , table(d$response) , ")" , sep = "" )
		# yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=1 ) 
		yLabels <- seq( round( min( d$gene , na.rm=TRUE ) , 1 ) ,  round( max( d$gene , na.rm=TRUE ) , 1 ) , by=.25 ) 
		boxplot( gene ~ response , data=d , ylab= paste( gene[ i ] , "mRNA expression" ) , xlab="ICB Response" , main="" , 
			col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .4), 	    	 
			boxlty = 1 ,outline=FALSE, axes=FALSE, ylim=c(min( d$gene ) ,  max( d$gene )))

		beeswarm( gene ~ response , data=d , 
		        pch = 19, corral="wrap", cex=.3, 
		        col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .6), 
		        add = TRUE)

		axis(side = 2, at=yLabels, labels=yLabels, las= 2,
		         cex.axis=1,tick=1,col="black")

		axis(side = 1, at=seq(1,length(xLabels),1) , padj=.5, labels=xLabels, las= 1,
		  cex.axis=1,tick=1,col="black")

		wil = wilcox.test( gene ~ response , data=d )
		mtext( paste( "Mann Whitney U P" , 
			ifelse( wil$p.value <= .001 , "≤ 0.001" , paste( "=" , format.pval( wil$p.value , 1 ) ) ) , sep="" ) , 
			col="#6D6D6D" )
	dev.off()

	################################
	################################

	d$response = as.numeric( ifelse( as.character( d$response ) %in% "R" , 0 , 1 ) )
	basicplot <- ggplot( d, aes( d = response , m = gene )) + geom_roc( n.cuts = 100, labels = FALSE )
	ROCPlot = basicplot + 
	style_roc() +
	theme( axis.text = element_text( colour = "black" ) ) +
	ggtitle( paste( gene[ i ] , "mRNA expression" ) ) + 
	annotate( "text", x = .75, y = .25 , 
	       label = paste( "AUC =" , round( calc_auc( basicplot )$AUC , 2 ) ) ) +
	scale_x_continuous( "1 - Specificity" , breaks = seq( 0 , 1 , by = .1 ) )

	auc = calc_auc( basicplot )$AUC

	ggsave( filename = paste( "../results/MouseModel/Gene/", gene[ i ] ,"/AUC_" , gene[ i ] , "_AB1.pdf" , sep = "" ) , plot = ROCPlot , device = "pdf" , height = 4 , width = 4 )

	ab1 = rbind( ab1 , c( gene[i] , meta_res$coef[ meta_res$gene %in% gene[i] ] , auc , wil$p.value ) )

	######################################################################################
	######################################################################################
	d = as.data.frame( cbind( annot[ colnames( tpm ) , ] ,  tpm[ gene[ i ] , ] ) )
	colnames(d) = c( colnames( annot ) , "gene" )
	d$gene = as.numeric( as.character( d$gene ) )
	d$patient = as.character( d$patient )
	d$tumorCells = as.character( d$tumorCells )
	d$tumorType = as.character( d$tumorType )
	d$response = as.character( d$response )

	d = d[ d$tumorCells %in% "Renca" , ]

	pdf( paste( "../results/MouseModel/Gene/", gene[ i ] ,"/Boxplot_" , gene[ i ] , "_Renca.pdf" , sep = "" ) , height=3.5,width=3,bg="transparent")
		
		xLabels <- paste( names(table(d$response)) , "\n(N=" , table(d$response) , ")" , sep = "" )
		# yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=1 ) 
		yLabels <- seq( round( min( d$gene , na.rm=TRUE ) , 1 ) ,  round( max( d$gene , na.rm=TRUE ) , 1 ) , by=.25 ) 
		boxplot( gene ~ response , data=d , ylab= paste( gene[ i ] , "mRNA expression" ) , xlab="ICB Response" , main="" , 
			col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .4), 	    	 
			boxlty = 1 ,outline=FALSE, axes=FALSE, ylim=c(min( d$gene ) ,  max( d$gene )))

		beeswarm( gene ~ response , data=d , 
		        pch = 19, corral="wrap", cex=.3, 
		        col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .6), 
		        add = TRUE)

		axis(side = 2, at=yLabels, labels=yLabels, las= 2,
		         cex.axis=1,tick=1,col="black")

		axis(side = 1, at=seq(1,length(xLabels),1) , padj=.5, labels=xLabels, las= 1,
		  cex.axis=1,tick=1,col="black")

		wil = wilcox.test( gene ~ response , data=d )
		mtext( paste( "Mann Whitney U P" , 
			ifelse( wil$p.value <= .001 , "≤ 0.001" , paste( "=" , format.pval( wil$p.value , 1 ) ) ) , sep="" ) , 
			col="#6D6D6D" )
	dev.off()

	################################
	################################

	d$response = as.numeric( ifelse( as.character( d$response ) %in% "R" , 0 , 1 ) )
	basicplot <- ggplot( d, aes( d = response , m = gene )) + geom_roc( n.cuts = 100, labels = FALSE )
	ROCPlot = basicplot + 
	style_roc() +
	theme( axis.text = element_text( colour = "black" ) ) +
	ggtitle( paste( gene[ i ] , "mRNA expression" ) ) + 
	annotate( "text", x = .75, y = .25 , 
	       label = paste( "AUC =" , round( calc_auc( basicplot )$AUC , 2 ) ) ) +
	scale_x_continuous( "1 - Specificity" , breaks = seq( 0 , 1 , by = .1 ) )

	auc = calc_auc( basicplot )$AUC

	ggsave( filename = paste( "../results/MouseModel/Gene/", gene[ i ] ,"/AUC_" , gene[ i ] , "_Renca.pdf" , sep = "" ) , plot = ROCPlot , device = "pdf" , height = 4 , width = 4 )


	renca = rbind( renca , c( gene[i] , meta_res$coef[ meta_res$gene %in% gene[i] ]  , auc , wil$p.value ) )
	
}
