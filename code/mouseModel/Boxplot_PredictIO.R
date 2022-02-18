

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
library(GSVA)	

get_Scale = function( x ){
	rid = rownames(x)
	cid = colnames(x)
	out = t( apply( x , 1 , scale ) )
	rownames(out) = rid
	colnames(out) = cid
	out
}

rownames( tpm ) = toupper( rownames( tpm ) )

sensitive = meta_res[ meta_res$coef < 0 & meta_res$include %in% 1 , ]$gene
resistance = meta_res[ meta_res$coef > 0 & meta_res$include %in% 1 , ]$gene

IO_resistance = NULL
if( ifelse( is.null( nrow( tpm[ rownames(tpm) %in% resistance , ]) ) , 1 , nrow( tpm[ rownames(tpm) %in% resistance , ] ) ) / length( resistance ) > 0.7 ){
	IO_resistance = as.numeric( gsva( get_Scale( x= tpm ) , list(resistance) , verbose=FALSE ) )
}

IO_sensitive = NULL
if( ifelse( is.null( nrow( tpm[ rownames(tpm) %in% sensitive , ]) ) , 1 , nrow( tpm[ rownames(tpm) %in% sensitive , ] ) ) / length( sensitive ) > 0.7 ){
	IO_sensitive = as.numeric( gsva( get_Scale( x= tpm ) , list(sensitive) , verbose=FALSE ) )
}

#######################
#######################
## Compte IO Meta-Score

PredictIO = NULL
if( !is.null( IO_resistance ) & !is.null( IO_sensitive ) ){
	PredictIO =  IO_sensitive - IO_resistance
	names(PredictIO) = colnames(tpm)
}

######################################################################################
######################################################################################

d = as.data.frame( cbind( annot[ colnames( tpm ) , ] ,  PredictIO ) )
colnames(d) = c( colnames( annot ) , "gene" )
d$gene = as.numeric( as.character( d$gene ) )
d$patient = as.character( d$patient )
d$tumorCells = as.character( d$tumorCells )
d$tumorType = as.character( d$tumorType )
d$response = as.character( d$response )

pdf(  "../results/MouseModel/PredictIO/Boxplot_PredictIO_ALL.pdf", height=3.5,width=3,bg="transparent")
	
	xLabels <- paste( names(table(d$response)) , "\n(N=" , table(d$response) , ")" , sep = "" )
	yLabels <- seq( round( min( d$gene , na.rm=TRUE ) , 1 ) ,  round( max( d$gene , na.rm=TRUE ) , 1 ) , by=.25 ) 
	boxplot( gene ~ response , data=d , ylab= "PredictIO" , xlab="ICB Response" , main="" , 
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
d$response = as.numeric( ifelse( as.character( d$response ) %in% "R" , 1 , 0 ) )

basicplot <- ggplot( d, aes( d = response , m = gene )) + geom_roc( n.cuts = 100, labels = FALSE )
ROCPlot = basicplot + 
style_roc() +
theme( axis.text = element_text( colour = "black" ) ) +
ggtitle( "PredictIO" ) + 
annotate( "text", x = .75, y = .25 , 
       label = paste( "AUC =" , round( calc_auc( basicplot )$AUC , 2 ) ) ) +
scale_x_continuous( "1 - Specificity" , breaks = seq( 0 , 1 , by = .1 ) )

auc = calc_auc( basicplot )$AUC

ggsave( filename = "../results/MouseModel/PredictIO/AUC_PredictIO_ALL.pdf" , plot = ROCPlot , device = "pdf" , height = 4 , width = 4 )


######################################################################################
######################################################################################
d = as.data.frame( cbind( annot[ colnames( tpm ) , ] , PredictIO ) )
colnames(d) = c( colnames( annot ) , "gene" )
d$gene = as.numeric( as.character( d$gene ) )
d$patient = as.character( d$patient )
d$tumorCells = as.character( d$tumorCells )
d$tumorType = as.character( d$tumorType )
d$response = as.character( d$response )

d = d[ d$tumorCells %in% "AB1" , ]

pdf( "../results/MouseModel/PredictIO/Boxplot_PredictIO_AB1.pdf" , height=3.5,width=3,bg="transparent")
	
	xLabels <- paste( names(table(d$response)) , "\n(N=" , table(d$response) , ")" , sep = "" )
	# yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=1 ) 
	yLabels <- seq( round( min( d$gene , na.rm=TRUE ) , 1 ) ,  round( max( d$gene , na.rm=TRUE ) , 1 ) , by=.25 ) 
	boxplot( gene ~ response , data=d , ylab= "PredictIO" , xlab="ICB Response" , main="" , 
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

d$response = as.numeric( ifelse( as.character( d$response ) %in% "R" , 1 , 0 ) )
basicplot <- ggplot( d, aes( d = response , m = gene )) + geom_roc( n.cuts = 100, labels = FALSE )
ROCPlot = basicplot + 
style_roc() +
theme( axis.text = element_text( colour = "black" ) ) +
ggtitle( "PredictIO" ) + 
annotate( "text", x = .75, y = .25 , 
       label = paste( "AUC =" , round( calc_auc( basicplot )$AUC , 2 ) ) ) +
scale_x_continuous( "1 - Specificity" , breaks = seq( 0 , 1 , by = .1 ) )

auc = calc_auc( basicplot )$AUC

ggsave( filename = "../results/MouseModel/PredictIO/AUC_PredictIO_AB1.pdf" , plot = ROCPlot , device = "pdf" , height = 4 , width = 4 )


######################################################################################
######################################################################################
d = as.data.frame( cbind( annot[ colnames( tpm ) , ] , PredictIO ) )
colnames(d) = c( colnames( annot ) , "gene" )
d$gene = as.numeric( as.character( d$gene ) )
d$patient = as.character( d$patient )
d$tumorCells = as.character( d$tumorCells )
d$tumorType = as.character( d$tumorType )
d$response = as.character( d$response )

d = d[ d$tumorCells %in% "Renca" , ]

pdf( "../results/MouseModel/PredictIO/Boxplot_PredictIO_Renca.pdf" , height=3.5,width=3,bg="transparent")
	
	xLabels <- paste( names(table(d$response)) , "\n(N=" , table(d$response) , ")" , sep = "" )
	yLabels <- seq( round( min( d$gene , na.rm=TRUE ) , 1 ) ,  round( max( d$gene , na.rm=TRUE ) , 1 ) , by=.25 ) 
	boxplot( gene ~ response , data=d , ylab= "PredictIO" , xlab="ICB Response" , main="" , 
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

d$response = as.numeric( ifelse( as.character( d$response ) %in% "R" , 1 , 0 ) )
basicplot <- ggplot( d, aes( d = response , m = gene )) + geom_roc( n.cuts = 100, labels = FALSE )
ROCPlot = basicplot + 
style_roc() +
theme( axis.text = element_text( colour = "black" ) ) +
ggtitle( "PredictIO" ) + 
annotate( "text", x = .75, y = .25 , 
       label = paste( "AUC =" , round( calc_auc( basicplot )$AUC , 2 ) ) ) +
scale_x_continuous( "1 - Specificity" , breaks = seq( 0 , 1 , by = .1 ) )

auc = calc_auc( basicplot )$AUC

ggsave( filename = "../results/MouseModel/PredictIO/AUC_PredictIO_Renca.pdf" , plot = ROCPlot , device = "pdf" , height = 4 , width = 4 )


