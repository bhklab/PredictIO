library(RColorBrewer)
library(beeswarm)
load( "../results/Summary_Figure/CorrPlot/CorrPlot.RData" )

signature = read.table( file= "../data/signatures/signature_INFO.txt" , sep="\t" , header=TRUE , stringsAsFactor=FALSE)
signature$Signature <- as.character( signature$Signature )
signature$method <- as.character( signature$method )
signature$association <- as.character( signature$association )
rownames( signature ) = signature$Signature

association = signature$association
names( association ) = signature$Signature

cor$association_1 = signature[ as.character( cor$var1 ) , ]$association
cor$association_2 = signature[ as.character( cor$var2 ) , ]$association

sensitive = cor[ cor$association_1 %in% "sensitive" & cor$association_2 %in% "sensitive" , ]
resistance = cor[ cor$association_1 %in% "resistance" & cor$association_2 %in% "resistance" , ]

sensitive$response = "sensitive"
resistance$response = "resistance"

d = rbind( sensitive , resistance )
d$value = as.numeric( as.character( d$value ) )
d$response = as.character( d$response )

id = sort( unique( d$id ))
data = NULL
for( i in 1:length( id) ){
	data = rbind( data , 
				cbind( 
					unique( d[ d$id %in% id[ i ] , ]$response ) , 
					mean( d[ d$id %in% id[ i ] , ]$value , na.rm = TRUE ) 
				)
			)
}
data = as.data.frame( data )
colnames( data ) = c( "response" , "value" )
data$value = as.numeric( as.character( data$value ) )
data$response = as.character( data$response )


pdf( "../results/Summary_Figure/CorrPlot/Boxplot_Correlation_Overall.pdf" , height=3.5,width=3.5,bg="transparent")
	
	xLabels <- paste( names(table(data$response)) , "\n(N=" , table(data$response) , ")" , sep = "" )
	yLabels <- seq( round( min( data$value , na.rm=TRUE ) , 1 ) ,  round( max( data$value , na.rm=TRUE ) , 1 ) , by=.25 ) 
	boxplot( value ~ response , data= data , ylab= "Spearman rho" , xlab="Signature Type" , main="" , 
		col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .4), 	    	 
		boxlty = 1 ,outline=FALSE, axes=FALSE, ylim = c( min( data$value , na.rm = TRUE ) ,  max( data$value , na.rm = TRUE ) ) )

	beeswarm( value ~ response , data= data , 
	        pch = 19, corral="wrap", cex=.2, 
	        col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .5), 
	        add = TRUE)

	axis(side = 2, at=yLabels, labels=yLabels, las= 2,
	         cex.axis=1,tick=1,col="black")

	axis(side = 1, at=seq(1,length(xLabels),1) , padj=.5, labels=xLabels, las= 1,
	  cex.axis=1,tick=1,col="black")

	wil = wilcox.test( value ~ response , data= data )
	mtext( paste( "Mann Whitney U P" , 
		ifelse( wil$p.value <= .001 , "≤ 0.001" , paste( "=" , format.pval( wil$p.value , 1 ) ) ) , sep="" ) , 
		col="#6D6D6D" )

dev.off()


################################################
################################################

tumor = c( "Melanoma" , "Lung" , "Kidney")

for( k in 1:length(tumor) ){

	data = NULL
	for( i in 1:length( id) ){
		data = rbind( data , 
					cbind( 
						unique( d[ d$id %in% id[ i ] & d$primary %in% tumor[k] , ]$response ) , 
						mean( d[ d$id %in% id[ i ] & d$primary %in% tumor[k] , ]$value , na.rm = TRUE ) 
					)
				)
	}
	data = as.data.frame( data )
	colnames( data ) = c( "response" , "value" )
	data$value = as.numeric( as.character( data$value ) )
	data$response = as.character( data$response )


	pdf( paste( "../results/Summary_Figure/CorrPlot/Boxplot_Correlation_", tumor[k] , ".pdf" , sep = "" ) , height=3.5 ,width=3.5 , bg="transparent")
		
		xLabels <- paste( names(table(data$response)) , "\n(N=" , table(data$response) , ")" , sep = "" )
		yLabels <- seq( round( min( data$value , na.rm=TRUE ) , 1 ) ,  round( max( data$value , na.rm=TRUE ) , 1 ) , by=.25 ) 
		boxplot( value ~ response , data= data , ylab= "Spearman rho" , xlab="Signature Type" , main="" , 
			col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .4), 	    	 
			boxlty = 1 ,outline=FALSE, axes=FALSE, ylim = c( min( data$value , na.rm = TRUE ) ,  max( data$value , na.rm = TRUE ) ) )

		beeswarm( value ~ response , data= data , 
		        pch = 19, corral="wrap", cex=.2, 
		        col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .5), 
		        add = TRUE)

		axis(side = 2, at=yLabels, labels=yLabels, las= 2,
		         cex.axis=1,tick=1,col="black")

		axis(side = 1, at=seq(1,length(xLabels),1) , padj=.5, labels=xLabels, las= 1,
		  cex.axis=1,tick=1,col="black")

		wil = wilcox.test( value ~ response , data= data )
		mtext( paste( "Mann Whitney U P" , 
			ifelse( wil$p.value <= .001 , "≤ 0.001" , paste( "=" , format.pval( wil$p.value , 1 ) ) ) , sep="" ) , 
			col="#6D6D6D" )

	dev.off()

}

################################################
################################################