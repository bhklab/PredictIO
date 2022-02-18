library(RColorBrewer)
library(beeswarm)
library(coin)

signature = read.table( file= "../data/signatures/signature_INFO.txt" , sep="\t" , header=TRUE , stringsAsFactor=FALSE)
signature$Signature <- as.character( signature$Signature )
signature$method <- as.character( signature$method )
signature$association <- as.character( signature$association )

association = signature$association
names( association ) = signature$Signature

studyName = as.data.frame( cbind( c( "Gide" , "INSPIRE" , "Kim" , "Puch" , "Shiuan" , "VanDenEnde" ) , 
									c( "Gide" , "INSPIRE" , "Kim" , "Puch" , "Shiuan" , "VanDenEnde" ) , 
									c( 34 , 60 , 31 , 49 , 13 , 32 ) ) )
colnames( studyName ) = c( "study" , "ID" , "N" )

studyName$study = as.character( studyName$study )
studyName$ID = as.character( studyName$ID )
studyName$N = as.numeric( as.character( studyName$N ) )
rownames( studyName ) = studyName$study

studyName = studyName[ order( studyName$study ) , ]


study_color = brewer.pal( n = nrow( studyName ) , name ="Dark2" )
names( study_color ) = studyName$study 

study_color = rev( sort( c("#d32f2f" , "#795548" , "#FFD700" , "#9370DB" , "#3CB371" , "#4682B4") ) )
names( study_color ) = studyName$study 

#####################################################################
#####################################################################

cohort = sort( c( "Kim" , "INSPIRE" , "Gide" , "Puch" , "VanDenEnde" , "Shiuan" ) ) 
metascore = "PredictIO"

auc_data = NULL

for(z in 1:length( cohort ) ){

	sig = NULL
	for(i in 1:length( association ) ){
		load( file= paste( "../results/signature_Validation/" , cohort[z] , "/" , names(association)[i] , "/" , names(association)[i] , "_AUC.RData" , sep="" ) ) 
		
		if( !is.null( auc ) ){
			sig = rbind( sig , c( names(association)[i] , association[i] , auc ) )
		}
	}

	meta = NULL
	for(i in 1:length( metascore ) ){
		load( file= paste( "../results/PredictIO_Validation/" , cohort[z] , "/" , metascore[i] , "/" , metascore[i] , "_AUC.RData" , sep="" ) ) 
		
		if( !is.null( auc ) ){
			meta = rbind( meta , c( metascore[i] , "sensitive" , auc ) )
		}
	}
	if( !is.null(meta) ){ meta[ 1 , 1 ] = "PredictIO" }

	sig = as.data.frame( sig )
	meta = as.data.frame( meta )
	colnames(sig) = colnames(meta) = c( "sig" , "association" , "auc" )

	sig$sig = as.character( sig$sig )
	sig$association = as.character( sig$association )
	sig$auc = as.numeric( as.character( sig$auc ) )

	meta$sig = as.character( meta$sig )
	meta$association = as.character( meta$association )
	meta$auc = as.numeric( as.character( meta$auc ) )

	output = rbind( 
			meta , 
			sig[ sig$association %in% "sensitive" , ][ rev( order( sig[ sig$association %in% "sensitive" , ]$auc ) ) , ] ,
			sig[ sig$association %in% "resistance" , ][ rev( order( sig[ sig$association %in% "resistance" , ]$auc ) ) , ] 
			)

	auc_data = rbind( auc_data , cbind( cohort[z] , output ) )

	bar = output$auc
	names( bar ) = output$sig
}


#####################################################################
#####################################################################

auc_data = as.data.frame( auc_data )
colnames(auc_data) = c( "cohort" , "sig" , "association" , "auc" )

auc_data$cohort = as.character( auc_data$cohort )
auc_data$sig = as.character( auc_data$sig )
auc_data$association = as.character( auc_data$association )
auc_data$auc = abs( as.numeric( as.character( auc_data$auc ) ) )

sig = as.character( unique( auc_data$sig ) ) 
median = sd = col = col_id = NULL
for( i in 1:length(sig) ){ 
	if( sig[i] %in% "PredictIO" ){
		col = c( col , "#f1c40f" )
		col_id = c( col_id , sig[i] )
	}

	c = ifelse( unique( auc_data[ auc_data$sig %in% sig[i] , ]$association ) %in% "sensitive" , "#03a9f4" , 
		ifelse( unique( auc_data[ auc_data$sig %in% sig[i] , ]$association ) %in% "resistance" , "#e53935" , NA ) )
	col = c( col , c )
	col_id = c( col_id , sig[i] )
	
	median = c( median , median( auc_data[ auc_data$sig %in% sig[i] , ]$auc , na.rm=TRUE ) ) 
	sd = c( sd , sd( auc_data[ auc_data$sig %in% sig[i] , ]$auc , na.rm=TRUE ) ) 
}
names(median) = names(sd) = sig

names( col ) = col_id


auc_data = auc_data[ !is.na( auc_data$auc ) , ]
auc_data = auc_data[ auc_data$cohort %in% auc_data[ auc_data$sig %in% "PredictIO", ]$cohort , ]

sigInData =  table( auc_data$sig )
auc_data = auc_data[ auc_data$sig %in% names( sigInData[ sigInData %in% max(sigInData) ] ) , ]

median = median[ names(median) %in% sort( unique( auc_data$sig ) ) ]
sd = sd[ names(median) ]

auc_data$sig = factor( as.character( auc_data$sig ) , levels= names( sort( median ) ) )


##########################################
##########################################

study = sort( unique( auc_data$cohort ) )

pdf( "../results/Summary_Figure/AUC/Boxplot_AUC_ALL.pdf" , height=5,width=10,bg="transparent")
	zones=matrix(c(1,2), ncol=2, byrow=TRUE)
	layout(zones, widths=c(3.8,1))


	################################################
	################################################
	## Plot
	par(mar=c(10,3,3,0))

	xLabels <- names( sort( median ) )
    yLabels <- seq( round( min( auc_data$auc , na.rm=TRUE ) , 1 ) ,  round( max( auc_data$auc , na.rm=TRUE ) , 1 ) , by=.1 ) 
    boxplot( auc ~ sig , data= auc_data , ylab= "AUC value" , xlab="" , main="" , 
    	col= "white" , 	    	 
    	boxlty = 1 , outline= FALSE , axes= FALSE , ylim= c( min( auc_data$auc ) ,  max( auc_data$auc )))

	beeswarm( auc ~ sig , data= auc_data , 
	        pwpch = ifelse( auc_data$cohort %in% "VanDenEnde" , 17 , 19 ) , corral= "wrap" , cex= .8 , 
	        pwcol = adjustcolor( study_color[ auc_data$cohort ] , alpha.f = .9 ) , 
	        add = TRUE)

    axis(side = 2, at=yLabels, labels= as.character( yLabels ), las= 2,
             cex.axis=1,tick=1 , col="black")

    mtext( "AUC Distribution across the validation cohorts" )

	Map(axis, side = 1, at=seq( 1 , length(xLabels) , 1 ) , 
		col.axis = ifelse( xLabels %in% "PredictIO" , "#f9a825" , 
					ifelse( association[ xLabels ] %in% "sensitive" , "#1976d2" , 
						ifelse( association[ xLabels ] %in% "resistance" , "#d32f2f" , NA ) ) ) , 
		padj=.5, labels=xLabels, las= 2,
      	cex.axis=1 , tick=1 , col= "black" )

	axis( 1 , at= seq( 1 , length(xLabels) , 1 ) , labels= FALSE )

	abline( h= .5 , lty= 2 , lwd= 1 )

	################################################
	################################################
	## Legend

	par(mar=c(0,0,1,0))
	plot.new()
	legend( x=0 , y=.95 ,title="Study",
			legend = paste( studyName[ sort( unique( auc_data$cohort ) ) , ]$ID , " (N=" , studyName[ sort( unique( auc_data$cohort ) ) , ]$N , ")" , sep = "" ) ,
			pt.bg= adjustcolor( study_color[ sort( unique( auc_data$cohort ) ) ] , alpha.f = .9 ) ,
			col= study_color[ sort( unique( auc_data$cohort ) ) ]  ,
			pt.cex= 1 ,
			y.intersp=1.2, cex=.7, pch= c( 21 , 21 , 21 , 21 , 21 , 24 ) , lwd=1, lty=0, bty='n', ncol=1)

	legend( x=0 , y= .57 ,title="Treatment",
			legend = c( "Adjuvant ICB" , "Neo-Adjuvant ICB" ) ,
			pt.bg= "white" ,
			col= "black"  ,
			pt.cex= 1 ,
			y.intersp=1.2, cex=.7, pch= c( 21 , 24 ) , lwd=1, lty=0, bty='n', ncol=1)

	legend( x=0 , y= .4 , title="Signature Type",
			legend = c( "Sensitive" , "Resistant" ,  "MetaScore") ,
			text.col = c( "#1976d2" , "#d32f2f" , "#f9a825" ) ,
			pt.bg= "white" ,
			col= "white" ,
			pt.cex= 1 ,
			y.intersp=1.2, cex=.7, pch= 21 , lwd=1, lty=0, bty='n', ncol=1)

dev.off()


################################################
################################################
set.seed( 1234567890 )

sig = sort( unique( as.character( auc_data[ !auc_data$sig %in% "PredictIO" , ]$sig ) ) )
ppw = pw = NULL
for(i in 1:length(sig)){
	d = auc_data[ auc_data$sig %in% c( "PredictIO" , sig[i] ) , ]

	y = as.numeric( as.character( auc_data[ auc_data$sig %in% "PredictIO" , ]$auc ) )
	x = as.numeric( as.character( auc_data[ auc_data$sig %in% sig[i] , ]$auc ) )

	## One-sided exact Wilcoxon signed-rank test
	wt <- wilcoxsign_test( y ~ x , distribution = approximate(nresample = 10000), 
									alternative = "greater" , zero.method = "Wilcoxon" )

	pw = c( pw , wilcox.test( y , x , paired = TRUE,
                     alternative = "greater", exact = FALSE)$p.value )

	ppw = c( ppw , midpvalue(wt) )

}
names( ppw ) = names( pw ) = sig

pvalue = as.data.frame( cbind( sig , round( pw , 4 ) ,  round( ppw , 4 ) , round( p.adjust( ppw , method='BH' ) , 4 ) ) )
colnames( pvalue ) = c( "signature" , "exact.P" , "permut.P" , "permut.FDR" )



