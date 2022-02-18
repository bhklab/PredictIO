##########################################################################################################################################################
##########################################################################################################################################################

library(survcomp)
library(genefu)
library(GSVA)
library(beeswarm)
library(plotROC)

#################################################
#################################################

get_Directory <- function( cohort ){

	file = paste( "../results/nsTMB_PredictIO/" , cohort , sep = "" )
	
	if( file.exists( file ) ) {
		unlink( file , recursive = TRUE )
	}

	dir.create( file )

}

#Remove rows if count is < zero in 50% of sample
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0) ))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}


get_Scale = function( x ){
	rid = rownames(x)
	cid = colnames(x)
	out = t( apply( x , 1 , scale ) )
	rownames(out) = rid
	colnames(out) = cid
	out
}

getIOscore <- function( data , meta_res ){
	
	remove <- rem(data)
	if( length(remove) ){
		data <- data[-remove,]
	}		
	
	sensitive = meta_res[ meta_res$coef < 0 & meta_res$include %in% 1 , ]$gene
	resistance = meta_res[ meta_res$coef > 0 & meta_res$include %in% 1 , ]$gene

	
	IO_resistance = NULL
	if( ifelse( is.null( nrow( data[ rownames(data) %in% resistance , ]) ) , 1 , nrow( data[ rownames(data) %in% resistance , ] ) ) / length( resistance ) > 0.8 ){
		IO_resistance = as.numeric( gsva( get_Scale( x= data ) , list(resistance) , verbose=FALSE ) )
	}
	
	IO_sensitive = NULL
	if( ifelse( is.null( nrow( data[ rownames(data) %in% sensitive , ]) ) , 1 , nrow( data[ rownames(data) %in% sensitive , ] ) ) / length( sensitive ) > 0.8 ){
		IO_sensitive = as.numeric( gsva( get_Scale( x= data ) , list(sensitive) , verbose=FALSE ) )
	}

	####################################################################################
	####################################################################################
	## Compte PredictIO
	signature = NULL
	if( !is.null( IO_resistance ) & !is.null( IO_sensitive ) ){
		signature = IO_sensitive - IO_resistance
		names(signature) = colnames(data)
	}
	signature
}

#################################################
#################################################

load( "../results/denovo_Single_Gene/Meta-analysis_Single_Gene_Response.RData" )

meta_res = as.data.frame(meta_res)
colnames(meta_res) = c( "gene" , "coef" , "se" , "pval")
meta_res$gene = as.character( meta_res$gene )
meta_res$coef = as.numeric( as.character( meta_res$coef ) )
meta_res$se = as.numeric( as.character( meta_res$se ) )
meta_res$pval = as.numeric( as.character( meta_res$pval ) ) 

m = meta_res[ meta_res$pval <= 0.05 , ]
m = m[ rev( order( abs( m$coef) ) ) , ]

for(i in 1:nrow(meta_res)){

	if( meta_res$pval[i] <= 0.05 ){
		meta_res$include[i] = ifelse( abs( meta_res$coef[i] ) >= abs( m$coef[ 100 ] ) , 1 , 0 )
	} else{
		meta_res$include[i] = 0
	}
}
meta_res$cutoff = abs( m$coef[ 100 ] )

#################################################
#################################################

source( "meta/Get_HR.R")

load("../results/nsTMB/nsTMB.cutoff.RData")

file = "../results/nsTMB_PredictIO/"
if( file.exists( file ) ) {	unlink( file , recursive = TRUE ) }
dir.create( file )

load("../data/process_data/ICB_exp_filtered.RData")

expr
study = names(expr)

for( i in 1:length(study)){

	data = exprs(expr[[i]])[ , phenoData(expr[[i]])$rna %in% c( "fpkm" , "tpm" ) & !is.na( phenoData(expr[[i]])$nsTMB_perMb ) ]
	
	nsTMB = NULL
	nsTMB = log10( phenoData(expr[[i]])$nsTMB_perMb[ phenoData(expr[[i]])$rna %in% c( "fpkm" , "tpm" ) & !is.na( phenoData(expr[[i]])$nsTMB_perMb ) ] + .5 )
	response = phenoData(expr[[i]])$response[ phenoData(expr[[i]])$rna %in% c( "fpkm" , "tpm" ) & !is.na( phenoData(expr[[i]])$nsTMB_perMb ) ]

	IOscore = NULL
	if( ncol(data)){ IOscore <- getIOscore( data = data , meta_res ) }

	if( !is.null( IOscore ) & !is.null( nsTMB ) ){

		get_Directory( cohort = study[i] )


		#######################################################################################################################
		#######################################################################################################################
		## Association with Response (Continous)
		x = ifelse( response %in% "R" , 1 , ifelse( response %in% "NR" , 0 , NA ) )

		fit = glm( x ~ IOscore , family = binomial( link="logit" ) )

		d = as.data.frame( cbind( x , IOscore ) )
  		colnames(d) = c( "response" , "sig" )
  		d = d[ !is.na( d$response ) , ]
  		d$response = as.numeric( as.character( d$response ) )
  		d$sig = as.numeric( as.character( d$sig ) )

		pdf( paste( "../results/nsTMB_PredictIO/" , study[i] , "/PredictIO_Response.pdf" , sep="" ) , height=3.5,width=3,bg="transparent")
	  		
	  		xLabels <- paste( c( "R" , "NR" ) , "\n(N=" , table(d$response) , ")" , sep = "" )
		    yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=.25 ) 
		    boxplot( sig ~ response , data=d , ylab= "PredictIO" , xlab="ICB Response" , main="" , 
		    	col=adjustcolor( c( "#00bcd4" , "#f44336" ), alpha.f = .4), 	    	 
		    	boxlty = 1 ,outline=FALSE, axes=FALSE, ylim=c(min( d$sig ) ,  max( d$sig )))

			beeswarm( sig ~ response , data=d , 
			        pch = 19, corral="wrap", cex=.3, 
			        col=adjustcolor( c( "#00bcd4" , "#f44336" ), alpha.f = .6), 
			        add = TRUE)

		    axis(side = 2, at=yLabels, labels=yLabels, las= 2,
		             cex.axis=1,tick=1,col="black")

	    	axis(side = 1, at=seq(1,length(xLabels),1) , padj=.5, labels=xLabels, las= 1,
		      cex.axis=1,tick=1,col="black")

			mtext( paste( "LogOR=" , round( summary(fit)$coefficients[ 2 , 1 ] , 2 ) , 
					" (95CI= [" , round( confint(fit)[ 2 , 1 ] , 2 ) , "; " , 
					round( confint(fit)[ 2 , 2 ] , 2 ) , "])\nP" , 
					ifelse( summary(fit)$coefficients[ 2 , 4 ] <= .001 , "≤ 0.001" , paste( "=" , format.pval( summary(fit)$coefficients[ 2 , 4 ] , 1 ) ) ) , sep="" ) , 
					col="#6D6D6D" )
		dev.off()

		################################################################
		################################################################

		basicplot <- ggplot(d, aes(d = response, m = sig )) + geom_roc(n.cuts = 100, labels = FALSE)
		ROCPlot = basicplot + 
		  style_roc() +
		  theme(axis.text = element_text(colour = "black")) +
		  ggtitle( "PredictIO" ) + 
		  annotate("text", x = .75, y = .25, 
		           label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2))) +
		  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))

		auc = calc_auc(basicplot)$AUC

		ggsave(filename= paste( "../results/nsTMB_PredictIO/" , study[i] , "/PredictIO_ROC.pdf" , sep="" ) , plot=ROCPlot, device="pdf", height=4 , width=4 )
		
		n_patient = length( IOscore )
		save( auc , n_patient , file = paste( "../results/nsTMB_PredictIO/" , study[i] , "/PredictIO_AUC.RData" , sep="" ) )

		#######################################################################################################################
		#######################################################################################################################
		## Association with Response (Continous)
		x = ifelse( response %in% "R" , 1 , ifelse( response %in% "NR" , 0 , NA ) )

		fit = glm( x ~ nsTMB , family = binomial( link="logit" ) )

  		d = as.data.frame( cbind( x , nsTMB ) )
  		colnames(d) = c( "response" , "sig" )
  		d = d[ !is.na( d$response ) , ]
  		d$response = as.numeric( as.character( d$response ) )
  		d$sig = as.numeric( as.character( d$sig ) )

		pdf( paste( "../results/nsTMB_PredictIO/" , study[i] , "/nsTMB_Response.pdf" , sep="" ) , height=3.5,width=3,bg="transparent")
	  		
	  		xLabels <- paste( c( "R" , "NR" ) , "\n(N=" , table(d$response) , ")" , sep = "" )
		    # yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=1 ) 
		    yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=.25 ) 
		    boxplot( sig ~ response , data=d , ylab= "nsTMB" , xlab="ICB Response" , main="" , 
		    	col=adjustcolor( c( "#00bcd4" , "#f44336" ), alpha.f = .4), 	    	 
		    	boxlty = 1 ,outline=FALSE, axes=FALSE, ylim=c(min( d$sig ) ,  max( d$sig )))

			beeswarm( sig ~ response , data=d , 
			        pch = 19, corral="wrap", cex=.3, 
			        col=adjustcolor( c( "#00bcd4" , "#f44336" ), alpha.f = .6), 
			        add = TRUE)

		    axis(side = 2, at=yLabels, labels=yLabels, las= 2,
		             cex.axis=1,tick=1,col="black")

	    	axis(side = 1, at=seq(1,length(xLabels),1) , padj=.5, labels=xLabels, las= 1,
		      cex.axis=1,tick=1,col="black")

			mtext( paste( "LogOR=" , round( summary(fit)$coefficients[ 2 , 1 ] , 2 ) , 
					" (95CI= [" , round( confint(fit)[ 2 , 1 ] , 2 ) , "; " , 
					round( confint(fit)[ 2 , 2 ] , 2 ) , "])\nP" , 
					ifelse( summary(fit)$coefficients[ 2 , 4 ] <= .001 , "≤ 0.001" , paste( "=" , format.pval( summary(fit)$coefficients[ 2 , 4 ] , 1 ) ) ) , sep="" ) , 
					col="#6D6D6D" )
		dev.off()

		################################################################
		################################################################

		basicplot <- ggplot(d, aes(d = response, m = sig )) + geom_roc(n.cuts = 100, labels = FALSE)
		ROCPlot = basicplot + 
		  style_roc() +
		  theme(axis.text = element_text(colour = "black")) +
		  ggtitle( "nsTMB" ) + 
		  annotate("text", x = .75, y = .25, 
		           label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2))) +
		  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))

		auc = calc_auc(basicplot)$AUC

		ggsave(filename= paste( "../results/nsTMB_PredictIO/" , study[i] , "/nsTMB_ROC.pdf" , sep="" ) , plot=ROCPlot, device="pdf", height=4 , width=4 )

		n_patient = length( nsTMB )
		save( auc , n_patient , file = paste( "../results/nsTMB_PredictIO/" , study[i] , "/nsTMB_AUC.RData" , sep="" ) )

	}
}

#######################################################################################################################
#######################################################################################################################

auc_data = NULL
study = list.files( "../results/nsTMB_PredictIO/" )
for(i in 1:length( study ) ){
	load( paste( "../results/nsTMB_PredictIO/" , study[i] , "/nsTMB_AUC.RData" , sep="" ) )
	auc_data = rbind( auc_data , c( study[i] , "nsTMB" , n_patient , auc ) )
	
	load( paste( "../results/nsTMB_PredictIO/" , study[i] , "/PredictIO_AUC.RData" , sep="" ) )
	auc_data = rbind( auc_data , c( study[i] , "PredictIO" , n_patient , auc ) )

}
auc_data = as.data.frame( auc_data )
colnames( auc_data ) = c( "cohort" , "sig" , "N" , "auc" )

auc_data$cohort = as.character( auc_data$cohort )
auc_data$sig = as.character( auc_data$sig )
auc_data$N = as.numeric( as.character( auc_data$N ) )
auc_data$auc = as.numeric( as.character( auc_data$auc ) )

#######################################################################################################################
#######################################################################################################################
library(RColorBrewer)

cohort = as.data.frame( unique( cbind( auc_data$cohort , auc_data$N ) ) )
colnames( cohort ) = c( "cohort" , "N" )

cohort$cohort = as.character( cohort$cohort )
cohort$N = as.numeric( as.character( cohort$N ) )

study_color = brewer.pal( n = nrow( cohort ) , name ="Paired" )
names( study_color ) = cohort$cohort 

wilcox.test( auc ~ sig , data= auc_data , paired = TRUE )

###################################################################
###################################################################

library(coin)
set.seed( 1234567890 )

y = as.numeric( as.character( auc_data[ auc_data$sig %in% "PredictIO" , ]$auc ) )
x = as.numeric( as.character( auc_data[ auc_data$sig %in% "nsTMB" , ]$auc ) )

## One-sided exact Wilcoxon signed-rank test
wt <- wilcoxsign_test( y ~ x , distribution = approximate(nresample = 10000), 
								alternative = "greater" , zero.method = "Wilcoxon" )



###################################################################
###################################################################



pdf( "../results/nsTMB_PredictIO/Boxplot_AUC.pdf" , height=5,width=4,bg="transparent")
	zones=matrix(c(1,2), ncol=2, byrow=TRUE)
	layout(zones, widths=c(2,2))


	################################################
	################################################
	## Plot
	par(mar=c(10,3,3,0))

	xLabels <- c( "nsTMB" , "PredictIO" )
    yLabels <- seq( round( min( c( auc_data$auc , .4 ) , na.rm=TRUE ) , 1 ) ,  round( max( auc_data$auc , na.rm=TRUE ) , 1 ) , by=.1 ) 
    boxplot( auc ~ sig , data= auc_data , ylab= "AUC value" , xlab="" , main="" , 
    	col= "white" , 	    	 
    	boxlty = 1 , outline= FALSE , axes= FALSE , ylim= c( min( c( auc_data$auc , .4 ) ) ,  max( auc_data$auc )))

	# beeswarm( auc ~ sig , data= auc_data , 
	#         pch = 19, corral= "wrap" , cex= 1 , 
	#         col=adjustcolor( col[ names( sort( median ) ) ] , alpha.f = .9 ) , 
	#         add = TRUE)

	beeswarm( auc ~ sig , data= auc_data , 
	        pch = 19 , corral= "wrap" , cex= 1 , 
	        pwcol = adjustcolor( study_color[ auc_data$cohort ] , alpha.f = .9 ) , 
	        add = TRUE)

    axis(side = 2, at=yLabels, labels= as.character( yLabels ), las= 2,
             cex.axis=1,tick=1 , col="black")
    axis(side = 1, at=seq( 1 , length(xLabels) , 1 ) , labels= xLabels, las= 2,
             cex.axis=1,tick=1 , col="black")

    mtext( "AUC Distribution" )

	abline( h= .5 , lty= 2 , lwd= 1 )

	################################################
	################################################
	## Legend

	par(mar=c(0,0,1,0))
	plot.new()
	legend( x=0 , y=.95 ,title="Study",
			legend = paste( cohort$cohort , " (N=" , cohort$N , ")" , sep = "" ) ,
			pt.bg= adjustcolor(  study_color[ cohort$cohort ] , alpha.f = 0.9 ) ,
			col= adjustcolor(  study_color[ cohort$cohort ] , alpha.f = 0.9 )  ,
			pt.cex= 1 ,
			y.intersp=1.2, cex=.7, pch= 21 , lwd=1, lty=0, bty='n', ncol=1)
dev.off()






