##########################################################################################################################################################
##########################################################################################################################################################

library(survcomp)
library(genefu)
library(GSVA)
library(beeswarm)
library(plotROC)
library(data.table)
library(pROC)
library( reshape2)

#################################################
#################################################

get_Scale = function( x ){
	rid = rownames(x)
	cid = colnames(x)
	out = t( apply( x , 1 , scale ) )
	rownames(out) = rid
	colnames(out) = cid
	out
}

#Remove rows if count is < zero in 50% of sample
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0) ))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
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
	## Compte IO Meta-Score
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

#######################################################################################################################
#######################################################################################################################

source( "meta/Get_HR.R")

snv = as.data.frame( fread( "../data/validation_cohort/snv/INSPIRE.csv.gz" , sep=";" , stringsAsFactors=FALSE ))
snv = snv[ snv$Effect %in% c( "Frame_Shift_Del" , "Missense_Mutation" , "Nonsense_Mutation" , "Splice_Site" ) ,]
tmb = table( snv$Sample ) / 50 

load("../data/validation_cohort/INSPIRE.RData")

file = "../results/nsTMB_PredictIO_Validation/"
dir.create( file )
dir.create( paste( file , "/INSPIRE" , sep="" ) )

response = clin$response
names(response) = rownames(clin)

IOscore = NULL
if( ncol(expr)){ IOscore <- getIOscore( data = expr , meta_res ) }
names(IOscore) = colnames( expr )

missing_patient = names(IOscore)[ !names(IOscore) %in% names(tmb) ]
id = c( names(tmb) , missing_patient )
tmb = c( tmb , rep( 0 , length( missing_patient ) ) )
names(tmb) = id

## Association with Response (Continous)

patient = intersect( intersect( names( response ) , names( IOscore ) ) , names( tmb ) )
response = response[ patient ]
IOscore = IOscore[ patient ]
tmb = tmb[ patient ]

x = ifelse( response %in% "R" , 1 , ifelse( response %in% "NR" , 0 , NA ) )

fit = glm( x ~ IOscore , family = binomial( link="logit" ) )

d = as.data.frame( cbind( x , IOscore ) )
	colnames(d) = c( "response" , "sig" )
	d = d[ !is.na( d$response ) , ]
	d$response = as.numeric( as.character( d$response ) )
	d$sig = as.numeric( as.character( d$sig ) )

pdf( "../results/nsTMB_PredictIO_Validation/INSPIRE/IO_MetaScore_Response.pdf" , height=3.5,width=3,bg="transparent")
		
		xLabels <- paste( c( "R" , "NR" ) , "\n(N=" , table(d$response) , ")" , sep = "" )
    # yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=1 ) 
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

ggsave(filename= "../results/nsTMB_PredictIO_Validation/INSPIRE/IO_MetaScore_ROC.pdf" , plot=ROCPlot, device="pdf", height=4 , width=4 )

n_patient = length( IOscore )
save( auc , n_patient , file = "../results/nsTMB_PredictIO_Validation/INSPIRE/IO_MetaScore_AUC.RData" )

########################################################
########################################################
## Association with Response (Continous)
x = ifelse( response %in% "R" , 1 , ifelse( response %in% "NR" , 0 , NA ) )

fit = glm( x ~ tmb , family = binomial( link="logit" ) )

	d = as.data.frame( cbind( x , tmb ) )
	colnames(d) = c( "response" , "sig" )
	d = d[ !is.na( d$response ) , ]
	d$response = as.numeric( as.character( d$response ) )
	d$sig = as.numeric( as.character( d$sig ) )

pdf( "../results/nsTMB_PredictIO_Validation/INSPIRE/tmb_Response.pdf" , height=3.5,width=3,bg="transparent")
		
		xLabels <- paste( c( "R" , "NR" ) , "\n(N=" , table(d$response) , ")" , sep = "" )
    # yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=1 ) 
    yLabels <- seq( round( min( d$sig , na.rm=TRUE ) , 1 ) ,  round( max( d$sig , na.rm=TRUE ) , 1 ) , by=.25 ) 
    boxplot( sig ~ response , data=d , ylab= "tmb" , xlab="ICB Response" , main="" , 
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
  ggtitle( "tmb" ) + 
  annotate("text", x = .75, y = .25, 
           label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2))) +
  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))

auc = calc_auc(basicplot)$AUC

ggsave(filename= "../results/nsTMB_PredictIO_Validation/INSPIRE/tmb_ROC.pdf" , plot=ROCPlot, device="pdf", height=4 , width=4 )

n_patient = length( tmb )
save( auc , n_patient , file = "../results/nsTMB_PredictIO_Validation/INSPIRE/tmb_AUC.RData" )


################################################################
################################################################

dat = as.data.frame( cbind( response , IOscore , tmb ) )
dat$response = as.numeric( as.character( ifelse( dat$response %in% "R" , 1 , ifelse( dat$response %in% "NR" , 0 , NA ) ) ) )
dat$IOscore = as.numeric( as.character( dat$IOscore ) )
dat$tmb = as.numeric( as.character( dat$tmb ) )

d <- melt(dat, id.vars = "response", measure.vars = colnames(dat)[-1] ) 
d$response = as.numeric( as.character( d$response ) )
d$value = as.numeric( as.character( d$value ) )
d$variable = as.character( d$variable )

logit = NULL
sig = sort( unique( d$variable ) )
for(j in 1:length(sig) ){
	fit = glm( d$response[ d$variable %in% sig[j] ] ~ d$value[ d$variable %in% sig[j] ] , family=binomial( link="logit" ) )	
	logit = rbind( logit , cbind( sig[j] , fit$y, fit$fitted.values ) )
}
logit = as.data.frame(logit)
colnames(logit) = c( "sig" , "y" , "fitted.values" )
logit$sig = as.character( logit$sig )
logit$y = as.numeric( as.character( logit$y ) )
logit$fitted.values = as.numeric( as.character( logit$fitted.values ) )


pdf( "../results/nsTMB_PredictIO_Validation/INSPIRE/ROCcurve_PredictIO_TMB.pdf" , height=4,width=4,bg="transparent")

	roc( logit$y[ logit$sig %in% "IOscore" ] , logit$fitted.values[ logit$sig %in% "IOscore" ] , plot=TRUE, legacy.axes=FALSE, percent=TRUE, col="#E53935", lwd=2, print.auc=TRUE)
	plot.roc( logit$y[ logit$sig %in% "tmb" ] , logit$fitted.values[ logit$sig %in% "tmb" ] , percent=TRUE, col="black", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y=40)

dev.off()

################################################################
################################################################
table(ifelse(tmb < median( tmb , na.rm=TRUE ) , 0 , 1 )  , response )
table(ifelse( IOscore < median( IOscore , na.rm=TRUE ) , 0 , 1 ) , response )

table(ifelse( tmb < 10 , 0 , 1 ) , response )
table(ifelse( IOscore < quantile( IOscore , na.rm=TRUE , probs=.66 ) , 0 , 1 ) , response )

