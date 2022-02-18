
##################################################################################################################################################################
##################################################################################################################################################################

library(meta)
library(metafor)
library(genefu)
library(forestplot)
library(dmetar)

Get_Cox_Forestplot = function( data , cancer , seq , label , prefix , dir , height_1 , height_2 , height_3 , height_4 , width_1 , width_2 , width_3 , width_4 ){
	data$study <- as.character( data$study )
	data$Sequencing <- as.character( data$Sequencing )
	data$Primary <- as.character( data$Primary )
	data$HR <- as.numeric(as.character( data$HR ))
	data$SE <- as.numeric(as.character( data$SE ))
	data$Pval <- as.numeric(as.character( data$Pval )) 

	data = data[ order( data$HR ) , ]
	cancer = cancer[ order( cancer$HR ) , ]
	seq = seq[ order( seq$HR ) , ]

	# Get xlim
	m <- c( min( c( 0 , data$HR ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(data$HR) ) , na.rm=TRUE) ) + .5 )

	data$study <- paste( data$study , ", " , data$Primary , ", n = " , data$N , sep= "" ) 
	
	meta <- metagen( TE = HR,
	                  seTE = SE ,
	                  data = data ,
	                  studlab = study ,
	                  fixed = FALSE ,
	                  random = TRUE ,
	                  control = list( maxiter = 10000 , stepadj=0.5 ) )

	######################################################################
	######################################################################
	## Save the merged HR and Pvalue
	meta_res <- c( unlist( strsplit( prefix , "/" , fixed=TRUE ) )[ length(unlist( strsplit( prefix , "/" , fixed=TRUE ) )) ] , 
					meta$TE.random ,  
					meta$seTE.random ,   
					meta$lower.random ,   
					meta$upper.random ,
					meta$pval.random , 
					meta$I2 ,
					meta$pval.Q )
	names(meta_res) <- c( "study" , "logHR" , "se_logHR" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" )
	save( meta_res , file=paste( prefix , "/Overall/" , dir , "/" , label , ".RData", sep="" ) )
	######################################################################
	######################################################################

	pdf( paste( prefix , "/Overall/" , dir , "/" , label , ".pdf", sep="" ), height= height_1, width= width_1 , bg="transparent" , onefile=FALSE )
		forest( meta , 
            leftcols = c("studlab", "effect.ci" , "Pval" ),
			leftlabs= c( "Study" , "logHR [95%CI]" , "P-value" ) , 
   			xlab = "logHR estimate",
			digits.se = 2 ,
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
	       	col.by = "#1565c0",
		    addrow.subgroups=TRUE 
	    )
	dev.off()


	######################################################################
	######################################################################
	## Small Sample Biais methods : Funnel Plots

	if( nrow(meta$data) >= 10 ){
		pdf( paste( prefix , "/Funnel/" , dir , "/" , label , ".pdf", sep="" ), height= height_4, width= width_4 , bg="transparent" , onefile=FALSE )
			egger = eggers.test(x = meta )
			funnel( meta , ylab="SE(logHR)" , xlab="Hedges' g", studlab = TRUE , comb.random= TRUE , comb.fixed= FALSE )
			mtext( paste( "Egger's test P" , ifelse( egger$p <0.01 , "<0.01" , paste( "=" , round( egger$p , 2 ) , sep="" ) ) , " (N=" , nrow(meta$data) , ")" ,sep="" ) , 
				font=2 , cex=1.2 , padj=-.5 )  
		dev.off()
	}

	#####################################
	#####################################
	# Per Cancer

	remove <- names( table( cancer$Primary)[ table(cancer$Primary) %in% c(1,2) ] )

	if( length( unique( cancer$Primary[ !cancer$Primary %in% remove ] ) ) > 1 ){
		
		# Get xlim
		m <- c( min( c( 0 , cancer$HR ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(cancer$HR) ) , na.rm=TRUE) ) + .5 )

		meta <- metagen( 	TE = HR ,
							seTE = SE ,
							data = cancer ,
							studlab = study ,
							fixed = FALSE ,
							random = TRUE ,
							control = list( maxiter = 10000 , stepadj=0.5 ) )

		if( length(remove) > 0 ){

			meta.subgroup <- update.meta(meta , 
									byvar=Primary , 
									exclude = cancer$Primary %in% remove ,
									fixed = FALSE ,
									random = TRUE ,
									control = list( maxiter = 10000 , stepadj=0.5 ) )
		} else{
			meta.subgroup <- update.meta(meta , 
									byvar=Primary , 
									comb.random = TRUE ,
									fixed = FALSE ,
									random = TRUE ,
									control = list( maxiter = 10000 , stepadj=0.5 ) )
		}

		pdf( paste( prefix , "/PerCancer/" , dir , "/" , label , "_byCancerType.pdf", sep="" ), height= height_2, width= width_2 , bg="transparent" , onefile=FALSE )


			forest( meta.subgroup , 
	   #          leftcols = c("studlab", "effect.ci", "Pval"),
				# leftlabs= c( "Study" , "logHR [95%CI]" , "P-value" ) , 
				digits.se = 2,
	      		colgap.forest=unit(10, "mm") ,
	      		plotwidth = unit( 30 , "mm") , 
				xlab = "logHR estimate",
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
		       	col.by = "#1565c0",
		       	addrow.subgroups=TRUE 
		    )


		dev.off()
	}

	######################################
	######################################
	## Per Sequencing

	remove <- names( table(seq$Sequencing)[ table(seq$Sequencing) %in% c(1,2) ] )
	
	if( length( unique( seq$Sequencing[ !seq$Sequencing %in% remove ] ) ) > 1 ){
	
		# Get xlim
		m <- c( min( c( 0 , seq$HR ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(seq$HR) ) , na.rm=TRUE) ) + .5 )

		meta <- metagen( 	TE = HR,
							seTE = SE,
							data = seq ,
							studlab = study ,
							fixed = FALSE ,
							random = TRUE ,
							control = list( maxiter = 10000 , stepadj=0.5 ) )

		if( length(remove) > 0 ){

			meta.subgroup <- update.meta(meta, 
										byvar=Sequencing, 
										exclude = seq$Sequencing %in% remove ,
										fixed = FALSE ,
										random = TRUE ,
										control = list( maxiter = 10000 , stepadj=0.5 ) )
		} else{
			meta.subgroup <- update.meta(meta, 
		                             	byvar=Sequencing ,
										fixed = FALSE ,
										random = TRUE ,
										control = list( maxiter = 10000 , stepadj=0.5 ) )
		}

		pdf( paste( prefix , "/PerSequencing/" , dir , "/" , label , "_bySequencing.pdf", sep="" ), height= height_3, width= width_3 , bg="transparent" , onefile=FALSE )
			forest( meta.subgroup , 
	   #          leftcols = c("studlab" , "effect.ci" , "Pval"),
				# leftlabs= c( "Study" , "logHR [95%CI]" , "P-value" ) , 
				digits.se = 2,
	      		colgap.forest=unit(10, "mm") , 
	      		plotwidth = unit( 30 , "mm") ,
				xlab = "logHR estimate",
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
		       	col.by = "#1565c0",
		       	addrow.subgroups=TRUE 
		    )

		dev.off()
	}

}


##################################################################################################################################################################
##################################################################################################################################################################

Get_DI_Forestplot = function( data , cancer , seq , label , prefix , dir , height_1 , height_2 , height_3 , height_4 , width_1 , width_2 , width_3 , width_4 ){
	
	data$study = as.character( data$study )
	data$Sequencing = as.character( data$Sequencing )
	data$Primary = as.character( data$Primary )
	data$DI = as.numeric(as.character( data$DI ))
	data$SE = as.numeric(as.character( data$SE ))
	data$Pval = as.numeric(as.character( data$Pval )) 


	data = data[ order( data$DI ) , ]
	cancer = cancer[ order( cancer$DI ) , ]
	seq = seq[ order( seq$DI ) , ]

	# Get xlim
	m <- c( min( c( data$DI , 1 ) , na.rm=TRUE) - .5 , ( max( c( 1 , abs(data$DI) ) , na.rm=TRUE) ) + .5 )

	data$study = paste( data$study , ", " , data$Primary , ", n = " , data$N , sep= "" ) 

	meta <- metagen( TE = DI ,
	                  seTE = SE ,
	                  data = data ,
	                  studlab = study ,
	                  fixed = FALSE ,
	                  random = TRUE ,
	                  control = list( maxiter = 10000 , stepadj=0.5 ) )

	######################################################################
	######################################################################
	## Save the merged DI and Pvalue
	meta_res = c( unlist( strsplit( prefix , "/" , fixed=TRUE ) )[ length(unlist( strsplit( prefix , "/" , fixed=TRUE ) )) ] , 
					meta$TE.random ,  
					meta$seTE.random ,   
					meta$lower.random ,   
					meta$upper.random ,
					meta$pval.random , 
					meta$I2 ,
					meta$pval.Q )
	names(meta_res) = c( "study" , "DI" , "se_DI" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
	save( meta_res , file=paste( prefix , "/Overall/" , dir , "/" , label , ".RData", sep="" ) )
	######################################################################
	######################################################################


	pdf( paste( prefix , "/Overall/" , dir , "/" , label , ".pdf", sep="" ), height= height_1, width= width_1 , bg="transparent" , onefile=FALSE )
		forest( meta , 
            leftcols = c("studlab", "effect.ci" , "Pval"),
			leftlabs= c( "Study" , "DI [95%CI]" , "P-value" ) , 
   			xlab = "D.Index estimate",
			xlim=m , ref=0 , 
			digits.se = 2,
			colgap.forest=unit(20, "mm") ,
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
	       	col.square= "black" ,  
	       	col.study= "black" ,  
	       	col.square.lines = "black" ,
	       	col.diamond.random  = "#1565c0"  ,
	       	col.diamond.lines.random  ="#1565c0" ,
	       	col.by = "#1565c0",
		    addrow.subgroups=TRUE
	    )

	dev.off()

	######################################################################
	######################################################################
	## Small Sample Bias methods : Funnel Plots

	if( nrow(meta$data) >= 10 ){
		pdf( paste( prefix , "/Funnel/" , dir , "/" , label , ".pdf", sep="" ), height= height_4, width= width_4 , bg="transparent" , onefile=FALSE )
			egger = eggers.test(x = meta )
			funnel( meta , ylab="SE(DI)" , xlab="Hedges' g", studlab = TRUE , comb.random= TRUE , comb.fixed= FALSE )
			mtext( paste( "Egger's test P" , ifelse( egger$p <0.01 , "<0.01" , paste( "=" , round( egger$p , 2 ) , sep="" ) ) , " (N=" , nrow(meta$data) , ")" ,sep="" ) , 
				font=2 , cex=1.2 , padj=-.5 )  
		dev.off()
	}

	#####################################
	#####################################
	# Per Cancer
	
	remove <- names( table( cancer$Primary)[ table(cancer$Primary) %in% c(1,2) ] )
	
	if( length( unique( cancer$Primary[ !cancer$Primary %in% remove ] ) ) > 1 ){

		meta <- metagen( 	TE = DI ,
							seTE = SE ,
							data = cancer ,
							studlab = study ,
							fixed = FALSE ,
							random = TRUE ,
							control = list( maxiter = 10000 , stepadj=0.5 ) )


		if( length(remove) > 0 ){

			meta.subgroup <- update.meta(meta , 
										byvar=Primary , 
										exclude = cancer$Primary %in% remove ,
										fixed = FALSE ,
										random = TRUE ,
										control = list( maxiter = 10000 , stepadj=0.5 ) )
						
		} else{
			meta.subgroup <- update.meta(meta , 
		                             	byvar=Primary ,
										fixed = FALSE ,
										random = TRUE ,
										control = list( maxiter = 10000 , stepadj=0.5 ) )
		}

		pdf( paste( prefix , "/PerCancer/" , dir , "/" , label , "_byCancerType.pdf", sep="" ), height= height_2, width= width_2 , bg="transparent" , onefile=FALSE )
			forest( meta.subgroup , 
	   #          leftcols = c("studlab", "effect.ci" , "Pval"),
				# leftlabs= c( "Study" , "DI [95%CI]" , "P-value" ) , 
				digits.se = 2,
	      		colgap.forest=unit(20, "mm") , 
	      		plotwidth = unit( 30 , "mm") ,
				xlab = "D.Index estimate",
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
				ref=0,
		       	addrow.subgroups=TRUE
			 )
		dev.off()
	}

	######################################
	######################################
	## Per Sequencing

	remove <- names( table(seq$Sequencing)[ table(seq$Sequencing) %in% c(1,2) ] )
	
	if( length(unique(seq$Sequencing[!seq$Sequencing%in%remove])) > 1 ){
		meta <- metagen( 	TE = DI ,
							seTE = SE ,
							data = seq ,
							studlab = study ,
							fixed = FALSE ,
							random = TRUE ,
							control = list( maxiter = 10000 , stepadj=0.5 ) )

		if( length(remove) > 0 ){

			meta.subgroup <- update.meta(meta , 
										byvar=Sequencing , 
										exclude = seq$Sequencing %in% remove ,
										fixed = FALSE ,
										random = TRUE ,
										control = list( maxiter = 10000 , stepadj=0.5 ) )
		} else{
			meta.subgroup <- update.meta(meta, 
		                             	byvar=Sequencing ,
										fixed = FALSE ,
										random = TRUE ,
										control = list( maxiter = 10000 , stepadj=0.5 ) )
		}

		pdf( paste( prefix , "/PerSequencing/" , dir , "/" , label , "_bySequencing.pdf", sep="" ), height= height_3, width= width_3 , bg="transparent" , onefile=FALSE )
			forest( meta.subgroup , 
	   #          leftcols = c("studlab" , "effect.ci" , "Pval"),
				# leftlabs= c( "Study" , "DI [95%CI]" , "P-value" ) , 
				digits.se = 2,
	      		colgap.forest=unit(20, "mm") ,
	      		plotwidth = unit( 30 , "mm") , 
				xlab = "D.Index estimate",
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
				ref=0,
		       	addrow.subgroups=TRUE
			)
		dev.off()
	}
}

##################################################################################################################################################################
##################################################################################################################################################################

Get_LogReg_Forestplot = function( data , cancer , seq , label , prefix , dir , height_1 , height_2 , height_3 , height_4 , width_1 , width_2 , width_3 , width_4 ){
	data$study = as.character( data$study )
	data$Sequencing = as.character( data$Sequencing )
	data$Primary = as.character( data$Primary )
	data$coef = as.numeric(as.character( data$coef ))
	data$SE = as.numeric(as.character( data$SE ))
	data$Pval = as.numeric(as.character( data$Pval )) 

	data = data[ order( data$coef ) , ]
	cancer = cancer[ order( cancer$coef ) , ]
	seq = seq[ order( seq$coef ) , ]

	# Get xlim
	m <- c( min( c( data$coef , 0 ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(data$coef) ) , na.rm=TRUE) ) + .5 )

	data$study = paste( data$study , ", " , data$Primary , ", n = " , data$N , sep= "" ) 

	meta <- metagen( TE = coef ,
	                  seTE = SE ,
	                  data = data ,
	                  studlab = study ,
	                  fixed = FALSE ,
	                  random = TRUE ,
	                  control = list( maxiter = 10000 , stepadj=0.5 ) )

	######################################################################
	######################################################################
	## Save the merged coef and Pvalue
	meta_res = c( unlist( strsplit( prefix , "/" , fixed=TRUE ) )[ length(unlist( strsplit( prefix , "/" , fixed=TRUE ) )) ] , 
					meta$TE.random ,  
					meta$seTE.random ,   
					meta$lower.random ,   
					meta$upper.random ,
					meta$pval.random , 
					meta$I2 ,
					meta$pval.Q )
	names(meta_res) = c( "study" , "coef" , "se_coef" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
	save( meta_res , file=paste( prefix , "/Overall/" , dir , "/" , label , ".RData", sep="" ) )
	######################################################################
	######################################################################

	pdf( paste( prefix , "/Overall/" , dir , "/" , label , ".pdf", sep="" ), height= height_1, width= width_1 , bg="transparent" , onefile=FALSE )
		forest( meta , 
            leftcols = c("studlab", "effect.ci" , "Pval"),
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
	
	######################################################################
	######################################################################
	## Small Sample Bias methods : Funnel Plots

	if( nrow(meta$data) >= 10 ){
		pdf( paste( prefix , "/Funnel/" , dir , "/" , label , ".pdf", sep="" ), height= height_4, width= width_4 , bg="transparent" , onefile=FALSE )
			egger = eggers.test(x = meta )
			funnel( meta , ylab="SE(logOR)" , xlab="Hedges' g", studlab = TRUE , comb.random= TRUE , comb.fixed= FALSE )
			mtext( paste( "Egger's test P" , ifelse( egger$p <0.01 , "<0.01" , paste( "=" , round( egger$p , 2 ) , sep="" ) ) , " (N=" , nrow(meta$data) , ")" ,sep="" ) , 
				font=2 , cex=1.2 , padj=-.5 )  
		dev.off()
	}

	#####################################
	#####################################
	# Per Cancer
	
	remove <- names( table( cancer$Primary)[ table(cancer$Primary) %in% c(1,2) ] )
	
	if( length( unique( cancer$Primary[ !cancer$Primary %in% remove ] ) ) > 1 ){

		# Get xlim
		m <- c( min( c( 0 , cancer$coef ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(cancer$coef) ) , na.rm=TRUE) ) + .5 )

		meta <- metagen( 	TE = coef ,
							seTE = SE ,
							data = cancer ,
							studlab = study ,
							fixed = FALSE ,
							random = TRUE ,
							control = list( maxiter = 10000 , stepadj=0.5 ) )


		if( length(remove) > 0 ){

			meta.subgroup <- update.meta(meta, 
										byvar=Primary, 
										exclude = cancer$Primary %in% remove ,
										fixed = FALSE ,
										random = TRUE ,
										control = list( maxiter = 10000 , stepadj=0.5 ) )
		} else{
			meta.subgroup <- update.meta(meta, 
		                             	byvar=Primary ,
										fixed = FALSE ,
										random = TRUE ,
										control = list( maxiter = 10000 , stepadj=0.5 ) )
		}

		pdf( paste( prefix , "/PerCancer/" , dir , "/" , label , "_byCancerType.pdf", sep="" ), height= height_2, width= width_2 , bg="transparent" , onefile=FALSE )
			forest( meta.subgroup , 
	   #          leftcols = c("studlab", "effect.ci" , "Pval"),
				# leftlabs= c( "Study" , "logOR [95%CI]" , "P-value" ) , 
				digits.se = 2,
	      		colgap.forest= unit(10, "mm") , 
	      		plotwidth = unit( 30 , "mm") ,
				xlab = "Estimated logOR",
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

	######################################
	######################################
	## Per Sequencing

	remove <- names( table(seq$Sequencing)[ table(seq$Sequencing) %in% c(1,2) ] )
	
	if( length(unique(seq$Sequencing[!seq$Sequencing%in%remove])) > 1 ){
	
		# Get xlim
		m <- c( min( c( 0 , seq$coef ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(seq$coef) ) , na.rm=TRUE) ) + .5 )

		meta <- metagen( 	TE = coef,
							seTE = SE,
							data = seq ,
							studlab = study ,
							fixed = FALSE ,
							random = TRUE ,
							control = list( maxiter = 10000 , stepadj=0.5 ) )

		if( length(remove) > 0 ){

			meta.subgroup <- update.meta(meta, 
										byvar=Sequencing, 
										exclude = seq$Sequencing %in% remove ,
										fixed = FALSE ,
										random = TRUE ,
										control = list( maxiter = 10000 , stepadj=0.5 ) )
		} else{
			meta.subgroup <- update.meta(meta, 
		                             	byvar=Sequencing ,
										fixed = FALSE ,
										random = TRUE ,
										control = list( maxiter = 10000 , stepadj=0.5 ) )
		}

		pdf( paste( prefix , "/PerSequencing/" , dir , "/" , label , "_bySequencing.pdf", sep="" ), height= height_3, width= width_3 , bg="transparent" , onefile=FALSE )
			forest( meta.subgroup , 
	   #          leftcols = c("studlab" , "effect.ci" , "Pval"),
				# leftlabs= c( "Study" , "logOR [95%CI]" , "P-value" ) , 
				digits.se = 2,
	      		colgap.forest=unit(10, "mm") , 
	      		plotwidth = unit( 30 , "mm") ,
				xlab = "Estimated logOR",
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
}
