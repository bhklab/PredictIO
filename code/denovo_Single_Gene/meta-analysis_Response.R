
library(meta)
library(metafor)
library(genefu)


load( "../results/denovo_Single_Gene/Single_Gene_LogReg_Response.RData" ) 

meta_res = NULL
for(i in 1:nrow(coef)){

	print( paste( i , ":" , nrow( coef ) ) )

	data = as.data.frame( cbind( colnames(coef) , coef[ i , ] , se[ i , ] , fdr[ i , ] ) )
	colnames(data) = c( "study" , "coef" , "se" , "fdr" )
	data$coef = as.numeric(as.character( data$coef ))
	data$se = as.numeric(as.character( data$se ))
	data$fdr = as.numeric(as.character( data$fdr )) 
	data = data[ !is.na(data$coef) , ]

	meta <- metagen( TE = coef,
					seTE = se,
					data = data,
					studlab = study ,
					fixed = FALSE ,
					random = TRUE ,
					control = list( maxiter = 10000 , stepadj=0.5 ) )

	meta_res <- rbind( meta_res , c( rownames(coef)[i] , 
				meta$TE.random ,  
				meta$seTE.random ,  
				meta$pval.random , 
				meta$I2 , 
				meta$pval.Q ) )
}


meta_res = as.data.frame(meta_res)
colnames(meta_res) = c( "gene" , "coef" , "se" , "pval" , "I2" , "I2_pval" )
meta_res$gene = as.character( meta_res$gene )
meta_res$coef = as.numeric(as.character( meta_res$coef ))
meta_res$se = as.numeric(as.character( meta_res$se ))
meta_res$pval = as.numeric(as.character( meta_res$pval )) 
meta_res$I2 = as.numeric(as.character( meta_res$I2 ))
meta_res$I2_pval = as.numeric(as.character( meta_res$I2_pval )) 

meta_res$fdr <-  p.adjust( meta_res$pval , method= "fdr" ) 

write.table(meta_res , file= "../results/denovo_Single_Gene/Meta-analysis_Single_Gene_Response.txt" , sep="\t" , quote=FALSE, row.names=FALSE )
save( meta_res , file= "../results/denovo_Single_Gene/Meta-analysis_Single_Gene_Response.RData" ) 
