
source("PredictIO/IO_Resistance.R")
source("PredictIO/IO_Sensitive.R")
source("PredictIO/PredictIO.R")
source("PredictIO/Meta-Analysis_PredictIO.R")

#################################################################################
#################################################################################

file = "../results/PredictIO"
	
metascore = c( "IO_Sensitive" , "IO_Resistance", "PredictIO" ) 

if( dir.exists( file ) ) {
	unlink( file , recursive = TRUE )
}

dir.create( file )

	
for( i in 1:length(metascore) ){

	dir.create( paste( file , "/" , metascore[i] , sep="" ) )

	dir.create( paste( file , "/" , metascore[i] , "/KMPlot" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/KMPlot/OS" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/KMPlot/PFS" , sep="" ) )

	dir.create( paste( file , "/" , metascore[i] , "/Funnel" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/Funnel/OS" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/Funnel/PFS" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/Funnel/Response" , sep="" ) )
	
	dir.create( paste( file , "/" , metascore[i] , "/Overall" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/Overall/OS" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/Overall/PFS" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/Overall/Response" , sep="" ) )
	
	dir.create( paste( file , "/" , metascore[i] , "/PerCancer" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/PerCancer/OS" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/PerCancer/PFS" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/PerCancer/Response" , sep="" ) )
	
	dir.create( paste( file , "/" , metascore[i] , "/PerSequencing" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/PerSequencing/OS" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/PerSequencing/PFS" , sep="" ) )
	dir.create( paste( file , "/" , metascore[i] , "/PerSequencing/Response" , sep="" ) )


}

## Per Tumor Type
metascore = c( "IO_Sensitive" , "IO_Resistance", "PredictIO" ) 
tumorID = c( "Melanoma" , "Lung" , "Kidney"  )

dir.create( "../results/Per_TumorType" )

file = "../results/Per_TumorType/PredictIO"

if( dir.exists( file ) ) {
	unlink( file , recursive = TRUE )
}

dir.create( file )

for( i in 1:length(metascore) ){
	dir.create( paste( file , "/" , metascore[i] , sep="" ) )
	
	for( j in 1:length(tumorID) ){
	
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , sep="" ) )

		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/Funnel" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/Funnel/OS" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/Funnel/PFS" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/Funnel/Response" , sep="" ) )
		
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/Overall" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/Overall/OS" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/Overall/PFS" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/Overall/Response" , sep="" ) )
		
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/PerCancer" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/PerCancer/OS" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/PerCancer/PFS" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/PerCancer/Response" , sep="" ) )
		
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/PerSequencing" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/PerSequencing/OS" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/PerSequencing/PFS" , sep="" ) )
		dir.create( paste( file , "/" , metascore[i] , "/" , tumorID[j] , "/PerSequencing/Response" , sep="" ) )
	}
}

#################################################################################
#################################################################################


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

print( "####### IO_Resistance #######" )
Get_IO_resistance( meta_res= meta_res )

print( "####### IO_Sensitive #######" )
Get_IO_sensitive( meta_res= meta_res )

print( "####### PredictIO #######" )
Get_PredictIO( meta_res= meta_res )

print( "####### Meta-Analysis #######" )
Get_Meta_Analysis_PredictIO( )
print( "####### Per Tumor-Type Meta-Analysis #######" )
Get_Meta_Analysis_PredictIO_PerTumor( )

#################################################################################
#################################################################################

