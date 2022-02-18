

##########################################################################################################################################################
##########################################################################################################################################################

library('survcomp')
library('Biobase')

load("../data/process_data/ICB_snv_filtered.RData")

snv
study = names(snv)
tumors = c( "Melanoma" , "Lung" , "Bladder" , "Kidney" , "HNC" , "Brain" , "Colon" , "Unknown" , "Esophagus" , "Ureteral" , "Stomach" , "Breast" , "Eye" )
melanoma = lung = bladder = kidney = HNC = brain = colon = unknown = esophagus = ureteral = stomach = breast = eye = NULL

for( i in 1:length(study)){

	# tumor = names( table( phenoData( snv[[i]] )$primary )[ table( phenoData( snv[[i]] )$primary ) >= 20 ] )
	tumor = names( table( phenoData( snv[[i]] )$primary ) )
	for( j in 1:length(tumor)){

			if(tumor[j] %in% "Melanoma" ){
				melanoma = c( melanoma , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "Lung" ){
				lung = c( lung , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "Bladder" ){
				bladder = c( bladder , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "Kidney" ){
				kidney = c( kidney , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "HNC" ){
				HNC = c( HNC , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "Brain" ){
				brain = c( brain , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "Colon" ){
				colon = c( colon , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "Unknown" ){
				unknown = c( unknown , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "Esophagus" ){
				esophagus = c( esophagus , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "Ureteral" ){
				ureteral = c( ureteral , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "Stomach" ){
				stomach = c( stomach , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "Breast" ){
				breast = c( breast , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
			if(tumor[j] %in% "Eye" ){
				eye = c( eye , phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] )
			}
	}
}

nsTMB.cutoff = c( median( melanoma , na.rm=TRUE) ,
					median( lung , na.rm=TRUE) ,
					median( bladder , na.rm=TRUE) ,
					median( kidney , na.rm=TRUE) ,
					median( HNC , na.rm=TRUE) ,
					median( brain , na.rm=TRUE) ,
					median( colon , na.rm=TRUE) ,
					median( unknown , na.rm=TRUE) ,
					median( esophagus , na.rm=TRUE) ,
					median( ureteral , na.rm=TRUE) ,
					median( stomach , na.rm=TRUE) ,
					median( breast , na.rm=TRUE) ,
					median( eye , na.rm=TRUE)  )
names(nsTMB.cutoff) = tumors

save( nsTMB.cutoff , file="../results/nsTMB/nsTMB.cutoff.RData" ) 

