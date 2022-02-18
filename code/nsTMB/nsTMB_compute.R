##########################################################################################################################################################
##########################################################################################################################################################

library('Biobase')
source('meta/Get_HR.R')
source('meta/Get_DI.R')

load("../data/process_data/ICB_snv_filtered.RData")
load("../results/nsTMB/nsTMB.cutoff.RData")

####################################################
####################################################

snv
study = names(snv)

cox_os = cox_pfs = cox_dicho_median_os = cox_dicho_median_pfs = cox_dicho_10TMB_os = cox_dicho_10TMB_pfs = NULL
di_os = di_pfs = NULL
log_response = log_dicho_median_response = log_dicho_10TMB_response = NULL

for( i in 1:length(study)){

	tumor = names( table( phenoData( snv[[i]] )$primary )[ table( phenoData( snv[[i]] )$primary ) >= 20 ] )

	if( length(tumor) > 0 ){
		for( j in 1:length(tumor)){

			tmb = log10( phenoData(snv[[i]])$nsTMB_perMb[ phenoData(snv[[i]])$primary %in% tumor[j] ] + .5 )

			if( length( tmb[ !is.na( tmb ) ] ) >= 20 ){
				if( length( tmb[ !is.na( phenoData(snv[[i]])$os[ phenoData(snv[[i]])$primary %in% tumor[j] ] ) ] ) >= 20 ){
					## Compute the association nsTMB with OS (36 month cutoff) using Cox Regression Model
					hr = Get_HR_continous( surv=phenoData(snv[[i]])$os[ phenoData(snv[[i]])$primary %in% tumor[j] ] , time=phenoData(snv[[i]])$t.os[ phenoData(snv[[i]])$primary %in% tumor[j] ] ,
											 time_censor=36 , variable= tmb )
					cox_os = rbind( cox_os , 
							c( study[i] , tumor[j] , toupper( phenoData(snv[[i]])$dna[1] ) , length( tmb[ !is.na( phenoData(snv[[i]])$os ) & phenoData(snv[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

					## Compute the association nsTMB with OS (36 month cutoff) using Concordence Index
					cci = Get_DI_continous( surv=phenoData(snv[[i]])$os[ phenoData(snv[[i]])$primary %in% tumor[j] ] , time=phenoData(snv[[i]])$t.os[ phenoData(snv[[i]])$primary %in% tumor[j] ] ,
											 time_censor=36 , variable= tmb )
					di_os = rbind( di_os , 
							c( study[i] , tumor[j] , toupper( phenoData(snv[[i]])$dna[1] ) , length( tmb[ !is.na( phenoData(snv[[i]])$os ) & phenoData(snv[[i]])$primary %in% tumor[j] ] ) ,  cci ) )

				
					## Compute the association nsTMB (median : HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
					hr = Get_HR_dicho( surv=phenoData(snv[[i]])$os[ phenoData(snv[[i]])$primary %in% tumor[j] ] , time=phenoData(snv[[i]])$t.os[ phenoData(snv[[i]])$primary %in% tumor[j] ] ,
											 time_censor=36 , variable= tmb , cutoff= log10( nsTMB.cutoff[ tumor[j] ] + .5 ) ,
											 title = paste( study[i] , tumor[j] , 36 , "OS_Overall" , sep="_" ) , ylab="Overall Survival" , dir="../results/nsTMB/KMPlot/OS" )
					cox_dicho_median_os = rbind( cox_dicho_median_os , 
							c( study[i] , tumor[j] , toupper( phenoData(snv[[i]])$dna[1] ) , length( tmb[ !is.na( phenoData(snv[[i]])$os ) & phenoData(snv[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

					## Compute the association nsTMB (10TMB : HIGH vs LOW) with OS (36 month cutoff) using Cox Regression Model
					hr = Get_HR_dicho( surv=phenoData(snv[[i]])$os[ phenoData(snv[[i]])$primary %in% tumor[j] ] , time=phenoData(snv[[i]])$t.os[ phenoData(snv[[i]])$primary %in% tumor[j] ] ,
											 time_censor=36 , variable= tmb , cutoff= log10(10.5) ,
											 title = paste( study[i] , tumor[j] , 36 , "OS_Overall" , sep="_" ) , ylab="Overall Survival" , dir="../results/nsTMB/KMPlot/OS" )
					cox_dicho_10TMB_os = rbind( cox_dicho_10TMB_os , 
							c( study[i] , tumor[j] , toupper( phenoData(snv[[i]])$dna[1] ) , length( tmb[ !is.na( phenoData(snv[[i]])$os ) & phenoData(snv[[i]])$primary %in% tumor[j] ] ) ,  hr ) )
				}

				if( length( tmb[ !is.na( phenoData(snv[[i]])$pfs[ phenoData(snv[[i]])$primary %in% tumor[j] ] ) ]  ) >= 20 ){
				
					## Compute the association of nsTMB with PFS using Cox Regression Model
					hr = Get_HR_continous( surv=phenoData(snv[[i]])$pfs[ phenoData(snv[[i]])$primary %in% tumor[j] ] , time=phenoData(snv[[i]])$t.pfs[ phenoData(snv[[i]])$primary %in% tumor[j] ] ,
											 time_censor=24 , variable= tmb )
					cox_pfs = rbind( cox_pfs , 
							c( study[i] , tumor[j] , toupper( phenoData(snv[[i]])$dna[1] ) , length( tmb[ !is.na( phenoData(snv[[i]])$pfs ) & phenoData(snv[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

					## Compute the association of nsTMB with PFS using Concordence Index
					cci = Get_DI_continous( surv=phenoData(snv[[i]])$pfs[ phenoData(snv[[i]])$primary %in% tumor[j] ] , time=phenoData(snv[[i]])$t.pfs[ phenoData(snv[[i]])$primary %in% tumor[j] ] ,
											 time_censor=24 , variable= tmb )
					di_pfs = rbind( di_pfs , 
							c( study[i] , tumor[j] , toupper( phenoData(snv[[i]])$dna[1] ) , length( tmb[ !is.na( phenoData(snv[[i]])$pfs ) & phenoData(snv[[i]])$primary %in% tumor[j] ] ) ,  cci ) )
				
					
					## Compute the association of nsTMB (median: HIGH vs LOW) with PFS using Cox Regression Model
					hr = Get_HR_dicho( surv=phenoData(snv[[i]])$pfs[ phenoData(snv[[i]])$primary %in% tumor[j] ] , time=phenoData(snv[[i]])$t.pfs[ phenoData(snv[[i]])$primary %in% tumor[j] ] ,
											 time_censor=24 , variable= tmb , cutoff= log10( nsTMB.cutoff[ tumor[j] ] + .5 ) ,
											 title = paste( study[i] , tumor[j] , 36 , "PFS_Overall" , sep="_" ) , ylab="Progression-Free Survival" , dir="../results/nsTMB/KMPlot/PFS" )
					cox_dicho_median_pfs = rbind( cox_dicho_median_pfs , 
							c( study[i] , tumor[j] , toupper( phenoData(snv[[i]])$dna[1] ) , length( tmb[ !is.na( phenoData(snv[[i]])$pfs ) & phenoData(snv[[i]])$primary %in% tumor[j] ] ) ,  hr ) )

					## Compute the association of nsTMB (10TMB: HIGH vs LOW) with PFS using Cox Regression Model
					hr = Get_HR_dicho( surv=phenoData(snv[[i]])$pfs[ phenoData(snv[[i]])$primary %in% tumor[j] ] , time=phenoData(snv[[i]])$t.pfs[ phenoData(snv[[i]])$primary %in% tumor[j] ] ,
											 time_censor=24 , variable= tmb , cutoff= log10(10.5) ,
											 title = paste( study[i] , tumor[j] , 36 , "PFS_Overall" , sep="_" ) , ylab="Progression-Free Survival" , dir="../results/nsTMB/KMPlot/PFS" )
					cox_dicho_10TMB_pfs = rbind( cox_dicho_10TMB_pfs , 
							c( study[i] , tumor[j] , toupper( phenoData(snv[[i]])$dna[1] ) , length( tmb[ !is.na( phenoData(snv[[i]])$pfs ) & phenoData(snv[[i]])$primary %in% tumor[j] ] ) ,  hr ) )
				}

				
				if( length( tmb[ !is.na( phenoData(snv[[i]])$response[ phenoData(snv[[i]])$primary %in% tumor[j] ] ) ]  ) >= 20 ){

					## Association with Response (Continous)
					x = ifelse( phenoData(snv[[i]])$response[ phenoData(snv[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(snv[[i]])$response[ phenoData(snv[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

					fit = glm( x ~ tmb , family=binomial( link="logit" ) )

					log_response = rbind( log_response , c( study[i] , 
								tumor[j] , 
								toupper( phenoData(snv[[i]])$dna[1] ) , 
								length( tmb[ !is.na( phenoData(snv[[i]])$response[ phenoData(snv[[i]])$primary %in% tumor[j] ] ) ]  ) , 
								round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
								round( confint(fit)[ 2 , ] , 2 ) , 
								summary(fit)$coefficients[ 2 , 4 ] ) )

					## Association with Response (median: HIGH vs Low)
					m = log10( nsTMB.cutoff[ tumor[j] ] + .5 )
					y = ifelse( tmb >= m , 1 ,0 )
					x = ifelse( phenoData(snv[[i]])$response[ phenoData(snv[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(snv[[i]])$response[ phenoData(snv[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

					
					if( sum(y) / length(y) >= .10 ){
						fit = glm( x ~ y , family=binomial( link="logit" ) )

						log_dicho_median_response = rbind( log_dicho_median_response , c( study[i] , 
									tumor[j] , 
									toupper( phenoData(snv[[i]])$dna[1] ) , 
									length( tmb[ !is.na( phenoData(snv[[i]])$response[ phenoData(snv[[i]])$primary %in% tumor[j] ] ) ]  ) , 
									round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
									round( confint(fit)[ 2 , ] , 2 ) , 
									summary(fit)$coefficients[ 2 , 4 ] ) )
					}

					## Association with Response (10TMB: HIGH vs Low)
					m = log10( 10.5 )
					y = ifelse( tmb >= m , 1 ,0 )
					x = ifelse( phenoData(snv[[i]])$response[ phenoData(snv[[i]])$primary %in% tumor[j] ] %in% "R" , 0 , ifelse( phenoData(snv[[i]])$response[ phenoData(snv[[i]])$primary %in% tumor[j] ] %in% "NR" , 1 , NA ) )

					if( sum(y) / length(y) >= .10 ){

						fit = glm( x ~ y , family=binomial( link="logit" ) )

						log_dicho_10TMB_response = rbind( log_dicho_10TMB_response , c( study[i] , 
									tumor[j] , 
									toupper( phenoData(snv[[i]])$dna[1] ) , 
									length( tmb[ !is.na( phenoData(snv[[i]])$response[ phenoData(snv[[i]])$primary %in% tumor[j] ] ) ]  ) , 
									round( summary(fit)$coefficients[ 2 , c( 1 , 2 ) ] , 2 ) , 
									round( confint(fit)[ 2 , ] , 2 ) , 
									summary(fit)$coefficients[ 2 , 4 ] ) )
					}
				}
			}
		}
	}
}

colnames(cox_os) = colnames(cox_pfs) = colnames(cox_dicho_median_os) = colnames(cox_dicho_median_pfs) = colnames(cox_dicho_10TMB_os) = colnames(cox_dicho_10TMB_pfs) = c( "study" , "Primary" , "Sequencing" , "N" , "HR" , "SE" , "95di_low" , "95di_high"  , "Pval" )
colnames(di_os) = colnames(di_pfs) = c( "study" , "Primary" , "Sequencing" , "N" , "DI" , "SE" , "95di_low" , "95di_high"  , "Pval" )
colnames(log_response) = colnames(log_dicho_median_response) = colnames(log_dicho_10TMB_response) = c( "study" , "Primary" , "Sequencing" , "N", "coef" , "SE" , "95di_low" , "95di_high" , "Pval" )


cox_os = as.data.frame( cox_os )
di_os = as.data.frame( di_os )

cox_pfs = as.data.frame( cox_pfs )
di_pfs = as.data.frame( di_pfs )

cox_dicho_median_os = as.data.frame( cox_dicho_median_os )
cox_dicho_10TMB_os = as.data.frame( cox_dicho_10TMB_os )

cox_dicho_median_pfs = as.data.frame( cox_dicho_median_pfs )
cox_dicho_10TMB_pfs = as.data.frame( cox_dicho_10TMB_pfs )

log_response = as.data.frame( log_response )
log_dicho_median_response = as.data.frame( log_dicho_median_response )
log_dicho_10TMB_response = as.data.frame( log_dicho_10TMB_response )

save( log_response , log_dicho_median_response , log_dicho_10TMB_response , file="../results/nsTMB/nsTMB_LogReg_result.RData" )
save( cox_os , cox_pfs , cox_dicho_median_os , cox_dicho_median_pfs , cox_dicho_10TMB_os , cox_dicho_10TMB_pfs , file="../results/nsTMB/nsTMB_COX_result.RData" ) 
save( di_os , di_pfs , file="../results/nsTMB/nsTMB_DI_result.RData" ) 
