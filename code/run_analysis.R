
#############################################
#############################################
## Run the nsTMB meta-analysis
print( "Run the nsTMB meta-analysis" )
source( "nsTMB/Run_nsTMB_analysis.R" )

#############################################
#############################################
## Run the signature meta-analysis in the discovery & validation cohorts
print( "Run the signature meta-analysis in the discovery & validation cohorts" )
source( "signatures/Compute_Signature.R" )
source( "signatures_Validation/Compute_Signature.R" )

#############################################
#############################################
## Run the de novo single gene meta-analysis
print( "Run the de novo single gene meta-analysis" )
source( "denovo_Single_Gene/run_denovo_SG.R" )

#############################################
#############################################
## Run the PredictIO meta-analysis in the discovery & validation cohorts
print( "Run the PredictIO meta-analysis in the discovery & validation cohorts" )
source( "PredictIO/run_PredictIO.R" )
source( "PredictIO_Validation/run_PredictIO.R" )


#############################################
#############################################
## Compare nsTMB and PredictIO prognostic value in the discovery & validation cohorts
print( "Compare nsTMB and PredictIO prognostic value in the discovery & validation cohorts" )
source( "nsTMB_PredictIO/AUC_nsTMB_PredictIO.R" )
source( "nsTMB_PredictIO_Validation/INSPIRE_AUC_nsTMB_PredictIO.R" )
source( "nsTMB_PredictIO_Validation/Kim_AUC_nsTMB_PredictIO.R" )

#############################################
#############################################
## Generate the figures of the Resistance gene of interest in the discovery & validation cohorts
print( "Generate the figures of the Resistance gene of interest in the discovery & validation cohorts" )
source( "Resistance_Gene_Validation/run_Resistance_Gene_Validation.R" )

#############################################
#############################################
## run the validation analysis on the mouse model
print( "run the validation analysis on the mouse model" )
source( "mouseModel/run_MouseModel.R" )

#############################################
#############################################
## Generate the summary figures
print( "Generate the summary figures" )
source( "Summary_Figures/create_Directory.R" )

print( "Generate the Tables" )
source( "Summary_Figures/Get_Summary_Table_nsTMB.R" )
source( "Summary_Figures/Get_Summary_Table_Signature.R" )

print( "Generate the nsTMB Heatmap" )
source( "Summary_Figures/Heatmap_nsTMB.R" )
source( "Summary_Figures/Forestplot_nsTMB.R" )

print( "Generate the Signature Network" )
source( "Summary_Figures/Signatures_Network_Affinity_Cluster.R" )
source( "Summary_Figures/Signatures_Network_EnrichR.R" )

print( "Generate the Signature Heatmap" )
source( "Summary_Figures/Heatmap_Signatures.R" )
source( "Summary_Figures/Heatmap_Signatures_PerTumorType.R" )

print( "Generate the AUC plot" )
source( "Summary_Figures/AUC_Plot.R" )

print( "Generate the Corrplot" )
source( "Summary_Figures/Corrplot.R" )
source( "Summary_Figures/Boxplot_Correlation_signature.R" )
