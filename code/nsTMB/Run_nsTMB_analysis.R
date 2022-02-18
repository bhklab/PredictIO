
####################################################
####################################################

file = "../results/nsTMB"

if( dir.exists( file ) ) {
	unlink( file , recursive = TRUE )
}

dir.create( file )

dir.create( paste( file , "/KMPlot" , sep="" ) )
dir.create( paste( file , "/KMPlot/OS" , sep="" ) )
dir.create( paste( file , "/KMPlot/PFS" , sep="" ) )
dir.create( paste( file , "/Funnel" , sep="" ) )
dir.create( paste( file , "/Funnel/OS" , sep="" ) )
dir.create( paste( file , "/Funnel/PFS" , sep="" ) )
dir.create( paste( file , "/Funnel/Response" , sep="" ) )
dir.create( paste( file , "/Overall" , sep="" ) )
dir.create( paste( file , "/Overall/OS" , sep="" ) )
dir.create( paste( file , "/Overall/PFS" , sep="" ) )
dir.create( paste( file , "/Overall/Response" , sep="" ) )
dir.create( paste( file , "/PerCancer" , sep="" ) )
dir.create( paste( file , "/PerCancer/OS" , sep="" ) )
dir.create( paste( file , "/PerCancer/PFS" , sep="" ) )
dir.create( paste( file , "/PerCancer/Response" , sep="" ) )
dir.create( paste( file , "/PerSequencing" , sep="" ) )
dir.create( paste( file , "/PerSequencing/OS" , sep="" ) )
dir.create( paste( file , "/PerSequencing/PFS" , sep="" ) )
dir.create( paste( file , "/PerSequencing/Response" , sep="" ) )

print( paste( "Created directory:" , file ) ) 


####################################################
####################################################

source( "nsTMB/nsTMB_cutoff.R")
source( "nsTMB/nsTMB_compute.R")
source( "nsTMB/nsTMB_meta-analysis.R")
