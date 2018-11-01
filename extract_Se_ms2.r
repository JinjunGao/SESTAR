args <- commandArgs(trailingOnly=T)
SESTAR_PATH <- Sys.getenv("SESTAR_PATH")
source(paste(SESTAR_PATH,"/function.r",sep=""))

e <- args[1]
ms2file <- args[2]
print(e)
print(ms2file)


#r2.cutoff <- 15
#ratio.Se.normal.cutoff <- 4
#rt.cut <- 10.0
isotope.mass.unit <- 1.0033548
mz.ppm.cut <- 0.00003 #30ppm
#charge.states <- seq(6,1)
Hplus <- 1.0072765

input.path <- getwd()
#ms2.names <- list.files(path="./",pattern="ms2$")

run_base <- strsplit(ms2file,'\\.')[[1]][1]
ms2_out <- paste(run_base,'_enriched','.ms2',sep='')
ms2 <- read.ms2.all(ms2file)
for (j in 1:length(ms2@name)) {
	count <- peak_extraction_judgement_ms2(ms2,j)
	if (count >= 2) {
		#break
		#Se_ms2 <- c(Se_ms2,ms2@name[j])
		#print(ms2@name[j])
		cat(paste("S\t",ms2@name[j],"\t",ms2@name[j],"\t",ms2@parent_ion[j],"\n",sep=""),file=ms2_out,append=TRUE)
		cat(paste("I\t","RetTime","\t",ms2@rt[j],"\n",sep=""),file=ms2_out,append=TRUE)
		cat(paste("I\t","PrecursorInt","\t",ms2@precursor_int[j],"\n",sep=""),file=ms2_out,append=TRUE)
		cat(paste("I\t","IonInjectionTime","\t",ms2@ion_injection_time[j],"\n",sep=""),file=ms2_out,append=TRUE)
		cat(paste("I\t","ActivationType","\t",ms2@activation_type[j],"\n",sep=""),file=ms2_out,append=TRUE)
		cat(paste("I\t","PrecursorFile","\t",ms2@precursor_file[j],"\n",sep=""),file=ms2_out,append=TRUE)
		cat(paste("I\t","PrecursorScan","\t",ms2@precursor_scan[j],"\n",sep=""),file=ms2_out,append=TRUE)
		cat(paste("I\t","InstrumentType","\t",ms2@instrument_type[j],"\n",sep=""),file=ms2_out,append=TRUE)
		cat(paste("Z\t",ms2@z[j],"\t",ms2@massplus[j],"\n",sep=""),file=ms2_out,append=TRUE)
		for (i in 1:length(ms2@peaks[j][[1]])) {
			cat(paste(ms2@peaks[j][[1]][i],"\t",ms2@intensity[j][[1]][i],"\t",ms2@peaks_z[j][[1]][i],"\n",sep=""),file=ms2_out,append=TRUE)
		}
	}
}




	
