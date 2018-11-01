args <- commandArgs(trailingOnly=T)
source("function.r")

e <- args[1]
print(e)
name <- args[2]

r2.cutoff <- 10
ratio.Se.normal.cutoff <- 6
rt.cut <- 10.0
isotope.mass.unit <- 1.0033548
mz.ppm.cut <- 0.00001 #10ppm
charge.states <- seq(6,1)
Hplus <- 1.0072765

cat(paste(name,"\n",sep=""))
run <- name
ms1 <- read.csv(name,stringsAsFactors = F)
cat( paste("detect",dim(ms1)[1], "ms1 spectra in total\n") )
envelop <- data.frame(run=numeric(1),rt=numeric(1),mz=numeric(1),charge=numeric(1),mass= numeric(1),
                     intensity=numeric(1),mass_mono=numeric(1),mz_mono=numeric(1),R2.Se=numeric(1),ratio.Se.normal=(1),max.index=(1),envelope_length=(1))
dir.create(paste(name,"PNG",sep="_"))
if (e == "Se_Pair") {
	envelop <- peak_extraction_judgement_HL(ms1)
	envelop <- envelop[-1,]
	out_matrix <- get_common_mass_single_run(envelop)
	out_matrix <- out_matrix[-1,]
} else {
	envelop <- peak_extraction_judgement(ms1)
	envelop <- envelop[-1,]
	out_matrix <- get_common_mass_single_run(envelop)
	out_matrix <- out_matrix[-1,]
}
colnames(envelope)[c(9,10) <- c("Score1","Score2")
cat("write csv file")
write.csv(envelop,paste(run,"envelop.csv",sep="_"))
write.csv(out_matrix,paste(run,"envelop_valid.csv",sep="_"))
file.remove(name)

