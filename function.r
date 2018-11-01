read.chem.table <- function( table.name ) {
  orig.table <- read.table(table.name, header=T, sep="\t", comment.char="!")
  named.table <- orig.table[,2:ncol(orig.table)]
  rownames(named.table) <- orig.table[,1]
  return(named.table)
}

init.atom.mass <- function() {
  atom.mass.vec <- c(12.000000, #C
                     1.007825,  #H
                     15.994915, #O
                     14.003074, #N
                     31.972072, #S
                     30.973763, #P
                     15.000109, #N15
                     2.014102,  #H2
                     13.003355, #C13
                     1.0072765,  #Hplus
                     34.96885, #Chlorine
                     78.91833, #Bromine
		     79.9165213 #Selenium
                     )
  names(atom.mass.vec) <- c("C","H","O","N","S","P","N15","H2","C13","Hplus","Cl","Br","Se")
  return(atom.mass.vec)
}

init.aa.mass <- function(atom.mass.vec, chem.table ) {
  aa.names <- rownames(chem.table)
  chem.names <- colnames(chem.table)
  aa.mass.vec <- rep(0, length(aa.names))
  names(aa.mass.vec) <- aa.names
  for (aa in aa.names ) {
    mass <- 0.0
    for (chem in chem.names) {
      mass <- mass + atom.mass.vec[chem]*chem.table[aa,chem] #C * 
    }
    aa.mass.vec[aa] <- mass
  }
  return(aa.mass.vec)
}

calc.peptide.mass <- function(sequence, aa.mass.vec) {
	mass <- 0
	if (sequence != '') {
		peptide.vec <- unlist(strsplit(sequence,"",fixed=T))
		defined.aas <- names(aa.mass.vec)
		for ( aa in peptide.vec ) {
			if ( aa %in% defined.aas ) {
			mass <- mass + aa.mass.vec[aa]
			}
		}
		mass <- mass + aa.mass.vec["NTERM"] + aa.mass.vec["CTERM"]
	}
	if (sequence == '') mass <- -611.30273
	return(mass)
}

averagine.count <- function(input.mass) {
	averagine.mass <- 111.1254
	elements <- c( "C", "H", "O", "N")
	averagine.comp <- c( 4.9348, 7.7583, 1.4773, 1.3577)
	names(averagine.comp) <- elements
	return( round(averagine.comp*(input.mass/averagine.mass)) )
}

isotope.dist <- function(elements.count, N15.enrichment=1.0) {
  elements <- c( "C", "H", "O", "N", "S", "P", "N15", "H2", "C13", "O18" , "Se", "Cl","Br")
  elements.count.local <- rep(0, length(elements))
  names(elements.count.local) <- elements
  for ( e in names(elements.count)) {
    elements.count.local[e] <- elements.count[e]
  }
  heavy <- c(1.10, 0.015, 0.20, 0.37, 4.21, 0, 100, 100, 100, 100 , 49.16, 24.24, 49.31)/100
  names(heavy) <- elements
  heavy["N15"] <- N15.enrichment

  light <- 1.00 - heavy
  light["Se"] <- 0.2377 # Special for selenium
  #names(elements.count) <- elements
  single.prob <- as.list( elements )
  names(single.prob) <- elements
  all.prob <- numeric(0)
  for ( e in elements ) {
    count <- elements.count.local[e]
    if (count == 0) next
    v <- seq(0,count)
    l <- light[e]
    h <- heavy[e]
    new.prob <- single.prob[[e]] <- round(choose(count,v)*(l^(count-v))*(h^v),4)
    if (e =="O" | e =="S" | e=="Cl" | e=="Br" | e== "O18") { # O, S, Cl and Br isotopes are 2 Da more
      new.prob <- rep(0,2*count+1)
      for( i in 1:(count+1)) {
        new.prob[(i-1)*2+1] <- single.prob[[e]][i]
      }
      single.prob[[e]] <- new.prob
    }
	if (e =="Se") {# Se has the special isotopic distribution
	  single.prob[[e]] <- c(0.0089,0.0000,0.0937,0.0763,0.2377,0.0000,0.4961,0.0000,0.0873)
	  new.prob <- single.prob[[e]]
	}
    #print(single.prob[[e]])
    if ( length(all.prob) > 0 ) {
      d <- length(all.prob)-length(new.prob)
      if ( d > 0) {
        new.prob <- c(new.prob,rep(0,d))
      } else if ( d< 0 ) {
        all.prob <- c( all.prob,rep(0,-d) )
      }
      all.prob <- round(convolve(all.prob, rev(new.prob),type="o"),4)
    } else {
      all.prob <- single.prob[[e]]
    }
  }
  return(all.prob)
##  return(single.prob)
}

ms2.checker <- function(envelop.mz,rt, ms2list) {
	rt <- as.numeric(rt)
	#mz_ppm_cut <- 0.000025
	#mzXML_files <- mzXML.files[[name]]
	#ms2_mz <- mzXML_files@msnPrecursorMz
	#ms2_id <- mzXML_files@msnAcquisitionNum
	ms2_mz <- as.numeric(ms2list@parent_ion)
	ms2_id <- ms2list@name
	ms2_rt <- as.numeric(ms2list@rt)
	## 6 minutes window
	rt_min_boundary <- rt/60 - 3
	rt_max_boundary <- rt/60 + 3
	scan_num_list <- c()
	for (mz in envelop.mz) {
		if (any(abs(ms2_mz - mz) < mz * 2 * mz.ppm.cut)) {
		scan_num <- ms2_id[which(abs(ms2_mz - mz) < mz * 2 * mz.ppm.cut)]
		scan_time <- ms2_rt[which(abs(ms2_mz - mz) < mz * 2 * mz.ppm.cut)]
			if (any(scan_time > rt_min_boundary & scan_time < rt_max_boundary)) {
				scan_num_temp <- scan_num[which(scan_time > rt_min_boundary & scan_time < rt_max_boundary)]
				scan_num_list <- c(scan_num_list,scan_num_temp)
			}
		}
	}
	return(scan_num_list)
}

read.ms1 <- function(ms1file) {
	setClass("ms1class", representation(name="list",rt = "list", peaks="list",intensity="list"))
	ms1file <- file(ms1file, "r")
	line=readLines(ms1file,n=1)
	ms1list <- new("ms1class", name=list(), rt=list(), peaks=list(), intensity=list())
	while( length(line) != 0 ) {
		line <- unlist(strsplit(line,"\\s"))
		##if this is a numeric line, conver them to numbers
		if (!is.na(as.numeric(line))) {
			line <- as.numeric(line)
		}
		if (line[1] == "S") {
			if (exists("spectrum_name")) {
				ms1list@name[[length(ms1list@name)+1]] <- spectrum_name
				ms1list@rt[[length(ms1list@rt)+1]] <- spectrum_rt
				ms1list@peaks[[length(ms1list@peaks)+1]] <- spectrum_peaks
				ms1list@intensity[[length(ms1list@intensity)+1]] <- spectrum_intensity
			}
			spectrum_name <- line[2]
			spectrum_peaks <- c()
			spectrum_intensity <- c()
		}
		if (line[2] == "RetTime") {
			spectrum_rt <- line[3]
		}
		if (is.numeric(line[1])) {
			spectrum_peaks <- append(spectrum_peaks, line[1])
			spectrum_intensity <- append(spectrum_intensity, line[2])
		} 
		line=readLines(ms1file,n=1)
	}
	close(ms1file)
	return(ms1list)
}

envelope.rating <- function(envelop.mz, envelop.charge, envelop.intensity) {
	if (any(envelop.mz == 0)) {
			zero.pos <- which(envelop.mz==0)
			for (i in zero.pos) {
			envelop.mz[i] <- (envelop.mz[i+1] + envelop.mz[i-1])/2
			}
		}
	envelop.mass <- max(envelop.mz) * envelop.charge
	max.index <- which.max(envelop.intensity)
	if (e=="Se") {
		elements.count <- averagine.count(envelop.mass)
		elements.count_Se <- averagine.count(envelop.mass-79.9165)
		#Se.prob <- c(0.0089,0.0000,0.0937,0.0763,0.2377,0.0000,0.4961,0.0000,0.0873)
		Se.prob <- c(0.0089,0.0000,0.0937,0.0763,0.2377,0.0000,0.4961,0.0000,0.0873)
		normal.prob <- isotope.dist(elements.count)[1:15]
		withSe.prob <- round(convolve(isotope.dist(elements.count_Se), rev(Se.prob),type="o"),4)[1:15]
		envelop_length <- length(envelop.intensity)
		envelop.length.diff <- 15 - envelop_length ## length diffence between experimental got envelop and simulated one
		Se.R2 <- normal.R2 <- rep(0,envelop.length.diff+1)
		for (i in 1:(envelop.length.diff+1))  {
			experimentint <- envelop.intensity
			simulatedint.Se <- withSe.prob[i:(length(envelop.intensity)-1+i)]
			simulatedint.Se[which(envelop.intensity == 0)] <- 0
			Se.R2[i] <- acos(sum(experimentint*simulatedint.Se) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint.Se*simulatedint.Se))))/pi*180
			simulatedint.normal <- normal.prob[i:(length(envelop.intensity)-1+i)]
			simulatedint.normal[which(envelop.intensity == 0)] <- 0
			normal.R2[i] <- acos(sum(experimentint*simulatedint.normal) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint.normal*simulatedint.normal))))/pi*180
		}
		R2.Se.index <- which.min(Se.R2)
		R2.Se <- Se.R2[R2.Se.index]
		R2.normal.index <- which.min(normal.R2)
		R2.normal <- normal.R2[R2.normal.index]
		ratio.Se.normal <- R2.Se/R2.normal
		results <- list(round(R2.Se,2),round(ratio.Se.normal,2),withSe.prob,normal.prob,max.index,envelop_length, R2.Se.index,R2.normal.index)
	}
	if (e == "Se74") {
		if (max.index != 7) {
			results <- list(0,0,0,0,0,0,0,0)
		} else {
			elements.count <- averagine.count(envelop.mass)
			elements.count_Se <- averagine.count(envelop.mass-79.9165)
			Se.prob <- c(0.0089,0.0000,0.0937,0.0763,0.2377,0.0000,0.4961,0.0000,0.0873)
			normal.prob <- isotope.dist(elements.count)[1:15]
			withSe.prob <- round(convolve(isotope.dist(elements.count_Se), rev(Se.prob),type="o"),4)[1:15]
			envelop_length <- length(envelop.intensity)
			experimentint <- envelop.intensity[3:envelop_length]
			simulatedint.Se <- withSe.prob[3:envelop_length]
			simulatedint.Se[which(envelop.intensity == 0)] <- 0
			R2.Se <- acos(sum(experimentint*simulatedint.Se) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint.Se*simulatedint.Se))))/pi*180
			simulatedint.normal <- normal.prob[3:envelop_length]
			simulatedint.normal[which(envelop.intensity == 0)] <- 0
			R2.normal <- acos(sum(experimentint*simulatedint.normal) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint.normal*simulatedint.normal))))/pi*180
			ratio.Se.normal <- R2.Se/R2.normal
			results <- list(round(R2.Se,2),round(ratio.Se.normal,2),withSe.prob,normal.prob,max.index,envelop_length,1,1)
		}
	}
	return(results)
}

envelope.rating.HL <- function(envelop.mz, envelop.charge, envelop.intensity) {
		envelop.mass <- max(envelop.mz) * envelop.charge
		max.index <- which.max(envelop.intensity)
		index_heavy <- 0
		if (max.index >= 9) {
			index_heavy <- max.index
			index_light <- max.index -6
		}
		if (max.index <= 7) {
			index_light <- max.index
			index_heavy <- max.index+6
		}
		if(index_heavy != 0 & index_heavy <= length(envelop.intensity)) {
			ratio_HL <- envelop.intensity[index_heavy]/envelop.intensity[index_light]
			elements.count <- averagine.count(envelop.mass)
			elements.count_H <- c(elements.count["C"]-6,elements.count[c("H","O","N")],6)
			names(elements.count_H)[5] <- "C13"
			elements.count_Se_L <- averagine.count(envelop.mass-79.9165)
			elements.count_Se_H <- c(elements.count_Se_L["C"]-6,elements.count_Se_L[c("H","O","N")],6)
			names(elements.count_Se_H)[5] <- "C13"
			normal.prob_L <- isotope.dist(elements.count)[1:25]
			normal.prob_H <- isotope.dist(elements.count_H)[1:25]*ratio_HL
			normal.prob <- normal.prob_L + normal.prob_H
			Se.prob <- c(0.0089,0.0000,0.0937,0.0763,0.2377,0.0000,0.4961,0.0000,0.0873)
			withSe.prob_L <- round(convolve(isotope.dist(elements.count_Se_L), rev(Se.prob),type="o"),4)[1:25]
			withSe.prob_H <- round(convolve(isotope.dist(elements.count_Se_H), rev(Se.prob),type="o"),4)[1:25]*ratio_HL
			withSe.prob <- withSe.prob_H + withSe.prob_L
			envelop_length <- length(envelop.intensity)
			envelop.length.diff <- 25 - envelop_length ## length diffence between experimental got envelop and simulated one
			Se.R2 <- normal.R2 <- rep(0,envelop.length.diff+1)
			for (i in 1:(envelop.length.diff+1))  {
				experimentint <- envelop.intensity
				simulatedint.Se <- withSe.prob[i:(length(envelop.intensity)-1+i)]
				simulatedint.Se[which(envelop.intensity == 0)] <- 0
				Se.R2[i] <- acos(sum(experimentint*simulatedint.Se) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint.Se*simulatedint.Se))))/pi*180
				
				experimentint <- envelop.intensity
				simulatedint.Se <- withSe.prob[i:(length(envelop.intensity)-1+i)]
				simulatedint.Se[which(envelop.intensity == 0)] <- 0
				Se.R2[i] <- acos(sum(experimentint*simulatedint.Se) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint.Se*simulatedint.Se))))/pi*180
				simulatedint.normal <- normal.prob[i:(length(envelop.intensity)-1+i)]
				simulatedint.normal[which(envelop.intensity == 0)] <- 0
				normal.R2[i] <- acos(sum(experimentint*simulatedint.normal) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint.normal*simulatedint.normal))))/pi*180
			}
			R2.Se.index <- which.min(Se.R2)
			R2.Se <- Se.R2[R2.Se.index]
			R2.normal.index <- which.min(normal.R2)
			R2.normal <- normal.R2[R2.normal.index]
			ratio.Se.normal <- R2.Se/R2.normal
			results <- list(round(R2.Se,2),round(ratio.Se.normal,2),withSe.prob,normal.prob,max.index,envelop_length, R2.Se.index,R2.normal.index,index_light)
		} else {
			results <- list(0,0,0,0,0,0,0,0,0)
		}
	return(results)
}

peak_extraction_judgement_HL <- function(ms1) {
	for (j in 1:length(ms1$Scan_No)) {
		#if (a == "T") break
		if (j %%100==0) cat(j)
		mz.this.rt.save <- mz.this.rt <- as.numeric(unlist(strsplit(ms1$ions_mz_list[j],";")))
		intensity.this.rt.save <- intensity.this.rt <- as.numeric(unlist(strsplit(ms1$ions_int_list[j],";")))
		rt <- ms1$RT[[j]]
		for (i in charge.states) {
			#if (a == "T") break
			mz.this.rt <- mz.this.rt.save
			intensity.this.rt <- intensity.this.rt.save
			mz.this.rt.charge <- mz.this.rt * i
			#at the begining, the candidate envelope only contains the first one of the mz list, so the envelop.end is 1
			#,and test wether the second is belong to this envelope, so the test.end is 2
			envelop.end <- 1
			envelop.mass <- c(mz.this.rt.charge[1])
			#set the candidate envelope member as NA, in order to delete it from the list after extracting the whole envelope
			#so that the computer could save its source and time.
			mz.this.rt.charge[1] <- NA
			mz.this.rt[1] <- NA
			envelop.intensity <- intensity.this.rt[1]
			test.end <- 2
			while(TRUE) {
				# if the tested mz index is out of the mz list or the testing mz index is too far away from the envelope end, finish doing this envelope extracting.
				if (test.end > length(mz.this.rt) | test.end - envelop.end > 10)  {
					#test whether it is a valid envelope
					if (length(envelop.mass) > 9 & length(envelop.mass) < 25) {
						##use match here to get the index of the envelop.mz in the mz.this.rt.save
						## if it is a real envleope, delete them from the saved list
						intensity.this.rt.save <- intensity.this.rt.save[-match(envelop.mass[which(envelop.mass!=0)],mz.this.rt.save * i)]
						mz.this.rt.save <- mz.this.rt.save[-match(envelop.mass[which(envelop.mass!=0)],mz.this.rt.save * i)]
						envelop.mz <- envelop.mass / i
						envelop.temp.mass <- paste(envelop.mass, collapse = ";")
						envelop.temp.mz <- paste(envelop.mz, collapse = ";")
						envelop.temp.intensity <- paste(envelop.intensity, collapse = ";")
						charge <- i
						results <- envelope.rating.HL(envelop.mz,charge,envelop.intensity)
						#print(results)
						#if (results == "F") break
						if (results[[1]] < r2.cutoff & results[[2]] < ratio.Se.normal.cutoff ) {
							mz_mono <- envelop.mz[results[[9]]]
							mass_mono <- envelop.mass[results[[9]]] - i *Hplus
							envelop.temp <- list(run,rt,envelop.temp.mz,charge,envelop.temp.mass,envelop.temp.intensity, mass_mono,mz_mono, results[[1]], results[[2]],results[[9]],results[[6]])
							envelop <- rbind(envelop,envelop.temp)
							out.picturename <- paste(results[[9]],results[[6]],mass_mono,mz_mono,rt,charge,results[[1]],results[[2]],".png",sep = "_")
							out.picturepath <- paste(name,"PNG",sep="_")
							png(paste(out.picturepath,out.picturename,sep ="/"),height = 400, width = 350, pointsize = 16)
							#layout(matrix(c(1,2),nr=1))
							ylimit <- c(0,max(envelop.intensity)*1.2)
							enlarge_factor_se <- max(envelop.intensity)/max(results[[3]][results[[7]]:(results[[7]]+length(envelop.mz)-1)])
							enlarge_factor_normal <- max(envelop.intensity)/max(results[[4]][results[[8]]:(results[[8]]+length(envelop.mz)-1)])
							plot(envelop.mz,envelop.intensity,type = "h", col = "black",xlab="m/z",ylab="Intensity",lwd=2,ylim=ylimit)
							par(new=T)
							plot(envelop.mz,(results[[3]][results[[7]]:(results[[7]]+length(envelop.mz)-1)])*enlarge_factor_se, col = "red",ylab="",xlab="",axes=F,pch=16,ylim=ylimit)
							par(new=T)
							plot(envelop.mz,(results[[4]][results[[8]]:(results[[8]]+length(envelop.mz)-1)])*enlarge_factor_normal, col = "springgreen2",ylab="",xlab="",axes=F,pch=16,ylim=ylimit)
							dev.off()
						}
					}
					#if (a == "T") break
					#delete the mz index whether it is a real or a pseudo envelope 
					intensity.this.rt <- intensity.this.rt[-which(is.na(mz.this.rt))]
					mz.this.rt <- mz.this.rt[-which(is.na(mz.this.rt))]
					mz.this.rt.charge <- mz.this.rt.charge[-which(is.na(mz.this.rt.charge))]
					## set the candidate envelope as the original status
					envelop.end <- 1
					envelop.mass <- c(mz.this.rt.charge[1])
					mz.this.rt.charge[1] <- NA
					mz.this.rt[1] <- NA
					envelop.intensity <- intensity.this.rt[1]
					test.end <- 2
				}
				#if (a == "T") break
				if (length(mz.this.rt) <= 5) break
				if (abs(abs(mz.this.rt.charge[test.end]-envelop.mass[length(envelop.mass)])-isotope.mass.unit) < 
					mz.ppm.cut * (mz.this.rt.charge[test.end] + envelop.mass[length(envelop.mass)]) / 2.0 ) {
					# if the differrence between the testing mass and the envelope end conforms with the iosotope mass unit,add it the the envelope and replace it with NA
					envelop.mass <- c(envelop.mass,mz.this.rt.charge[test.end])
					#print(envelop.mass)
					envelop.intensity <- c(envelop.intensity, intensity.this.rt[test.end])
					mz.this.rt.charge[test.end] <- NA
					mz.this.rt[test.end] <- NA
					envelop.end <- test.end
					test.end <- test.end + 1
					next
				} 
				test.end <- test.end + 1
			}
		}
	}
	return(envelop)
}

peak_extraction_judgement <- function(ms1) {
	for (j in 1:length(ms1$Scan_No)) {
		#if (a == "T") break
		if (j %%100==0) cat(j)
		mz.this.rt.save <- mz.this.rt <- as.numeric(unlist(strsplit(ms1$ions_mz_list[j],";")))
		intensity.this.rt.save <- intensity.this.rt <- as.numeric(unlist(strsplit(ms1$ions_int_list[j],";")))
		rt <- ms1$RT[[j]]
		for (i in charge.states) {
			#if (a == "T") break
			mz.this.rt <- mz.this.rt.save
			intensity.this.rt <- intensity.this.rt.save
			mz.this.rt.charge <- mz.this.rt * i
			#at the begining, the candidate envelope only contains the first one of the mz list, so the envelop.end is 1
			#,and test wether the second is belong to this envelope, so the test.end is 2
			envelop.end <- 1
			envelop.mass <- c(mz.this.rt.charge[1])
			#set the candidate envelope member as NA, in order to delete it from the list after extracting the whole envelope
			#so that the computer could save its source and time.
			mz.this.rt.charge[1] <- NA
			mz.this.rt[1] <- NA
			envelop.intensity <- intensity.this.rt[1]
			test.end <- 2
			while(TRUE) {
				# if the tested mz index is out of the mz list or the testing mz index is too far away from the envelope end, finish doing this envelope extracting.
				if (test.end > length(mz.this.rt) | test.end - envelop.end > 10)  {
					#test whether it is a valid envelope
					if (length(envelop.mass) > 5 & length(envelop.mass) < 15) {
						##use match here to get the index of the envelop.mz in the mz.this.rt.save
						## if it is a real envleope, delete them from the saved list
						intensity.this.rt.save <- intensity.this.rt.save[-match(envelop.mass[which(envelop.mass!=0)],mz.this.rt.save * i)]
						mz.this.rt.save <- mz.this.rt.save[-match(envelop.mass[which(envelop.mass!=0)],mz.this.rt.save * i)]
						envelop.mz <- envelop.mass / i
						mz_mono <- envelop.mz[which.max(envelop.intensity)]
						mass_mono <- envelop.mass[which.max(envelop.intensity)] - i*Hplus
						envelop.temp.mass <- paste(envelop.mass, collapse = ";")
						envelop.temp.mz <- paste(envelop.mz, collapse = ";")
						envelop.temp.intensity <- paste(envelop.intensity, collapse = ";")
						charge <- i
						results <- envelope.rating(envelop.mz,charge,envelop.intensity)
						#print(results)
						#if (results == "F") break
						if (results[[1]] < r2.cutoff & results[[2]] < ratio.Se.normal.cutoff & results[[5]] > 3 & results[[6]] > 5 & (results[[6]]-results[[5]]) >1 ) {
							#a <- "T"
							envelop.temp <- list(run,rt,envelop.temp.mz,charge,envelop.temp.mass,envelop.temp.intensity, mass_mono,mz_mono, results[[1]], results[[2]],results[[5]],results[[6]])
							envelop <- rbind(envelop,envelop.temp)
							out.picturename <- paste(results[[6]],results[[5]],mass_mono,mz_mono,rt,charge,results[[1]],results[[2]],".png",sep = "_")
							out.picturepath <- paste(name,"PNG",sep="_")
							png(paste(out.picturepath,out.picturename,sep ="/"),height = 400, width = 350, pointsize = 16)
							#layout(matrix(c(1,2),nr=1))
							ylimit <- c(0,max(envelop.intensity)*1.2)
							enlarge_factor_se <- max(envelop.intensity)/max(results[[3]][results[[7]]:(results[[7]]+length(envelop.mz)-1)])
							enlarge_factor_normal <- max(envelop.intensity)/max(results[[4]][results[[8]]:(results[[8]]+length(envelop.mz)-1)])
							plot(envelop.mz,envelop.intensity,type = "h", col = "black",xlab="m/z",ylab="Intensity",lwd=2,ylim=ylimit)
							par(new=T)
							plot(envelop.mz,(results[[3]][results[[7]]:(results[[7]]+length(envelop.mz)-1)])*enlarge_factor_se, col = "red",ylab="",xlab="",axes=F,pch=16,ylim=ylimit)
							par(new=T)
							plot(envelop.mz,(results[[4]][results[[8]]:(results[[8]]+length(envelop.mz)-1)])*enlarge_factor_normal, col = "springgreen2",ylab="",xlab="",axes=F,pch=16,ylim=ylimit)
							dev.off()
						}
					}
					#if (a == "T") break
					#delete the mz index whether it is a real or a pseudo envelope 
					intensity.this.rt <- intensity.this.rt[-which(is.na(mz.this.rt))]
					mz.this.rt <- mz.this.rt[-which(is.na(mz.this.rt))]
					mz.this.rt.charge <- mz.this.rt.charge[-which(is.na(mz.this.rt.charge))]
					## set the candidate envelope as the original status
					envelop.end <- 1
					envelop.mass <- c(mz.this.rt.charge[1])
					mz.this.rt.charge[1] <- NA
					mz.this.rt[1] <- NA
					envelop.intensity <- intensity.this.rt[1]
					test.end <- 2
				}
				#if (a == "T") break
				if (length(mz.this.rt) <= 5) break
				if (abs(abs(mz.this.rt.charge[test.end]-envelop.mass[length(envelop.mass)])-isotope.mass.unit) < 
					mz.ppm.cut * (mz.this.rt.charge[test.end] + envelop.mass[length(envelop.mass)]) / 2.0 ) {
					# if the differrence between the testing mass and the envelope end conforms with the iosotope mass unit,add it the the envelope and replace it with NA
					envelop.mass <- c(envelop.mass,mz.this.rt.charge[test.end])
					#print(envelop.mass)
					envelop.intensity <- c(envelop.intensity, intensity.this.rt[test.end])
					mz.this.rt.charge[test.end] <- NA
					mz.this.rt[test.end] <- NA
					envelop.end <- test.end
					test.end <- test.end + 1
					next
				} 
				test.end <- test.end + 1
			}
		}
	}
	return(envelop)
}

get_common_mass_single_run <- function(envelop) {
	out_matrix <- data.frame(run=numeric(1),count=numeric(1),rt=numeric(1),rt_mean=numeric(1),rt_variance=numeric(1),mass=numeric(1),mass_mean=numeric(1),mass_variance=numeric(1),Score1=numeric(1),
				Score1_mean=numeric(1),Score1_variance=numeric(1),Score2 = numeric(1),Score2_mean=numeric(1),Score2_variance=numeric(1),charge=numeric(1),charge_max=numeric(1),mz=numeric(1),envelope_length=numeric(1),max_index=numeric(1))
	z_candidate <- envelop[,'charge']
	rt_candidate <- envelop[,'rt']
	mass_candidate <- envelop[,'mass_mono']
	R2_candidate <- envelop[,'R2.Se']
	RR2_candidate <- envelop[,'ratio.Se.normal']
	mz_candidate <- envelop[,'mz_mono']
	envelope_length_candidate <- envelop[,'envelope_length']
	max_index_candidate <- envelop[,'max.index']
	while (TRUE) {
		if (length(mass_candidate) <= 1) break
		mass_this_i <- mass_candidate[1]
		rt_this_i <- rt_candidate[1]
		if (any(abs(mass_candidate[-1] - mass_this_i)/mass_this_i <mz.ppm.cut)) {
			mass_match_pos <- which(abs(mass_candidate - mass_this_i)/mass_this_i <mz.ppm.cut)
			if (any(abs(rt_candidate[mass_match_pos[-1]]-rt_this_i) < 2)) {
				rt_mass_match_pos <- mass_match_pos[which(abs(rt_candidate[mass_match_pos]-rt_this_i) < 2)]
				count <- length(rt_mass_match_pos)
				mass_pick <- paste(mass_candidate[rt_mass_match_pos],collapse=';')
				mass_pick_mean <- mean(mass_candidate[rt_mass_match_pos])
				mass_pick_variance <-sd(mass_candidate[rt_mass_match_pos])
				rt_pick <- paste(rt_candidate[rt_mass_match_pos],collapse=';')
				rt_pick_mean <- mean(rt_candidate[rt_mass_match_pos])
				rt_pick_variance <- sd(rt_candidate[rt_mass_match_pos])
				z_pick <- paste(z_candidate[rt_mass_match_pos],collapse=';')
				z_pick_max <- as.numeric(names(table(z_candidate[rt_mass_match_pos]))[table(z_candidate[rt_mass_match_pos]) == max(table(z_candidate[rt_mass_match_pos]))])[1]
				R2_pick <- paste(R2_candidate[rt_mass_match_pos], collapse=';')
				R2_pick_mean <- mean(R2_candidate[rt_mass_match_pos])
				R2_pick_variance <- sd(R2_candidate[rt_mass_match_pos])
				RR2_pick <- paste(RR2_candidate[rt_mass_match_pos],collapse=';')
				RR2_pick_mean <- mean(RR2_candidate[rt_mass_match_pos])
				RR2_pick_variance <- sd(RR2_candidate[rt_mass_match_pos])
				mz_pick <- paste(mz_candidate[rt_mass_match_pos],collapse=';')
				envelope_length_pick <- paste(envelope_length_candidate[rt_mass_match_pos],collapse=';')
				max_index_pick <- paste(max_index_candidate[rt_mass_match_pos],collapse=';')
				new_list <- list(run,count,rt_pick,rt_pick_mean,rt_pick_variance,mass_pick,mass_pick_mean,mass_pick_variance,
								R2_pick,R2_pick_mean,R2_pick_variance,RR2_pick,RR2_pick_mean,RR2_pick_variance,z_pick,z_pick_max,mz_pick,envelope_length_pick,max_index_pick)
				out_matrix <- rbind(out_matrix,new_list)
				z_candidate <- z_candidate[-rt_mass_match_pos]
				rt_candidate <- rt_candidate[-rt_mass_match_pos]
				mass_candidate <- mass_candidate[-rt_mass_match_pos]
				R2_candidate <- R2_candidate[-rt_mass_match_pos]
				RR2_candidate <- RR2_candidate[-rt_mass_match_pos]
				mz_candidate <- mz_candidate[-rt_mass_match_pos]
				envelope_length_candidate <- envelope_length_candidate[-rt_mass_match_pos]
				max_index_candidate <- max_index_candidate[-rt_mass_match_pos]
			} else {
				z_candidate <- z_candidate[-1]
				rt_candidate <- rt_candidate[-1]
				mass_candidate <- mass_candidate[-1]
				R2_candidate <- R2_candidate[-1]
				RR2_candidate <- RR2_candidate[-1]
				mz_candidate <- mz_candidate[-1]
				envelope_length_candidate <- envelope_length_candidate[-1]
				max_index_candidate <- max_index_candidate[-1]
			}
		} else {
			z_candidate <- z_candidate[-1]
			rt_candidate <- rt_candidate[-1]
			mass_candidate <- mass_candidate[-1]
			R2_candidate <- R2_candidate[-1]
			RR2_candidate <- RR2_candidate[-1]
			mz_candidate <- mz_candidate[-1]
			envelope_length_candidate <- envelope_length_candidate[-1]
			max_index_candidate <- max_index_candidate[-1]
		}
	}
	return(out_matrix)
}

get_envelope_from_one_peak <- function(mz,mz_list,int,int_list,z) {
	envelop_mz <- c(mz)
	envelop_int <- c(int)
	while (TRUE) {
		if (any(abs(z*(mz_list - mz) - isotope.mass.unit) < mz.ppm.cut * mz*z )) {
			envelop_end <- which(abs(z*(mz_list - mz) - isotope.mass.unit) < mz.ppm.cut * mz*z )
			envelop_mz <- c(envelop_mz,mz_list[envelop_end])
			envelop_int <- c(envelop_int,int_list[envelop_end])
			mz <- mz_list[envelop_end]
		} else break
	}
	envelop <- list(envelop_mz,envelop_int)
	return(envelop)
}

pair_finder <- function (ms1) {
	envelop_pair <- data.frame(mz_L=numeric(1),mz_H=numeric(1),rt=numeric(1),ms1_num=numeric(1),envelop_mz=numeric(1),envelop_int=numeric(1),charge= numeric(1), score=numeric(1))
	for (j in 1:length(ms1@name)) {
		if (j %%100==0) cat(j)
		mz.this.rt.save <- mz.this.rt <- ms1@peaks[[j]]
		intensity.this.rt.save <- intensity.this.rt <- ms1@intensity[[j]]
		rt <- ms1@rt[[j]]
		ms1_num <- ms1@name[[j]]
		for (z in charge.states) {
			mz.this.rt <- mz.this.rt.save
			intensity.this.rt <- intensity.this.rt.save
			while(TRUE) {
				#if (mz.this.rt[1]==a) break
				#print(length(mz.this.rt))
				if (length(mz.this.rt) < 6) break
				if (any(abs(z*(mz.this.rt - mz.this.rt[1])-pair_delta) < mz.ppm.cut * mz.this.rt[1]*z )) {
					pos <- which(abs(z*(mz.this.rt - mz.this.rt[1])-pair_delta) < mz.ppm.cut * mz.this.rt[1]*z )
					if (length(pos) >1) {
						pos <- pos[1]
					}
					charge <- z
					mz_L <- mz.this.rt[1]
					#if (mz_L == a) break
					mz_H <- mz.this.rt[pos]
					int_L <- intensity.this.rt[1]
					int_H <- intensity.this.rt[pos]
					index_L <- which(mz.this.rt.save==mz_L)
					index_H <- which(mz.this.rt.save==mz_H)
					pair_mz_vector <- c(mz_L,mz_H)
					pair_int_vector <- c(int_L,int_H)
					if (int_L/int_H<0.7 | int_L/int_H>1.3 | index_H-index_L < 6) {
						mz.this.rt <- mz.this.rt[-match(mz.this.rt[1],mz.this.rt)]
						intensity.this.rt <- intensity.this.rt[2:length(intensity.this.rt)]
						next
					}
					envelop_L <- get_envelope_from_one_peak(pair_mz_vector[1],mz.this.rt,pair_int_vector[1],intensity.this.rt,z)
					#if (any(abs(envelop_L[[1]]-a)<0.001)) break
					if (mz_H %in% envelop_L[[1]]) {
						mz.this.rt <- mz.this.rt[-match(envelop_L[[1]],mz.this.rt)]
						intensity.this.rt <- intensity.this.rt[-match(envelop_L[[2]],intensity.this.rt)]
						if (length(envelop_L[[1]]) <= 9) next
						if (length(envelop_L[[1]])>13) {
							envelop_L[[1]] <- envelop_L[[1]][1:13]
							envelop_L[[2]] <- envelop_L[[2]][1:13]
						}
						mz.this.rt.save <- mz.this.rt.save[-match(envelop_L[[1]],mz.this.rt.save)]
						intensity.this.rt.save <- intensity.this.rt.save[-match(envelop_L[[2]],intensity.this.rt.save)]
						mass_L <- mz_L*z
						elements.count <- averagine.count(mass_L)
						elements.count_H <- c(elements.count["C"]-6,elements.count[c("H","O","N")],6)
						names(elements.count_H)[5] <- "C13"
						normal.prob_L <- isotope.dist(elements.count)[1:15]
						normal.prob_H <- isotope.dist(elements.count_H)[1:15]
						normal.prob <- normal.prob_L + normal.prob_H
						envelop.length.diff <- 15 - length(envelop_L[[2]]) ## length diffence between experimental got envelop and simulated one
						normal.R2 <- rep(0,envelop.length.diff+1)
						for (i in 1:(envelop.length.diff+1))  {
							experimentint <- envelop_L[[2]]
							simulatedint <- normal.prob[i:(length(envelop_L[[2]])-1+i)]
							normal.R2[i] <- acos(sum(experimentint*simulatedint) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint*simulatedint))))/pi*180
						}
						R2.normal.index <- which.min(normal.R2)
						R2.normal <- normal.R2[R2.normal.index]
						output_mz <- paste(envelop_L[[1]],collapse=";")
						output_int <- paste(envelop_L[[2]],collapse=";")
						if (R2.normal < r2.cutoff) {
							temp <- c(mz_L,mz_H,rt,ms1_num,output_mz,output_int,z,R2.normal)
							envelop_pair <- rbind(envelop_pair,temp)
							out.picturename <- paste(mz_L,rt,ms1_num,z,R2.normal,".png",sep = "_")
							out.picturepath <- paste(name,"PNG",sep="_")
							png(paste(out.picturepath,out.picturename,sep ="/"),height = 400, width = 350*3, pointsize = 16)
							layout(matrix(c(1,2),nr=1))
							plot(envelop_L[[1]],envelop_L[[2]],type = "h", col = "black",xlab="m/z",ylab="Intensity")
							#par(new=T)
							plot(normal.prob, type = "h", col = "black", xlab = "Difference to mono mass", ylab = "Probobility" )
							dev.off()
						}
					} else {
						envelop_H <- get_envelope_from_one_peak(pair_mz_vector[2],mz.this.rt,pair_int_vector[2],intensity.this.rt,z)
						if (length(envelop_L[[1]])>6) {
							envelop_L[[1]] <- envelop_L[[1]][1:6]
							envelop_L[[2]] <- envelop_L[[2]][1:6]
						}
						if (length(envelop_H[[1]])>6) {
							envelop_H[[1]] <- envelop_H[[1]][1:6]
							envelop_H[[2]] <- envelop_H[[2]][1:6]
						}
						#if (any(abs(envelop_H[[1]]-a)<0.001)) break
						#if(length(envelop_L[[1]])>4) break
						mz.this.rt <- mz.this.rt[-match(c(envelop_L[[1]],envelop_H[[1]]),mz.this.rt)]
						intensity.this.rt <- intensity.this.rt[-match(c(envelop_L[[2]],envelop_H[[2]]),intensity.this.rt)]
						if (length(envelop_L[[1]]) < 4 | length(envelop_H[[1]]) < 4) next
						mass_L <- mz_L*z
						elements.count <- averagine.count(mass_L)
						elements.count_H <- c(elements.count["C"]-6,elements.count[c("H","O","N")],6)
						names(elements.count_H)[5] <- "C13"
						normal.prob_L <- isotope.dist(elements.count)[1:8]
						normal.prob_H <- isotope.dist(elements.count_H)[7:14]
						normal.prob <- normal.prob_L + isotope.dist(elements.count_H)[1:14]
						envelop.length.diff_L <- 6 - length(envelop_L[[2]]) ## length diffence between experimental got envelop and simulated one
						envelop.length.diff_H <- 6 - length(envelop_H[[2]])
						L.R2 <- rep(0,envelop.length.diff_L+1)
						H.R2 <- rep(0,envelop.length.diff_H+1)
						for (i in 1:(envelop.length.diff_L+1))  {
							experimentint <- envelop_L[[2]]
							simulatedint <- normal.prob_L[i:(length(envelop_L[[2]])-1+i)]
							L.R2[i] <- acos(sum(experimentint*simulatedint) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint*simulatedint))))/pi*180
						}
						for (i in 1:(envelop.length.diff_H+1))  {
							experimentint <- envelop_H[[2]]
							simulatedint <- normal.prob_H[i:(length(envelop_H[[2]])-1+i)]
							H.R2[i] <- acos(sum(experimentint*simulatedint) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint*simulatedint))))/pi*180
						}
						R2.L.index <- which.min(L.R2)
						R2.L <- L.R2[R2.L.index]
						R2.H.index <- which.min(H.R2)
						R2.H <- H.R2[R2.H.index]
						score <- mean(c(R2.L,R2.H))
						output_mz <-paste(c(envelop_L[[1]],envelop_H[[1]]),collapse=";")
						output_int <-paste(c(envelop_L[[2]],envelop_H[[2]]),collapse=";")
						if (R2.L < r2.cutoff & R2.H < r2.cutoff) {
							mz.this.rt.save <- mz.this.rt.save[-match(c(envelop_L[[1]],envelop_H[[1]]),mz.this.rt.save)]
							intensity.this.rt.save <- intensity.this.rt.save[-match(c(envelop_L[[2]],envelop_H[[2]]),intensity.this.rt.save)]
							temp <- c(mz_L,mz_H,rt,ms1_num,output_mz,output_int,z,score)
							envelop_pair <- rbind(envelop_pair,temp)
							out.picturename <- paste(mz_L,rt,ms1_num,z,score,".png",sep = "_")
							out.picturepath <- paste(name,"PNG",sep="_")
							png(paste(out.picturepath,out.picturename,sep ="/"),height = 400, width = 350*3, pointsize = 16)
							layout(matrix(c(1,2),nr=1))
							plot(c(envelop_L[[1]],envelop_H[[1]]),c(envelop_L[[2]],envelop_H[[2]]),type = "h", col = "black",xlab="m/z",ylab="Intensity")
							#par(new=T)
							plot(normal.prob, type = "h", col = "black", xlab = "Difference to mono mass", ylab = "Probobility" )
							dev.off()
						}
					}
				} else {
					mz.this.rt <- mz.this.rt[-match(mz.this.rt[1],mz.this.rt)]
					intensity.this.rt <- intensity.this.rt[2:length(intensity.this.rt)]
				}
			}
				
			}
		}
	return(envelop_pair)
}

get_common_mass_pair <- function(envelop) {
	out_matrix <- data.frame(run=numeric(1),count=numeric(1),rt=numeric(1),rt_mean=numeric(1),rt_variance=numeric(1),ms1_num=numeric(1),mass=numeric(1),mass_mean=numeric(1),mass_variance=numeric(1),R2=numeric(1),
				R2_mean=numeric(1),R2_variance=numeric(1),z=numeric(1))
	envelop <- envelop[order(envelop[,"mz_L"]), ]
	mass_candidate <- as.numeric(envelop[,'mz_L'])*as.numeric(envelop[,'charge'])
	mass_delta <- as.numeric(mass_candidate[2:length(mass_candidate)])-as.numeric(mass_candidate[1:length(mass_candidate)-1])
	mass_delta_ppm <- mz.ppm.cut * ( as.numeric(mass_candidate[2:length(mass_candidate)]) + as.numeric(mass_candidate[1:length(mass_candidate)-1] )) / 2.0
	#check whether the difference is less than max ppm
	#mass_delta_pass <- mass_delta < mass_delta_ppm
	mass_delta_pass <- abs(mass_delta) < mass_delta_ppm
	pass_index <- which(mass_delta_pass == TRUE)
	pass_index_new <- as.list(pass_index)
	begin <- 1
	end <- 1
	# get consecutive index of the pass index list and retrive the corresponding information to the output
	while (TRUE) {
		if (end %% 100 == 0) cat(end)
		if (begin > length(mass_delta_pass)[1]) break 
		if (mass_delta_pass[end] == 1) {
			end = end+1
		}
		##cat(end)
		if (end >= dim(envelop)[1] | mass_delta_pass[end] != 1) {
			candidate_begin <- begin
			candidate_end <- end
			candidate_list <- seq(candidate_begin, candidate_end, 1)
			count <- length(candidate_list)

			if (count > 1) {
				#print(c(begin,end))
				candidate_z <- envelop$charge[candidate_list[1]]
				candidate_mass <- paste(envelop$mz_L[candidate_list],collapse = ';')
				candidate_mass_mean <- mean(as.numeric(envelop$mz_L[candidate_list]))
				candidate_mass_variance <- sd(as.numeric(envelop$mz_L[candidate_list]))
				candidate_rt <- paste(envelop$rt[candidate_list],collapse = ';')
				candidate_rt_mean <- mean(as.numeric(envelop$rt[candidate_list]))
				candidate_rt_variacne <- sd(as.numeric(envelop$rt[candidate_list]))
				candidate_ms1_num <- paste(envelop$ms1_num[candidate_list],collapse = ';')
				candidate_R2 <- paste(envelop$score[candidate_list],collapse = ';')
				candidate_R2_mean <- mean(as.numeric(envelop$score[candidate_list]))
				candidate_R2_variance <- sd(as.numeric(envelop$score[candidate_list]))
				new_list <- list(run,count,candidate_rt,candidate_rt_mean,candidate_rt_variacne,candidate_ms1_num,candidate_mass,candidate_mass_mean,candidate_mass_variance,
								candidate_R2,candidate_R2_mean,candidate_R2_variance,candidate_z)
				out_matrix <- rbind(out_matrix,new_list)
			}
			begin <- end+1
			end <- begin
		}
	}
	return(out_matrix)
}

calc.b.y.ions <- function(sequence, aa.mass.vec) {
	peptide.vec <- unlist(strsplit(sequence,"",fixed=T))
	b_y_matrix <- data.frame(b_index=numeric(1),b_ions=numeric(1),aa=numeric(1),y_ions=numeric(1),y_index=numeric(1))
	mod_pos <- which(peptide.vec == "*") - 1
	peptide.vec <- peptide.vec[-which(peptide.vec=="*")]
	for (i in 1:(length(peptide.vec)-1)) {
		if ( i < mod_pos) {
			b_index <- i
			y_index <- length(peptide.vec) - b_index 
			b_ions <- sum(aa.mass.vec[peptide.vec[1:i]]) + aa.mass.vec["NTERM"]
			y_ions <- sum(aa.mass.vec[peptide.vec[(i+1):length(peptide.vec)]]) +aa.mass.vec["CTERM"] + 2*1.0072765 +aa.mass.vec["*"]
		}
		if ( i >= mod_pos) {
			b_index <- i 
			y_idnex <- length(peptide.vec) - b_index 
			b_ions <- sum(aa.mass.vec[peptide.vec[1:i]]) + aa.mass.vec["NTERM"] + aa.mass.vec["*"]
			y_ions <- sum(aa.mass.vec[peptide.vec[(i+1):length(peptide.vec)]]) +aa.mass.vec["CTERM"] + 2*1.0072765
		}
		temp <- c(b_index,b_ions,peptide.vec[i],y_ions,y_index)
		b_y_matrix <- rbind(b_y_matrix,temp)
	}
	return(b_y_matrix)
}

ms2_comparison <- function(mass_list, ms2) {
	for (i in 1:length(ms2@name)) {
		count <- 0
		candidate_peaks <- c()
		id <- ms2@name[[i]]
		parent_ion <- ms2@parent_ion[[i]]
		rt <- ms2@rt[[i]]
		peaks <- ms2@peaks[[i]]
		for (j in mass_list) {
			if (any(abs(peaks - j) < mz.ppm.cut * j )) {
				pos <- which(abs(peaks - j) < mz.ppm.cut * j )
				count <- count +1
				candidate_peaks <- c(candidate_peaks,peaks[pos] )
			}
		}
		if ( count > length(mass_list)/2 ) {
			candidate_peaks_temp <- paste(candidate_peaks, collapse=";")
			output_temp <- c(id, parent_ion, rt, count, candidate_peaks_temp)
			output_matrix <- rbind(output_matrix,output_temp)
		}
	}
	return(output_matrix)
}

pair_finder_2 <- function (ms1) {
	envelop_pair <- data.frame(mz_L=numeric(1),mz_H=numeric(1),envelop_mz=numeric(1),envelop_int=numeric(1),charge= numeric(1), score=numeric(1))
	for (j in 1:length(ms1@name)) {
		if (j %%100==0) cat(j)
		mz.this.rt.save <- mz.this.rt <- ms1@peaks[[j]]
		intensity.this.rt.save <- intensity.this.rt <- ms1@intensity[[j]]
		rt <- ms1@rt[[j]]
		for (z in charge.states) {
			mz.this.rt <- mz.this.rt.save
			intensity.this.rt <- intensity.this.rt.save
			while(TRUE) {
				print(length(mz.this.rt))
				if (length(mz.this.rt) < 6) break
				if (any(z*(mz.this.rt - mz.this.rt[1])-pair_delta < mz.ppm.cut * mz.this.rt[1]*z & z*(mz.this.rt - mz.this.rt[1])-pair_delta>0)) {
					pos <- which(z*(mz.this.rt - mz.this.rt[1])-pair_delta < mz.ppm.cut * mz.this.rt[1]*z & z*(mz.this.rt - mz.this.rt[1])-pair_delta>0)
					charge <- z
					mz_L <- mz.this.rt[1]
					mz_H <- mz.this.rt[pos]
					mz_L_2_s <- mz_L+isotope.mass.unit/6
					mz_H_2_s <- mz_H+isotope.mass.unit/6
					mz_L_3_s <- mz_L+isotope.mass.unit*2/6
					mz_H_3_s <- mz_H+isotope.mass.unit*2/6
					if(any(mz.this.rt-mz_L_2_s < mz_L_2_s*mz.ppm.cut) & any(mz.this.rt-mz_L_3_s < mz_L_3_s*mz.ppm.cut)
						& any(mz.this.rt-mz_H_2_s < mz_H_2_s*mz.ppm.cut) & any(mz.this.rt-mz_H_3_s < mz_H_3_s*mz.ppm.cut)) {
						mz_L_2 <- mz.this.rt(which(mz.this.rt-mz_L_2_s < mz_L_2_s*mz.ppm.cut))
						mz_H_2 <- mz.this.rt(which(mz.this.rt-mz_L_3_s < mz_L_3_s*mz.ppm.cut))
						mz_L_3 <- mz.this.rt(which(mz.this.rt-mz_H_2_s < mz_H_2_s*mz.ppm.cut))
						mz_H_3 <- mz.this.rt(which(mz.this.rt-mz_H_3_s < mz_H_3_s*mz.ppm.cut))
						int_L_2 <- intensity.this.rt(which(mz.this.rt-mz_L_2_s < mz_L_2_s*mz.ppm.cut))
						int_H_2 <- intensity.this.rt(which(mz.this.rt-mz_L_3_s < mz_L_3_s*mz.ppm.cut))
						int_L_3 <- intensity.this.rt(which(mz.this.rt-mz_H_2_s < mz_H_2_s*mz.ppm.cut))
						int_H_3 <- intensity.this.rt(which(mz.this.rt-mz_H_3_s < mz_H_3_s*mz.ppm.cut))
					}
					
					
				} else {
					mz.this.rt <- mz.this.rt[-match(mz.this.rt[1],mz.this.rt)]
				}
			}
				
			}
		}
}

read.ms2.all <- function(ms2file) {
	setClass("ms2class", representation(name="list",parent_ion="list",rt = "list", precursor_int = "list", ion_injection_time = "list", activation_type = "list", precursor_file = "list", precursor_scan = "list",
				instrument_type = "list", z = "list", massplus = "list", peaks="list",intensity="list",peaks_z = "list"))
	ms2file <- file(ms2file, "r")
	line=readLines(ms2file,n=1)
	ms2list <- new("ms2class", name=list(), parent_ion=list(), rt=list(), precursor_int = list(), ion_injection_time = list(), activation_type = list(), precursor_file = list(), precursor_scan = list(),
				instrument_type = list(), z = list(), massplus = list(), peaks=list(), intensity=list(),peaks_z = list())
	while( length(line) != 0 ) {
		line <- unlist(strsplit(line,"\\s"))
		##if this is a numeric line, conver them to numbers
		if (!is.na(as.numeric(line))) {
			line <- as.numeric(line)
		}
		if (line[1] == "S") {
			if (exists("spectrum_name")) {
				ms2list@name[[length(ms2list@name)+1]] <- spectrum_name
				ms2list@rt[[length(ms2list@rt)+1]] <- spectrum_rt
				ms2list@parent_ion[[length(ms2list@parent_ion)+1]] <- spectrum_parent_ion
				ms2list@precursor_int[[length(ms2list@precursor_int)+1]] <- spectrum_PI
				ms2list@ion_injection_time[[length(ms2list@ion_injection_time)+1]] <- spectrum_IIT
				ms2list@activation_type[[length(ms2list@activation_type)+1]] <- spectrum_AT
				ms2list@precursor_file[[length(ms2list@precursor_file)+1]] <- spectrum_PF
				ms2list@precursor_scan[[length(ms2list@precursor_scan)+1]]<- spectrum_PS
				ms2list@instrument_type[[length(ms2list@instrument_type)+1]] <- spectrum_IT
				ms2list@z[[length(ms2list@z)+1]] <- spectrum_z
				ms2list@massplus[[length(ms2list@massplus)+1]] <- spectrum_massplus
				ms2list@peaks[[length(ms2list@peaks)+1]] <- spectrum_peaks
				ms2list@intensity[[length(ms2list@intensity)+1]] <- spectrum_intensity
				ms2list@peaks_z[[length(ms2list@peaks_z)+1]] <- spectrum_peaks_z
			}
			spectrum_name <- line[2]
			if (as.numeric(spectrum_name) %% 1000==0) cat(spectrum_name)
			spectrum_parent_ion <- line[4]
			spectrum_peaks <- c()
			spectrum_intensity <- c()
			spectrum_peaks_z <- c()
		}
		if (line[2] == "RetTime") {
			spectrum_rt <- line[3]
		}
		if (line[2] == "PrecursorInt") {
			spectrum_PI <- line[3]
		}
		if (line[2] == "IonInjectionTime") {
			spectrum_IIT <- line[3]
		}
		if (line[2] == "ActivationType") {
			spectrum_AT <- line[3]
		}
		if (line[2] == "PrecursorFile") {
			spectrum_PF <- line[3]
		}
		if (line[2] == "PrecursorScan") {
			spectrum_PS <- line[3]
		}
		if (line[2] == "InstrumentType") {
			spectrum_IT <- line[3]
		}
		if (line[1] == "Z") {
			spectrum_z <- line[2]
			spectrum_massplus <- line[3]
		}
		if (is.numeric(line[1])) {
			spectrum_peaks <- append(spectrum_peaks, line[1])
			spectrum_intensity <- append(spectrum_intensity, line[2])
			spectrum_peaks_z <- append(spectrum_peaks_z, line[3])
		} 
		line=readLines(ms2file,n=1)
		if (length(line) == 0) {
			#break
			ms2list@name[[length(ms2list@name)+1]] <- spectrum_name
			ms2list@rt[[length(ms2list@rt)+1]] <- spectrum_rt
			ms2list@parent_ion[[length(ms2list@parent_ion)+1]] <- spectrum_parent_ion
			ms2list@precursor_int[[length(ms2list@precursor_int)+1]] <- spectrum_PI
			ms2list@ion_injection_time[[length(ms2list@ion_injection_time)+1]] <- spectrum_IIT
			ms2list@activation_type[[length(ms2list@activation_type)+1]] <- spectrum_AT
			ms2list@precursor_file[[length(ms2list@precursor_file)+1]] <- spectrum_PF
			ms2list@precursor_scan[[length(ms2list@precursor_scan)+1]]<- spectrum_PS
			ms2list@instrument_type[[length(ms2list@instrument_type)+1]] <- spectrum_IT
			ms2list@z[[length(ms2list@z)+1]] <- spectrum_z
			ms2list@massplus[[length(ms2list@massplus)+1]] <- spectrum_massplus
			ms2list@peaks[[length(ms2list@peaks)+1]] <- spectrum_peaks
			ms2list@intensity[[length(ms2list@intensity)+1]] <- spectrum_intensity
			ms2list@peaks_z[[length(ms2list@peaks_z)+1]] <- spectrum_peaks_z
		}
	}
	close(ms2file)
	return(ms2list)
}

peak_extraction_judgement_ms2 <- function(ms2,index) {
		count <- 0
		mz.this.rt.save <- mz.this.rt <- ms2@peaks[[index]]
		intensity.this.rt.save <- intensity.this.rt <- ms2@intensity[[index]]
		max_int_this_scan <- max(intensity.this.rt)
		charge.states <- seq((as.numeric(ms2@z[index])-1),1)
		rt <- ms2@rt[[index]]
		for (i in charge.states) {
			mz.this.rt <- mz.this.rt.save
			intensity.this.rt <- intensity.this.rt.save
			mz.this.rt.charge <- mz.this.rt * i
			#at the begining, the candidate envelope only contains the first one of the mz list, so the envelop.end is 1
			#,and test wether the second is belong to this envelope, so the test.end is 2
			envelop.end <- 1
			envelop.mass <- c(mz.this.rt.charge[1])
			#set the candidate envelope member as NA, in order to delete it from the list after extracting the whole envelope
			#so that the computer could save its source and time.
			mz.this.rt.charge[1] <- NA
			mz.this.rt[1] <- NA
			envelop.intensity <- intensity.this.rt[1]
			test.end <- 2
			while(TRUE) {
				# if the tested mz index is out of the mz list or the testing mz index is too far away from the envelope end, finish doing this envelope extracting.
				if (test.end > length(mz.this.rt) | test.end - envelop.end > 5)  {
					#if (envelop.mass[1] == a) break
					max_int_this_envelop <- max(envelop.intensity)
					#test whether it is a valid envelope
					if (length(envelop.mass) > 2 & max_int_this_envelop/max_int_this_scan > 0.1 & length(envelop.mass) < 15) {
						if (length(envelop.mass) > 3) {
							r2.cutoff <- 20
							ratio.Se.normal.cutoff <- 1/2.5
						} else {
							r2.cutoff <- 15
							ratio.Se.normal.cutoff <- 1/4
						}
						##use match here to get the index of the envelop.mz in the mz.this.rt.save
						## if it is a real envleope, delete them from the saved list
						intensity.this.rt.save <- intensity.this.rt.save[-match(envelop.mass[which(envelop.mass!=0)],mz.this.rt.save * i)]
						mz.this.rt.save <- mz.this.rt.save[-match(envelop.mass[which(envelop.mass!=0)],mz.this.rt.save * i)]
						envelop.mz <- envelop.mass / i
						mass_mono <- envelop.mass[which.max(envelop.intensity)] - i*Hplus
						envelop.temp.mass <- paste(envelop.mass, collapse = ";")
						envelop.temp.mz <- paste(envelop.mz, collapse = ";")
						envelop.temp.intensity <- paste(envelop.intensity, collapse = ";")
						charge <- i
						results <- envelope.rating(envelop.mz,charge,envelop.intensity)
						#print(results)
						if (results[[1]] < r2.cutoff & results[[2]] < ratio.Se.normal.cutoff) {
							#break
							#print(envelop.mass)
							#print(envelop.intensity)
							#print(c(results[[1]],results[[2]]))
							count <- count + 1
						}
					}
					#delete the mz index whether it is a real or a pseudo envelope 
					intensity.this.rt <- intensity.this.rt[-which(is.na(mz.this.rt))]
					mz.this.rt <- mz.this.rt[-which(is.na(mz.this.rt))]
					mz.this.rt.charge <- mz.this.rt.charge[-which(is.na(mz.this.rt.charge))]
					## set the candidate envelope as the original status
					envelop.end <- 1
					envelop.mass <- c(mz.this.rt.charge[1])
					mz.this.rt.charge[1] <- NA
					mz.this.rt[1] <- NA
					envelop.intensity <- intensity.this.rt[1]
					test.end <- 2
				}
				if (length(mz.this.rt) <= 5) break
				if (abs(abs(mz.this.rt.charge[test.end]-envelop.mass[length(envelop.mass)])-isotope.mass.unit) < 
					mz.ppm.cut * (mz.this.rt.charge[test.end] + envelop.mass[length(envelop.mass)]) / 2.0 ) {
					# if the differrence between the testing mass and the envelope end conforms with the iosotope mass unit,add it the the envelope and replace it with NA
					envelop.mass <- c(envelop.mass,mz.this.rt.charge[test.end])
					#print(envelop.mass)
					envelop.intensity <- c(envelop.intensity, intensity.this.rt[test.end])
					mz.this.rt.charge[test.end] <- NA
					mz.this.rt[test.end] <- NA
					envelop.end <- test.end
					test.end <- test.end + 1
					next
				}
				test.end <- test.end + 1
			}
		}
	return(count)
}

envelope.rating.ms2 <- function(envelop.mz, envelop.charge, envelop.intensity) {
	if (any(envelop.mz == 0)) {
			zero.pos <- which(envelop.mz==0)
			for (i in zero.pos) {
			envelop.mz[i] <- (envelop.mz[i+1] + envelop.mz[i-1])/2
			}
		}
	envelop.mass <- max(envelop.mz) * envelop.charge
	max.index <- which.max(envelop.intensity)
	if (e=="Se") {
		elements.count <- averagine.count(envelop.mass)
		elements.count_Se <- averagine.count(envelop.mass-79.9165)
		#Se.prob <- c(0.0089,0.0000,0.0937,0.0763,0.2377,0.0000,0.4961,0.0000,0.0873)
		Se.prob <- c(0.0089,0.0000,0.0937,0.0763,0.2377,0.0000,0.4961,0.0000,0.0873)
		normal.prob <- isotope.dist(elements.count)[1:6]
		withSe.prob <- round(convolve(isotope.dist(elements.count_Se), rev(Se.prob),type="o"),4)[4:10]
		envelop.length.diff <- 15 - length(envelop.intensity) ## length diffence between experimental got envelop and simulated one
		Se.R2 <- normal.R2 <- rep(0,envelop.length.diff+1)
		for (i in 1:(envelop.length.diff+1))  {
			experimentint <- envelop.intensity
			simulatedint.Se <- withSe.prob[i:(length(envelop.intensity)-1+i)]
			simulatedint.Se[which(envelop.intensity == 0)] <- 0
			Se.R2[i] <- acos(sum(experimentint*simulatedint.Se) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint.Se*simulatedint.Se))))/pi*180
			simulatedint.normal <- normal.prob[i:(length(envelop.intensity)-1+i)]
			simulatedint.normal[which(envelop.intensity == 0)] <- 0
			normal.R2[i] <- acos(sum(experimentint*simulatedint.normal) / (sqrt(sum(experimentint*experimentint))*sqrt(sum(simulatedint.normal*simulatedint.normal))))/pi*180
		}
		R2.Se.index <- which.min(Se.R2)
		R2.Se <- Se.R2[R2.Se.index]
		R2.normal.index <- which.min(normal.R2)
		R2.normal <- normal.R2[R2.normal.index]
		ratio.Se.normal <- R2.Se/R2.normal
		results <- list(R2.Se,ratio.Se.normal,withSe.prob,normal.prob,max.index)
	}
	return(results)
}

Runtime <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
