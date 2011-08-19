# This function reads a set of MultiDimensionalHistograms written to a file
# using the << operator. It returns a list of matricies.
#
# intervalsigfigs: the number of significant figures used for interval labels
# intervalscale: the scaling factor applied to the histogram range (e.g. 180/pi)
read.mdhists <- function(filename, intervalsigfigs=3, intervalscale=1) {

	filelines <- readLines(filename)
	filelines <- filelines[nchar(filelines) != 0]

	histlist <- list()
	num <- 1

	while (num < length(filelines)) {

		linesplit <- strsplit(filelines[num], " ")[[1]]
		num <- num+1
		numdims <- as.integer(linesplit[1])
		label <- paste(linesplit[-1], collapse=" ")

		dnames <- list()
		dranges <- list()
		for (i in seq_len(numdims)) {
			linesplit <- strsplit(filelines[num], " ")[[1]]
			num <- num+1
			numbins <- as.integer(linesplit[1])
			start <- as.numeric(linesplit[2])*intervalscale
			end <- as.numeric(linesplit[3])*intervalscale
			breaks <- seq(start, end, length.out=numbins+1)
			if (length(breaks) > 2) {
				breaktext <- paste("[",
				                   signif(breaks[1:(length(breaks)-2)], intervalsigfigs),
				                   ",",
				                   signif(breaks[2:(length(breaks)-1)], intervalsigfigs),
				                   ")", sep="")
			} else {
				breaktext <- character()
			}
			breaktext <- c(breaktext, paste("[",
			                                signif(breaks[length(breaks)-1], intervalsigfigs),
			                                ",",
			                                signif(breaks[length(breaks)], intervalsigfigs),"]",
			                                sep=""))

			newdnames <- list(c("underflow", breaktext, "overflow"))
			names(newdnames) <- paste(linesplit[-(1:3)], collapse=" ")
			dnames <- c(dnames, newdnames)

			newdranges <- list(c(start, end))
			names(newdranges) <- paste(linesplit[-(1:3)], collapse=" ")
			dranges <- c(dranges, newdranges)
		}

		newhist <- list(array(as.integer(strsplit(filelines[num], " ")[[1]]), sapply(dnames, length), dnames))
		num <- num+1
		if (nchar(label)) names(newhist) <- label
		attr(newhist[[1]], "dimranges") <- dranges

		histlist <- c(histlist, newhist)
	}

	histlist
}
