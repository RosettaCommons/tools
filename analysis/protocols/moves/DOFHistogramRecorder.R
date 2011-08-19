source(file.path(rosetta_analysis_dir, "numeric", "MultiDimensionalHistogram.R"))

plot_dof_histograms <- function(filename) {

	mdhists <- read.mdhists(filename, 8, 180/pi)
	
	filename_noext <- gsub("\\..+$", "", filename)
	
	pdf(paste(filename_noext, ".pdf", sep=""))
	on.exit(dev.off())
	
	for (histnum in seq_along(mdhists)) {
	
		mdhist <- mdhists[[histnum]]
	
		for (i in seq_along(dimnames(mdhist))) {
		
			dimname <- names(dimnames(mdhist))[[i]]
			
			counts <- apply(mdhist, i, sum)
			if (sum(counts[c(1, length(counts))])) warning("Underflow/Overflow Detected")
			counts <- counts[-c(1, length(counts))]
			
			breaks <- as.numeric(unique(strsplit(paste(gsub("[^-0-9.,]", "", names(counts)), collapse=","), ",")[[1]]))
			
			counts <- unname(counts)
			
			dof_type <- as.integer(gsub(".*type= ([0-9]+).*", "\\1", dimname))
			
			if (dof_type == 1) {
			
				expected_counts <- sum(counts)*uniform_dof_distribution(dof_type, breaks)
			
			} else if (dof_type == 2) {
				
				counts <- rev(counts)
				breaks <- rev(180-breaks)
				
				expected_counts <- sum(counts)*uniform_dof_distribution(dof_type, breaks)
			
			} else {
				
				break
			}
			
			mse <- mean((counts/sum(counts) - expected_counts/sum(expected_counts))^2)
			
			plot(0, 0, type="n", xlim=range(breaks), ylim=range(0, counts, expected_counts),
			     main=paste(names(mdhists)[histnum], ": ", dimname, " (MSE: ", signif(mse, 4), ")", sep=""), xlab="Degrees", ylab="Counts")
			
			rect(breaks[-length(breaks)], numeric(length(counts)), breaks[-1], counts, col="gray", border="black")
			
			points(breaks[-1] - diff(breaks)/2, expected_counts, type="p", col="blue")
			points(breaks[-1] - diff(breaks)/2, expected_counts, type="l", col="blue")
			
			cat(dimname, " : ", signif(mse), "\n", sep="")
		}
	}
}

uniform_dof_distribution <- function(dof_type, breaks) {

	if (dof_type == 1) {
	
		frequencies <- rep(1, length(breaks)-1)
	}
	
	if (dof_type == 2) {
	
		frequencies <- cos(pi/180*breaks[-1]) - cos(pi/180*breaks[-length(breaks)])
	}
	
	frequencies/sum(frequencies)
}
