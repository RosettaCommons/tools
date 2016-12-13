# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# Define vectors for mapping between 1 and 3 character residue names

aa3 <- c("TRP", "PHE", "TYR", "MET", "LEU", "ILE", "VAL", "ALA", "GLY", "SER", 
         "THR", "ARG", "LYS", "HIS", "ASN", "GLN", "ASP", "GLU", "PRO", "CYS")

aa1 <- c("W", "F", "Y", "M", "L", "I", "V", "A", "G", "S", 
         "T", "R", "K", "H", "N", "Q", "D", "E", "P", "C")

names(aa3) <- aa1
names(aa1) <- aa3

# This function reads *.ga.entities checkpoint files. It assumes that all 
# entities have the same composition. It reads checkpointed Entity and 
# MultiStateEntity objects. It returns a data.frame object with a row for each 
# entity. The traits at each sequence position are given first, then the 
# overall fitness, then the state fitnesses and metric values.

read_ga_entities <- function(filename, restypes = NULL) {

	metric_types <- list(Real = numeric(), Int = integer(), Size = integer(), Bool = logical())
	
	file_con <- if (length(grep("\\.gz$", filename))) gzfile(filename) else file(filename)
	tokens <- scan(file_con, character(), 300, quiet = TRUE)
	close(file_con)
	
	oldformat <- length(grep("AA:", tokens[2])) == 0
	
	if (tokens[1] != "traits") stop(paste(filename, "does not appear to be an entities file"))
	
	scan_what <- list(NULL)
	
	# set up the traits input
	num_traits <- grep("fitness", tokens)[1]-2
	if (oldformat) {
		res_nums <- lapply(strsplit(tokens[seq_len(num_traits)+1], "\\."), "[[", 1)
	} else {
		res_nums <- lapply(strsplit(tokens[seq_len(num_traits)+1], ":"), "[[", 2)
	}
	traits_what <- rep(list(character()), num_traits)
	names(traits_what) <- paste("AA", res_nums, sep = "")
	
	scan_what <- c(scan_what, traits_what, list(NULL, fitness = numeric()))
	
	# handle additional MultiStateEntity data
	if (tokens[num_traits+4] == "states") {
	
		scan_what <- c(scan_what, list(NULL, NULL))
		
		# iterate over the number of states
		num_states <- as.integer(tokens[num_traits+5])
		token_offset <- num_traits+6
		for (i in seq_len(num_states)) {
			state_what <- list(NULL, numeric(), NULL, NULL)
			names(state_what) <- c("", paste("state", i, "_fitness", sep = ""), "", "")
			scan_what <- c(scan_what, state_what)
			
			# iterate over the number of metrics
			num_metrics <- as.integer(tokens[token_offset+3])
			token_offset <- token_offset+4
			for (j in seq_len(num_metrics)) {
			
				metric_name <- paste("state", i, "_", tokens[token_offset], sep = "")
				metric_type <- tokens[token_offset+1]
				
				scan_what <- c(scan_what, list(NULL, NULL))
				token_offset <- token_offset + 2
				
				if (length(grep("\\[$", metric_type))) {
					# handle vector metrics
					metric_length <- which(tokens[-seq_len(token_offset-1)] == "]")[1] - 1
					metric_type <- substr(metric_type, 1, nchar(metric_type)-1)
					metric_what <- rep(list(metric_types[[metric_type]]), metric_length)
					names(metric_what) <- paste(metric_name, seq_len(metric_length), sep = "")
					scan_what <- c(scan_what, metric_what, list(NULL))
					token_offset <- token_offset + metric_length + 1
				} else {
					# handle scalar metrics
					metric_what <- list(metric_types[[metric_type]])
					names(metric_what) <- metric_name
					scan_what <- c(scan_what, metric_what)
					token_offset <- token_offset + 1
				}
			}
		}
	}
	
	file_con <- if (length(grep("\\.gz$", filename))) gzfile(filename) else file(filename)
	result <- scan(file_con, scan_what, quiet = TRUE)[names(scan_what) != ""]
	close(file_con)
	
	if (is.null(restypes)) {
		for (i in seq_len(num_traits)) {
			if (oldformat) {
				result[[i]] <- unname(aa1[sub(".+\\.", "", result[[i]])])
			} else {
				result[[i]] <- sub(".+:", "", result[[i]])
			}
		}
	} else {
		for (i in seq_len(num_traits)) {
			if (oldformat) {
				result[[i]] <- factor(unname(aa1[sub(".+\\.", "", result[[i]])]), restypes)
			} else {
				result[[i]] <- factor(sub(".+:", "", result[[i]]), restypes)
			}
		}
	}
	
	as.data.frame(result)
}

# Read a list of *.ga.entities checkpoint files out of a directory. By default,
# the parsed data is saved in the R format

read_ga_entities_list <- function(dirpath, filepattern = NULL, recompute = FALSE, savedata = FALSE, readgen = FALSE) {

	filename <- file.path(dirpath, paste("entities", filepattern, ".Rda", sep = ""))
	
	if (file.exists(filename) && !recompute) {
		load(filename)
	} else {
		simpattern <- if (is.null(filepattern)) "" else filepattern
		
		entitiesfiles <- list.files(dirpath, pattern = paste(simpattern, ".*\\.ga\\.entities", sep = ""),
		                            full.names = TRUE, recursive = TRUE)
		entitiesfiles <- entitiesfiles[file.info(entitiesfiles)$size > 0]
		
		entitieslist <- vector("list", length(entitiesfiles))
		generationslist <- vector("list", length(entitiesfiles))
		
		for (i in seq(along = entitiesfiles)) {
			print(entitiesfiles[i])
			entitieslist[[i]] <- read_ga_entities(entitiesfiles[i], unname(aa1))
			if (readgen) generationslist[[i]] <- read_ga_generations(sub("entities", "generations", entitiesfiles[i]), entitieslist[[i]])
		}
		names(entitieslist) <- gsub(".+/", "", entitiesfiles)
		
		if (readgen) attr(entitieslist, "generations") <- generationslist
		
		if (savedata == TRUE) save(entitieslist, file = filename)
	}
	
	entitieslist
}

# This function reads *.ga.generations checkpoint files. It requires that the 
# output of the read_ga_entities() function also be provided. It returns a list 
# of integer vectors with the indices to the entites in each generation.

read_ga_generations <- function(filename, entities) {

	file_con <- if (length(grep("\\.gz$", filename))) gzfile(filename) else file(filename)
	gen_lines <- readLines(file_con)
	close(file_con)
	
	oldformat <- length(grep("AA:", gen_lines[2])) == 0
	
	if (oldformat) {
		genenerations_traits <- gsub("[^ ]+\\.", "", gen_lines)
		replacefunc <- function(x) {
			paste(if (x[1] == "generation") x else unname(aa1[x]), collapse=" ")
		}
		genenerations_traits <- sapply(strsplit(genenerations_traits, " "), replacefunc)
	} else {
		genenerations_traits <- gsub("[^ ]+:", "", gen_lines)
	}
	entities_traits <- do.call("paste", entities[,grep("^AA", names(entities))])
	
	trait_matches <- match(genenerations_traits, entities_traits)
	
	gen_line_indices <- grep("^generation ", gen_lines)
	gen_line_numbers <- as.integer(sub("^generation ", "", gen_lines[gen_line_indices]))
	gen_lengths <- diff(c(gen_line_indices, length(gen_lines)+1)) - 1
	
	result <- vector("list", length(max(gen_line_numbers)))
	
	for (i in seq_along(gen_line_indices)) {
		result[[i]] <- trait_matches[seq_len(gen_lengths[i])+gen_line_indices[i]]
	}
	
	result
}

# Thes function writes *.ga.generations checkpoint files. It takes a filename
# and a data.frame, matrix, or list thereof.

write_ga_generations <- function(filename, traits, oldformat=FALSE) {

	file_con <- if (length(grep("\\.gz$", filename))) gzfile(filename, "w") else file(filename, "w")
	on.exit(close(file_con))

	if (is.data.frame(traits) || is.matrix(traits)) {
		traits <- list(traits)
	}

	for (i in seq_along(traits)) {
	
		poscols <- grep("^AA[0-9]+$", colnames(traits[[i]]))
		posnums <- as.integer(sub("AA", "", colnames(traits[[i]])[poscols]))
		posmat <- as.matrix(traits[[i]][,poscols])
		for (j in seq_len(ncol(posmat))) {
			if (oldformat) {
				posmat[,j] <- paste(posnums[j], aa3[posmat[,j]], sep = ".")
			} else {
				posmat[,j] <- paste("AA", posnums[j], posmat[,j], sep = ":")
			}
		}
		
		cat("generation ", i, "\n", sep="", file=file_con)
		cat(apply(posmat, 1, paste, collapse=" "), sep = "\n", file=file_con)
	}
}

# This function returns the fitness of a given set of entities

entities_fitness <- function(entities, fitness_coef = NULL) {

	if (is.null(fitness_coef)) {
		fitness_coef <- c(fitness = 1)
	}
	
	fitness_matrix <- as.matrix(entities[,names(fitness_coef),drop=FALSE])
	
	fitness_matrix %*% fitness_coef
}

# This function takes an entities data frame as read by read_ga_entities and 
# determines a position weight matrix for the sampled sequence positions. It 
# uses either a fitness cutoff above the sequence with the best fitness, or
# Boltzmann weighting of the those energies. The total fitness is calculated
# by weighting numeric data read along with the sequences. A list of data 
# frames can also be provided, in which case the PWM will be based on merging
# the sequence from all data frames. The minimum score from each individual
# data fram will be used as the reference.
# WARNING: This function assumes the levels of all factors are identical!

entities_pwm <- function(entities, temp_or_thresh, fitness_coef = NULL, 
                         type = c("boltzmann", "cutoff")) {

	type <- match.arg(type)
	
	entities_list <- is.list(entities) && ! is.data.frame(entities)
	
	if (entities_list) {
		nrows <- sapply(entities, nrow)
		offsets <- c(0,cumsum(nrows))
		entities <- do.call(rbind, entities)
	}
	
	fitness <- entities_fitness(entities, fitness_coef)
	
	if (entities_list) {
		min_fitness <- numeric(length(fitness))
		for (i in seq_along(nrows)) min_fitness[seq_len(nrows[i])+offsets[i]] <- min(fitness[seq_len(nrows[i])+offsets[i]])
	} else {
		min_fitness <- min(fitness)
	}
	
	if (type == "cutoff") {
		weight <- fitness <= min_fitness+temp_or_thresh
	} else {
		if (temp_or_thresh != 0) {
			weight <- exp(-(fitness-min_fitness)/temp_or_thresh)
		} else {
			weight <- fitness == min_fitness
		}
	}
	
	pos_columns <- grep("^AA", colnames(entities))
	freqmat <- matrix(nrow = length(levels(entities[,1])), ncol = length(pos_columns))
	
	weight_sum <- sum(weight)
	for (i in seq_along(pos_columns)) {
		freqmat[,i] <- tapply(weight, entities[,pos_columns[i]], sum)/weight_sum
		freqmat[is.na(freqmat[,i]),i] <- 0
	}
	rownames(freqmat) <- levels(entities[,1])
	
	#print(freqmat)
	
	freqmat
}

# This function takes a list of entities data frames and returns a list of
# position weight matrices, where each matrix correspons to a single sequence
# position. The matrices will have one column for every input data frame.
# If combine is used, the PWMs will be combined by weighting all sequences
# together.
# WARNING: This function assumes the levels of all factors are identical!

entities_pwms <- function(entitieslist, temp_or_thresh, fitness_coef = NULL, 
                          type = c("boltzmann", "cutoff"), combine=FALSE) {

	type <- match.arg(type)
	naa <- length(levels(entitieslist[[1]][,1]))
	
	pwmlist <- rep(list(matrix(nrow = naa, ncol = 0)), length(grep("^AA", colnames(entitieslist[[1]]))))
	
	if (combine) {
	
		freqmat <- entities_pwm(entitieslist, temp_or_thresh, fitness_coef, type)
		
		for (j in seq(along = pwmlist)) {
			
			pwmlist[[j]] <- cbind(pwmlist[[j]], freqmat[,j])
		}
	
	} else {
		
		for (i in seq(along = entitieslist) ){
		
			freqmat <- entities_pwm(entitieslist[[i]], temp_or_thresh, fitness_coef, type)
			
			for (j in seq(along = pwmlist)) {
			
				pwmlist[[j]] <- cbind(pwmlist[[j]], freqmat[,j])
			}
		}
	}
	
	#print(pwmlist)
	
	pwmlist
}

# This function takes a list of entities data frames and returns an array of
# position weight matrices. The first dimension has the amino acids types, 
# the second dimension has the replicate, and the third dimension has the
# sequence position

entities_list_pwms <- function(entities, fitness_coef = c(1/2.5, 1/2.5, 1/2.5, 1),
                               temp_or_thresh = 0.228, 
                               type = c("boltzmann", "cutoff")) {

	type <- match.arg(type)
	if (is.null(names(fitness_coef))) {
		names(fitness_coef) <- paste("state1_fitness_comp", seq_along(fitness_coef), sep="")
	}

	pwms <- entities_pwms(entities, temp_or_thresh, fitness_coef, type)
	
	posnames <- colnames(entities[[1]])[seq_along(pwms)]
	posnames <- gsub("AA", "", posnames)
	
	pwmsdimnames <- list(aa=rownames(pwms[[1]]), rep=NULL, pos=posnames)
	
	pwms <- array(do.call("c", pwms), dim = c(dim(pwms[[1]]), length(pwms)))
	dimnames(pwms) <- pwmsdimnames
	
	pwms
}

# This function takes an array of pwms as returned by entities_pwms_array and
# collapses the arry into a single PWM, using a provided percentile cutoff.

collapse_pwms <- function(pwms, percentile = .5) {

	pwm <- apply(pwms, c(1,3), quantile, percentile)
	
	minnotzero <- function(x) {
		x <- x[x!=0]
		if (length(x)) return(min(x))
		NA
	}
	plastmin <- apply(pwms, c(1,3), minnotzero)
	correcteddist <- apply(plastmin, 2, function(x) as.numeric(!is.na(x) & x==min(x, na.rm = TRUE)))
	
	for (i in which(colSums(pwm) == 0)) {
		#print(paste("Correcting", i))
		pwm[,i] <- correcteddist[,i]
	}
	
	pwm <- apply(pwm, 2, function(x) x/sum(x))
	
	pwm
}

# This function extracts the sequence from a PDB file. The residue IDs
# (<chainID><resSeq><iCode>) are given in as the names.

pdb_sequence <- function(pdbpath) {

	if (length(grep(".gz$", pdbpath))) {
		# if a gzip was passed in, unzip and read the lines into an array
		pdbcon <- gzfile(pdbpath)
		pdblines <- readLines(pdbcon)
		close(pdbcon)
	} else {
		# otherwise, just read the lines into an array directly
		pdblines <- readLines(pdbpath)
	}
	
	# Create an array of the lines in the file starting with ATOM
	atomlines <- grep("^ATOM", pdblines, value=TRUE)
	
	# Create arrays of the residue names e.g. GLU, and IDs e.g. 'A 318 '
	resName <- substr(atomlines, 18, 20)
	resID <- substr(atomlines, 22, 27)
	
	resID <- gsub("^ ", "_", resID)
	resID <- gsub(" ", "", resID)
	
	# Assign the corresponding residue ID as a name to each resName
	names(resName) <- resID
	
	# Get a boolean array determining which lines are duplicates
	# Pointwise negate this array to get an array of unique lines
	# Then return an array of residue names whose residue IDs are unique
	resName[!duplicated(resID)]
}

# The function converts a position weight matrix to a matrix of sequences with 
# the same approximate distribution as the original PWM. 

pwm_to_seqmat <- function(pwm, numseq=100) {

	seqmat <- matrix(character(), nrow=numseq, ncol=ncol(pwm))
	for (i in seq_len(ncol(pwm))) {
		colfun <- stepfun(c(0,cumsum(pwm[,i])), c(1,seq_along(pwm[,i]),length(pwm[,i])))
		funx <- seq(0, 1, length.out=numseq+1)
		funx <- funx[-1] - mean(diff(funx))/2
		seqmat[,i] <- names(pwm[,i])[colfun(funx)]
	}
	
	colnames(seqmat) <- colnames(pwm)

	seqmat
}

# This function plots a character matrix, with optional foreground and
# background color matrices.

plot_charmat <- function(charmat, col = NULL, bg = NULL, cex=1, xlim=par("usr")[1:2], ylim=par("usr")[3:4]) {
	
	xlines <- seq(xlim[1], xlim[2], length.out=ncol(charmat)+1)
	xleft <- rep(xlines[-length(xlines)], each=nrow(charmat))
	xright <- rep(xlines[-1], each=nrow(charmat))
	
	ylines <- seq(ylim[1], ylim[2], length.out=nrow(charmat)+1)
	ybottom <- rep(ylines[-length(ylines)], nrow(charmat))
	ytop <- rep(ylines[-1], nrow(charmat))
	
	xcenter <- (xleft+xright)/2
	ycenter <- (ybottom+ytop)/2
	
	if (!is.null(bg)) {
		rect(xleft, ybottom, xright, ytop, col=bg, border=NA)
	}
	
	text(xcenter, ycenter, charmat, col = col, cex = cex)
}

# This function plots a ranked table of amino acid types for each position
# It takes a PWM, an optional experimental PWM, and an optional wild type
# sequence.

plot_seqrank <- function(freq_mat, exp_freq_mat = NULL, wt_seq = NULL, star_mat = NULL, rank_line = 0, wt_col = "red", other_col = "black") {

	# Create three matrices of the dimensions of the pwm
	# char_mat is an matrix of residue names (as 1-character codes), ordered by increasing rank of the residue in the pwm
	# bg_freq_mat is a corresponding matrix of rank values
	# col_mat is a matrix of color names, defaulting to "black"
	char_mat <- matrix(nrow=nrow(freq_mat), ncol=ncol(freq_mat))
	col_mat <- matrix(other_col, nrow=nrow(freq_mat), ncol=ncol(freq_mat))
	bg_freq_mat <- matrix(nrow=nrow(freq_mat), ncol=ncol(freq_mat))
	
	# Loop over all columns (residue positions)
	for (i in seq_len(ncol(freq_mat))) {

		char_mat[,i] <- rownames(freq_mat)[order(freq_mat[,i])]
		if (!is.null(star_mat)) {
			star_mat[,i] <- rev(star_mat[char_mat[,i],i])
		}
		# Loop over all amino acids
		for (j in seq_len(nrow(freq_mat))) {
		
			if (is.null(exp_freq_mat)) {
				bg_freq_mat[j,i] <- freq_mat[char_mat[j,i],i]
			} else {
				bg_freq_mat[j,i] <- exp_freq_mat[char_mat[j,i],i]
			}
			
			if (!is.null(wt_seq)) {
				if (char_mat[j,i] == wt_seq[i]) col_mat[j,i] <- wt_col
			}
		}
	}
	
	# Color shading
	col_levels <- seq(0,1,by=.1)
	col_levels <- seq(0,ceiling(max(bg_freq_mat)/.1)*.1,by=.1)
	#cols <- gray(seq(1,0,length.out=length(col_levels)-1))
	cols <- rev(c(topo.colors(length(col_levels)-2), "white"))
	
	bg_mat <- matrix(cols[pmin(floor((bg_freq_mat)*length(cols)/max(col_levels))+1, length(cols))], nrow=nrow(bg_freq_mat))
	
	op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    mar1 <- c(0.6, 2.7, 2.7, 0.4)
    mar2 <- c(0.6, 0.4, 2.7, 3.4)
    
    devwidth <- par("din")[1]*2.54
    charheight <- par("cin")[2]*2.54
    width1 <- (mar1[2]+mar1[4])*charheight
    width2 <- (mar2[2]+mar2[4])*charheight
    boxwidth <- (devwidth - sum(width1+width2))/(ncol(freq_mat)+1)
    layout(matrix(1:2, nrow=1,ncol=2), widths=c(width1+boxwidth*ncol(freq_mat),width2+boxwidth))
	
	par(mar=mar1, mgp=c(1.5, .25, 0), cex=1)
	
	plot(0, 0, type="n", xlim=c(0.5,0.5+ncol(freq_mat)), ylim=c(20.5,0.5), xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab="")
	plot_charmat(char_mat, col_mat, bg_mat)
	mtext("Predicted Rank", 2, par("mgp")[1])
	axis(2, 1:20, tick=FALSE, las=2)
	mtext("Residue", 3, par("mgp")[1])
	residx <- seq(1, ncol(freq_mat), by=2)
	axis(3, residx, colnames(freq_mat)[residx], tick=FALSE)
	if (ncol(freq_mat) >= 2) {
		residx <- seq(2, ncol(freq_mat), by=2)
		axis(3, residx, colnames(freq_mat)[residx], tick=FALSE)
	}
	box(lwd=.5)
    
    if (!is.null(star_mat)) {
    	points(t(t(which(t(star_mat), arr.ind=TRUE))+c(.3,0)), pch="*")
    }
    
    if (rank_line) {
    	abline(h=rank_line+0.5, lty="dashed")
    }
	
	maradj <- (1-length(cols)/nrow(col_mat))*0.5*par("pin")[2]/par("cin")[2]
	mar2[1] <- mar2[1]+maradj
	mar2[3] <- mar2[3]+maradj
	
	par(mar=mar2, mgp=c(2.2, .25, 0), cex=1)
	
	plot.new()
    plot.window(xlim = c(0, 1), ylim = range(col_levels), xaxs = "i", yaxs = "i")
    rect(0, col_levels[-length(col_levels)], 1, col_levels[-1L], col = cols, lwd=.5)
    axis(4, col_levels[seq(1,length(col_levels),1)], paste(round(col_levels[seq(1,length(col_levels),1)]*100), "%", sep=""), tick=FALSE, las=2)
    bg_title <- "Predicted Frequency"
    if (!is.null(exp_freq_mat)) bg_title <- "Experimental Frequency"
    mtext(bg_title, 4, par("mgp")[1])
    box(lwd=.5)
    
    invisible(bg_freq_mat)
}

# This function produces boxplots showing the amount each generation contributes to
# the the position weight matrix.

plot_gen_contrib <- function(entitieslist, generationslist, 
                             fitness_coef = c(1/2.5, 1/2.5, 1/2.5, 1),
                             temp_or_thresh = 0.228, 
                             type = c("boltzmann", "cutoff"),
                             main = "") {

	type <- match.arg(type)

	names(fitness_coef) <- paste("state1_fitness_comp", seq_along(fitness_coef), sep="")

	gen_contrib <- matrix(nrow=length(entitieslist), ncol=length(generationslist[[1]]))

	for (i in seq_along(entitieslist)) {

		fitness <- entities_fitness(entitieslist[[i]], fitness_coef)
		min_fitness <- min(fitness)
		
		if (type == "cutoff") {
			weight <- fitness <= min_fitness+temp_or_thresh
		} else {
			if (temp_or_thresh != 0) {
				weight <- exp(-(fitness-min_fitness)/temp_or_thresh)
			} else {
				weight <- fitness == min_fitness
			}
		}
		
		weight <- weight/sum(weight)
		first_gen <- integer(length(fitness))
		for (j in rev(seq_along(generationslist[[i]]))) {
			first_gen[generationslist[[i]][[j]]] <- j
		}
		for (j in seq_along(generationslist[[i]])) {
			gen_contrib[i,j] <- sum(weight[first_gen == j])
		}
	}
	
	colnames(gen_contrib) <- seq_along(generationslist[[1]])
	
	boxplot(as.data.frame(gen_contrib), xlab="Generation", ylab="Sequence Contribution", main=main, ylim=c(0,1), yaxt="n")
	axis(2, seq(0, 1, by=.25), labels=FALSE)
	axis(2, seq(0, 1, by=.5), seq(0, 1, by=.5), tick=FALSE, line=FALSE)
}

# This function plots the fitnesses for those sequences that contribute
# significantly to the position weight matrix.

plot_seq_contrib <- function(entitieslist, 
                             fitness_coef = c(1/2.5, 1/2.5, 1/2.5, 1),
                             temp_or_thresh = 0.228, 
                             type = c("boltzmann", "cutoff"),
                             main = "") {

	type <- match.arg(type)
	
	if (temp_or_thresh == 0) {
		maxfit <- 1
	} else if (type == "boltzmann") {
		maxfit <- -log(.01)*temp_or_thresh*2
	} else {
		maxfit <- 3*temp_or_thresh
	}
	
	layout(matrix(1:2, nrow=2), heights=c(0.2, 0.8))
	
	mar1 <- mar2 <- par("mar")
	
	mar1[1] <- 0.1
	mar2[3] <- 0.1
	
	par(mar=mar1)
	
	plot(0, 0, xlim=c(0, maxfit), type="n", ylim=c(0, 1), xaxt="n", yaxt="n", xlab="", ylab="Weight", main=main)
	
	if (type == "cutoff") {
		segments(0, 1, temp_or_thresh, 1)
		segments(temp_or_thresh, 1, temp_or_thresh, 0)
		segments(temp_or_thresh, 0, maxfit*1.1, 0)
	} else {
		fitval <- seq(0, maxfit*1.1, length.out=100)
		points(fitval, exp(-fitval/temp_or_thresh), type="l")
	}
	
	axis(2, labels=FALSE)
	axis(2, c(0,1), tick=FALSE)
	
	par(mar=mar2)
	
	plot(0, 0, xlim=c(0, maxfit), type="n", ylim=c(length(entitieslist), 1), xlab="Normalized Fitness", ylab="Backbone")

	for (i in seq_along(entitieslist)) {

		fitness <- entities_fitness(entitieslist[[i]], fitness_coef)
		min_fitness <- min(fitness)
		
		norm_fit <- fitness-min_fitness
		plot_idx <- which(norm_fit < maxfit*1.1)
		
		points(norm_fit[plot_idx], rep(i, length(plot_idx)), pch=20, cex=.5)
	}
	
	if (type == "cutoff") abline(v=temp_or_thresh, lty="dashed")
}

# This function is the main data processing procedure. It takes a directory 
# path which contains *.ga.entities files. It reads all those files and
# produces a set of boxplots in several different file formats. It also 
# generates a position weight matrix and FASTA file for producing a sequence 
# logo. By specifying plotgen=TRUE, it will produce a plot similar to 
# Figure 5 in the PLoS One manuscript.

process_seqtol <- function(dirpath = ".", fitness_coef = c(1/2.5, 1/2.5, 1/2.5, 1),
                           temp_or_thresh = 0.228, 
                           type = c("boltzmann", "cutoff"),
                           percentile = .5, prefix = "seqtol",
                           plotgen = FALSE, plotseq = TRUE) {

	type <- match.arg(type)
	names(fitness_coef) <- paste("state1_fitness_comp", seq_along(fitness_coef), sep="")
	
	entities <- read_ga_entities_list(dirpath, readgen=plotgen)
	pwms <- entities_list_pwms(entities, fitness_coef, temp_or_thresh, type)
	pwm <- collapse_pwms(pwms, percentile)
	posnames <- colnames(pwm)
	
	inputseq <- NULL
	
	seqtoloutput <- file.path(dirpath, "seqtol_1_stdout.txt")
	if (!file.exists(seqtoloutput)) {
		seqtoloutput <- file.path(dirpath, sub(".ga.entities.*", "_seqtol.out", names(entities)[1]))
	}
	if (file.exists(seqtoloutput)) {
		seqtolcmd <- readLines(seqtoloutput, 2)
		seqtolcmd <- grep("core.init: command", seqtolcmd, value=TRUE)
		if (length(seqtolcmd)) {
			startpdbfile <- gsub("^.+ -s ([^ ]+) .+$", "\\1", seqtolcmd)
			if (!file.exists(startpdbfile)) {
				startpdbfile <- paste(startpdbfile, ".gz", sep="")
			}
			if (file.exists(startpdbfile)) {
				# Parse the file into a set of unique residue names/position ID pairs (GLU with name A318, ...)
				pdbseq <- pdb_sequence(file.path(dirpath, startpdbfile))
				# Index into pdbseq (residue names e.g. GLU named by position ID e.g. A318) with posnames
				# Store the one-character corresponding codes e.g. E into inputseq
				inputseq <- aa1[pdbseq[as.integer(posnames)]]
				colnames(pwm) <- posnames <- names(pdbseq)[as.integer(posnames)]
			}
		}
	}
	
	# Write pwm to a tab-separated file (do not add quotation marks)
	# row.names defaults to TRUE so we add a blank column name for the first column
	write.table(pwm, paste(prefix, "_pwm.txt", sep=""), quote=FALSE, sep="\t", col.names=NA)
	
	seqmat <- pwm_to_seqmat(pwm)
	
	cat(paste(">", seq_len(nrow(seqmat)), "\n", apply(seqmat, 1, paste, collapse=""), sep=""), file=paste(prefix, "_sequences.fasta", sep=""), sep="\n")
	
	plotwidth <- 7
	plotheight <- 3
	pointsize <- 12
	
	pdf(paste(prefix, "_boxplot.pdf", sep=""), width=plotwidth, height=plotheight, pointsize=pointsize)
	pdfdev <- dev.cur()
	png(paste(prefix, "_boxplot.png", sep=""), width=plotwidth*72, height=plotheight*72*length(posnames), pointsize=3/2*pointsize)
	pngdev <- dev.cur()
	par(mfrow=c(length(posnames), 1))
	
	for (i in seq_along(posnames)) {
		
		for (imgtype in c("pdf", "png", "pngsep")) {
			
			if (imgtype == "pdf") dev.set(pdfdev)
			if (imgtype == "png") dev.set(pngdev)
			if (imgtype == "pngsep") png(paste(paste(prefix, "_boxplot_", sep=""), posnames[i],".png", sep=""), width=plotwidth*72, height=plotheight*72, pointsize=pointsize)
			
			par(mar = c(2.8, 2.8, 1.5, 0.1), mgp = c(1.7, 0.6, 0))
			main <- paste("Residue", posnames[i], "Sequence Tolerance Boxplot")
			plot(0, 0, type="n", xlim=c(1,20), ylim=c(0,1), main=main, xlab="Amino Acid", ylab="Predicted Frequency", axes=FALSE)
			abline(h=seq(0, 1, by=.2), col="gray")
			boxplot(as.data.frame(t(pwms[,,i])), col="white", add=TRUE)
			points(1:20, pwm[,i], pch=4, col="blue")
			
			if (imgtype == "pngsep") dev.off()
		}
	}
	dev.off(pdfdev)
	dev.off(pngdev)
	
	seqrank_width <- 2.921+(1+ncol(pwm))*.2
	seqrank_height <- 4
	png_scale <- 1.5
	pdf(paste(prefix, "_seqrank.pdf", sep=""), width=seqrank_width, height=seqrank_height, pointsize=10)
	
	# inputseq is an array of one-character residue names e.g. E corresponding to the design
	plot_seqrank(pwm, wt_seq=inputseq, rank_line=5)
	dev.off()
	png(paste(prefix, "_seqrank.png", sep=""), width=seqrank_width*72*png_scale, height=seqrank_height*72*png_scale, pointsize=10*png_scale)
	plot_seqrank(pwm, wt_seq=inputseq, rank_line=5)
	dev.off()
	
	if (plotgen) {
		pdf(paste(prefix, "_gencontrib.pdf", sep=""), width=7, height=3, pointsize=pointsize)
		par(mar=c(2.7,2.7,1.5,0.2), mgp=c(1.7, 0.6, 0))
		plot_gen_contrib(entities, attr(entities, "generations"), fitness_coef, temp_or_thresh, type, "Generation Contributions")
		dev.off()
	}
	
	if (plotseq) {
		pdf(paste(prefix, "_seqcontrib.pdf", sep=""), width=6, height=6, pointsize=pointsize)
		par(mar=c(2.7,2.7,1.5,0.2), mgp=c(1.7, 0.6, 0))
		plot_seq_contrib(entities, fitness_coef, temp_or_thresh, type, "Sequence Contributions")
		dev.off()
	}
}

# This function is the main data processing procedure. It takes a directory 
# path which contains *.ga.entities files. It reads all those files and
# produces a set of boxplots in several different file formats. It also 
# generates a position weight matrix and FASTA file for producing a sequence 
# logo.
# THIS IS DEPRECATED AND WILL BE REMOVED IN THE FUTURE

process_specificity <- function(dirpath = ".", fitness_coef = c(1/2.5, 1/2.5, 1/2.5, 1),
                                temp_or_thresh = 0.228, 
                                type = c("boltzmann", "cutoff"),
                                percentile = .5) {

	warning("process_specificity is deprecated, please switch to process_seqtol")

	process_seqtol(dirpath, fitness_coef, temp_or_thresh, type, percentile, "specificity")
}
