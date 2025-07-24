#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)

for (arg in seq(1, length(argv), 2)){
	print(argv[arg])
	#Initialize header row and column
	species <- list.files(path = ".",pattern = argv[arg], full.names = F, recursive = T)
	genes <- read.csv(species[1], header = T, stringsAsFactors = F)
	result <- as.data.frame(genes[,1])
	names(result) <- "Gene"
	result[2:(length(species) + 1)] <- NA

	#Extract the information in each file and combine
	column <- if (argv[arg] == "depth_of_coverage_.*csv") 6 else 2
	for (i in 1:length(species)){
			gene_data <- read.csv(species[i], header = T, stringsAsFactors = F)
			result[,i+1] <- gene_data[,column]
			names(result)[i+1] <- strsplit(species[i], "/")[[1]][1]
	}

	#Calculate averages for each gene
	gene_avg <- data.frame("Gene Average" = rowMeans(result[,-1]))
	result[,ncol(result) + 1] <- NA
	result <- cbind(result, gene_avg)
	
	#Calculate averages for each species
	species_avg <- matrix(c("Average", colMeans(result[2:(length(species) + 3)], na.rm = T, dims = 1)), nrow = 1)
	species_avg <- data.frame(species_avg)
	names(species_avg) <- names(result)
	result[nrow(result) + 1,] <- NA
	result <- rbind(result, species_avg)
	
	#Pretty up the matrix
	names(result)[length(species) + 2] <- ""
	result[nrow(result), length(species) + 2] = result[nrow(result), length(species) + 3] = ""
	
	#Export
	write.csv(result, file = argv[arg + 1], row.names = F, na = "")
}
