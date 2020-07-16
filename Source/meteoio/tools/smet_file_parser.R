#####################################################################################
#  Copyright 2019 Michael Reisecker                                                 #
#####################################################################################
#  This is free software: you can redistribute and/or modify it under the terms of  #
#  the GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
#####################################################################################

#SYNOPSIS: smet <- read_smet("<smetfile>"); print(smet@data$TA[[1]])

#Methods to parse a SMET file. Note that there exists a more powerful and complicated
#R package for this called RSMET (CRAN archive).

library(methods)

########################################
#  SMET CLASS DEFINITION               #
########################################

setClass("SMET", slots = list(header = "list", data = "data.frame", time = "POSIXlt",
    fields = "character", station_id = "character"))

setMethod("initialize", signature = "SMET",
	definition = function(.Object, infile) {
		.Object <- read.smet(.Object, infile)
	  	return(.Object)
	}
)

setGeneric("read.smet",
	function(smet, infile) {
		standardGeneric("read.smet")
	}
)

setMethod("read.smet", signature = "SMET",
	function(smet, infile) { #read and process a smet file
		raw <- get_raw_smet(infile)

		smet@header <- raw[[1]]
		smet@data <- as.data.frame(raw[[2]])
		smet@station_id <- get_header_value(smet, "station_id")
		smet@fields <- get_header_value(smet, "fields")

		colnames(smet@data) <- smet@fields

		#convert time strings to time list and save to dedicated slot:
		smet@time <- as.POSIXlt(smet@data$time, format = "%Y-%m-%dT%H:%M:%OS")
		smet@data$time <- NULL #remove time from dataframe
		return(smet)
	}
)

########################################
#  SMET METHODS                        #
########################################

setGeneric("get.timeframe", #timeframe of dataset in days
	function(smet) {
		standardGeneric("get.timeframe")
	}
)

setMethod("get.timeframe", signature = "SMET",
	function(smet) { #return difference in days:
		return(difftime(smet@time[length(smet@time)], smet@time[1], units = "days"))
  }
)

########################################
#  FILE READING                        #
########################################

get_raw_smet <- function(infile) { #read header and data lines from file system
	idx.datasection <- grep("\\[DATA\\]", readLines(infile), value = FALSE) #return index of line "[DATA]"
	data <- read.csv(file = infile, sep = "", skip = idx.datasection, head = FALSE, na.strings = "-999",
	    stringsAsFactors = FALSE, strip.white = TRUE)

	idx.headersection <- grep("\\[HEADER\\]", readLines(infile), value = FALSE) #return index of line "[HEADER]"
	header <- read.table(file = infile, sep = "=", skip = idx.headersection, nrows = idx.datasection - idx.headersection - 1,
	    stringsAsFactors = FALSE, strip.white = TRUE)

	return(list(header, data))
}

get_header_value <- function(smet, key) { #find a key in the header and return its value ("key = value")
	idx <- grep(key, smet@header[[1]])
	str <- smet@header[[2]][idx]
	return(strsplit(str, " ", fixed = TRUE)[[1]]) #unlist
}

########################################
#  WRAPPER                             #
########################################

read_smet <- function(inputfile) { #wrapper for class initialization
	smet <- new("SMET", infile = inputfile)
	return(smet)
}

