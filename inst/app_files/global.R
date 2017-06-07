load("preparedData.RData")

fileBaseNames <- function(stringlist) {
  stringlist <- as.character(stringlist)
  for (i in 1:length(stringlist)) {
    # take the first filename, split by directory separator.
    splitstring <- strsplit(stringlist[i], '/')[[1]]
    # keep the last three elements (corresponds to the _processed directory and onwards)
    nstrings <- length(splitstring)
    stringlist[i] <- paste(splitstring[(nstrings-3):nstrings], collapse="/")
  }
  return(stringlist)
}

# replace the filenames in each statsList entry with the truncated, relative version
for (i in 1:length(statsList)) {
  statsList[[i]]$filenames <- fileBaseNames(statsList[[i]]$filenames)
}



