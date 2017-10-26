#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("tools", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))

# Get command line options & arguments

arguments <- parse_args(OptionParser(usage = "%prog mtbr.folder output.mtbr.csv"), positional_arguments = 2)
kMtbrPath <- arguments$args[1]
kOutputCsv <- arguments$args[2]

if(!file.exists(kMtbrPath)){
  stop("mtbr path \"", kMtbrPath ,"\" does not exist.")
}

output.filename <- kOutputCsv
if (file.exists(output.filename)) {
	unlink(output.filename)
}

cat("chrom,posi,rC_n,rC_p,rT_n,rT_p\n", file = output.filename)
mtbr.files <- list.files(kMtbrPath, full.names = TRUE)
for (mtbr.file in mtbr.files) {
	message("converting ", mtbr.file)
	load(mtbr.file)
	cg.mtbr$posi <- as.integer(cg.mtbr$posi)
	cg.mtbr$rC_n <- as.integer(cg.mtbr$rC_n)
	cg.mtbr$rC_p <- as.integer(cg.mtbr$rC_p)
	cg.mtbr$rT_n <- as.integer(cg.mtbr$rT_n)
	cg.mtbr$rT_p <- as.integer(cg.mtbr$rT_p)

	write.table(cg.mtbr, output.filename, row.names = FALSE, quote = FALSE, append = TRUE, col.names = FALSE, sep = ",")
}

