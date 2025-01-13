library(DiffBind)
library(GenomicRanges)
library(limma)
library(ChIPpeakAnno)
library(Rsamtools)
library(GenomicAlignments)
library(data.table) 
library(foreach) 
library(doMC) 

# Set the number of cores for parallel processing
N.CORES <- 30
registerDoMC(cores=N.CORES)

# Read the meta.table
meta.table <- read.table("meta.table.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")

# Update paths in meta.table if necessary
meta.table$bamReads <- meta.table$Bam_Filtered_Path
meta.table$Peaks <- meta.table$Peak_narrow_Path
str(meta.table)

# Ensure all files exist
stopifnot(all(file.exists(meta.table$bamReads)))
stopifnot(all(file.exists(meta.table$Peaks)))

# Adjusted processPeaks function
processPeaks <- function(peakPath) {
  if (is.na(peakPath) || length(readLines(peakPath)) == 0) {
    return(NULL) # Handle empty or non-existent files
  }
  peakData <- fread(peakPath, header = FALSE) # Read without expecting a header
  setnames(peakData, c("chromosome", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
  significantPeaks <- peakData[, .(chromosome, start, end)]
  makeGRangesFromDataFrame(significantPeaks, keep.extra.columns = TRUE)
}

# Process Peaks
peaksList <- lapply(meta.table$Peaks, processPeaks)
peaksList <- Filter(Negate(is.null), peaksList) # Remove NULL elements

# Convert to GRangesList and reduce
peaksGRangesList <- GRangesList(peaksList)
consensusPeaks <- reduce(unlist(peaksGRangesList))


# GET READS ---------------------------------------------------------------
# Process each BAM file
processBamFile <- function(bamPath, consensusPeaks) {
  bamFile <- BamFile(bamPath, index = paste0(bamPath, ".bai"), yieldSize = 2000000L, asMates = FALSE)
  
  result <- tryCatch({
    countData <- summarizeOverlaps(features=consensusPeaks, reads=bamFile, ignore.strand=TRUE, singleEnd=TRUE, mode="Union", param=ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE, isSecondaryAlignment=FALSE)))
    assays(countData)$counts
  }, error=function(e) {
    message("Failed to process bam file: ", bamPath, "\nError: ", e$message)
    NULL
  })
  
  return(result)
}

# Iterate over each row in meta.table to process BAM files
res.peaks <- foreach(bamPath = meta.table$bamReads, .packages = c("Rsamtools", "GenomicAlignments", "GenomicRanges")) %dopar% {
  processBamFile(bamPath, consensusPeaks)
}

# Combine results and clean up column names
#countMT <- do.call(cbind, res.peaks)
#colnames(countMT) <- gsub(".filtered.bam$", "", meta.table$sample_name)

# Filter out NULL values from the list of results
res.peaks.non.null <- res.peaks[sapply(res.peaks, function(x) !is.null(x))]

if (length(res.peaks.non.null) > 0) {
  # Combine the non-NULL results. This will give a matrix if elements are vectors or matrices.
  combinedResults <- do.call(cbind, res.peaks.non.null)
  
  # Convert the combined results into a dataframe
  countMT <- as.data.frame(combinedResults)
  
  # Update column names. This assumes the sample names are correctly aligned with the non-NULL entries.
  colnames(countMT) <- gsub(".filtered.bam$", "", meta.table$sample_name[sapply(res.peaks, function(x) !is.null(x))])
} else {
  message("All results were NULL. Please check your input data and processing steps.")
}


# Ensure only valid results are included
countMT <- countMT[, !sapply(countMT, is.null)]

#Verification of Sample Names in Outputs
expectedSampleNames <- meta.table$sample_name

# Ensure that all processed sample names are in the expected sample names list
allProcessedNamesMatch <- all(colnames(countMT) %in% expectedSampleNames)

# Ensure that all expected sample names have been processed
allExpectedNamesProcessed <- all(expectedSampleNames %in% colnames(countMT))

# Check if both conditions are TRUE
if(allProcessedNamesMatch && allExpectedNamesProcessed) {
  message("All sample names match correctly.")
} else {
  message("There is a mismatch in sample names.")
  # Optional: Identify which names are mismatched
  mismatchedNames <- setdiff(union(expectedSampleNames, colnames(countMT)), intersect(expectedSampleNames, colnames(countMT)))
  message("Mismatched names: ", paste(mismatchedNames, collapse = ", "))
}


# Annotate Peaks -----------------------------------------------------------

# Region to Gene(s) map
gr <- consensusPeaks
stopifnot(length(gr) == length(unique(gr)))
data(TSS.mouse.GRCm38)
x <- annotatePeakInBatch(gr, AnnotationData=TSS.mouse.GRCm38,output="overlapping", maxgap=5000)
x <- addGeneIDs(x, orgAnn="org.Mm.eg.db", IDs2Add=c("symbol"))
gmap <- data.table(as.data.frame(x))
gmap[,region:=paste(seqnames, start,end, sep="_")]
save(gmap, file=("./OutATAC-seq_eosinophils/GeneMap.Full.RData"))
gmap[,middle := start+width/2]
gmap[,distance := ifelse(feature_strand == "+", middle - start_position, end_position - middle)]
gmap[,distance := abs(distance)]
gmap <- gmap[,c("region", "width", "feature", "symbol", "distance"),with=F]
colnames(gmap)[1:5] <- c("probe", "width", "ensG", "gene", "distance")

# Filter Blacklisted Regions
# Note: Adjust the path to the blacklist file as per your reference genome.
#testblacklist <- blacklist <- fread("./mm10-blacklist.v2.bed")
#head(testblacklist)

blacklist <- fread("./mm10-blacklist.v2.bed")[,1:3, with=F]
colnames(blacklist) <- c("CHR", "START", "END")
blacklist <- as(blacklist, "GRanges")
ol <- findOverlaps(blacklist, gr)
goodRegions <- data.table(as.data.frame(gr[-subjectHits(ol)]))
goodRegions.str <- with(goodRegions, paste(seqnames, start, end, sep="_"))
goodRegions.str <- goodRegions.str[goodRegions.str %in% gmap$probe]
gmap[,blacklisted := "yes"]
gmap[probe %in% goodRegions.str, blacklisted :="no"]
table(gmap$blacklisted)

# OUTPUT ------------------------------------------------------------------
row.names(countMT) <- with(consensusPeaks, paste(seqnames, start, end, sep="_"))
stopifnot(all(colnames(countMT) %in% meta.table$sample_name))
write.table(countMT, sep=",", quote=F, row.names=T, file=(paste0( "./OutATAC-seq_eosinophils/Counts.csv")))
write.table(gmap, (paste0( "./OutATAC-seq_eosinophils/Probes.csv")), sep="\t")
write.table(data.table(as.character(seqnames(consensusPeaks)), start(consensusPeaks), end(consensusPeaks)), col.names=F, (paste0("./OutATAC-seq_eosinophils/Peaks.bed")), sep="\t")
save(gmap, countMT, file=(paste0("./OutATAC-seq_eosinophils/.RData")))

