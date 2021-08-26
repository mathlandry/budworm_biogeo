library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(picante); packageVersion("picante")
library(reshape2)
library(seqTools)
# set file paths
path <- "/mnt/storage/mathieu/runfeb2019_DADA2" # CHANGE ME to the directory containing the fastq files after unzipping.
#path <- "/data/lab_data_archive/2017-06-miseq-run/dada2/idemp-demultiplexed" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="reads2.", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="reads1.", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
sample.names <- gsub('.{9}$', '', sample.names)
sample.names

# Plot quality profile of forward reads
plotQualityProfile(fnFs[1:4])
# Plot quality profile of reverse reads
plotQualityProfile(fnRs[1:4])

# need to trim 20 from start of both reads (gets rid of primer)
# Forward reads quality crashes around 175
# Reverse reads quality crashes around 150

# set filtered file folder path
filt_path <- file.path(path, "dada2-filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter reads at quality crashpoints identified earlier
# multithread argument = number of cores to use
# note need to trim out primers at the start of the reads
# and trim the end of the reads when the quality crashes
# original - try to merge by overlap (did not work well)
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(20,120), truncLen=c(220,185),
#                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#                     compress=TRUE, multithread=8, verbose=TRUE) # On Windows set multithread=FALSE
# current - trim to avoid overlap, concatenate reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(20,20), truncLen=c(220,220),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=4, verbose=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn error rates
# multithread argument = number of cores to use
errF <- learnErrors(filtFs, multithread=4)
errR <- learnErrors(filtRs, multithread=4)

# sanity check - visualize error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE) #here
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=4)
dadaRs <- dada(derepRs, err=errR, multithread=4)
# e.g. inspect results
dadaFs[[1]]

# Merge the denoised forward and reverse reads:
# TODO for now there is a problem - many many mismatches - why?
## tried - remoev first 20, trim reads shorter, try removing first 120 of reverse? nothing changes.
## currently - trimming reads shorter and allowing mismatches! re-do with trim&concat approach...

# original - try to overlap (does not work well)
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=10, returnRejects=TRUE, verbose=TRUE)
# current - concatenate only
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=4, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# summary - track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

# identify taxonomy
# note had to add tryRC=TRUE option for it work properly
taxa <- assignTaxonomy(seqtab.nochim, "/mnt/storage/mathieu/runfeb2019_DADA2/SILVA128-training-dada2/silva_nr_v128_train_set.fa.gz", multithread=4, tryRC=TRUE)
# exact species matching
# won't work if concatenating sequences
#taxa <- addSpecies(taxa,  "/data/shared/SILVA128-training-dada2/silva_species_assignment_v128.fa.gz")

# inspect the taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#write output tables 
write.csv(seqtab.nochim, "community_DADA2.csv")
write.csv(taxa, "taxonomy_DADA2.csv")
