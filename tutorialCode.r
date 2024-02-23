#seems to work on HPC using R version 4.2.1
#srun --nodes=1 --ntasks=50 --partition=normal --time=47:0:0 --mem=10G --pty bash -i
#module load R/4.2.1
#change directoory to where the tutorial pack as been unzipped
#then start R

##########################################################################33
#tutorial0_prerequisites
######################################################################33

#make sure all of these libraries load properly

#bioconductor
library("GenomicRanges")
library("DESeq2")
library("ACME")
library("GEOquery")
library("EnhancedVolcano")

#R 
library("stringr")
library("dplyr")
library("parallel")
library("rmarkdown")
library("knitr")
library("ggfortify")
library("data.table")
library("ggrepel")



#if there were any errors, use this code to install missing libraries
#bioconductor packages
bioc.packages = c(
    "GenomicRanges",
    "DESeq2",
    "ACME",
    "GEOquery",
    "EnhancedVolcano"
)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(bioc.packages)


#R packages
r.packages = c(
    "stringr",
    "dplyr",
    "parallel",
    "rmarkdown",
    "knitr",
    "ggfortify",
    "data.table",
    "ggrepel"
)
ix.installed = r.packages %in% installed.packages()[,1]
install.packages(r.packages[!ix.installed])



#################################################
#tutorial1_rna
#########################################33333

#########
#RNA section 2.6

#get the count files
pfiles = list.files(path = "data/counts_rna",
    pat = ".*star.*.tab", full = TRUE, recursive = TRUE)

#load the count files and just use the first 2 columns
pcounts = lapply(pfiles, function(f){
    print(f)
    try({
        x = read.csv(f, sep="\t", header=F)
        #remove bad rows
        x = x[grepl(x[,1], pat="ENSG"),]
        y = data.frame(ensg=x[,1], counts=x[,2])
        y
    }, silent=T)
})

#remove NAs
names(pcounts) = pfiles
ix = unlist(lapply(pcounts, function(a){class(a) != "try-error"}))
pcounts = pcounts[ix]
is.na(pcounts)

#merge into a data.frame
x = do.call(cbind, lapply(pcounts, function(a){
	a$counts
}))
#make the names easier to read
colnames(x) = gsub(names(pcounts), pat=".*/(.*)_ReadsPerGene.out.tab", rep="\\1")
rownames(x) = pcounts[[1]][,1]


#get rid of NA counts
ix2 = which(apply(x, 1, function(a){
    !any(is.na(a))
}))
x=x[ix2,]

#take out the rownames that have __ in them
#only use the ENSG rows
x = x[grepl(rownames(x), pat="ENSG"),]

#take out the dot?
rownames(x) = gsub(rownames(x), pat="\\..*", rep="")

head(x)


##################
#RNA section 3.1

#create a folder/directory for where geo downloads will be cached
if (!dir.exists("geo")){
    dir.create("geo")
}

#download geo data with id GSE132040 and save to "geo" dir
dat <- getGEO("GSE132040", destdir = "geo")
head(dat)


#load the phenotype data from GEO, obtain phenodata
#some GSE's have more than 1 of dataset, getGEO returns them in a list
phenodata <- pData(dat[[1]])
head(phenodata)
#Alternatively, you can do the following 
#phenodata <- dat[["GSE132040_series_matrix.txt.gz"]]@phenoData@data


###############
#RNA section 3.2

#Load the expression data
rnadata.geo <- exprs(dat[[1]])
dim(rnadata.geo) #0 rows?!

#the columns of rnadata should match up to the rows of phenodata)
#cbind(colnames(rnadata),rownames(phenodata))
all(colnames(rnadata.geo) == rownames(phenodata))



#if you want to use read.csv you'll need to unzip the file first
#rnadata <- read.csv("GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv")

#here's a trick to load large files using the data.table function.  It's much faster than read.csv and even works on zipped files
rnadata <- data.frame(fread("data/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv.gz"))
dim(rnadata)

#look at the first 5 rows and first 10 columns to get an idea of the data
#View(rnadata[1:5, 1:10])
#print(rnadata[1:5, 1:10])
rnadata[1:5, 1:10]


table(duplicated(rnadata$gene)) # check if any duplication gene names 


rownames(rnadata) <- rnadata$gene  # Assign the first column to the rownames
rnadata$gene <- NULL #remove the column from the rest of the matrix

class(rnadata)

rnadata[1:5, 1:10]

################
#RNA section 3.3

cbind(head(colnames(rnadata)))

cbind(head(phenodata$title))

# we can use gsub to manipulate the string. in this case, replace a part of a string.
colnames(rnadata) <- gsub(colnames(rnadata), pattern=".gencode.vM19", replace="")
                     
cbind(head(colnames(rnadata))) # after


# create a new column in phenodata called "tmp_id" 
# that has the sample identifier extracted from the title
# we can use replace again noticed I used \\ in front of ( ) and [ ], because they are special characters 
# a more general approach is to use regular expression 
phenodata$tmp_id <- gsub(phenodata$title, pattern="Tabula Muris Senis \\(bulk RNA-seq\\) \\[", replace="")
phenodata$tmp_id <- gsub(phenodata$tmp_id, pattern="\\]", replace="")

cbind(head(phenodata$tmp_id))

# lets also get rid of :chr1 in some colnames
colnames(phenodata) <- gsub(colnames(phenodata), pattern=":ch1", replace="")


#is every sample accounted for? yes
all(phenodata$tmp_id %in% colnames(rnadata))

#but they don't match
head(cbind(phenodata$tmp_id, colnames(rnadata)))

#we need to reorder the rows of phenodata to match the columns of rnadata

#better to work with a temporary variable so you can check your work
phenodata_reordered <- phenodata[match(colnames(rnadata), phenodata$tmp_id),]

#check it again (now they match)
head(cbind(phenodata_reordered$tmp_id, colnames(rnadata)))

#from now on we will replace phenodata with phenodata_reordered
phenodata <- phenodata_reordered
rm(phenodata_reordered)


#################
#RNA section 3.4

#run deseq on a subset of rnadata for demonstration purposes
table(phenodata$tissue, phenodata$age)

#lets just look at bone vs brain
#in 1 months postnatal mice

#extract the relevant rows from phenotype data
phenodata_small <- phenodata[phenodata$tissue %in% c("Bone", "Brain"),]
phenodata_small <- phenodata_small[phenodata_small$age %in% c("1 months postnatal"),]


#extract the relevant columns from count data
rnadata_small <- rnadata[,colnames(rnadata) %in% phenodata_small$tmp_id]

dim(phenodata_small)
dim(rnadata_small)

#make sure the IDs match up
table(phenodata_small$tmp_id == colnames(rnadata_small))

#use the new ids instead of the GSM ids
#deseq requires the rownames from the info match the colnames of the counts
rownames(phenodata_small) <- phenodata_small$tmp_id


#make sure the column names don't have any special characters like ":"
colnames(phenodata_small) <- make.names(colnames(phenodata_small))


#################
#RNA section 4.1

#it's always a good idea to explicitly set a categorical variable as a factor
#so you can control which is the baseline reference.  Otherwise, it sorts 
#alphabetically and chooses the first as the ref.
phenodata_small$tissue <- factor(phenodata_small$tissue, levels = c("Bone", "Brain"))

dds <- DESeqDataSetFromMatrix(countData=rnadata_small, colDat=phenodata_small, design = ~tissue)
dds <- DESeq(dds) #takes a minute
res <- results(dds)
res <- as.data.frame(res) #easier to work with data.frames

#take a look at the results in excel
write.table(res, file="res.csv",sep = ",",quote = F,row.names = T, col.names = NA)

summary(res)

head(res, 30)


###############
#RNA section 4.2

#remove results with no counts
res <- res[res$baseMean > 0,]
#sort by fold change and 
res <- res[order(res$log2FoldChange),] # turn decreasing = TRUE if want from the largest to the smallest 
#write the output to a spreadsheet

write.table(res,file="de_boneVsbrain_age1month.csv",sep = ",",quote = F,row.names = T,col.names = NA)

dim(res)

#how many genes show significant DE after adjusting for multiple comparisons?
sum(res$padj < .01, na.rm=T)

#what are the genes that have  a log2FC > 2 and an adjusted pvalue < 0.01?
res_up <- res[!is.na(res$log2FoldChange),]
res_up <- res_up[res_up$log2FoldChange > 2,]
res_up <- res_up[!is.na(res_up$padj),]
res_up <- res_up[res_up$padj < 0.01,]
res_up <- res_up[order(res_up$log2FoldChange,decreasing = TRUE),]

res_down <- res[!is.na(res$log2FoldChange),]
res_down <- res_down[res_down$log2FoldChange < -2,]
res_down <- res_down[!is.na(res_down$padj),]
res_down <- res_down[res_down$padj < 0.01,]


nrow(res_up)
nrow(res_down)

#look at the top 20 in each
head(rownames(res_up), 20)
head(rownames(res_down), 20)

##################
#RNA section 4.3

# when the expression data is alinged to the phenotype data
# it's easy to run all the usual informatics procedures 
dim(rnadata_small)


# the prcomp function complains when there are columns with 0 variance
# to fix it, remove the genes that have no variance
gene_var <- apply(rnadata_small, 1, var) > 0 
table(gene_var)

rnadata_filtered <- rnadata_small[gene_var,]
dim(rnadata_filtered)
pca1 <- prcomp(t(rnadata_filtered), scale=T)


#plot the first 2 principal components, colored by tissue type
#pdf("test1.pdf")
autoplot(pca1, data=phenodata_small, colour='tissue', shape="Sex")
#dev.off()



#################
#RNA section 4.4

vlnPlot1 <- EnhancedVolcano(res,
    lab = rownames(res),
    #subtitle = NULL, #get rid of the subtitle
    #colCustom = keyvals.colour, #give customized color
    x = 'log2FoldChange',
    y = 'padj',
    #xlim = c(-30,30), #x axis range
    #ylim = c(0,300), # y axis range
    title = 'DE_Bone_vs_Brain: 1month', # title label
    xlab= bquote(~Log[2]~ 'fold change'), # x axis label
    ylab= bquote(~-Log[10]~ 'Padj'), # y axis label
    caption = NULL,
    pCutoff = 0.05,
    FCcutoff = 2,
    pointSize = 2.0)

#pdf("test2.pdf")
vlnPlot1
#dev.off()



#################################################
#tutorial2_atac
#########################################33333


##############
#ATAC section 2.6.1

#load in the peak files
pfiles = list.files(path="data/macsResults", recursive=T, pat="^.*xls$", full=T)
peaks = lapply(pfiles, read.table, header=T)

#contruct genomic ranges of the peak summits +/- 200bp
ranges = lapply(peaks, function(peak){
	GRanges(seqnames = peak[,1],
        strand=rep("*", nrow(peak)),
            IRanges(start = peak[,2], end = peak[,3]))
})
peaks.reduced = reduce(do.call(c, ranges))

#save genomic range as a bed file
#https://www.biostars.org/p/89341/
grange2bed <- function(gr, bfile="foo.bed"){

	df <- data.frame(seqnames=seqnames(gr),
		starts=start(gr)-1,
		ends=end(gr),
		names=c(rep(".", length(gr))),
		scores=c(rep(".", length(gr))),
		strands=strand(gr))

	write.table(df, file=bfile, quote=F, sep="\t", row.names=F, col.names=F)
}

grange2bed(peaks.reduced, bfile="peaks_stat.bed")



##############
#ATAC section 2.7
pfiles = list.files(path="data/counts_atac", pat="*counts", full=T)

#load in 1 file to have a look at the format
x = read.csv(pfiles[1], sep="\t", header=F)
head(x)

#process all the count files
pcounts = lapply(pfiles, function(file){
    #print(file)
	x = read.csv(file, sep="\t", header=F)
    #extract the number that appears after the *| in the 6th column of the input
	hits = as.numeric(gsub(x[,6], pat="\\*\\|", rep=""))
	y = data.frame(chr=x[,1], start=x[,2], end=x[,3], counts=hits)
	y
})
#bind the counts together as columns in a matrix
peakdata = do.call(cbind, lapply(pcounts, function(a){
	a$counts
}))
#give the columns better looking names. 
#grab the strings between the last "/" and ".counts" 
colnames(peakdata) = gsub(pfiles, pat=".*/(.*).counts", rep="\\1")


#peak locations
peaks = pcounts[[1]][,1:3]
head(peaks)
#make a descriptive string ID for the peak 
peaknames = sapply(1:nrow(peaks), function(i){
	paste0(peaks[i,1:3], collapse="_")
})
rownames(peakdata) = peaknames
cbind(head(rownames(peakdata)))



#get rid of NA counts
ix2 = which(apply(peakdata, 1, function(a){
    !any(is.na(a))
}))
peakdata=peakdata[ix2,]

dim(peakdata)
head(peakdata)



##########
#ATAC section 3

info <- colnames(peakdata)
comp1 <- rep(NA, length(info))
comp1[1:2] <- "ctrl"
comp1[3:4] <- "drug_treated"
coldata <- data.frame(info, comp1)

#always make sanity checks with your own data
#the row order of coldata needs to match the column order of peakdata
all(coldata$info == colnames(peakdata))
head(coldata)


#run the analysis
dds <- DESeqDataSetFromMatrix(countData = peakdata,
								  colData = coldata,
								  design = ~ comp1)
dds <- DESeq(dds)
res <- results(dds)
res <- as.data.frame(res)

head(res, 50)


##########
#ATAC section 3.1


#run DESeq, add peak info, sort and write as a csv.
#countdata: a matrix of integer counts, (nPeak x nSamples)
#coldata: a data.frame with info on the samples (nSamples x nVariables)
#ix: the indices of the samples you want to use in the analysis.
#formula: shows which columns of coldata you want to use in the analysis
#resFile: is where to save the result file
runDE <- function(countdata, coldata, ix, formula = formula("~ ep300"), resFile){

	dds = DESeqDataSetFromMatrix(countData = countdata[,ix],
								  colData = coldata[ix,],
								  design = formula)
	dds <- updateObject(dds)
	dds = DESeq(dds)
	res = results(dds)
	y = data.frame(res)

    #break the peak identifier into parts
    pos = data.frame(pid = rownames(y), str_split_fixed(rownames(y), "_", 3))
    colnames(pos)[2:4] = c("chrom", "start", "stop")
    pos$start = as.numeric(pos$start)
    pos$stop = as.numeric(pos$stop)
    pos$mid = round((pos$stop+pos$start)/2)
	y = cbind(pos, y)

    #add peak counts
	diffexp = merge(y, countdata[,ix], by.x=1, by.y=0, all=T)

	#append the normalized counts to the matrix
	norm.counts = counts(dds, normalized=T)
	colnames(norm.counts) = paste0("norm_",dds$info)
	y = merge(diffexp, norm.counts, by.x=1, by.y=0, all=T)

    #sort results by fold change
	y = y[order(y[,"log2FoldChange"]),]

    #add a chromosomal position to make IGV navigation easier
    y$chrpos = paste0(y$chrom, ":", y$start, "-", y$stop)

	write.table(file=resFile, x=y, row.names=F, sep="\t")

    #also return the results
    y
}


##########
#ATAC section 3.2

#use the non-NA samples in comp1
ix = which(!is.na(coldata$comp1))
#or you can do it like
#ix = 1:4

#run the analysis
res = runDE(peakdata, coldata, ix, formula("~ comp1"), "de_CtrlVsDrug.txt")

#view some of the results
head(res, 50)


##########
#ATAC section 3.3

appendClosestGeneInfo <- function(x, chr=1, pos = 2, genome="hg19", numCores = 1){
	res = mclapply(1:nrow(x), function(i){
		if (i %% 1000 == 0){
			print(i)
		}
		genes = findClosestGene(x[i,chr],x[i, pos],genome)
		#if there's more than one, just take the first
		data.frame(genes[1,])
	}, mc.cores = numCores)

	x.res = bind_rows(res)
	cbind(x.res, x)
}

#how many peaks show significant DE?
ix.sig = which(res$padj < 0.05)
length(ix.sig)

#show off the function using a small portion of the data
res2 = appendClosestGeneInfo(x=res[ix.sig,], chr=2, pos=5, genome="hg38") 

#if you get an error, try running on a single thread
#res2 = appendClosestGeneInfo(x=res[ix.sig,], chr=2, pos=5, genome="hg38", numCores=1) 

#if you still get an error, there may be issues with the web host(see below)

head(res2, 50)
write.csv(file="deAnno_padj0.05_CtrlVsDrug.txt", x=res2, row.names=F)



#################################
#I've run into intermittent problems with the code above.
#the original site that hosts the reference database occasionally goes offline
#If you still want to use the code above, you can hack some of the functions inside
#the ACME package to make it point to a different host
mygetRefflat <- function (genome = "hg38"){
    tmpfile <- tempfile()
    download.file(paste("http://hgdownload.soe.ucsc.edu/goldenPath/", 
        genome, "/database/refFlat.txt.gz", sep = ""), tmpfile, 
        mode = "wb")
    rf <- read.delim(conn <- gzfile(tmpfile, open = "rt"), header = FALSE, 
        sep = "\t")
    close(conn)
    colnames(rf) <- c("geneName", "name", "chrom", "strand", 
        "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", 
        "exonStarts", "exonEnds")
    txEndNeg <- rf$txStart
    txStartNeg <- rf$txEnd
    cdsStartNeg <- rf$cdsEnd
    cdsEndNeg <- rf$cdsStart
    NegStrand <- rf$strand == "-"
    rf[NegStrand, "cdsEnd"] <- cdsEndNeg[NegStrand]
    rf[NegStrand, "cdsStart"] <- cdsStartNeg[NegStrand]
    rf[NegStrand, "txEnd"] <- txEndNeg[NegStrand]
    rf[NegStrand, "txStart"] <- txStartNeg[NegStrand]
    return(rf)
}


myfindClosestGene <- function (chrom, pos, genome = "hg38", position = "txStart"){
    if (!exists("refflat")) {
        reftmp <- list()
        reftmp[[genome]] <- mygetRefflat(genome)
        assign("refflat", reftmp, .GlobalEnv)
    }
    else if (!(genome %in% names(refflat))) {
        refflat[[genome]] <<- mygetRefflat(genome)
    }
    rf <- refflat[[genome]]
    chromsub <- rf$chrom == chrom
    diffdist <- rf[chromsub, position] - pos
    sub <- which(abs(diffdist) == min(abs(diffdist)))
    rf <- rf[chromsub, 1:9][sub, ]
    return(data.frame(rf, Distance = diffdist[sub]))
}

#now you can try again
appendClosestGeneInfo2 <- function(x, chr=1, pos = 2, genome="hg38", numCores = 4){
	res = mclapply(1:nrow(x), function(i){
		if (i %% 1000 == 0){
			print(i)
		}
		genes = myfindClosestGene(x[i,chr],x[i, pos],genome)
		#if there's more than one, just take the first
		data.frame(genes[1,])
	}, mc.cores = numCores)

	x.res = bind_rows(res)
	cbind(x.res, x)
}


res2 = appendClosestGeneInfo2(x=res[ix.sig,], chr=2, pos=5, genome="hg38") 
head(res2, 50)
write.csv(file="deAnno_padj0.05_CtrlVsDrug.txt", x=res2, row.names=F)




