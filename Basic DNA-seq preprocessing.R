library(tidyverse)
library(edgeR)


#read in data
  counts <- read.csv(url("https://research.stowers.org/compbio/course/r_intro/data/star_count.csv"),row.names = 1)

#overwrite old column names with something nicer:
colnames(counts) <- gsub("^s_slj866\\d_*","",colnames(counts))
samples <- colnames(counts)

head(counts)

#creating a logical index vector that will help us remove the non genes from df
non_genes.iv <- grepl("^N_",rownames(counts))

#save these values in case they become important
align_stats <- counts[non_genes.iv,]

#using the LIV to remove the non genes from the count df
counts <- counts[!non_genes.iv,]

#more data cleaning, remove _repX from end of column names, these are now groups rather than replicates
#These groups will be used to set up our differential expression test
groups <- gsub("_rep\\d$","",colnames(counts))


#create a DGEList object using counts df and groups vector
y <- DGEList(counts=counts, group=groups)

#Now we will normalize counts
normCounts <- cpm(y,normalized.lib.sizes=T)
colnames(normCounts) <- paste("cpm.",samples,sep='')
tempnorm <- data.frame(normCounts)
tempnorm$gene <-rownames(tempnorm)

#making the data tidy for processing
normCounts_tidy <- pivot_longer(data = tempnorm, cols= -gene, names_prefix="cpm.",names_to="sample",values_to = "normCount")

#we are removing the lowly expressed genes, by default
keep.iv<-filterByExpr(y)

#checking how many we kept
length(which(keep.iv))/length(y$counts[,1])


design <- model.matrix(~0+y$samples$group)
colnames(design) <- levels(y$samples$group)
y <-estimateDisp(y,design)
fit <- glmQLFit(y, design)
mut_v_wt <- c(-1,1)

#perform differential expression test
qlf.mut_v_wt <- glmQLFTest(fit, contrast = mut_v_wt)

#output
result <- qlf.mut_v_wt$table

#adjust p-values for multiple hypothesis testing
result$padj <-p.adjust(result$PVal, method ="BH")

#sort by p-value
result.sort <-result[order(result$padj),]

#MA Plot from the result df. x axis is logCPM, y is logFC. general amount of variation in data
ggplot(result,aes(x=logCPM,y=logFC))+geom_point()+ ggtitle("MA Plot") + xlab("logCPM") + ylab("logFC")

#filter how many genes have a logFC greater than 1 and a PVaue less than .05
query.results <- filter(result,PValue<0.05 &logFC>1)
dim(query.results)

#tidy up the counts dataframe
tempcounts <- counts
#make a column from the rownames
tempcounts$gene <- rownames(tempcounts)
counts_tidy <- pivot_longer(tempcounts, cols = -gene, names_to="Sample",values_to="count")

#barplot of gene YBR171W (filter to get here)

counts_gene <- filter(counts_tidy,gene=="YBR171W")
ggplot(counts_gene,aes(x=Sample,y=count))+geom_col()


#boxplot log2(count)

counts_tidy$log2Counts <- log2(counts_tidy$count)

ggplot(counts_tidy,aes(y=log2Counts,x=Sample))+geom_boxplot()


#boxplot of log2(normCount)
normCounts_tidy$log2Counts <- log2(normCounts_tidy$normCount)
ggplot(normCounts_tidy,aes(y=log2Counts,x=sample))+geom_boxplot()


download.file("http://research.stowers.org/compbio/course/r_intro/data/membrane_genes.txt","data/membrane_genes.txt")
download.file("http://research.stowers.org/compbio/course/r_intro/data/allgenes.txt","data/allgenes.txt")
download.file("http://research.stowers.org/compbio/course/r_intro/data/degenes.txt","data/degenes.txt")

memb<-read.table("data/membrane_genes.txt",sep='\t',header=F)
allgenes<-read.table("data/allgenes.txt",sep='\t',header=T,quote="")
degenes<-read.table("data/degenes.txt",sep='\t',header=T)
#install clusterProfiler for GO analysis and eulerr to make venn diagrams. post to teams if you have issues with these, they worked on my computer...
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version='devel')
BiocManager::install("clusterProfiler")
BiocManager::install("org.Sc.sgd.db")
library("clusterProfiler")
library("org.Sc.sgd.db")
library("eulerr")
library("ggplot2")


#pulls out the genes in degenes & memb
degenes1 <- degenes$sys_id
memb1 <- memb$V1
#venn diagram of two gene datasets
sets <- list(memb2 = memb1, degenes2 = degenes1)
vcombo <- euler(sets)
plot(vcombo, quantities = T, fills = list(fill = c("red", "steelblue4"), alpha = 0.5))

#SI gene YLR152C in the intesection?

intersection <- intersect(memb1,degenes1)
grep("YLR152C",intersection)


#Go Enrichment Analysis


