########################################################
## This is the final Pathway Analysis ##
## Need to push this up to github
## October 7 ####
######################################################
rm(list = ls())

######################
# Pathway Analysis 
#####################
####################
## Gene Set Analysis
#setwd("/Users/asigdel/Documents/Conception Per Insemination/Geneset_Analysis")
setwd("/Users/asigdel/Documents/Documents - Anil’s MacBook Air/Conception Per Insemination/Geneset_Analysis")

## FILE freqdata.count.after.clean
filter = read.table("freqdata.count.after.clean", header = F)
colnames(filter) = c("SNP", "Frequency", "Filter")
f = filter$Filter == 0; table(f)


# Prepare data for Gene Set Analysis 
total_genes <- read.table("unique_all_slopes_Total_genes.txt", header= T)
head(total_genes)
total.genes <- total_genes$x;length(total.genes) # 19305

Sig_genes <- read.table("genes_repeat_2_3times.txt",header = T)
head(Sig_genes);dim(Sig_genes) # 928
sig.genes <- unique(Sig_genes$Var1);length(sig.genes) # 928

# Unique genes that are repeated 2 times 
genes_2_times <- Sig_genes[Sig_genes$Freq == 2,]
length(unique(genes_2_times$Var1)) # 765unique genes are repeated 2 times

# Unique genes that are repeated 3 times 
genes_3_times <- Sig_genes[Sig_genes$Freq==3,]
length(unique(genes_3_times$Var1)) # 163 unique genes are repeated 3 times

########################################################################
# KEGG Analysis
# 2019-08-23
#############################################################################
library(goseq)
library(org.Bt.eg.db)
library(GO.db)
library(KEGG.db)
library(GOSemSim)
library(corrplot)

assayed.genes = array(total.genes)
de.genes = array(sig.genes)
gene.vector = as.integer(assayed.genes %in% de.genes)
names(gene.vector) = assayed.genes
table = table(gene.vector)

cat(paste("Significant Genes: ", table[2], 
          " Backgroung Genes:", table[1], "\n")) 

pwf = nullp(gene.vector, "bosTau4", "ensGene", plot.fit = FALSE)
head(pwf)
############################################################################
KEGG.hiper <- goseq(pwf, "bosTau4", "ensGene", 
                    method = "Hypergeometric", test.cats = "KEGG", use_genes_without_cat = TRUE)

nKEGG.hiper = KEGG.hiper[(KEGG.hiper$numInCat <= 500 & KEGG.hiper$numInCat >= 5),] 
enriched.KEGG = nKEGG.hiper$category[nKEGG.hiper$over_represented_pvalue <= 0.05]
length(enriched.KEGG)

kegg = as.list(KEGGPATHID2NAME)
for(j in 1:length(enriched.KEGG)){
  for (i in 1:length(names(kegg))) 
  {
    if(names(kegg[i]) == enriched.KEGG[j]){
      
      cat("#############################################\n")	
      cat(paste("KEGG ID: ", enriched.KEGG[j], "\n"))
      cat(paste("KEGG Term Name: ", kegg[[i]], "\n"))
      pvalue = nKEGG.hiper[nKEGG.hiper$category == enriched.KEGG[j],]$over_represented_pvalue
      pvalue = round(pvalue,4)
      nGenes = nKEGG.hiper[nKEGG.hiper$category == enriched.KEGG[j],]$numInCat
      nDEG = nKEGG.hiper[nKEGG.hiper$category == enriched.KEGG[j],]$numDEInCat
      cat(paste("Total Genes:", nGenes, " Number DEG:", nDEG, " P-value:", pvalue))
      cat("\n")
    }
  }
}


## Revealing Significant Genes in KEGG Terms

KEGGterm = "04114"    #change the term
kegg = getgo(names(gene.vector), "bosTau4", "ensGene", fetch.cats = c("KEGG"))
genes = as.character(); m = 1

for(i in 1:length(kegg)){
  if(length(kegg[[i]]) > 0){
    for(j in 1:length(kegg[[i]])){
      if(kegg[[i]][j] == KEGGterm){
        genes[m] = names(kegg[i])
        m = m + 1
      }
    }
  }
}

table(genes%in%sig.genes)
sort(genes[genes%in%sig.genes])

##################################################################
## GO Analysis
library(goseq)
library(org.Bt.eg.db)
library(GO.db)
library(KEGG.db)
library(GOSemSim)
library(corrplot)

assayed.genes = array(total.genes)
de.genes = array(sig.genes)
gene.vector = as.integer(assayed.genes %in% de.genes)
names(gene.vector) = assayed.genes
table = table(gene.vector)

cat(paste("Significant Genes: ", table[2], 
          " Backgroung Genes:", table[1], "\n")) 

pwf = nullp(gene.vector, "bosTau4", "ensGene", plot.fit = FALSE)
head(pwf)

## GO Biological Processes
GO.hiper = goseq(pwf, "bosTau4", "ensGene", 
                 method = "Hypergeometric", test.cats="GO:BP",use_genes_without_cat = TRUE)

nGO.hiper = GO.hiper[(GO.hiper$numInCat <= 500 & GO.hiper$numInCat >= 5),] 
enriched.GO = nGO.hiper$category[nGO.hiper$over_represented_pvalue <= 0.05]
length(enriched.GO)

for (go in enriched.GO[1:length(enriched.GO)]) {	
  
  cat("#############################################\n")
  cat(paste("GOID: ", GOTERM[[go]]@GOID, "\n"))
  cat(paste("GO Term Name: ", GOTERM[[go]]@Term, "\n"))
  cat(paste("Definition: ", GOTERM[[go]]@Definition, "\n"))
  
  pvalue = nGO.hiper[nGO.hiper$category == go,]$over_represented_pvalue
  pvalue = round(pvalue,5)
  nGenes = nGO.hiper[nGO.hiper$category == go,]$numInCat
  nDEG = nGO.hiper[nGO.hiper$category == go,]$numDEInCat
  cat(paste("Total Genes:", nGenes, " Number DEG:", nDEG, " P-value:", pvalue))
  cat("\n")
}

## GO Molecular Function
GO.hiper = goseq(pwf, "bosTau4", "ensGene", 
                 method = "Hypergeometric", test.cats="GO:MF",use_genes_without_cat = TRUE)

nGO.hiper = GO.hiper[(GO.hiper$numInCat <= 500 & GO.hiper$numInCat >= 5),] 
enriched.GO = nGO.hiper$category[nGO.hiper$over_represented_pvalue <= 0.05]
length(enriched.GO)

for (go in enriched.GO[1:length(enriched.GO)]) {	
  
  cat("#############################################\n")
  cat(paste("GOID: ", GOTERM[[go]]@GOID, "\n"))
  cat(paste("GO Term Name: ", GOTERM[[go]]@Term, "\n"))
  cat(paste("Definition: ", GOTERM[[go]]@Definition, "\n"))
  
  pvalue = nGO.hiper[nGO.hiper$category == go,]$over_represented_pvalue
  pvalue = round(pvalue,5)
  nGenes = nGO.hiper[nGO.hiper$category == go,]$numInCat
  nDEG = nGO.hiper[nGO.hiper$category == go,]$numDEInCat
  cat(paste("Total Genes:", nGenes, " Number DEG:", nDEG, " P-value:", pvalue))
  cat("\n")
}

## GO Cellular Components ##
GO.hiper = goseq(pwf, "bosTau4", "ensGene", 
                 method = "Hypergeometric", test.cats="GO:CC",use_genes_without_cat = TRUE)

nGO.hiper = GO.hiper[(GO.hiper$numInCat <= 500 & GO.hiper$numInCat >= 5),] 
enriched.GO = nGO.hiper$category[nGO.hiper$over_represented_pvalue <= 0.05]
length(enriched.GO)

for (go in enriched.GO[1:length(enriched.GO)]) {	
  
  cat("#############################################\n")
  cat(paste("GOID: ", GOTERM[[go]]@GOID, "\n"))
  cat(paste("GO Term Name: ", GOTERM[[go]]@Term, "\n"))
  cat(paste("Definition: ", GOTERM[[go]]@Definition, "\n"))
  
  pvalue = nGO.hiper[nGO.hiper$category == go,]$over_represented_pvalue
  pvalue = round(pvalue,5)
  nGenes = nGO.hiper[nGO.hiper$category == go,]$numInCat
  nDEG = nGO.hiper[nGO.hiper$category == go,]$numDEInCat
  cat(paste("Total Genes:", nGenes, " Number DEG:", nDEG, " P-value:", pvalue))
  cat("\n")
}

## Revealing Significant Genes in GO Terms

GOterm = "GO:0034332"
go = getgo(names(gene.vector), "bosTau4", "ensGene")
genes = as.character(); m = 1

for(i in 1:length(go)){
  if(length(go[[i]]) > 0){
    for(j in 1:length(go[[i]])){
      if(go[[i]][j] == GOterm){
        genes[m] = names(go[i])
        m = m + 1
      }
    }
  }
}

table(genes %in% sig.genes)
sort(genes[genes %in% sig.genes])

# Revealing Significant genes in GO TERM

GOterm = "GO:0006886"
go = getgo(names(gene.vector), "bosTau4", "ensGene")
genes = as.character(); m = 1

for(i in 1:length(go)){
  if(length(go[[i]]) > 0){
    for(j in 1:length(go[[i]])){
      if(go[[i]][j] == GOterm){
        genes[m] = names(go[i])
        m = m + 1
      }
    }
  }
}

table(genes%in%sig.genes)
sort(genes[genes%in%sig.genes])

### MeSH Analysis ####


################
## MeSH Analysis

source("https://bioconductor.org/biocLite.R")

biocLite("org.Bt.eg.db")
biocLite("meshr")
biocLite("MeSH.db")
biocLite("MeSH.Bta.eg.db")

library(org.Bt.eg.db)
library(meshr) # They use this program in the package
library(MeSH.db)
library(MeSH.Bta.eg.db)

key.symbol = keys(org.Bt.eg.db,  keytype = c("SYMBOL"))
entrezUniverse = select(org.Bt.eg.db, as.character(key.symbol), 
                        columns = c("ENTREZID", "ENSEMBL"),keytype = "SYMBOL") # They use this program in the file
entrezUniverse2 <- entrezUniverse[!duplicated(entrezUniverse[,2]),]
entrezUniverse3 <- entrezUniverse2[!duplicated(entrezUniverse2[,1]),]

## FULL GENES
genes.back = data.frame(total.genes)
colnames(genes.back) <- "ENSEMBL"
geneID.back <- merge(genes.back, entrezUniverse3, by ="ENSEMBL")
geneID2.back <- geneID.back[ !duplicated(geneID.back[,2]),]

## SIGNIFICANT GENES
genes.sig = data.frame(sig.genes)
colnames(genes.sig) <- "ENSEMBL"
geneID.sig <- merge(genes.sig, entrezUniverse3, by ="ENSEMBL")
geneID2.sig <- geneID.sig[ !duplicated(geneID.sig[,2]),]

## Total Genes

ns = length(geneID2.sig[,1])
nt = length(geneID2.back[,1])
cat(paste("Significant Genes:", ns, " and Backgroung Genes:", nt - ns), "\n")

### MeSH Phenomena and Processes

meshParams <- new("MeSHHyperGParams", geneIds = geneID2.sig[,3], 
                  universeGeneIds = geneID2.back[,3], 
                  annotation = "MeSH.Bta.eg.db", category = "G", database = "gene2pubmed", 
                  pvalueCutoff = 0.05, pAdjust = "none")
meshR <- meshHyperGTest(meshParams)
out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM, 
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH Term ID", "MeSH Term Name", 
                  "Total Genes", "DE Genes", "P-value")
print(unique(out), row.names = F)


### MeSH Chemicals and Drugs

meshParams <- new("MeSHHyperGParams", geneIds = geneID2.sig[,3], 
                  universeGeneIds = geneID2.back[,3], 
                  annotation = "MeSH.Bta.eg.db", category = "D", database = "gene2pubmed", 
                  pvalueCutoff = 0.05, pAdjust = "none")
meshR <- meshHyperGTest(meshParams)
out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM, 
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH Term ID", "MeSH Term Name", 
                  "NT.Genes", "DE Genes", "P-value")
print(unique(out), row.names = F)

### MeSH Diseases

meshParams <- new("MeSHHyperGParams", geneIds = geneID2.sig[,3], 
                  universeGeneIds = geneID2.back[,3], 
                  annotation = "MeSH.Bta.eg.db", category = "A", database = "gene2pubmed", 
                  pvalueCutoff = 0.05, pAdjust = "none")
meshR <- meshHyperGTest(meshParams)
out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM, 
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH Term ID", "MeSH Term Name", 
                  "Total Genes", "DE Genes", "P-value")
print(unique(out), row.names = F)

