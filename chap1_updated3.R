library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(picante); packageVersion("picante")
library(reshape2)
library(seqTools)
library(adespatial)
library(ade4)
library(adegraphics)
library(spdep)
library(maptools)
library(akima)
library(gridExtra)
library(ggpubr)


######Community tables preparation
##INPUTS
meta <- read.table(file="organized2018_extraction_metadata.txt", sep="\t", header=TRUE, row.names="ID",stringsAsFactors = FALSE)
taxo <- read.csv(file="taxonomy_DADA2.csv", header=TRUE, row.names=1)
comm <- read.csv(file="community_DADA2.csv", header=TRUE, row.names=1)

##PREPROCESSING OF THE TABLES
par(mfrow=c(1,1))
all.comm <- comm
all.meta <- meta
all.taxo <- taxo

#replace all NAs by unassigned if we want the following part to work
for(i in colnames(all.taxo)){
  levels(all.taxo[,i])<-c(levels(all.taxo[,i]), "unassigned", NA)
  for(j in rownames(all.taxo)){
    if(is.na(all.taxo[j,i])){
      all.taxo[j,i] <- "unassigned"
    }
  }
}
all.trash <- all.taxo[c(all.taxo$Kingdom == "unassigned" |all.taxo$Kingdom == "Eukaryota" | all.taxo$Family == "Mitochondria"| all.taxo$Class == "Chloroplast"  ),]
all.taxo <- subset(all.taxo, !(rownames(all.taxo) %in% rownames(all.trash)))

#match the metadata and community tables
all.comm <- all.comm[rownames(all.meta),]
all.comm <- all.comm[,rownames(all.taxo)]

all.taxo <- data.frame(lapply(all.taxo, as.character), stringsAsFactors=FALSE)
rownames(all.taxo) <- colnames(all.comm)
#replace unassigned by the lowest  assigned taxonomic level
for(i in c(1:length(colnames(all.taxo)))){
  for(j in rownames(all.taxo)){
    if(all.taxo[j,i]=="unassigned"){
      if(substr(all.taxo[j,i-1], start=1,stop=2)=="UA"){
        all.taxo[j,i] <- all.taxo[j,i-1]
      }
      else{
        all.taxo[j,i] <- paste0("UA_",all.taxo[j,i-1])
      }
    }
  }
}


#change the names of our ASVs to something more practical
n <- 0
sequences <- character()
for(i in colnames(all.comm)){
  n <- n +1
  sequences[n] <- colnames(all.comm)[n]
  colnames(all.comm)[n] <- paste("ASV",n, sep="")
  names(sequences)[n] <- colnames(all.comm)[n]
}
rm(n,i,j)
rownames(all.taxo) <- colnames(all.comm)

#change headers in the metadata
n <- 0
for(colname in colnames(all.meta)){
  n <- n + 1
  if(colname == "Dev._Stage"){
    colnames(all.meta)[n] <-"sample_type"
  }
  else if(colname == "Tree_spec."){
    colnames(all.meta)[n] <-"tree_species"
  }
  if(colname == "Defoliation_level"){
    colnames(all.meta)[n] <-"host_tree_damage"
  }
}
rm(n)

#remove the controls
all.meta <- all.meta[!is.na(all.meta$sample_type),]


#add more metadata to our sites
sitesSpecs <- data.frame(read.csv(file="site_specs_2.csv", header=TRUE, row.names="ID1", stringsAsFactors = FALSE))
tempMeta <- data.frame()
n <- 0
for (rowname in rownames(all.meta)){
  for (siteName in rownames(sitesSpecs)){
    if (siteName == all.meta[rowname,"Site"]){
      tempMeta <- rbind(tempMeta,sitesSpecs[siteName,c(5:11)])
      n <- n + 1
      rownames(tempMeta)[n] <- rowname
    }
  }
}
all.meta <- cbind(all.meta[rownames(tempMeta),], tempMeta)
rm(n,sitesSpecs,tempMeta,all.trash, colname, rowname, siteName)
all.meta$site_specie <- paste(all.meta$Site, all.meta$tree_species)

#make sure all the numbers are numeric
list_of_numerics <- c("host_tree_damage","Latitude","Longitude","Elevation","Mean_temperature","Total_precipitations","Solar_radiation","Reference_evaporation","Climatic_moisture_deficit","Relative_humidity")
for (number in list_of_numerics){
  all.meta[,number] <- as.numeric(all.meta[,number])
}
rm(number,list_of_numerics)

#keep only larvae with all the metadata
all.meta <- all.meta[all.meta$sample_type == "L6"&!is.na(all.meta$host_tree_damage),]
all.comm <- all.comm[rownames(all.meta),]

#add a region factor
all.meta$region <- NA 
for (i in c(1:length(rownames(all.meta)))){
  x <- substr(as.character(all.meta[i,"Site"]),1,1)
  #print(x)
  if(x=="T"|x=="E"){
    all.meta[i,"region"] <- "Côte-Nord"
  }
  else if(x=="Q"){
    all.meta[i,"region"] <- "Gaspésie"
  }
  else if(x=="N"){
    all.meta[i,"region"] <- "Labrador"
  }
}

##rarefy our dataset
quantile(apply(all.comm, 1,sum))
all.ab.comm <- rrarefy(all.comm, 4000)
all.ab.comm <- all.ab.comm[apply(all.ab.comm,1,sum)>= 4000,]
rarecurve(all.ab.comm, step=100, label=FALSE)
all.ab.meta <- all.meta[rownames(all.ab.comm),]
all.ab.taxo <- all.taxo[colnames(all.ab.comm),]

##MEMS
coords <- cbind(as.numeric(unique(all.ab.meta$Latitude)), as.numeric(unique(all.ab.meta$Longitude)))
rownames(coords) <- unique(all.ab.meta$Site)
knn1 <- knearneigh(coords, k = 1)
nbknn1 <- knn2nb(knn1, sym = TRUE)
knn2 <- knearneigh(coords, k = 2)
nbknn2 <- knn2nb(knn2, sym = TRUE)
g1 <- s.label(coords, nb = nbknn1, pnb.edge.col = "red", main = "Nearest neighbors (k=1)", plot = FALSE)
g2 <- s.label(coords, nb = nbknn2, pnb.edge.col = "red", main = "Nearest neighbors (k=2)", plot = FALSE)
cbindADEg(g1, g2, plot = TRUE)
n.comp.nb(nbknn1)
nbgab <- graph2nb(gabrielneigh(coords), sym = TRUE)
g2 <- s.label(coords, nb = nbgab, pnb.edge.col = "red", main = "Gabriel", plot = FALSE)
distgab <- nbdists(nbgab, coords)
fdist <- lapply(distgab, function(x) 1 - x/max(dist(coords)))
listwgab <- nb2listw(nbgab, glist = fdist, style = "B")
mem.gab <- mem(listwgab)
barplot(attr(mem.gab, "values"), 
        main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7)
moranI <- moran.randtest(mem.gab, listwgab, 99)
signi <- which(moranI$pvalue < 0.05)
plot(mem.gab[,signi])
plot(mem.gab, SpORcoords = coords)
plot(mem.gab[,signi])
ADEgS(list(g2))
plot(mem.gab[,signi], SpORcoords = coords)
mem.signi <- mem.gab[,signi]
rownames(mem.signi) <- rownames(coords)

rm(coords,knn1,nbknn1,knn2,nbknn2,g1,g2,nbgab,distgab,)
#add the MeMs to the metadata
tempMeta <- data.frame()
n <- 0
for (rowname in rownames(all.ab.meta)){
  for (siteName in rownames(mem.signi)){
    if (siteName == all.ab.meta[rowname,"Site"]){
      tempMeta <- rbind(tempMeta,mem.signi[siteName,])
      n <- n + 1
      rownames(tempMeta)[n] <- rowname
    }
  }
}
all.ab.meta <- cbind(all.ab.meta[rownames(tempMeta),], tempMeta)
rm(n,tempMeta)
MEMs <- colnames(mem.signi)

###LARVAE VS FOLIAGE
###LARVAE VS FOLIAGE
##INPUTS
#meta <- read.table(file="organized2018_extraction_metadata.txt", sep="\t", header=TRUE, row.names="ID",stringsAsFactors = FALSE)
#taxo <- read.csv(file="taxonomy_DADA2.csv", header=TRUE, row.names=1)
#comm <- read.csv(file="community_DADA2.csv", header=TRUE, row.names=1)

##PREPROCESSING OF THE TABLES
compare.comm <- comm
compare.meta <- meta
compare.taxo <- taxo

#replace all NAs by unassigned if we want the following part to work
for(i in colnames(compare.taxo)){
  levels(compare.taxo[,i])<-c(levels(compare.taxo[,i]), "unassigned", NA)
  for(j in rownames(compare.taxo)){
    if(is.na(compare.taxo[j,i])){
      compare.taxo[j,i] <- "unassigned"
    }
  }
}
compare.trash <- compare.taxo[c(compare.taxo$Kingdom == "unassigned" |compare.taxo$Kingdom == "Eukaryota" | compare.taxo$Family == "Mitochondria"| compare.taxo$Class == "Chloroplast"  ),]
compare.taxo <- subset(compare.taxo, !(rownames(compare.taxo) %in% rownames(compare.trash)))

#match the metadata and community tables
compare.comm <- compare.comm[rownames(compare.meta),]
compare.comm <- compare.comm[,rownames(compare.taxo)]


compare.taxo <- data.frame(lapply(compare.taxo, as.character), stringsAsFactors=FALSE)
rownames(compare.taxo) <- colnames(compare.comm)
#replace unassigned by the lowest  assigned taxonomic level
for(i in c(1:length(colnames(compare.taxo)))){
  for(j in rownames(compare.taxo)){
    if(compare.taxo[j,i]=="unassigned"){
      if(substr(compare.taxo[j,i-1], start=1,stop=2)=="UA"){
        compare.taxo[j,i] <- compare.taxo[j,i-1]
      }
      else{
        compare.taxo[j,i] <- paste0("UA_",compare.taxo[j,i-1])
      }
    }
  }
}
#change the names of our ASVs to something more practical
n <- 0
compare.sequences <- character()
for(i in colnames(compare.comm)){
  n <- n +1
  sequences[n] <- colnames(compare.comm)[n]
  colnames(compare.comm)[n] <- paste("ASV",n, sep="")
  names(sequences)[n] <- colnames(compare.comm)[n]
}
rm(n,i,j,compare.trash)
rownames(compare.taxo) <- colnames(compare.comm)

#change headers in the metadata
n <- 0
for(colname in colnames(compare.meta)){
  n <- n + 1
  if(colname == "Dev._Stage"){
    colnames(compare.meta)[n] <-"sample_type"
  }
  else if(colname == "Tree_spec."){
    colnames(compare.meta)[n] <-"tree_species"
  }
  if(colname == "Defoliation_level"){
    colnames(compare.meta)[n] <-"host_tree_damage"
  }
}
rm(n)

#remove the controls
compare.meta <- compare.meta[!is.na(compare.meta$sample_type),]
compare.meta <- compare.meta[!is.na(compare.meta$host_tree_damage),]
compare.comm <- compare.comm[rownames(compare.meta),]

#add more metadata to our sites
sitesSpecs <- data.frame(read.csv(file="site_specs_2.csv", header=TRUE, row.names="ID1", stringsAsFactors = FALSE))
tempMeta <- data.frame()
n <- 0
for (rowname in rownames(compare.meta)){
  for (siteName in rownames(sitesSpecs)){
    if (siteName == compare.meta[rowname,"Site"]){
      tempMeta <- rbind(tempMeta,sitesSpecs[siteName,c(5:11)])
      n <- n + 1
      rownames(tempMeta)[n] <- rowname
    }
  }
}
compare.meta <- cbind(compare.meta[rownames(tempMeta),], tempMeta)
rm(n,sitesSpecs,tempMeta, colname, rowname, siteName)
compare.meta$site_specie <- paste(compare.meta$Site, compare.meta$tree_species)

#add the MeMs to the metadata
tempMeta <- data.frame()
n <- 0
for (rowname in rownames(compare.meta)){
  for (siteName in rownames(mem.signi)){
    if (siteName == compare.meta[rowname,"Site"]){
      tempMeta <- rbind(tempMeta,mem.signi[siteName,])
      n <- n + 1
      rownames(tempMeta)[n] <- rowname
    }
  }
}
compare.meta <- cbind(compare.meta[rownames(tempMeta),], tempMeta)
rm(n,tempMeta)
MEMs <- colnames(mem.signi)

#make sure all the numbers are numeric
list_of_numerics <- c("host_tree_damage","Latitude","Longitude")
for (number in list_of_numerics){
  compare.meta[,number] <- as.numeric(compare.meta[,number])
}
rm(number,list_of_numerics)

tempMeta <- data.frame()
n <- 0


##rarefy our dataset
quantile(apply(compare.comm, 1,sum))
compare.abund.comm <- rrarefy(compare.comm, 4000)
compare.abund.comm <- compare.abund.comm[apply(compare.abund.comm,1,sum)>= 4000,]
rarecurve(compare.abund.comm, step=100, label=FALSE)
compare.abund.meta <- compare.meta[rownames(compare.abund.comm),]
compare.abund.taxo <- compare.taxo[colnames(compare.abund.comm),]
split.comm <- compare.abund.comm
split.taxo <- compare.abund.taxo
split.meta <- compare.abund.meta

#function to identify samples for which foliage and larvae are present
pairs <- table(compare.abund.meta$Tree_ID) == 2
compareList <- c()
for (element in names(pairs)){
  if (pairs[element] == TRUE){
    compareList <- c(compareList,element) 
  }
}
rm(element)
rm(pairs)
compareListOut <- c()
for (tree_ID in compareList){
  for (rowname in rownames(compare.abund.comm)){
    if (tree_ID == compare.abund.meta[rowname,"Tree_ID"]){
      compareListOut <- c(compareListOut,rowname)
    }
  }
}
rm(rowname,tree_ID,compareList)

#make tables
compare.abund.comm <- compare.abund.comm[compareListOut,]
compare.abund.meta <- compare.abund.meta[rownames(compare.abund.comm),]
compare.abund.taxo <- compare.abund.taxo[colnames(compare.abund.comm),]


















#########Community analysis
####all samples community analysis
##PERMANOVA
#MARKDOWN
adonis(all.ab.comm~ Site/tree_species,data=all.ab.meta)
adonis(all.ab.comm~ region/Site/tree_species,data=all.ab.meta)

### try aggregating genus for PERMANOVA
all.genus <- aggregate(t(all.ab.comm), by=list(all.ab.taxo$Genus), sum)
rownames(all.genus) <- all.genus[,1]
all.genus <- all.genus[,-1]
all.genus <- t(all.genus)
adonis(all.genus~ region/Site/tree_species,data=all.ab.meta)

### try aggregating family for PERMANOVA
all.family <- aggregate(t(all.ab.comm), by=list(all.ab.taxo$Family), sum)
rownames(all.family) <- all.family[,1]
all.family <- all.family[,-1]
all.family <- t(all.family)
adonis(all.family~ region/Site/tree_species,data=all.ab.meta)


##NMDS with all variables significatively correlated with ordination axes
#MARKDOWN
par(mfrow=c(1,1))
all.ab.mds <- metaMDS(decostand(all.ab.comm, method="hellinger"),try=200,autotransform=FALSE,k=3)
ordiplot(all.ab.mds, display="sites")
all.variables <- c("host_tree_damage","Latitude","Longitude","Elevation","Mean_temperature","Total_precipitations","Solar_radiation","Reference_evaporation","Climatic_moisture_deficit","Relative_humidity",MEMs)
plot(envfit(ord=all.ab.mds, env=all.ab.meta[,all.variables], na.rm =TRUE, labels=TRUE))
all.envfit <-envfit(ord=all.ab.mds, env=all.ab.meta[,c("host_tree_damage","Elevation","Mean_temperature","Total_precipitations","Solar_radiation","Reference_evaporation","Climatic_moisture_deficit","Relative_humidity", MEMs)], na.rm =TRUE, labels=TRUE)
all.explanatories <- data.frame(cbind(all.envfit$vectors$r,all.envfit$vectors$pvals))
colnames(all.explanatories) <- c("R2", "P_Value")
all.explanatories <- all.explanatories[order(all.explanatories$P_Value),]
all.envfit.sig <- all.explanatories[all.explanatories$P_Value<=0.05,]
NMDS.env <- all.envfit.sig

ordiplot(all.ab.mds, display="sites")
plot(envfit(ord=all.ab.mds, env=all.ab.meta[,rownames(all.envfit.sig)], na.rm =TRUE, labels=TRUE))
for (row in c(1:length(rownames(all.ab.meta)))){
  if (all.ab.meta$tree_species[row] == "Balsam_fir"){
    all.ab.meta$HST_color[row] <- "red"
  } 
  else if (all.ab.meta$tree_species[row] == "Spruce"){
    all.ab.meta$HST_color[row] <- "green"
  } 
}
orditorp(all.ab.mds,display="sites", pcol = all.ab.meta$HST_color,labels = FALSE, pch = 19)
plot(envfit(ord=all.ab.mds, env=all.ab.meta[,rownames(all.envfit.sig)], na.rm =TRUE, labels=TRUE))

library(ggordiplots)
par(mfrow=c(1,1))
all.ab.mds <- metaMDS(decostand(all.ab.comm, method="hellinger"),try=200,autotransform=FALSE,k=3)
my.plot <- gg_ordiplot(all.ab.mds, groups = all.ab.meta$tree_species, conf=0.68) 
ord.data <- my.plot$df_ord
ord.stress<- all.ab.mds$stress
ellipse.data <- my.plot$df_ellipse
en <- envfit(ord=all.ab.mds, env=all.ab.meta[,rownames(all.envfit.sig)], na.rm =TRUE, labels=TRUE)
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

rownames(en_coord_cont) <- gsub(rownames(en_coord_cont),pattern = "_",replacement = " ")

ord.data$region <- all.ab.meta[rownames(ord.data),"region"]
library(viridis)
library(ggrepel)
NMDS <- ggplot(data = ord.data, aes(x = x, y = y, color = Group, shape=region)) + geom_point(size = 4) + 
  xlab("NMDS 1") + ylab("NMDS 2") + labs(color="Larvae found on:",shape="Region",caption = paste0("NMDS stress=",as.character(round(ord.stress,5)))) + scale_color_manual(values=c(viridis(5)[2],viridis(5)[4]),labels=c("Balsam fir", "Black spruce"))
NMDS

NMDS.envfit <- NMDS + geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, shape=NULL), data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30", arrow = arrow(length = unit(0.5, "cm"))) + geom_text_repel(data = en_coord_cont, aes(x = NMDS1, y = NMDS2,shape=NULL), colour = "grey30", fontface = "bold", label = row.names(en_coord_cont),segment.color = 'grey50') 
NMDS.envfit

tiff("figure4.tiff", units="in", width=15, height=10, res=600)

NMDS.envfit

dev.off()

NMDS.ellipse <- NMDS + geom_polygon(aes(x = x, y = y, group = Group, shape=NULL, colour=Group), alpha=0, data=ellipse.data )
NMDS.ellipse

NMDS.nolegend <- NMDS + guides(shape=FALSE,color=FALSE)





####try with PCOA
library(ggordiplots)
par(mfrow=c(1,1))
all.ab.mds <- cmdscale(vegdist(decostand(all.ab.comm, method="hellinger")))
all.ab.mds.GOF <- cmdscale(vegdist(decostand(all.ab.comm, method="hellinger")),eig=TRUE)$GOF
all.ab.mds.eig <- cmdscale(vegdist(decostand(all.ab.comm, method="hellinger")),eig=TRUE)$eig
variance.explained <- (100*all.ab.mds.eig/sum(all.ab.mds.eig))[1:2]
variance.explained.x <- format(round(variance.explained[1], 2), nsmall = 2)
variance.explained.y <- format(round(variance.explained[2], 2), nsmall = 2)
all.envfit <-envfit(ord=all.ab.mds, env=all.ab.meta[,c("host_tree_damage","Elevation","Mean_temperature","Total_precipitations","Solar_radiation","Reference_evaporation","Climatic_moisture_deficit","Relative_humidity", MEMs)], na.rm =TRUE, labels=TRUE)
all.explanatories <- data.frame(cbind(all.envfit$vectors$r,all.envfit$vectors$pvals))
colnames(all.explanatories) <- c("R2", "P_Value")
all.explanatories <- all.explanatories[order(all.explanatories[,"P_Value"]),]
all.envfit.sig <- all.explanatories[all.explanatories$P_Value<=0.05,]
all.env.table <- as.data.frame(cbind(all.envfit$vectors$arrows,all.envfit$vectors$r,all.envfit$vectors$pvals))
colnames(all.env.table) <- c("Correlation with PCoA 1","Correlation with PCoA 2","R2", "P-Value")
all.env.table <-all.env.table[order(all.env.table[,"P-Value"]),]
write.csv(x = all.env.table,"supp2.csv")
PCOA.env <- all.envfit.sig
my.plot <- gg_ordiplot(all.ab.mds, groups = all.ab.meta$tree_species, conf=0.68) 
ord.data <- my.plot$df_ord
ellipse.data <- my.plot$df_ellipse
en <- envfit(ord=all.ab.mds, env=all.ab.meta[,rownames(all.envfit.sig)], na.rm =TRUE, labels=TRUE)
en_coord_cont = as.data.frame(scores(en, "vectors")) 

rownames(en_coord_cont) <- gsub(rownames(en_coord_cont),pattern = "_",replacement = " ")
rownames(en_coord_cont) <- gsub(rownames(en_coord_cont),pattern = "precipitations",replacement = "precipitation")

ord.data$region <- all.ab.meta[rownames(ord.data),"region"]
ord.data$Site <- all.ab.meta[rownames(ord.data),"Site"]
ord.data$ID <- all.ab.meta[rownames(ord.data),"Tree-ID"]
ord.data$date <- all.ab.meta[rownames(ord.data),"Date_of_extraction"]
library(viridis)
library(ggrepel)
PCOA <- ggplot(data = ord.data, aes(x = x, y = y, color = Group, shape=region)) + geom_point(size = 4,) + 
  xlab(paste0("PCOA 1 (",variance.explained.x,"%)")) + ylab(paste0("PCOA 2 (",variance.explained.y,"%)")) + labs(color="Larvae found on:",shape="Region") + scale_color_manual(values=c(viridis(5)[2],viridis(5)[4]),labels=c("Balsam fir", "Black spruce"))

PCOA 
#PCOA + geom_text(aes(label=Site),hjust=0, vjust=0)


PCOA.envfit <- PCOA + geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2, shape=NULL), data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30", arrow = arrow(length = unit(0.5, "cm"))) + geom_text_repel(data = en_coord_cont, aes(x = Dim1, y = Dim2,shape=NULL), colour = "grey30", fontface = "bold", label = row.names(en_coord_cont),segment.color = 'grey50') 
PCOA.envfit

tiff("PCOA4.tiff", units="in", width=8, height=5, res=320)

PCOA.envfit

dev.off()

PCOA.ellipse <- PCOA + geom_polygon(aes(x = x, y = y, group = Group, shape=NULL, colour=Group), alpha=0, data=ellipse.data )
PCOA.ellipse



PCOA.nolegend <- PCOA + guides(shape=FALSE,color=FALSE)

grid.arrange(NMDS.nolegend,PCOA.nolegend,ncol=2)


### try aggregating genus for PERMANOVA
all.genus <- aggregate(t(all.ab.comm), by=list(all.ab.taxo$Genus), sum)
rownames(all.genus) <- all.genus[,1]
all.genus <- all.genus[,-1]
all.genus <- t(all.genus)
adonis(all.genus~ Site/tree_species,data=all.ab.meta)

apply(all.genus,2,sum)
listofUA <- rownames(all.ab.taxo)[grepl("UA",all.ab.taxo$Genus)]
table(grepl("UA",all.ab.taxo$Genus))
sum(apply(all.ab.comm[,listofUA],2,"sum"))/sum(apply(all.ab.comm,2,"sum"))

### try aggregating family for PERMANOVA
all.family <- aggregate(t(all.ab.comm), by=list(all.ab.taxo$Family), sum)
rownames(all.family) <- all.family[,1]
all.family <- all.family[,-1]
all.family <- t(all.family)
adonis(all.family~ Site/tree_species,data=all.ab.meta)

apply(all.family,2,sum)
listofUA <- rownames(all.ab.taxo)[grepl("UA",all.ab.taxo$Family)]
table(grepl("UA",all.ab.taxo$Family))
sum(apply(all.ab.comm[,listofUA],2,"sum"))/sum(apply(all.ab.comm,2,"sum"))







####scatterplot of pairwise geographic distance (x axis) versus pairwise beta-diversity (y – axis)
#matrix of betadiversity
beta.dist <- vegdist(all.ab.comm)
#install.packages("geodist")
library(geodist)
coords <- all.ab.meta[,c("Latitude","Longitude")]
geo.dist <- dist(coords)
table(names(geo.dist)==names(beta.dist))

mantel(geo.dist,beta.dist)
mantel(beta.dist,geo.dist)
both.dist <- cbind(geo.dist,beta.dist)

tiff("scatterplot.tiff", units="in", width=15, height=10, res=300)

ggplot(data = as.data.frame(both.dist), aes(x = geo.dist.km, y = beta.dist)) + geom_point(size = 2) + xlab("Geographic distance (km)") + ylab("Bray-Curtis dissimilarity") + theme(axis.text = element_text(size=14), axis.title = element_text(size=16) ) + xlim(c(0,1200))

dev.off()


####scatterplot of pairwise geographic distance (x axis) versus pairwise beta-diversity (y – axis) while aggregated to the site
table(rownames(all.ab.comm)==rownames(all.ab.meta))
siteagg.comm <- aggregate(all.ab.comm, by=list(all.ab.meta$Site), mean)
rownames(siteagg.comm) <- siteagg.comm$Group.1
siteagg.comm <- siteagg.comm[,-1]
#matrix of betadiversity
beta.dist <- vegdist(siteagg.comm)
#install.packages("geodist")
library(geodist)
siteagg.meta <- aggregate(all.ab.meta, by=list(all.ab.meta$Site), mean)
rownames(siteagg.meta) <- siteagg.meta$Group.1
siteagg.meta <- siteagg.meta[,-1]
siteagg.meta <- siteagg.meta[rownames(siteagg.comm),]

coords <- siteagg.meta[,c("Latitude","Longitude")]
geo.dist <- dist(coords)
geo.dist.km <- geo.dist*111
table(names(geo.dist.km)==names(beta.dist))

both.dist <- cbind(geo.dist.km,beta.dist)

tiff("scatterplot_agg.tiff", units="in", width=16, height=8, res=300)

ggplot(data = as.data.frame(both.dist), aes(x = geo.dist.km, y = beta.dist)) + geom_point(size = 2) + xlab("Geographic distance (km)") + ylab("Bray-Curtis dissimilarity") + theme(axis.text = element_text(size=14), axis.title = element_text(size=16) ) + xlim(c(0,1200))

dev.off()










###Ancom-BC analysis
library(ANCOMBC)

all.meta$tree_species <- factor(all.meta$tree_species)
all.phylo<-phyloseq(otu_table(all.comm,taxa_are_rows = FALSE), tax_table(as.matrix(all.taxo)), sample_data(all.meta))

ancom.prep <- function(phylo,factor){
  
  comm <- otu_table(phylo)
  taxo <- data.frame(tax_table(phylo))
  meta <- data.frame(sample_data(phylo))
  
  
  vector.factor <- factor(meta[,factor])
  groups <- levels(vector.factor)
  #set up the prevalence percentage so ASVs in at least 2 samples
  perc <- 1-(2/dim(comm)[1])
  ancom <- ancombc(phylo,formula = factor,group=factor,zero_cut = perc,struc_zero=TRUE)
  ancom.data <- as.data.frame(ancom$res)
  colnames(ancom.data) <- names(ancom$res)
  #add the taxonomy
  ancom.data <- cbind(ancom.data,taxo[rownames(ancom.data),])
  #order with adj p-values
  ancom.data <- ancom.data[order(ancom.data$q_val),]
  #keep track which group has which beta values
  under0 <- groups[1]
  over0 <- groups[2]
  #keep track which is the positively-associated group
  for (i in c(1:length(rownames(ancom.data)))) {
    if (ancom.data$beta[i]<0){
      ancom.data$Associated_with[i] <- under0
    }
    else if (ancom.data$beta[i]>0){
      ancom.data$Associated_with[i] <- over0
    }
  }
  
  #check which ASVs are common to both chapters
  list1 <- colnames(comm)[apply(comm[meta[,factor]==under0,],2,sum) > 0]
  list2 <- colnames(comm)[apply(comm[meta[,factor]==over0,],2,sum) > 0]
  list.total <- colnames(comm)[apply(comm,2,sum) > 0]
  list.common <- list1[list1%in%list2]
  table(list.total%in%list.common)
  
  #add if they are common to both groups or not
  for (rowname in rownames(ancom.data)){
    if (rowname %in% list.common){
      ancom.data[rowname,"Common"] <- "Common"
    }
    else {
      ancom.data[rowname,"Common"] <- "Unique"
    }
  }
  
  ancom.sigtab <- ancom.data[ancom.data$q_val<0.05,]
  library("ggplot2")
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  
  # Phylum order
  x = tapply(ancom.sigtab$beta, ancom.sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  ancom.sigtab$Phylum = factor(as.character(ancom.sigtab$Phylum), levels=names(x))
  # Genus order
  x = tapply(ancom.sigtab$beta, ancom.sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  ancom.sigtab$Genus = factor(as.character(ancom.sigtab$Genus), levels=names(x))
  ancom.ggplot <- ggplot(ancom.sigtab, aes(x=Genus, y=beta, color=Phylum, shape=Associated_with)) + theme_grey()+ geom_point(size=4) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),axis.title.y.left = element_text(hjust = 0.5)) + geom_hline(yintercept=0, linetype="dashed", 
                                                                                                                                    color = "black", size=0.5) + ylab("ANCOM-BC log-linear coefficient") 
  return(list(ancom.sigtab,ancom.ggplot))
}


all.ancom <- ancom.prep(all.phylo,"tree_species")
all.ancom.ggplot <- all.ancom[[2]]
all.ancom.table <- all.ancom[[1]]

all.ancom.ggplot.complete <- all.ancom.ggplot + scale_shape_manual(labels= c("Balsam fir","Black spruce"),values = c(16,18)) + labs(shape="Associated with larvae feeding on:")
all.ancom.ggplot.complete

all.ancom.table <- all.ancom.table[order(decreasing = TRUE,abs(all.ancom.table$beta)),]

tiff("all.ancom.tiff", units="in", width=10, height=6, res=300)

all.ancom.ggplot.complete

dev.off()

write.csv(all.ancom.table, "all.ancom.csv")


#make a figure of all the most differentially abundant genera
all.ancom.betas <- all.ancom.table[,c("beta","Genus")]
all.ancom.betas$beta <- abs(all.ancom.betas$beta)
all.ancom.betas <- all.ancom.betas[order(all.ancom.betas$beta,decreasing = TRUE),]
all.ancom.top12 <- all.ancom.betas[c(1:12),]
top12.comm <- data.frame(all.comm[,rownames(all.ancom.top12)])
colnames(top12.comm) <- all.ancom.top12$Genus
top12.comm$group <- all.meta[rownames(top12.comm),"tree_species"]
top12.melt <- melt(top12.comm,id.vars="group")
all.top12 <- ggplot(data=top12.melt,aes(x=variable,y=value,fill=group,color=group)) + geom_jitter(position=position_jitterdodge(dodge.width=0.9)) + scale_color_manual(labels= c("Balsam fir","Black spruce"),values = c(viridis(5)[2],viridis(5)[4])) + labs(color="Associated with larvae feeding on:") + guides(fill=FALSE) + ylab("Abundance") + xlab("Bacterial genus")

tiff("all.top12.tiff", units="in", width=15, height=6, res=300)

all.top12

dev.off()



















##Comparison between host species
#split between host species
fir.comm <- all.ab.comm[all.ab.meta$tree_species == "Balsam_fir",]
spruce.comm <- all.ab.comm[all.ab.meta$tree_species == "Spruce",]
fir.meta <- all.ab.meta[all.ab.meta$tree_species == "Balsam_fir",]
spruce.meta <- all.ab.meta[all.ab.meta$tree_species == "Spruce",]

#ordinations to compare spruce and fir
fir.meta <- fir.meta[order(fir.meta$Site),]
spruce.meta <- spruce.meta[order(spruce.meta$Site),]

#make sure both spruce and fir samples are comparable
fir.comp.meta <- data.frame()
for (site in unique(spruce.meta$Site)){
  for(n in 1:length(spruce.meta$Site[spruce.meta$Site == site])){
    if (n <=   length(fir.meta$Site[fir.meta$Site == site])){
      fir.comp.meta <- rbind(fir.comp.meta, fir.meta[fir.meta$Site == site,][n,])
    }
  }
}
spruce.comp.meta <- data.frame()
for (site in unique(fir.comp.meta$Site)){
  for(n in 1:length(fir.meta$Site[fir.meta$Site == site])){
    if (n <=   length(spruce.meta$Site[spruce.meta$Site == site])){
      spruce.comp.meta <- rbind(spruce.comp.meta, spruce.meta[spruce.meta$Site == site,][n,])
    }
  }
}
rm(site,n)

spruce.comp.comm <- spruce.comm[rownames(spruce.comp.meta),]
fir.comp.comm <- fir.comm[rownames(fir.comp.meta),]
fir.comp.meta <- fir.meta[rownames(fir.comp.comm),]
spruce.comp.meta <- spruce.meta[rownames(spruce.comp.comm),]

#ordination of larvae on both host-tree species

#fir
fir.mds <- cmdscale(vegdist(decostand(fir.comp.comm, method="hellinger")))
fir.mds.eig <- cmdscale(vegdist(decostand(fir.comp.comm, method="hellinger")),eig=TRUE)$eig
fir.variance.explained <- (100*fir.mds.eig/sum(fir.mds.eig))[1:2]
fir.variance.explained.x <- format(round(fir.variance.explained[1], 2), nsmall = 2)
fir.variance.explained.y <- format(round(fir.variance.explained[2], 2), nsmall = 2)
fir.envfit <-envfit(ord=fir.mds, env=fir.comp.meta[,c("host_tree_damage","Elevation","Mean_temperature","Total_precipitations","Solar_radiation","Reference_evaporation","Climatic_moisture_deficit","Relative_humidity", MEMs)], na.rm =TRUE, labels=TRUE)
fir.explanatories <- data.frame(cbind(fir.envfit$vectors$r,fir.envfit$vectors$pvals))
colnames(fir.explanatories) <- c("R2", "P_Value")
fir.explanatories <- fir.explanatories[order(fir.explanatories[,"P_Value"]),]
fir.envfit.sig <- fir.explanatories[fir.explanatories$P_Value<0.05,]
fir.env.table <- as.data.frame(cbind(fir.envfit$vectors$arrows,fir.envfit$vectors$r,fir.envfit$vectors$pvals))
colnames(fir.env.table) <- c("Correlation with PCoA 1","Correlation with PCoA 2","R2", "P-Value")
fir.env.table <-fir.env.table[order(fir.env.table[,"P-Value"]),]
write.csv(x = fir.env.table,"fir.supp2.csv")
PCOA.env <- fir.envfit.sig
my.plot <- gg_ordiplot(fir.mds, groups = fir.comp.meta$tree_species, conf=0.68) 
fir.ord.data <- my.plot$df_ord
fir.ord.data$facet <- "Balsam fir"
fir.ellipse.data <- my.plot$df_ellipse
fir.ellipse.data$facet <- "Balsam fir"
en <- envfit(ord=fir.mds, env=fir.comp.meta[,rownames(fir.envfit.sig)], na.rm =TRUE, labels=TRUE)
fir.en_coord_cont = as.data.frame(scores(en, "vectors")) 
fir.en_coord_cont$facet <- "Balsam fir"
rownames(en_coord_cont) <- gsub(rownames(en_coord_cont),pattern = "_",replacement = " ")
fir.ord.data$region <- fir.comp.meta[rownames(fir.ord.data),"region"]
rownames(fir.en_coord_cont) <- gsub(rownames(fir.en_coord_cont),pattern = "_",replacement = " ")
rownames(fir.en_coord_cont) <- gsub(rownames(fir.en_coord_cont),pattern = "precipitations",replacement = "precipitation")


#spruce
spruce.mds <- cmdscale(vegdist(decostand(spruce.comp.comm, method="hellinger")))
spruce.mds.eig <- cmdscale(vegdist(decostand(spruce.comp.comm, method="hellinger")),eig=TRUE)$eig
spruce.variance.explained <- (100*spruce.mds.eig/sum(spruce.mds.eig))[1:2]
spruce.variance.explained.x <- format(round(spruce.variance.explained[1], 2), nsmall = 2)
spruce.variance.explained.y <- format(round(spruce.variance.explained[2], 2), nsmall = 2)
spruce.envfit <-envfit(ord=spruce.mds, env=spruce.comp.meta[,c("host_tree_damage","Elevation","Mean_temperature","Total_precipitations","Solar_radiation","Reference_evaporation","Climatic_moisture_deficit","Relative_humidity", MEMs)], na.rm =TRUE, labels=TRUE)
spruce.explanatories <- data.frame(cbind(spruce.envfit$vectors$r,spruce.envfit$vectors$pvals))
colnames(spruce.explanatories) <- c("R2", "P_Value")
spruce.explanatories <- spruce.explanatories[order(spruce.explanatories[,"P_Value"]),]
spruce.envfit.sig <- spruce.explanatories[spruce.explanatories$P_Value<0.05,]
spruce.env.table <- as.data.frame(cbind(spruce.envfit$vectors$arrows,spruce.envfit$vectors$r,spruce.envfit$vectors$pvals))
colnames(spruce.env.table) <- c("Correlation with PCoA 1","Correlation with PCoA 2","R2", "P-Value")
spruce.env.table <-spruce.env.table[order(spruce.env.table[,"P-Value"]),]
write.csv(x = spruce.env.table,"spruce.supp2.csv")
PCOA.env <- spruce.envfit.sig
my.plot <- gg_ordiplot(spruce.mds, groups = spruce.comp.meta$tree_species, conf=0.68) 
spruce.ord.data <- my.plot$df_ord
spruce.ord.data$facet <- "Black spruce"
spruce.ellipse.data <- my.plot$df_ellipse
spruce.ellipse.data$facet <- "Black spruce"
en <- envfit(ord=spruce.mds, env=spruce.comp.meta[,rownames(spruce.envfit.sig)], na.rm =TRUE, labels=TRUE)
spruce.en_coord_cont = as.data.frame(scores(en, "vectors")) 
spruce.en_coord_cont$facet <- "Black spruce"
rownames(en_coord_cont) <- gsub(rownames(en_coord_cont),pattern = "_",replacement = " ")
spruce.ord.data$region <- spruce.comp.meta[rownames(spruce.ord.data),"region"]
rownames(spruce.en_coord_cont) <- gsub(rownames(spruce.en_coord_cont),pattern = "_",replacement = " ")
rownames(spruce.en_coord_cont) <- gsub(rownames(spruce.en_coord_cont),pattern = "precipitations",replacement = "precipitation")


HTS.ord.data <- rbind(fir.ord.data,spruce.ord.data)
HTS.ellipse.data <- rbind(fir.ellipse.data,spruce.ellipse.data)
HTS.en_coord_cont <- rbind(fir.en_coord_cont,spruce.en_coord_cont)
rownames(HTS.en_coord_cont) <- gsub(rownames(HTS.en_coord_cont),pattern = "_",replacement = " ")
rownames(HTS.en_coord_cont) <- gsub(rownames(HTS.en_coord_cont),pattern = "precipitations",replacement = "precipitation")

library(viridis)
library(ggrepel)


fir.PCOA <- ggplot(data = fir.ord.data, aes(x = x, y = y, color = Group, shape=region)) + geom_point(size = 4) + 
  xlab(paste0("PCOA 1 (",fir.variance.explained.x,"%)")) + ylab(paste0("PCOA 2 (",fir.variance.explained.y,"%)")) + labs(color="Larvae found on:",shape="Region",title = "Balsam fir") + scale_color_manual(values=c(viridis(5)[2]),labels=c("Balsam fir")) 
fir.PCOA 

fir.PCOA.envfit <- fir.PCOA + geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2, shape=NULL), data = fir.en_coord_cont, size =1, alpha = 0.5, colour = "grey30", arrow = arrow(length = unit(0.5, "cm"))) + geom_text_repel(data = fir.en_coord_cont, aes(x = Dim1, y = Dim2,shape=NULL), colour = "grey30", fontface = "bold", label = row.names(fir.en_coord_cont),segment.color = 'grey50') + theme(legend.position = "none")
fir.PCOA.envfit

spruce.PCOA <- ggplot(data = spruce.ord.data, aes(x = x, y = y, color = Group, shape=region)) + geom_point(size = 4) + 
  xlab(paste0("PCOA 1 (",spruce.variance.explained.x,"%)")) + ylab(paste0("PCOA 2 (",spruce.variance.explained.y,"%)")) + labs(color="Larvae found on:",shape="Region",title = "Black spruce") + scale_color_manual(values=viridis(5)[4],labels=c("Black spruce")) 
spruce.PCOA 

spruce.PCOA.envfit <- spruce.PCOA + geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2, shape=NULL), data = spruce.en_coord_cont, size =1, alpha = 0.5, colour = "grey30", arrow = arrow(length = unit(0.5, "cm"))) + geom_text_repel(data = spruce.en_coord_cont, aes(x = Dim1, y = Dim2,shape=NULL), colour = "grey30", fontface = "bold", label = row.names(spruce.en_coord_cont),segment.color = 'grey50') + theme(legend.position = "none")
spruce.PCOA.envfit



#facet way(only to keep the legend)
PCOA <- ggplot(data = HTS.ord.data, aes(x = x, y = y, color = Group, shape=region)) + geom_point(size = 4) + 
  xlab("PCOA 1") + ylab("PCOA 2") + labs(color="Larvae found on:",shape="Region") + scale_color_manual(values=c(viridis(5)[2],viridis(5)[4]),labels=c("Balsam fir", "Black spruce")) + facet_wrap(facets=~facet)
PCOA 
#PCOA + geom_text(aes(label=Site),hjust=0, vjust=0)
PCOA.envfit <- PCOA + geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2, shape=NULL), data = HTS.en_coord_cont, size =1, alpha = 0.5, colour = "grey30", arrow = arrow(length = unit(0.5, "cm"))) + geom_text_repel(data = HTS.en_coord_cont, aes(x = Dim1, y = Dim2,shape=NULL), colour = "grey30", fontface = "bold", label = c(row.names(HTS.en_coord_cont)[HTS.en_coord_cont$facet=="Balsam fir"],gsub(pattern = "1",replacement = "",row.names(HTS.en_coord_cont)[HTS.en_coord_cont$facet=="Black spruce"])),segment.color = 'grey50')+ facet_wrap(facets=~facet,) 
PCOA.envfit
leg <- get_legend(PCOA.envfit)

grid.arrange(fir.PCOA.envfit,spruce.PCOA.envfit,leg, ncol=3, widths=c(5,5,1.6))

tiff("HTS.PCOA4.tiff", units="in", width=10, height=6, res=360)

grid.arrange(fir.PCOA.envfit,spruce.PCOA.envfit,leg, ncol=3, widths=c(5,5,1.6))

dev.off()


#Procrustes test between fir and spruce
par(mfrow=c(1,1))
hosts.proc <- procrustes(fir.mds, spruce.mds)
hosts.test <- protest(fir.mds, spruce.mds)
plot(hosts.proc)
hosts.test

















#####Community analysis of Foliage vs larvae

#PERMANOVAS ON FOLIAGE AND LARVAE
adonis(compare.abund.comm[compare.abund.meta$sample_type == "L6",]~ Site/tree_species,data=compare.abund.meta[compare.abund.meta$sample_type == "L6",])
adonis(compare.abund.comm[compare.abund.meta$sample_type == "Foliage",]~ Site/tree_species,data=compare.abund.meta[compare.abund.meta$sample_type == "Foliage",])



#with pcoa
compare.larvae.mds <- cmdscale(vegdist(decostand(compare.abund.comm[compare.abund.meta$sample_type == "L6",], method="hellinger")))
compare.foliage.mds <- cmdscale(vegdist(decostand(compare.abund.comm[compare.abund.meta$sample_type == "Foliage",], method="hellinger")))
#procrustes of foliage vs larvae
#MARKDOWN
par(mfrow=c(1,1))
comp.proc <- procrustes(compare.larvae.mds, compare.foliage.mds)
comp.test <- protest(compare.larvae.mds, compare.foliage.mds)
plot(comp.proc)
comp.test


###diversity stats
foliage.comm <- split.comm[split.meta$sample_type=="Foliage",]
larvae.comm <- split.comm[split.meta$sample_type=="L6",]
foliage.comm <- foliage.comm[,apply(foliage.comm,2,sum)> 0]
larvae.comm <- larvae.comm[,apply(larvae.comm,2,sum)> 0]

#richness
mean(apply(foliage.comm[]>0,1,sum))
sd(apply(foliage.comm[]>0,1,sum))
mean(apply(larvae.comm[]>0,1,sum))
sd(apply(larvae.comm[]>0,1,sum))
rich.foliage <- data.frame(apply(foliage.comm[]>0,1,sum))
rich.foliage$type <- "foliage"
rownames(rich.foliage) <- paste(rownames(rich.foliage),rich.foliage$type )
colnames(rich.foliage)[1] <- "richness"
rich.larvae <- data.frame(apply(larvae.comm[]>0,1,sum))
rich.larvae$type <- "larvae"
rownames(rich.larvae) <- paste(rownames(rich.larvae),rich.larvae$type )
colnames(rich.larvae)[1] <- "richness"
rich.data <- rbind(rich.larvae,rich.foliage)
rich.aov <- aov(richness~type,data=rich.data)
summary(rich.aov)

#adapt taxonomic files
foliage.taxo <- split.taxo[colnames(foliage.comm),]
larvae.taxo <- split.taxo[colnames(larvae.comm),]

#most abundant genera
foliage.taxocomm.genus <- aggregate(t(foliage.comm), by=list(foliage.taxo$Genus), sum)
rownames(foliage.taxocomm.genus)<-foliage.taxocomm.genus[,1]
foliage.taxocomm.genus <- foliage.taxocomm.genus[,-1]
foliage.taxocomm.genus$Total <- apply(foliage.taxocomm.genus,1,sum)
foliage.taxocomm.genus <- foliage.taxocomm.genus[order(foliage.taxocomm.genus$Total),]

larvae.taxocomm.genus <- aggregate(t(larvae.comm), by=list(larvae.taxo$Genus), sum)
rownames(larvae.taxocomm.genus)<-larvae.taxocomm.genus[,1]
larvae.taxocomm.genus <- larvae.taxocomm.genus[,-1]
larvae.taxocomm.genus$Total <- apply(larvae.taxocomm.genus,1,sum)
larvae.taxocomm.genus <- larvae.taxocomm.genus[order(larvae.taxocomm.genus$Total),]

#most abundant phyla
foliage.taxocomm.phylum <- aggregate(t(foliage.comm), by=list(foliage.taxo$Phylum), sum)
rownames(foliage.taxocomm.phylum)<-foliage.taxocomm.phylum[,1]
foliage.taxocomm.phylum <- foliage.taxocomm.phylum[,-1]
foliage.taxocomm.phylum$Total <- apply(foliage.taxocomm.phylum,1,sum)
foliage.taxocomm.phylum <- foliage.taxocomm.phylum[order(foliage.taxocomm.phylum$Total),]
foliage.taxocomm.phylum$Prop <- round(foliage.taxocomm.phylum$Total/sum(foliage.taxocomm.phylum$Total)*100,digits = 1)
foliage.taxocomm.phylum$Std <- round(apply(foliage.taxocomm.phylum[,c(1:22)],1,sd)/sum(foliage.taxocomm.phylum$Total)*100,digits = 1)

larvae.taxocomm.phylum <- aggregate(t(larvae.comm), by=list(larvae.taxo$Phylum), sum)
rownames(larvae.taxocomm.phylum)<-larvae.taxocomm.phylum[,1]
larvae.taxocomm.phylum <- larvae.taxocomm.phylum[,-1]
larvae.taxocomm.phylum$Total <- apply(larvae.taxocomm.phylum,1,sum)
larvae.taxocomm.phylum <- larvae.taxocomm.phylum[order(larvae.taxocomm.phylum$Total),]
larvae.taxocomm.phylum$Prop <- round(larvae.taxocomm.phylum$Total/sum(larvae.taxocomm.phylum$Total)*100,digits = 1)
larvae.taxocomm.phylum$Std <- round(apply(larvae.taxocomm.phylum[,c(1:99)],1,sd)/sum(larvae.taxocomm.phylum$Total)*100,digits = 1)



## Steve's updated analyses - model simplification to reduce collinearity
## to reduce collinearity - check for each matrix separately?
## or subset matrix of all variables?
#community and metadata matrix preparation
larva.hell <- decostand(compare.abund.comm[compare.abund.meta$sample_type == "L6",],method = "hellinger")
larva.meta <- compare.abund.meta[rownames(larva.hell),]
#make tables
# env, spatial, MEM
environmentals <- scale(larva.meta[,c("host_tree_damage","Mean_temperature","Total_precipitations","Solar_radiation","Reference_evaporation","Climatic_moisture_deficit","Relative_humidity")])
#pairs(larva.meta[,c("host_tree_damage","Mean_temperature","Total_precipitations","Solar_radiation","Reference_evaporation","Climatic_moisture_deficit","Relative_humidity")])
environmentals <- scale(larva.meta[,c("host_tree_damage","Mean_temperature","Total_precipitations")])
spatials <- scale(larva.meta[,c("Latitude","Longitude")])
MEM <- larva.meta[,MEMs]
HST <- larva.meta[,"tree_species"]
HST.num <- c()
for(item in HST){
  if(item=="Spruce"){
    HST.num<-c(HST.num,0)
  }
  else{
    HST.num<-c(HST.num,1)
  }
}
HST <- larva.meta[,"tree_species"]
#diet
#PCOA
compare.foliage.mds <- cmdscale(vegdist(decostand(compare.abund.comm[compare.abund.meta$sample_type == "Foliage",], method="hellinger")))
foliage.scores <- scores(compare.foliage.mds)
#NMDS
#compare.foliage.mds <- metaMDS(decostand(compare.abund.comm[compare.abund.meta$sample_type == "Foliage",],method="hellinger"), k=5)
#foliage.scores <- scores(compare.foliage.mds)

#colnames(foliage.scores) <- sub(pattern = "NMDS",replacement = "foliage_NMDS", x=colnames(foliage.scores))
larva.meta <- cbind(larva.meta,foliage.scores)
diet <- larva.meta[,colnames(foliage.scores)]
# combine all variables
everyVariable <- cbind(environmentals,spatials,MEM,diet)
# CCA/cca model - can try both here
# create subset model of full model using ordistep
mod<- cca(larva.hell ~ ., everyVariable)
#mod0 <- cca(larva.hell ~ 1, everyVariable)
mod.step <- ordistep(mod)
#mod.stepR2 <- ordiR2step(mod0, mod)
plot(mod)
plot(mod.step)
#plot(mod.stepR2)
vif.cca(mod.step)
#vif.cca(mod.stepR2)
barplot(vif.cca(mod.step))
pairs(everyVariable[,names(vif.cca(mod.step))])
#subset tables individually
#environmentals
mod.env <- cca(larva.hell ~ ., data.frame(environmentals))
mod.env.step <- ordistep(mod.env)
vif.cca(mod.env)
vif.cca(mod.env.step)
environmentals.sub <- environmentals[,names(vif.cca(mod.env.step)),drop=FALSE]
#spatials
mod.spatial <- cca(larva.hell ~ ., data.frame(spatials))
mod.spatial.step <- ordistep(mod.spatial)
vif.cca(mod.spatial)
vif.cca(mod.spatial.step)
spatial.sub <- spatials[,names(vif.cca(mod.spatial)),drop=FALSE]
#MEM
mod.MEM <- cca(larva.hell ~ ., data.frame(MEM))
mod.MEM.step <- ordistep(mod.MEM)
vif.cca(mod.MEM)
vif.cca(mod.MEM.step)
MEM.sub <- MEM[,names(vif.cca(mod.MEM.step)),drop=FALSE]
#diet
mod.diet <- cca(larva.hell ~ ., data.frame(diet))
mod.diet.step <- ordistep(mod.diet)
vif.cca(mod.diet)
vif.cca(mod.diet.step)
diet.sub <- diet[,names(vif.cca(mod.diet.step)),drop=FALSE]
#combine subsetted variables and calculate varpart/ordination

pairs(everyVariable.sub)
# model on subsetted variables
mod.sub<-cca(larva.hell~., everyVariable.sub)
plot(mod.sub)
vif.cca(mod.sub)
mod.sub.step <- ordistep(mod.sub)
plot(mod.sub.step)
# variance partitioning on subsetted matrices
larva.vp.sub <- varpart(larva.hell, environmentals.sub,MEM.sub,diet.sub)
plot(larva.vp.sub)
# try with HST
larva.vp.sub.HST <- varpart(larva.hell, environmentals.sub,MEM.sub,diet.sub,HST.num)
plot(larva.vp.sub.HST,Xnames = (c("Environment","Spatial autocorrelation","Foliage communities","Host-tree species")),cex=1.5,cex.sub=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.5)


larva.vp.sub.HST <- varpart(vegdist(larva.hell), environmentals.sub,MEM.sub,diet.sub,HST.num)
plot(larva.vp.sub.HST,Xnames = (c("Environment","Spatial autocorrelation","Foliage communities","Host-tree species")),cex=1.5,cex.sub=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.5,bg=c("blue","red","green","yellow"))
everyVariable.sub <- cbind(environmentals.sub,MEM.sub,diet.sub,HST.num)

#larva.vp.sub.HST <- varpart(larva.hell, environmentals.sub,MEM.sub,diet.sub,HST.num)
#plot(larva.vp.sub.HST,Xnames = (c("Environment","Spatial autocorrelation","Foliage communities","Host-tree species")),cex=1.5,cex.sub=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.5,bg=c("blue","red","green","yellow"))

tiff("varpartPCOA.tiff", units="in", width=9, height=7.5, res=300)

plot(larva.vp.sub.HST,Xnames = (c("Environment","Spatial structure","Foliar communities","Host-tree species")),cex=1.5,cex.sub=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.5,bg=c("blue","red","green","yellow"))

dev.off()


corr.env <- cor.table((environmentals))
write.csv(corr.env,file="correnv.csv")

corr.mem <- cor.table((MEM))
write.csv(corr.mem,file="corrmem.csv")

#use the unrarefied dataset, but only samples who have reads in both foliage in larvae samples
comp.comm <- compare.comm[rownames(compare.abund.comm),]
comp.taxo <- compare.taxo[colnames(comp.comm),]
comp.meta <- compare.meta[rownames(comp.comm),]


###Ancom-BC analysis
library(ANCOMBC)

comp.meta$sample_type <- factor(comp.meta$sample_type)
comp.phylo<-phyloseq(otu_table(comp.comm,taxa_are_rows = FALSE), tax_table(as.matrix(comp.taxo)), sample_data(comp.meta))

comp.ancom <- ancom.prep(comp.phylo,"sample_type")
comp.ancom.ggplot <- comp.ancom[[2]]
comp.ancom.table <- comp.ancom[[1]]

comp.ancom.ggplot.complete <- comp.ancom.ggplot + scale_shape_manual(labels= c("Foliage","Larvae"),values = c(16,18)) + labs(shape="Associated with type of sample:")
comp.ancom.ggplot.complete


tiff("comp.ancom.tiff", units="in", width=10, height=6, res=300)

comp.ancom.ggplot.complete

dev.off()


comp.ancom.table <- comp.ancom.table[order(decreasing = TRUE,abs(comp.ancom.table$beta)),]
write.csv(comp.ancom.table, "comp.ancom.csv")


#make a figure of all the most differentially abundant genera
comp.ancom.betas <- comp.ancom.table[,c("beta","Genus")]
comp.ancom.betas$beta <- abs(comp.ancom.betas$beta)
comp.ancom.betas <- comp.ancom.betas[order(comp.ancom.betas$beta,decreasing = TRUE),]
comp.ancom.top12 <- comp.ancom.betas[c(1:12),]
top12.comm <- data.frame(comp.comm[,rownames(comp.ancom.top12)])
colnames(top12.comm) <- comp.ancom.top12$Genus
top12.comm$group <- comp.meta[rownames(top12.comm),"sample_type"]
top12.melt <- melt(top12.comm,id.vars="group")
comp.top12 <- ggplot(data=top12.melt,aes(x=variable,y=value,fill=group,color=group)) + geom_jitter(position=position_jitterdodge(dodge.width=0.9)) + scale_color_manual(label=c("Foliage","Larvae"),values = c(viridis(5)[1],viridis(5)[3])) + labs(color="Associated with type of sample:") + guides(fill=FALSE) + ylab("Abundance") + xlab("Bacterial genus")

tiff("comp.top12.tiff", units="in", width=15, height=6, res=300)

comp.top12

dev.off()


tiff("MEMs.leg.tiff", units="in", width=30, height=18, res=300)

plot(mem.gab[,signi])

dev.off()

tiff("MEMs.tiff", units="in", width=15, height=6, res=300)

plot(mem.gab[,signi], SpORcoords = coords)

dev.off()

