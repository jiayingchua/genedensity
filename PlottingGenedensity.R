rm(list=ls())
graphics.off()

######### LOAD GFF #########
##Ridaeus_Ras1
require(ape)

gfffile <- "C:/Users/JiaYing/Group Project/GFF/Ridaeus_Ras1_v1.genes.gff3"
gfffile <- "C:/Users/JiaYing/Group Project/GFF/AnitraGenome1161_2022_01_08.gff"
gfffile <- "C:/Users/JiaYing/Group Project/GFF/Hillquist_release2.gff"

gff<- read.gff(gfffile, na.strings = c(".","?"), GFF3 = TRUE)

## scale the references
BP_query <- 31454125 #rid1
#BP_reference <- 28685549 #anitra
BP_reference <- 31951077 #Hillquist
scale <- BP_query/BP_reference

## user to give number of chromosomes (noc) in genome, eg. 7
noc <- 7

genes <- data.frame()
seqid_list <- c()
for (i in 1:nrow(gff)){
  if (gff$type[i] == "gene"){
    seqid <- as.character(gff$seqid[i])
    if (!(seqid %in% seqid_list)){
      seqid_list <- append(seqid_list, seqid)
    }
    seqid <- gsub("[^0-9]", "", seqid)
    if (seqid == noc+1){
      break
    }
    start <- gff$start[i]*scale
    end <- gff$end[i]*scale
    genes <- rbind(genes, list(seqid, start, end))
  }
}
colnames(genes) <- c("seqid", "start", "end")
#print(seqid_list)

##gene density
# genes$seqid <- gsub("[^0-9]", "", genes$seqid)
for (i in 1:length(genes$seqid)){
  if (startsWith(genes$seqid[i], "0")){
    genes$seqid[i] <- gsub("[0]", "", genes$seqid[i])
  }
}
hist(genes$start, breaks = 1000)

## plot gene density
# choose specific chromosome, c
c=1
positions <- data.frame()
for (i in 1:nrow(genes)){
  if (genes$seqid[i] == c) {
    positions <- rbind(positions, genes$start[i])
  }
}
window <- 100000
count = 0
lengthcounter = 0
density <- c()
start = 0
end = start + window
loops = 0

alloutputs <- c()
ends <- c()

while (lengthcounter < nrow(positions)) {
  temp = end
  ends <- append(ends,temp)
  for (i in 1:nrow(positions)) {
    if ((positions[i,1] >= start) & (positions[i,1] < end)) {
      count=count+1;
      lengthcounter=lengthcounter+1;
    }
  }
  density <- append(density, count)
  count = 0;
  start = temp
  end = start + window
}

x11()
#barplot(density, xlab = 'Bases', ylab = 'Gene density', names.arg = ends, main = paste("Variant density across base window of", window))
plot(density, type = "l", xlab = 'Bp (e+05)', ylab = 'Gene density',
     ylim = c(0,max(density)+10))


# save unique seqid


##Hillquist (gff and assembly matches)
gfffile_hq <- "C:/Users/JiaYing/Group Project/Rubus_argutus/Hillquist_release2.gff/Hillquist_release2.gff"
gff_hq<- read.gff(gfffile_hq, na.strings = c(".","?"), GFF3 = TRUE)

genes_hq <- data.frame()
seqid_list_hq <- c()
for (i in 1:nrow(gff_hq)){
  if (gff_hq$type[i] == "gene"){
    seqid <- as.character(gff_hq$seqid[i])
    if (!(seqid %in% seqid_list_hq)){
      seqid_list_hq <- append(seqid_list_hq, seqid)
    }
    start <- gff_hq$start[i]
    end <- gff_hq$end[i]
    genes_hq <- rbind(genes_hq, list(seqid, start, end))
  }
}
colnames(genes_hq) <- c("seqid", "start", "end")
print(seqid_list_hq)


#lines(density2)

png("chr1_genedensity.png", width = 30, height = 5, units = "in", res = 100)
barplot(density, main = paste("Gene density of Chr", c, "with base window size:", window))
dev.off()


##### HIDE #####
##Rubus_Idaeus (gff and assembly doesn't match)
gfffile_RI <- "C:/Users/JiaYing/Group Project/Rubus_idaeus/AnitraGenome1161_2022_01_08_masked.nohints.gff"
gff_RI<- read.gff(gfffile_RI, na.strings = c(".","?"), GFF3 = TRUE)

genes_RI <- data.frame()
seqid_list_RI <- c()
for (i in 1:nrow(gff_RI)){
  if (gff_RI$type[i] == "gene"){
    seqid <- as.character(gff_RI$seqid[i])
    if (!(seqid %in% seqid_list_RI)){
      seqid_list_RI <- append(seqid_list_RI, seqid)
    }
    start <- gff_RI$start[i]
    end <- gff_RI$end[i]
    genes_RI <- rbind(genes_RI, list(seqid, start, end))
  }
}
colnames(genes_RI) <- c("seqid", "start", "end")
print(seqid_list_RI)