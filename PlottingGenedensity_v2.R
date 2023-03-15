## PLOTS GENE DENSITY OF MULTIPLE GFF FILES
## INPUTS - .gff files, length of chromosomes, noc

rm(list = ls())
graphics.off()

# install.packages("ape")
# install.packages("stringr")
# install.packages("RColorBrewer")
# install.packages("seqinr")
require(ape)
require(stringr)
require(RColorBrewer)
require(seqinr)

# For plot later
cols <- brewer.pal(9, 'Set1')
legend <- c()

# Vector of file paths
gfffiles <- c() ## first index = query, second index = ref

gfffile <-
  "C:/Users/JiaYing/Group Project/GFF/Ridaeus_Ras1_v1.genes.gff3"
legend <- append(legend, tail(str_split(gfffile, "/")[[1]], 1))
gfffiles <- append(gfffiles, gfffile)

gfffile <-
  "C:/Users/JiaYing/Group Project/GFF/Hillquist_release2.gff"
legend <- append(legend, tail(str_split(gfffile, "/")[[1]], 1))
gfffiles <- append(gfffiles, gfffile)

gfffile <- "C:/Users/JiaYing/Group Project/GFF/Rubus_occ_V3.all.gff"
legend <- append(legend, tail(str_split(gfffile, "/")[[1]], 1))
gfffiles <- append(gfffiles, gfffile)

gfffile <- "C:/Users/JiaYing/Group Project/GFF/Rubus_chingii_Hu.gff"
legend <- append(legend, tail(str_split(gfffile, "/")[[1]], 1))
gfffiles <- append(gfffiles, gfffile)

# Require vector of chromosome basepair lengths

# user to give number of chromosomes (noc) in genome, eg. 7
noc <- 7

# Option 1: Load FASTA and count from FASTA
# Read FASTA
bplen <- c()
fastafiles <- c() ##IMPORTANT! FIRST FILEPATH MUST BE QUERY!!

fastafile <-
  "C:/Users/JiaYing/Group Project/FASTA/Ridaeus_Ras1_v1.fasta.gz"
fastafiles <- append(fastafiles, fastafile)
fastafile <-
  "C:/Users/JiaYing/Group Project/Rubus_argutus/Hillquist_genome_v1_purged_primary_contigs_HiC.fasta.gz"
fastafiles <- append(fastafiles, fastafile)
fastafile <-
  "C:/Users/JiaYing/Group Project/FASTA/Rubus_occ_V3_10-12-17.fasta.gz"
fastafiles <- append(fastafiles, fastafile)
fastafile <-
  "C:/Users/JiaYing/Group Project/FASTA/Rubus_chingii_Hu_fina.fa.gz"
fastafiles <- append(fastafiles, fastafile)

# fastafile <-
#   "C:/Users/JiaYing/Group Project/FASTA/Burbank_genome_v1_purged_primary_contigs_HiC.fasta.gz"
# fastafiles <- append(fastafiles, fastafile)
# fastafile <-
#   "C:/Users/JiaYing/Group Project/FASTA/Rubus_Idaeus.genome.fa.gz"
# fastafiles <- append(fastafiles, fastafile)
# fastafile <-
#   "C:/Users/JiaYing/Group Project/FASTA/Rubus_occidentalis_masked_v1.1.fasta.gz"
# fastafiles <- append(fastafiles, fastafile)
# fastafile <-
#   "C:/Users/JiaYing/Group Project/FASTA/Rubus_occidentalis_v1.1.fasta.gz"
# fastafiles <- append(fastafiles, fastafile)

header_length_table <- data.frame()

for (file in fastafiles) {
  lengths <- c()
  headers <- c()
  fasta <- read.fasta(file, seqtype = 'DNA')
  lengths <- getLength(fasta)
  header <- gsub("[^0-9A-Za-z]", "", getAnnot(fasta))
  for (h in 1:length(header)) {
    if (getLength(getAnnot(fasta)[h]) < 8) {
      header_length_table <-
        rbind(header_length_table, c(header[[h]], lengths[h]))
    }
  }
}

# Line up all chromosomes (row) of each fastafile (col)
headpoints <- nrow(header_length_table) / length(fastafiles)
header_length_table2 <- data.frame()
for (a in 1:noc) {
  chrlen <- c()
  for (i in 1:nrow(header_length_table)) {
    if (grepl(a, header_length_table[i, 1])) {
      chrlen <- append(chrlen, header_length_table[i, 2])
    }
  }
  header_length_table2 <- rbind(header_length_table2, chrlen)
}


# # Option 2: User input
# bplen <- c() ## first index = query, second index  = ref
# bp_input <- 31454125 #Rid1
# bplen <- append(bplen, bp_input)
# bp_input <- 31951077 #Hillquist
# bplen <- append(bplen, bp_input)
# bp_input <- 34159342 #OccV3 (random number)
# bplen <- append(bplen, bp_input)
# bp_input <- 25609021 #Hu (random number)
# bplen <- append(bplen, bp_input)

## remove fasta information
rm(fasta, file, fastafile, fastafiles, headpoints, header_length_table, chrlen, headers, lengths)


# Set window size
window <- 100000

if (window < 1000) {
  windowlabel = paste(window, "b")
}
if ((window < 1000000) && (window >= 1000)) {
  windowlabel = paste(window / 1000, "kb")
}
if (window >= 1000000) {
  windowlabel = paste(window / 1000000, "Mb")
}

# Plots all chromosomes (noc = 7)
for (c in 3:noc) {
  # declare density table
  density_table <- data.frame()
  density_table2 <- data.frame()
  
  for (f in 1:length(gfffiles)) {
    # read gff
    gff <- read.gff(gfffiles[f],
                    na.strings = c(".", "?"),
                    GFF3 = TRUE)
    # scale the references
    scale <- as.numeric(header_length_table2[c,1]) / as.numeric(header_length_table2[c,f])
    
    genes <- data.frame()
    seqid_list <- c()
    for (r in 1:nrow(gff)) {
      if (gff$type[r] == "gene") {
        seqid <- as.character(gff$seqid[r])
        if (!(seqid %in% seqid_list)) {
          seqid_list <- append(seqid_list, seqid)
        }
        seqid <- gsub("[^0-9]", "", seqid)
        if (seqid == noc + 1) {
          break
        }
        start <- gff$start[r] * scale
        end <- gff$end[r] * scale
        genes <- rbind(genes, list(seqid, start, end))
      }
    }
    colnames(genes) <- c("seqid", "start", "end")
    
    # gene density
    for (s in 1:length(genes$seqid)) {
      if (startsWith(genes$seqid[s], "0")) {
        genes$seqid[s] <- gsub("[0]", "", genes$seqid[s])
      }
    }
    
    # plot gene density
    positions <- data.frame()
    for (g in 1:nrow(genes)) {
      if (genes$seqid[g] == c) {
        positions <- rbind(positions, genes$start[g])
      }
    }
    
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
      ends <- append(ends, temp)
      for (i in 1:nrow(positions)) {
        if ((positions[i, 1] >= start) & (positions[i, 1] < end)) {
          count = count + 1
          
          lengthcounter = lengthcounter + 1
          
        }
      }
      density <- append(density, count)
      count = 0
      
      start = temp
      end = start + window
    }
    
    density_table <- data.frame(density)
    density_table2 <- rbind(density_table2, density_table)
  }
  
  densitypoints <- nrow(density_table2) / length(gfffiles)
  
  ##PLOT
  # for 4 files, 1:315, 316:630, 631:945, 946:1260
  #x11()
  #barplot(density, xlab = 'Bases', ylab = 'Gene density', names.arg = ends, main = paste("Variant density across base window of", window))
  png(
    paste("genedensity_chr", c, ".png", sep = ""),
    width = 1200,
    height = 300,
    units = "px"
  )
  plot(
    density_table2[1:densitypoints, ],
    type = "l",
    xlab = "",
    ylab = "",
    main = paste("Chr", c, "\n Window size:", windowlabel),
    cex.main = 1,
    cex.axis = 0.8,
    ylim = c(0, max(density) + 5),
    col = cols[1]
  )
  for (i in 2:length(gfffiles)) {
    startindex <- (densitypoints * (i - 1)) + 1
    endindex <- densitypoints * i
    lines(density_table2[startindex:endindex, ],
          col = cols[i])
  }
  legend(
    "topright",
    legend = legend,
    col = cols[1:length(gfffiles)],
    lty = 1,
    cex = 0.7
  )
  dev.off()
}
