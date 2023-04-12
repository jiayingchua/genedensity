## PLOTS GENE DENSITY OF MULTIPLE GFF FILES
## INPUTS - .gff files, length of chromosomes, noc

rm(list = ls())
graphics.off()

args = commandArgs(trailingOnly = TRUE)
# argument 1 = gff paths text file
# argument 2 = fasta paths text path
# argument 3 = output file path
# argument 4 = window size
# argument 5 = no. of chromosomes
gffpaths <- args[1]
#gffpaths <- "C:\\Users\\JiaYing\\GP\\gfffiles.txt"
#gffpaths <- "C:\\Users\\JiaYing\\OneDrive - Cranfield University\\Documents\\Group Project\\jy1\\gfffiles.txt"
fastapaths <- args[2]
#fastapaths <- "C:\\Users\\JiaYing\\GP\\fastafiles.txt"
#fastapaths <- "C:\\Users\\JiaYing\\OneDrive - Cranfield University\\Documents\\Group Project\\jy1\\fastafiles.txt"
outputfile <- args[3]
#outputfile <- "C:\\Users\\JiaYing\\GP\\OUT"
#outputfile <- "C:\\Users\\JiaYing\\OneDrive - Cranfield University\\Documents\\Group Project\\jy1\\output"
window <- as.numeric(args[4])
#window <- 1000000
noc <- as.numeric(args[5])
#noc <- 2

# List of packages for session
.packages = c("ape", "stringr", "RColorBrewer", "seqinr")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if (length(.packages[!.inst]) > 0)
  install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only = TRUE)

# For plot later
cols <- brewer.pal(9, 'Set1')

# Read gfffile.txt
readGFF = function(filepath) {
  gfffiles <- c()
  con = file(filepath, "r")
  while (TRUE) {
    gfffile = readLines(con, n = 1)
    if (length(gfffile) == 0) {
      break
    }
    gfffiles <- append(gfffiles, gfffile)
  }
  return(gfffiles)
  close(con)
}

# Read fastafiles.txt
readFASTA = function(filepath) {
  fastafiles <- c()
  con = file(filepath, "r")
  while (TRUE) {
    fastafile = readLines(con, n = 1)
    if (length(fastafile) == 0) {
      break
    }
    fastafiles <- append(fastafiles, fastafile)
  }
  return(fastafiles)
  close(con)
}

gfffiles <- readGFF(gffpaths)
fastafiles <- readFASTA(fastapaths)

# Set keys for legend based on name of gff file
legend <- c()
for (i in 1:length(gfffiles)){
  gfffile <- gfffiles[i]
  legend <- append(legend, tail(str_split(gfffile, "[^-_A-Za-z0-9.]")[[1]],1))
}

###
header_length_table_v3 <- data.frame()

for (file in fastafiles) {
  lengths <- c()
  headers <- c()
  # read the fasta file
  fasta <- read.fasta(file, seqtype = 'DNA')
  # extract lengths and headers of each sequence
  lengths <- getLength(fasta)
  header <- gsub("[^_0-9A-Za-z]", "", getAnnot(fasta))
  for (h in 1:length(header)) {
    seqnumber <- gsub("[^0-9]", "", header[[h]])
    # header length table consists of id, sequence number and length of each fastafile
  header_length_table_v3 <-
    rbind(header_length_table_v3, c(header[[h]], seqnumber, lengths[h]))
  }
  
}
  colnames(header_length_table_v3)[1] <- "seqid"
  colnames(header_length_table_v3) [2] <- "seqnum"
  colnames(header_length_table_v3) [3] <- "seqlen"

# Line up all chromosomes (row) of each fastafile (col)
header_length_table2_v3 <- data.frame()
for (a in 1:noc) {
  chrlen <- c()
  for (i in 1:nrow(header_length_table_v3)) {
    if (header_length_table_v3[i, 2] == a) {
      print(i)
      print(header_length_table_v3[i, 2])
      chrlen <- append(chrlen, header_length_table_v3[i, 3])
    }
  }
  header_length_table2_v3 <- rbind(header_length_table2_v3, chrlen)
}

## remove fasta information
rm(a, i, h, fasta, file, fastafiles, gfffile, gffpaths, fastapaths, header_length_table_v3, chrlen, headers, lengths)


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
for (c in 1:noc) {
  # declare density table
  density_table <- data.frame()
  density_table2 <- data.frame()
  
  for (f in 1:length(gfffiles)) {
    # read gff
    gff <- read.gff(gfffiles[f],
                    na.strings = c(".", "?"),
                    GFF3 = TRUE)
    # scale the references
    scale <-
      as.numeric(header_length_table2_v3[c, 1]) / as.numeric(header_length_table2_v3[c, f])
    
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
  mindensity <- min(density_table2)
  maxdensity <- max(density_table2)
  densitypoints <- nrow(density_table2) / length(gfffiles)
  
  ##PLOT
  # for 4 files, 1:315, 316:630, 631:945, 946:1260
  #x11()
  #barplot(density, xlab = 'Bases', ylab = 'Gene density', names.arg = ends, main = paste("Variant density across base window of", window))
  png(
    paste(outputfile, "_genedensity_chr", c, ".png", sep = ""),
    width = 1200,
    height = 300,
    units = "px"
  )
  plot(
    density_table2[1:densitypoints,],
    type = "l",
    xlab = "",
    ylab = "",
    main = paste("Chr", c, "\n Window size:", windowlabel),
    cex.main = 1,
    cex.axis = 0.8,
    ylim = c(mindensity - 5, maxdensity + 5),
    col = cols[1]
  )
  for (i in 2:length(gfffiles)) {
    startindex <- (densitypoints * (i - 1)) + 1
    endindex <- densitypoints * i
    lines(density_table2[startindex:endindex,],
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
