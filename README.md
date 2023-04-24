# genedensity

This repository was created for Cranfield BIX Group Project, Group 2.
It contains the R script that performs Gene Density analysis on the server. The code is run via command line.
The script that is currently on the server, and that the application runs, is PlottingGenedensity_v3.R

The code goes through a text file of gff file paths (gfffiles.txt) and calculates how many genes appear within each specified window.
A line graph is plotted for each gff file listed in the text file, and exported as an image (.png) to a specifed folder.
The code also reads through fastafiles.txt to determine the scaling of the plots in relation to the first file listed.

The order of the gff files in gfffiles.txt must correlate with that of their associated fasta file in fastafiles.txt.

Jia Ying Chua
April 2023
