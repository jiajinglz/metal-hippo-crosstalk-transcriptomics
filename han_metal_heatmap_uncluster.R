#########################################
###import counts data and format matrix##
countsTable <- read.csv(file="C:/Users/JJ/Desktop/han_metal/counts summary.csv",header=T,stringsAsFactors =FALSE)
countsMatrix <- data.matrix(countsTable[,2:ncol(countsTable)])
rownames(countsMatrix) <- countsTable[,1] #set 1st column as row names for matrix
######################
###heat map and gene cluster###
library("RColorBrewer")
library('gplots') ##can be installed in packages##
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),               # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green
# creates a 5 x 5 inch image
png("C:/Users/JJ/Desktop/han_metal/heatmaps_in_r.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
heatmap.2(countsMatrix,
          cellnote = countsMatrix,  # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering