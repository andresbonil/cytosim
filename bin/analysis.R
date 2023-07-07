# Analysis of cytosim output results
# Written by A.G., 20230618

# Cytosim Run script integration 20230705

library(ggplot2)
library(stringr)

# This script will be run directly after each simulation, outputting a graph to the directory

# Argument will contain path to folder/output file to open, this is to be specified after calling script
args <- commandArgs(trailingOnly = TRUE)

# This code takes the folder name with parameters/values to use as the filename
split <- unlist(strsplit(args[1], "/"))
location <- args[2]
fileName <- split[length(split)]

# Initial paramter before simulation
MT_num = 40
MT_length = 40
segmentation = 0.1
nb_frames = 200

# Creating path for the output file and reading it 
data <- read.csv(paste(args[1], "/output.txt", sep=""))

sum = 0
k = 0
for (i in nrow(data):1) {
  sum = sum + data$X400[[i]]
  k = k + 1
  if (sum >= (MT_num * MT_length/segmentation)) {
      break
    }
  }
  
lastFrameIndex = seq(nrow(data)-k, nrow(data))
lastFrameData = data.frame(MT_length = data[lastFrameIndex, ]*segmentation)

# Embedded paste() used to create path to output folder and add .pdf extension
# Change the path in the first paste to the folder's name on your machine 
location <- paste(location, "/", sep = "")
pdf(file = paste(location, paste(fileName, ".pdf", sep=""), sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4)
print(
  ggplot(lastFrameData, aes(x = MT_length)) +
    geom_histogram(bins = 10) +
    scale_x_continuous(limits = c(0, 45))+
    scale_y_continuous(limits = c(0, 45))+
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18)
    )+
    xlab("MT length")+
    ylab("Count")
)
  
dev.off()

  


