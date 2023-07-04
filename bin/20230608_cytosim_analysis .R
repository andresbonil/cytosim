# Analysis of cytosim output results
# Written by A.G., 20230618

library(ggplot2)

# Run script for each output file once simulation is done

# Argument will contain path to folder/output file to open
args <- commandArgs(trailingonly = TRUE)
setwd("/Users/qgeng/Desktop/cytosim_outputs")

# TODO: make sure script works for an output.txt, and saves graphs to a directory within the cytosim folder.
# ---------------------------------------------------------------------------------------------------------------------

# fileList <- list.files(path=".", pattern="*.txt", full.names=TRUE)
# fileNameList = tools::file_path_sans_ext(basename(fileList))

# Initial paramter before simulation
MT_num = 40
MT_length = 40
segmentation = 0.1
nb_frames = 200


for (f in 1:length(fileList)) {
  data <- read.csv(fileList[[f]])
  # # For test, MTs in the first frame
  # sum = 0
  # k = 0
  # for (i in 1:nrow(data)) {
  #   sum = sum + data$X400[[i]]
  #   k = k + 1
  #   if (sum >= (MT_num * MT_length/segmentation)) {
  #     print(sum)
  #     print(k)
  #     break
  #   }
  # }
  
  # extract MTs in the last frame; 
  # add up MT length and store in "sum" until it reaches MT_num * MT_length
  # while "k" stores the count of the MT number in the last frame
  sum = 0
  k = 0
  for (i in nrow(data):1) {
    sum = sum + data$X400[[i]]
    k = k + 1
    if (sum >= (MT_num * MT_length/segmentation)) {
      # print(sum)
      # print(k)
      break
    }
  }
  
  lastFrameIndex = seq(nrow(data)-k, nrow(data))
  lastFrameData = data.frame(MT_length = data[lastFrameIndex, ]*segmentation)
  
  fileName = fileNameList[f]
  pdf(file = paste(fileName, ".pdf", sep=""),   # The directory you want to save the file in
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
}
  


