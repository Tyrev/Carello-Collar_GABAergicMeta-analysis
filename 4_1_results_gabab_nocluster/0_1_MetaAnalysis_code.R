#############################
# Data
#############################
files <- dir()
files <- files[grep(pattern = "Giovanna_MetaAnalysis.txt", x = files)]
data_original <- read.delim(file = files, header = TRUE, sep = "\t", dec = ",")
source(file = paste(getwd(), "0_2_MetaAnalysis_code.R", sep = "/"))
