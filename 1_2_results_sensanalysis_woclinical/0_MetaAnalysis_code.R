#############################
# Data
#############################
files <- dir()
files <- files[grep(pattern = "Giovanna_MetaAnalysis.txt", x = files)]
data_original <- read.delim(file = files, header = TRUE, sep = "\t", dec = ",")
source(file = paste(dirname(getwd()), "MetaA_GregorAdapt_202206.R", sep = "/"))
