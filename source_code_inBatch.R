res_folders <- dir()
res_folders <- res_folders[grep(pattern = "results", x = res_folders)]
res_folders <- res_folders[-length(res_folders)]
original_folder <- getwd()
for(d in res_folders) {
     setwd(d)
     code_file <- dir()
     code_file <- code_file[grep(pattern = "_code.R", x = code_file)]
     source(file = code_file)
     setwd(original_folder)
}
