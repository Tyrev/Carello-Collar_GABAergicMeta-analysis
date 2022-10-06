#############################
# Package Install/Load
#############################
## DOI: 10.1016/j.brs.2021.05.014
library(metafor)
library(dplyr)
library(metaviz)
library(ggplot2)

#############################
# Data
#############################
## Data.frame ##
### Columns description ###
#### - study: list of different studies
#### - brain_region: list of brain regions in each study
#### - n_c: number of samples in reference group for each study
#### - m_c: mean of samples in reference group for each study
#### - sd_c: standard deviation in reference group for each study
#### - n_t: number of samples in test group for each study
#### - m_t: mean of samples in test group for each study
#### - sd_t: standard deviation in test group for each study
#### - n_s: number of subjects in study
colnames(data_original) <- tolower(colnames(data_original))
data_original$id <- as.character(data_original$id)

data_split <- split(x = data_original, f = data_original$table.id)

resList <- lapply(data_split, function(data) {
     #############################
     # Calculate Effect Sizes
     #############################
     data_escalc <- escalc(n1i = n_t, m1i = m_t, sd1i = sd_t, n2i = n_c, m2i = m_c, sd2i = sd_c, measure = "SMD",
                           data = data, vtype = "UB", append = F) 
     # I use vtype = UB in my publications, it's more conservative
     
     data <- cbind(data, data_escalc)
     
     dup <- data$study[duplicated(data$study)]
     
     data_dup <- data[data$study%in%dup,]
     data_unique <- data[!data$study%in%dup,]
     row.names(data_unique) <- data_unique$study
     
     if(nrow(data_dup)>0) {
          n_samples <- split(x = data_dup, f = data_dup$study)
          n_samples <- lapply(n_samples, function(x) c(sum(x$n_c), sum(x$n_t)))
          n_samples <- do.call("rbind", n_samples)
          colnames(n_samples) <- c("n_c", "n_t")
          n_samples <- rbind(n_samples, data_unique[,c("n_c","n_t")])
          colnames(n_samples) <- c("n_samples_c", "n_samples_t")
     } else {
          n_samples <- data_unique[,c("n_c","n_t")]
          colnames(n_samples) <- c("n_samples_c", "n_samples_t")
     }
     
     #############################
     # Combined SMD (adjusted standardized mean difference by fixed-effect model)
     #############################
     studs <- unique(data_dup$study)
     res <- lapply(studs, function(x, dt) {
          ourdata <- dt[dt$study==x,]
          z <- rma.uni(yi = yi, vi = vi, measure = "SMD",
                       data = ourdata, method = "FE")
          c(z[1],x)
     }, dt = data_dup)
     data_comb <- data_dup[!duplicated(data_dup$study),]
     row.names(data_comb) <- data_comb$study
     for(s in studs) {
          data_comb[s, "n_c"] <- max(data_comb[s, "n_c"])
          data_comb[s, "n_t"] <- max(data_comb[s, "n_t"])
     }
     data_comb$yi <- as.numeric(lapply(res, '[[',1))
     data_comb$vi <- ((data_comb$n_c + data_comb$n_t)/(data_comb$n_c * data_comb$n_t) + 
                           data_comb$yi^2)/(2*(data_comb$n_c + data_comb$n_t))
     # here you have to be consistent with the vtype above 
     # (you can find the formula corresponding to your outcome measure in the escalc function script, 
     # for SMD and vtype UB  it is: (n1+n2)/(n1*n2) +y^2 / (2*(n1+n2)) )
     
     #############################
     #
     #############################
     data_final <- rbind(data_unique, data_comb)
     data_final$ci_low <- data_final$yi - (stats::qnorm(1 - (1 - 0.95)/2) * as.numeric(sqrt(data_final$vi)))
     data_final$ci_high <- data_final$yi + (stats::qnorm(1 - (1 - 0.95)/2) * as.numeric(sqrt(data_final$vi)))
     row.names(data_final) <- data_final$id
     data_final <- data_final[unique(data$id),]
     
     return(data_final)
})
resData <- do.call("rbind",resList)

#############################
# Save
#############################
write.table(x = resData, file = "4_1_results_gabab_nocluster.txt", sep = "\t", row.names = FALSE)
