############################
# Power Analysis for Random effects Model
#############################
library(readxl)
tbl <- read_excel("Statistical_power.xlsx")
tbl <- as.data.frame(tbl)

for(r in seq(row.names(tbl))) {
     ### https://osf.io/w4xrs/
     ### https://www.linkedin.com/pulse/how-calculate-statistical-power-meta-analysis-jakob-tiebel/
     es <- abs(tbl[r,"Effect size"]) # Overall Effect Size
     k <- tbl[r,"Number of studies"] # Number of included studies
     # hg <- seq(from = 0, to = 3, by = 0.25) # Heterogeneity (0 no, .33 small, 1 moderate, 3 high)
     hg <- tbl[r,"Heterogeneity (I2)"] # Heterogeneity (0 no, .33 small, 1 moderate, 3 high)
     nt <- tbl[r,"Sample size (AD)"] # Number of participants in treatment group
     nc <- tbl[r,"Sample size (HC)"] # Number of participants in control group
     eq1 <- ((nt+nc)/((nt)*(nc))) + ((es^2)/(2*(nt+nc)))
     eq2 <- hg*(eq1)
     eq3 <- eq2+eq1
     eq4 <- eq3/k
     eq5 <- (es/sqrt(eq4))
     tbl[r,"Statistical power"] <- (1-pnorm(1.96-eq5)) # when alpha .05 -> 1.64 for one-tailed and 1.96 for two tailed
}
write.table(x = tbl, file = "Statistical_power.txt", row.names = FALSE, sep = "\t")