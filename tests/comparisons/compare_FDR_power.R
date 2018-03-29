lee = readRDS(file="lee.rds")
liu = readRDS(file="liu.rds")
kn = readRDS(file="knockoff.rds")
kn_plus = readRDS(file="knockoff+.rds")
randomized_full=readRDS(file="randomized_full.rds")
randomized_partial=readRDS(file="randomized_partial.rds")


table = matrix(c(mean(lee$FDR_sample), mean(lee$power_sample), 
                 mean(liu$FDR_sample), mean(liu$power_sample),
                 mean(kn$FDR_sample), mean(kn$power_sample),
                 mean(kn_plus$FDR_sample), mean(kn_plus$power_sample),
                 mean(randomized_full$FDR_sample), mean(randomized_full$power_sample),
                 mean(randomized_partial$FDR_sample), mean(randomized_partial$power_sample)), ncol=2, byrow=TRUE)
colnames(table)=c("FDR", "power")
rownames(table)=c("Lee", "Liu", "knockoff", "knockoff+", "randomized full", "randomized partial")
table  
