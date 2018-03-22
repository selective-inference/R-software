lee = readRDS(file="lee_full.rds")
liu = readRDS(file="liu_full.rds")
kn = readRDS(file="knockoff.rds")
kn_plus = readRDS(file="knockoff+.rds")


table = matrix(c(mean(lee$FDR_sample), mean(lee$power_sample), 
                 mean(liu$FDR_sample), mean(liu$power_sample),
                 mean(kn$FDR_sample), mean(kn$power_sample),
                 mean(kn_plus$FDR_sample), mean(kn_plus$power_sample)), ncol=2, byrow=TRUE)
colnames(table)=c("FDR", "power")
rownames(table)=c("Lee", "Liu", "knockoff", "knockoff+")
table  
