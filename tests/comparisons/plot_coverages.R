
lee = readRDS(file="lee_full.rds")
liu = readRDS(file="liu_full.rds")

n=lee$n
p=lee$p
s=lee$s
snr=lee$snr
rhp=lee$rho

lee_cov = mean(lee$sel_coverages)
lee_len = median(lee$sel_lengths)

liu_cov=mean(liu$sel_coverages)
liu_len=median(liu$sel_lengths)

naive_cov = mean(liu$naive_coverages)
naive_len=median(liu$naive_lengths)

labels = rep("",3)
labels[1] = paste("len: ", round(lee_len,2), ", cov: ", round(lee_cov,2), sep = "")
labels[2] = paste("len: ", round(liu_len,2), ", cov: ", round(liu_cov,2), sep = "")
labels[3] = paste("len: ", round(naive_len,2), ", cov: ", round(naive_cov,2), sep = "")

title = paste("n=", n,", p=", p, ", s=", s, ", signal=", snr, ", rho=", rho, sep = "")

boxplot(lee$sel_lengths[which(lee$sel_lengths!=Inf)], liu$sel_lengths, liu$naive_lengths,
        at = seq(1,11,4) - 0.25, col = "pink", las = 2, outline = F,
        names = c("Lee", "Liu", "naive"), boxwex = 0.6, range = 0.5,
        log="x", horizontal = T, xlim = c(0,12.5), cex.axis = 1)

title(title, cex.main = 1, font.main = 1, line = 0.2)
text(labels = labels, x = c(lee_len, liu_len, naive_len)+0.3, y = seq(2,11,4) - 0.25, cex = 1)

ninf_lee = length(which(lee$sel_lengths==Inf))/length(lee$sel_lengths)
ninf_liu = length(which(liu$sel_lengths==Inf))/length(liu$sel_lengths)
c(ninf_lee, ninf_liu)

