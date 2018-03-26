
lee = readRDS(file="lee.rds")
liu = readRDS(file="liu.rds")
rand_full=readRDS(file="randomized_full.rds")
rand_partial=readRDS(file="randomized_partial.rds")

n=lee$n
p=lee$p
s=lee$s
snr=lee$snr
rho=lee$rho

lee_cov = mean(lee$sel_coverages)
lee_len = median(lee$sel_lengths)

liu_cov=mean(liu$sel_coverages)
liu_len=median(liu$sel_lengths)

naive_cov = mean(liu$naive_coverages)
naive_len=median(liu$naive_lengths)

rand_full_cov = mean(rand_full$sel_coverages)
rand_full_len = median(rand_full$sel_lengths)

rand_partial_cov=mean(rand_partial$sel_coverages)
rand_partial_len=median(rand_partial$sel_lengths)

ninf_lee = length(which(lee$sel_lengths==Inf))/length(lee$sel_lengths)
ninf_liu = length(which(liu$sel_lengths==Inf))/length(liu$sel_lengths)
c(ninf_lee, ninf_liu)

labels = rep("",5)
labels[1] = paste("len: ", round(lee_len,2), ", cov: ", round(lee_cov,2), ", ninf: ", round(ninf_lee,2), sep = "")
labels[2] = paste("len: ", round(liu_len,2), ", cov: ", round(liu_cov,2), ", ninf: ", round(ninf_liu,2), sep = "")
labels[3] = paste("len: ", round(naive_len,2), ", cov: ", round(naive_cov,2), sep = "")
labels[4] = paste("len: ", round(rand_full_len,2), ", cov: ", round(rand_full_cov,2), sep = "")
labels[5] = paste("len: ", round(rand_partial_len,2), ", cov: ", round(rand_partial_cov,2), sep = "")


title = paste("n=", n,", p=", p, ", s=", s, ", signal=", round(snr,2), ", rho=", rho, sep = "")

boxplot(lee$sel_lengths[which(lee$sel_lengths!=Inf)], liu$sel_lengths, liu$naive_lengths,
        rand_full$sel_lengths, rand_partial$sel_lengths, 
        at = seq(1,11,2.5) - 0.25, col = "pink", las = 2, outline = F,
        names = c("Lee \n (full)", "Liu \n (full)", "naive", "rand \n (full)", "rand \n (partial)"), 
        boxwex = 0.6, range = 0.5,
        log="x", horizontal = T, xlim = c(0,12.5), cex.axis = 1)

title(title, cex.main = 1, font.main = 1, line = 0.2)
text(labels = labels, x = c(lee_len, liu_len, naive_len, rand_full_len,rand_partial_len)+0.05, y = seq(1,11,2.5) + 0.5, cex = 1)




plot(ecdf(lee$pvalues), xlab="values", ylab="Empirical cdf", main="Selective pvalues", 
     col="red")
lines(ecdf(liu$pvalues), col="blue")
abline(c(0,1))
legend("bottomright", legend=c("Lee et al", "Liu et al"),
       col=c("red", "blue"), lty=1:1, cex=1.5)



