result.path = "/Users/Jelena/Documents/collect_cluster_results/full/"
setwd(result.path)

#method="liu"
#method="lee"
#method="knockoff_result"
#method="knockoff+_result"
#method="randomized_full"
#method="randomized_partial"

txt_files = list.files(pattern = paste(method, "*", sep=""))

i = 0

pvalues = NULL
sel_intervals=NULL
sel_coverages=NULL
sel_lengths=NULL
naive_pvalues = NULL
naive_intervals=NULL
naive_coverages=NULL
naive_lengths=NULL
FDR_sample = NULL
power_sample=NULL

for (file in txt_files){
  result=readRDS(file)
  pvalues=c(pvalues, result$pvalues)
  sel_intervals=c(sel_intervals, result$sel_intervals)
  sel_coverages=c(sel_coverages, result$sel_coverages)
  sel_lengths=c(sel_lengths, result$sel_lengths)
  naive_pvalues=c(naive_pvalues, result$naive_pvalues)
  naive_intervals=c(naive_intervals, result$naive_intervals)
  naive_coverages=c(naive_coverages, result$naive_coverages)
  naive_lengths=c(naive_lengths, result$naive_lengths)
  FDR_sample=c(FDR_sample, result$FDR_sample)
  power_sample=c(power_sample, result$power_sample)
  n=result$n
  p=result$p
  s=result$s
  rho=result$rho
  snr=result$snr
  i=i+1
}

print(i)

setwd("/Users/Jelena/GitHub Jelena/R-selective/tests/comparisons/")

saveRDS(list(pvalues=pvalues, 
             sel_intervals=sel_intervals, sel_coverages=sel_coverages, sel_lengths=sel_lengths,
             naive_pvalues=naive_pvalues, 
             naive_intervals=naive_intervals, naive_coverages=naive_coverages, naive_lengths=naive_lengths,
             n=n, p=p, s=s, rho=rho, snr=snr,
             FDR_sample=FDR_sample, power_sample=power_sample), file=paste(method,".rds", sep=""))




