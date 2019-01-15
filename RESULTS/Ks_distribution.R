


############################ Ks distribution #############################

MCL_paml_result <- read.delim("~/Desktop/MCL_paml_result.txt", header=FALSE)
View(MCL_paml_result)

res= MCL_paml_result
res=res[,-c(2,3)]
colnames(res)=c('family','omega','ka','ks')
lengh(unique(res$family))

hist(res$ks, n=2000, xlim=c(0,10), xlab='ks values', ylab='frequency', col='red')

length(res[which(res$ks <1),1]) / nrow(res)
 