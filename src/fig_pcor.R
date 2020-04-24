tbl <- read.table('cache/glasso_sassa_cmp.tsv')

pdf('fig/pcor.pdf', width=3, height=3)
par(mfrow=c(1, 2), mar=c(6, 4, 4, 2), cex=0.7)
corr <- list()
corr[['GL']] <- tbl$cor_rmse_glasso
corr[['SAGL']] <- tbl$cor_rmse_sassa
boxplot(corr, main='A', ylab='', vertical=T, las=2, outpch=NA, yaxt='n')
axis(2, at=c(0.022, 0.023, 0.024))
title(ylab='RMSE', line=2)
stripchart(corr, vertical=T, pch=20, add=T, method='jitter')
pcor <- list()
pcor[['GL']] <- tbl$pcor_rmse_glasso
pcor[['SAGL']] <- tbl$pcor_rmse_sassa
boxplot(pcor, main='B', ylab='', vertical=T, las=2, outpch=NA, yaxt='n')
axis(2, at=c(0.016, 0.017, 0.018))
title(ylab='RMSE', line=2)
stripchart(pcor, vertical=T, pch=20, add=T, method='jitter')
dev.off()
