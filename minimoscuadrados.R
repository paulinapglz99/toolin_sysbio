plot(log10(table[,9]),table[,1],pch=16,cex=2,col=rgb(0,0,0,0.5),xlab="Log10 Dilution",ylab="Signal",ylim=c(0,13),main="Fitting")


x<-10^seq(-7,-1,length=100)
lines(log10(x),f(x,2.2,9.5,0.0001))

residuals<-table[,1]-f(table[,9],2.2,9.5,0.0001)

plot(log10(table[,9]),residuals,pch=16,cex=2,col=rgb(0,0,0,0.5),xlab="Log10 Dilution",ylab="Residuals",main="Residuals")
abline(h=0)

residuals_mean=mean(residuals)
residuals_sd=sd(residuals)

residuals_regress<-lm(residuals ~ log10(table[,9])  )
residuals_regress_intersect<-coef(residuals_regress)[[1]]
residuals_regress_slope<-coef(residuals_regress)[[2]]
residuals_regress_intersect
residuals_regress_slope

text(-6,-1,paste("mean=",  round(residuals_mean,3),"\n", "sd=",round(residuals_sd,3)))
lines(log10(x), residuals_regress_intersect+residuals_regress_slope*log10(x),col="red")
