pst3 = pst

## Choose a set of soft-thresholding powers
powers = c(c(1:22))
cex1 = 0.9

E<-list()
# Set the loop for 1-n number of permutations ran
for(i in 1:5){
  E[[i]]<-new.env()
  load(paste0('~/Downloads/',i,'.pst.RData'),env=E[[i]])
}

# Set the number of rows and columns for the plots to fill in
par(mfrow=c(2,3))

## Plot
for(i in 1:5){
  pst<-E[[i]]$pst
# Scale free topology fit
plot(pst$fitIndices[,1],
     -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     xlab="Soft pst (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(pst$fitIndices[,1],
     
     -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-psting power
#plot(pst$fitIndices[,1], pst$fitIndices[,5],
#     xlab="Soft pst (power)",ylab="Mean Connectivity", type="n",
#     main = paste("Mean connectivity"))
#text(pst$fitIndices[,1], pst$fitIndices[,5], labels=powers, cex=cex1,col="red")
}