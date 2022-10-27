#breast cancer data file
Breastdata=read.csv("C:\\Users\\Acer\\OneDrive\\Desktop\\cancer genomic\\GSE214507.csv",header=T,row.names = 1)
print(Breastdata)

matrix=Breastdata
for(i in 1:ncol(Breastdata)){
  matrix[,i]=(Breastdata[,i]/sum(Breastdata[,i]))*1000000
}
print(matrix)

log_matrix=log2(matrix+1)
saveRDS(log_matrix,file="log_matrix.rds")
summary(log_matrix)
print(log_matrix)

library(matrixStats)
z_score = (log_matrix - rowMeans(log_matrix))/rowSds(as.matrix(log_matrix))[row(log_matrix)]
print(z_score)

variance = apply(log_matrix,1,var)
variance = sort(variance,decreasing = T)
top50 = variance[1:50]
mat = z_score[names(top50),]


library(ComplexHeatmap)
Heatmap(mat)









