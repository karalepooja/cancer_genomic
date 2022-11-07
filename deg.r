data = source("Heatmap.r")
# converting data of log_matrix values into matrix
Bmatrix = matrix(NA,ncol = 4,nrow = nrow(log_matrix))
print(Bmatrix)
rownames(Bmatrix)=rownames(log_matrix)
colnames(Bmatrix)=c('LiverTumor','Control','pvalue','log2FC')    #log2FC = log2foldchange

for (i in 1:nrow(log_matrix)) {
  vector1 = as.numeric(log_matrix[i,1:3])
  vector2 = as.numeric(log_matrix[i,4:6])
  
  test1 = t.test(vector1,vector2,paired = F,alternative = "two.sided")
  Bmatrix[i,1]=test1$estimate[[1]]
  Bmatrix[i,2]=test1$estimate[[2]]
  Bmatrix[i,3]=test1$p.value
  Bmatrix[i,4]=Bmatrix[i,1]-Bmatrix[i,2]
  
}

Bmatrix = as.data.frame(Bmatrix)
num=which(is.nan(Bmatrix$pvalue))
Bmatrix[num,'pvalue'] = 1

# save Bmatrix to rds format
saveRDS(Bmatrix,file = "Bmatrix.rds")


# Visualizing of differentially expressed genes using Volcano plot  using package 'EnhancedVolcano'

library(EnhancedVolcano)
vplot = EnhancedVolcano(Bmatrix,lab = rownames(Bmatrix),x = 'log2FC', y = 'pvalue' )
print(vplot)