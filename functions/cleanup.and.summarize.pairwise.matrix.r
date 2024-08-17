cleanup.and.summarize.pairwise.matrix <- function(matrix) {
  matlist <- subset(reshape(melt(as.matrix(matrix)), value!=0))
  
  listReverse <- matlist[,c(2,1,3)]
  colnames(listReverse) <- c("Var1", "Var2", "value")
  fulllist <- rbind(matlist,listReverse)
  
  listbypop <- as.data.frame(fulllist) %>%
    dplyr::group_by(Var1) %>%
    dplyr::summarise(matlist = mean(value))
}