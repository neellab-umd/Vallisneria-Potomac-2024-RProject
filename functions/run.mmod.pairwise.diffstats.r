
genind<-list.of.genind.objects$genind.allpotomac

test.diff.list <- map(list.of.genind.objects, ~run.mmod.pairwise.diffstats(.x))

run.mmod.pairwise.diffstats(genind)
run.mmod.pairwise.diffstats <- function(genind) {
  matrixD.Jost <- mmod::pairwise_D(genind)
  
  #generate pairwise Hedrick's G''st
  matrixGst.Hedrick <- mmod::pairwise_Gst_Hedrick(genind)
  
  matrixFst <- genet.dist(genind,method = "WC84")
 
  diff.list<-list(matrixD.Jost,matrixGst.Hedrick,matrixFst)
}

list.of.distance.matrices<-list("matrixD.Jost", "matrixGst.Hedrick","matrixFst")

cleanup.and.summarize.pairwise.matrix(matrixD.Jost)



diveRsity_diff <- cbind(WC84bypopdiv, as.data.frame(Gprimebypop$GprimeSt),JostsDbypop$JostsD)
colnames(diveRsity_diff) <- c("Pop", "ThetaDiv", "G'ST", "JostsD")

