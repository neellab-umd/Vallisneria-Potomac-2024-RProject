## This script uses poppr to generate minimum spanning networks for MLGs.  

## Load Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               here,
               poppr,
               igraph,
               RColorBrewer,
               update = FALSE)
## Read in rds file with genind objects

list.of.genind.objects<-readRDS(here("processed.data", "list.of.genind.objects.rds"))

## If the list.of.genind.objects.rds file is not in the processed.data.folder, you need to create it with the script 01.2018_PR_Manuscript_adegenet_create.R.

## Set color Ramp
pcpal <- colorRampPalette(c("blue", "gold"))
set.seed(9001)

## Create minimum spanning network ----

## sortorder<-sort(as.numeric(popNames(allpotomac.genind)))

allpotomac.gc <- as.genclone(list.of.genind.objects[[1]])
setPop(allpotomac.gc) <- ~(OrderPop)

allpotomac.gc.nomissing<-missingno(allpotomac.gc, type = "geno", cutoff = 0.0, quiet = FALSE, freq = FALSE)

## replen is the number of repeats for our microsats
allpotomac.msn <- bruvo.msn(allpotomac.gc.nomissing, replen = rep(3, nLoc(allpotomac.gc.nomissing)), showplot = FALSE)

## ordering does not work - Still need to figure it out.
allordered<- as.data.frame(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "23", "24","25", "26","27", "28", "29","30","31", "32", "33", "34"))

set.seed(120)

png(file= "./figures/msn.all.potomac.png",width=600, height=350)
plot_poppr_msn(allpotomac.gc.nomissing, allpotomac.msn, inds = "none", palette = pcpal, nodescale = 7.5)
dev.off()


## NonTidal Only MSN ----

subnontidalordered<- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "23", "24")

nontidal.msn <- bruvo.msn(allpotomac.gc.nomissing, replen = rep(3, nLoc(allpotomac.gc.nomissing)), sublist = subnontidalordered, showplot = FALSE)

set.seed(120)
png(file= "./figures/msn.nontidal.potomac.png",width=700, height=350)
plot_poppr_msn(allpotomac.gc.nomissing, nontidal.msn, inds = "none", palette = pcpal, nodescale = 7.5)
dev.off()

## Tidal Only MSN ----

subtidalordered<- c("25", "26","27", ".28", "29","30","31", "32", "33", "34", "35","36")
tidal.msn <- bruvo.msn(allpotomac.gc.nomissing, replen = rep(3, nLoc(allpotomac.gc.nomissing)), sublist =as.vector(subtidalordered), showplot = FALSE)

set.seed(120)
png(file= "./figures/msn.tidal.potomac.png",width=500, height=350)
plot_poppr_msn(allpotomac.gc.nomissing, tidal.msn, inds = "none", palette = pcpal, nodescale = 7.5)
dev.off()
