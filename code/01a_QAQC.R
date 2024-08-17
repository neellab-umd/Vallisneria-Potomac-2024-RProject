## Author: Maile Neel

## This code does QA/QC on the dataframes for the Potomac sites that are used in analyses of population genetic diversity. These genind objects were created in the script 01_Create_adegenet_geninds.R

## The tests follow procedures in Helen Wagner's Landscape Genetics course materials.  https://bookdown.org/hhwagner1/LandGenCourse_book/WE_3.html


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Install & load packages ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
## install missing packages that are found on CRAN, and then load the packages
pacman::p_load(devtools,
               adegenet,
               poppr,
               tidyverse,
               PopGenReport,
               here,
               ade4)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Load genind and genpop objects ######
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

list.of.genind.objects <- readRDS(here::here("processed.data","list.of.genind.objects.rds"))


list.of.genpop.objects <- readRDS(here::here("processed.data","list.of.genpop.objects.rds"))

adegenet::summary(list.of.genind.objects$genind.noreps.potomac)$NA.perc
adegenet::summary(list.of.genind.objects$genind.noreps.potomac.NonTidal)$NA.perc
adegenet::summary(list.of.genind.objects$genind.noreps.potomac.Tidal)$NA.perc
adegenet::summary(list.of.genind.objects$genind.noreps.atall)$NA.perc


gendata.reallynoreps.potomac<-gendata.allpotomac %>% 
  dplyr::distinct(Clone.ID.2018, .keep_all=TRUE) 
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#####    CODE FOR ALL SUMMARIES   ### 
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

##To get QA/QC summaries for the different data sets, change which genind is assigned to the object my.genind.

#my.genind <- list.of.genind.objects$genind.noreps.potomac
#my.genind <- list.of.genind.objects$genind.noreps.potomac.NonTidal
#my.genind <- list.of.genind.objects$genind.noreps.potomac.Tidal
#my.genind <- list.of.genind.objects$genind.noreps.atall


## Calculate allele frequencies for one rep of each MLG by population in adegenet

no.reps.allele.freqs.by.pop <- adegenet::makefreq(my.genind, quiet = FALSE, missing = NA,truenames = TRUE)

#export to csv.

readr::write_csv(as.data.frame(no.reps.allele.freqs.by.pop),"./processed.data/allele.freqs.by.site.no.reps.csv")

## Calculate HWE overall ----

## The values in the manuscript are based on the list.of.genind.objects$genind.noreps.potomac dataset.

overall.HWE <- round(pegas::hw.test(my.genind, B = 999), digits = 3) %>% 
  tibble::rownames_to_column("Locus")


### Write results to csv ----

readr::write_csv(as.data.frame(overall.HWE),"./processed.data/overall.hwe.csv")

## Results indicated that all but aagx012 and M16, and possibly aagx030 are out of HWE in the full data set.

## Check HWE by site ----

#seppop() splits the genind object by population. sapply() then applies the function 'hw.test' from package 'pegas' to each population.  'B=0' to specifies that we do no permutations. The function t() transposes the resulting matrix. To make this work, we have to use data.matrix() to temporarily interpret the data frame as a matrix. 
## Chi-squared test: p-value

### Chi Square Uncorrected ----

HWE.test.chisq <- sapply(seppop(my.genind), 
                         function(ls) pegas::hw.test(ls, B=10)[,3])

HWE.test.chisq.pvalues <- t(data.matrix(HWE.test.chisq)) %>% 
  round(.,3) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Pop") %>% 
  mutate(x2.Out.of.HWE := rowSums(across(is.numeric) < 0.01)/10)

### Chi Square Corrected with FDR ----
Chisq.fdr <- matrix(p.adjust(HWE.test.chisq,method="bonferroni"), nrow = nrow(HWE.test.chisq)) 
colnames(Chisq.fdr)<-colnames(HWE.test.chisq)
rownames(Chisq.fdr)<-rownames(HWE.test.chisq)

Chisq.fdr.pvalues <- t(data.matrix(Chisq.fdr)) %>% 
  round(.,3) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Pop") %>% 
  mutate(x2.fdr.Out.of.HWE := rowSums(across(is.numeric) < 0.01)/10)

### Monte Carlo Test Uncorrected----

HWE.test.MC <- sapply(seppop(my.genind), 
                      function(ls) pegas::hw.test(ls, B=999)[,4])
HWE.test.MC.pvalues <- t(data.matrix(HWE.test.MC)) %>% 
  round(.,3) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Pop") %>% 
  mutate(MC.Out.of.HWE := rowSums(across(is.numeric) < 0.01)/10)

### Monte Carlo Test Corrected with FDR ----

MC.fdr <- matrix(p.adjust(HWE.test.MC,method="bonferroni"), nrow = nrow(HWE.test.MC)) 
colnames(MC.fdr)<-colnames(HWE.test.MC)
rownames(MC.fdr)<-rownames(HWE.test.MC)

MC.fdr.pvalues <- t(data.matrix(MC.fdr)) %>% 
  round(.,3) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Pop") %>% 
  mutate(MC.fdr.Out.of.HWE := rowSums(across(is.numeric) < 0.01)/10)

### Total loci out of HWE by site and test ----

## Summarize the proportion of loci out of HWE in each site based on chi-square and MC values tests, with and without correction

Total.Prop.loci.out.of.HWE.by.Pop <- HWE.test.chisq.pvalues %>% 
  select(Pop, contains("Out.of.HWE")) %>%
  left_join(Chisq.fdr.pvalues %>% select(Pop, contains("Out.of.HWE")), by = "Pop") %>%
  left_join(HWE.test.MC.pvalues %>% select(Pop, contains("Out.of.HWE")), by = "Pop") %>%
  left_join(MC.fdr.pvalues %>% select(Pop, contains("Out.of.HWE")), by = "Pop")

Total.Prop.loci.out.of.HWE.by.Pop %>% 
  summarise(across(where(is.numeric), ~ sum(. > 0)))

### Write results to csv---- 

readr::write_csv(Total.Prop.loci.out.of.HWE.by.Pop, "./processed.data/proportion.of.loci.out.of.HWE.by.site.csv")

## Check by Locus ----

## For each locus, calculate the number of sites where it was out of HWE using alpha = 0.01 for each test, with and with out correction done above.

Number.sites.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq < alpha, 
                                                  1, 
                                                  sum), 
                                      MC = apply(HWE.test.MC < alpha, 
                                                 1, 
                                                 sum),
                                      Chisq.fdr = apply(Chisq.fdr < alpha, 
                                                        1, sum),
                                      MC.fdr = apply(MC.fdr < alpha, 
                                                     1, 
                                                     sum)) %>% 
  tibble::rownames_to_column("Locus")


### Write results to csv---- 

readr::write_csv(Number.sites.out.of.HWE, "./processed.data/number.of.sites.out.of.HWE.by.locus.csv")

  
### Test for linkage disequilibrium ---- 

## This test of LD follows https://bookdown.org/hhwagner1/LandGenCourse_book/WE_3.html
  
poppr::ia(my.genind, sample=199)
  
LD.pair <- poppr::pair.ia(my.genind) %>% 
    as.data.frame() %>% 
    mutate(rbarD.squared = rbarD^2) %>% 
    mutate(across(where(is.numeric), \(x) round(x, 3)))
  
##Effect size:  rbarD can be interpreted similarly to a linear correlation coefficient r. As in regular regression, square the rbarD to figure out how much variance is due to linkage for each locus. 
  
## Test for Null Alleles  ----

## from https://bookdown.org/hhwagner1/LandGenCourse_book/WE_3.html

## null.all() returns a lists with two components (‘homozygotes’ and ‘null.allele.freq’), and each of these is again a list. See ‘?null.all’ for details and choice of method. 

Null.alleles <- PopGenReport::null.all(my.genind)

## get the probability of homozygotes out of the list

Null.alleles$homozygotes$probability.obs
  
## get the null allele frequency summary out of the list  Each summary table contains a summary with observed, median, 2.5th percentile and 97.5the percentile. The percentiles form a 95% confidence interval. From the help file: “If the 95% confidence interval includes zero, it indicates that the frequency of null alleles at a locus does not significantly differ from zero.”

## Chakraborty et al. (1994)’s method (e.g. summary1) should be used if there are individuals with no bands at a locus seen, but they are discounted as possible artefacts. If all individuals have one or more bands at a locus then Brookfield (1996)’s method (e.g. summary2) should be used.

cat(" summary1 (Chakraborty et al. 1994):", "\n") 
    round(Null.alleles$null.allele.freq$summary1,2) 

## This summary option is better for our data - the other one has many distributions that are all negative. There is no discussion of what that means.

cat("summary2 (Brookfield et al. 1996):", "\n")
    round(Null.alleles$null.allele.freq$summary2,2)


