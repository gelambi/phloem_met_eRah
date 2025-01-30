# ======================================================================
# GOAL 1 
# Investigating Metabolite Differences in Tree of Heaven Plants
# Attractors vs Repellents of Spotted Lanternfly
# ======================================================================

# Clear the workspace
rm(list=ls())

# Load required libraries

#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("mzR")

library(erah)
library(Rcpp)
library(mzR)

# ======================================================================
# Load samples, all chromatograms without the Grob
# ======================================================================

setwd("~/Desktop/phloem_met_eRah/samples")
createdt("~/Desktop/phloem_met_eRah/samples") # This step creates the two .csv files 
ex_alk <- newExp(instrumental="~/Desktop/phloem_met_eRah/samples/samples_inst.csv",
                 phenotype = "~/Desktop/phloem_met_eRah/samples/samples_pheno.csv") # create a new experiment
metaData(ex_alk) # check that everything looks good
phenoData(ex_alk)

# ======================================================================
# DECONVOLUTION
# ======================================================================

ex.dec.par <- setDecPar(min.peak.width=1, min.peak.height=2500, ### Default parameters
                        noise.threshold=500, avoid.processing.mz=c(73:75, 147:149))
ex <- deconvolveComp(ex_alk, ex.dec.par)

ex.dec.par_0.7_5000_1000 <- setDecPar(min.peak.width=0.7, min.peak.height=5000, ### I think most of my peaks might be <1 second
                        noise.threshold=1000, avoid.processing.mz=c(73:75, 147:149))
ex <- deconvolveComp(ex_alk, ex.dec.par_0.7_5000_1000) # Warning messages: 1: In processSample(Experiment, .x, plotting, down.sample, virtualScansPerSecond) : Unable to extract factors from low/Q_102623_027.cdf. Data may be corrupted. 2: In processSample(Experiment, .x, plotting, down.sample, virtualScansPerSecond) :Unable to extract factors from standard/Q_102623_026.cdf. Data may be corrupted.

# Save
save(ex, file = "deconvolution_ALLsamples_NOgrob.rda")
# Load
load("deconvolution_ALLsamples_NOgrob.rda")

# ======================================================================
# ALIGNMENT, dafault dec paramenters
# ======================================================================

ex.al.par <- setAlPar(min.spectra.cor = 0.80,
                      max.time.dist = 3)
ex <- alignComp(ex, alParameters = ex.al.par, blocks.size = 30)
alignment <- alignList(ex, by.area=TRUE)
alignment ### 1905 obs 119 variables
write.csv(alignment, file = "alignment_eRah_norecovered.csv")

ex.al.par <- setAlPar(min.spectra.cor = 0.70,
                      max.time.dist = 5)
ex <- alignComp(ex, alParameters = ex.al.par, blocks.size = 30)
alignment <- alignList(ex, by.area=TRUE)
alignment ### 1792 obs 119 variables
write.csv(alignment, file = "alignment_eRah_norecovered.csv")

ex_alk <- recMissComp(ex, min.samples = 25, free.model = FALSE) # got a warning message. 
alignment_recovered_25 <- alignList(ex_alk, by.area=TRUE)

ex_alk <- recMissComp(ex, min.samples = 30, free.model = FALSE) # did not get the error, seems the same though 
alignment_recovered_30 <- alignList(ex_alk, by.area=TRUE)
write.csv(alignment_recovered_30, file = "alignment_ALLsamples_recovered.csv")

# ======================================================================
# LIBRARY SEARCH
# ======================================================================

ex <- identifyComp(ex_alk)
id.list <- idList(ex)
export2MSP(
  ex,
  export.id = NULL,
  id.database = mslib,
  store.path = getwd(),
  alg.version = 2
) ### In cbind(SpectNames.3, met.name) : number of rows of result is not a multiple of vector length (arg 2)

ex <- identifyComp(ex_alk)
id.list <- idList(ex)

### Using the Golm Metabolome Database

g.info <- "
GOLM Metabolome Database
------------------------
Kopka, J., Schauer, N., Krueger, S., Birkemeyer, C., Usadel, B., Bergmuller, E., Dor-
mann, P., Weckwerth, W., Gibon, Y., Stitt, M., Willmitzer, L., Fernie, A.R. and Stein-
hauser, D. (2005) GMD.CSB.DB: the Golm Metabolome Database, Bioinformatics, 21, 1635-
1638."

golm.database <- importGMD(filename="GMD_20111121_VAR5_ALK_MSP.txt", 
                           DB.name="GMD", 
                           DB.version="GMD_20111121", 
                           DB.info= g.info,type="VAR5.ALK")

# The library in R format can now be stored for a posterior faster loading
save(golm.database, file= "golmdatabase.rda")

load("golmdatabase.rda")
mslib <- golm.database
ex <- identifyComp(ex_alk)
id.list <- idList(ex) ### seems wrong, I am confident the big peak at 25 is sucrose. I might be using the wrong library

export2MSP(
  ex,
  export.id = NULL,
  id.database = mslib,
  store.path = getwd(),
  alg.version = 2
)
# compare spectra with NIST through the MassHunter workstation

# ======================================================================
# ALIGNMENT, modified dec paramenters
# ======================================================================

### new parameters: peak.width=0.7, min.peak.height=5000, noise.threshold=1000

save(ex, file = "deconvolution_ALLsamples_NOgrob_0,7.rda")
# Load
load("deconvolution_ALLsamples_NOgrob_0,7.rda")

# Alignment
ex.al.par <- setAlPar(min.spectra.cor = 0.90,
                      max.time.dist = 1)
ex <- alignComp(ex, alParameters = ex.al.par, blocks.size = 50) # I changed the 'block.size' parameter until I got the alignment. 
### Error in matrix(0, nrow = length(x$ID)) : data is too long

alignment <- alignList(ex, by.area=TRUE)
alignment ### 1905 obs 119 variables
write.csv(alignment, file = "alignment_eRah_norecovered.csv")

ex.al.par <- setAlPar(min.spectra.cor = 0.70,
                      max.time.dist = 5)
ex <- alignComp(ex, alParameters = ex.al.par, blocks.size = 30)
alignment <- alignList(ex, by.area=TRUE)
alignment ### 1792 obs 119 variables
write.csv(alignment, file = "alignment_eRah_norecovered.csv")

ex_alk <- recMissComp(ex, min.samples = 25, free.model = FALSE) # got a warning message. 
alignment_recovered_25 <- alignList(ex_alk, by.area=TRUE)

ex_alk <- recMissComp(ex, min.samples = 30, free.model = FALSE) # did not get the error, seems the same though 
alignment_recovered_30 <- alignList(ex_alk, by.area=TRUE)
write.csv(alignment_recovered_30, file = "alignment_ALLsamples_recovered.csv")

# ======================================================================
# SECOND LIBRARY SEARCH
# ======================================================================

ex <- identifyComp(ex_alk)
id.list <- idList(ex)
export2MSP(
  ex,
  export.id = NULL,
  id.database = mslib,
  store.path = getwd(),
  alg.version = 2
) ### In cbind(SpectNames.3, met.name) : number of rows of result is not a multiple of vector length (arg 2)

ex <- identifyComp(ex_alk)
id.list <- idList(ex)
