#!/mnt/software/bin/Rscript-3.1.2
rm (list = ls()) 
base <- './'
parameters <- read.table(paste0(base,"BVS.results.parameters.txt"), head = T)
## filter out values with Bayes factor < 10
parameters[parameters$significance < 10, 4] <- 0

## exploratory: paramter distribution?
### growth rate
library(fitdistrplus)
gr <- parameters[parameters$parameter_type=='growth_rate',]$value
self.int <- with(parameters, parameters[as.character(source_taxon)==as.character(target_taxon) & parameter_type == 'interaction',]$value )
int <- with(parameters, parameters[as.character(source_taxon)!=as.character(target_taxon) & parameter_type == 'interaction',]$value )

## initial abundances getting from mdsine output


save.image(file="dat.gen.model.RData")

#############################################################################
############################ End of preprocessing############################
#############################################################################

