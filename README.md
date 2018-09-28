# BEEM-{Static}

## Description

BEEM-{Static} is an extension of the original [BEEM](https://github.com/lch14forever/BEEM) algorithm to handle cross-sectional microbial profiling data with an assumption that large fraction of the samples are at their steady states. As similar to BEEM, it simultaneously estimates the biomass for each samples and the generalized Lotka-Volterra model parameters for all the species in the data.

## Installation

```r
devtools::install('lch14forever/beem_static')
```

## Example usage

### Simulated data

The demo dataset is a simulated community of 20 species and 500 samples. All of the samples are at the equilibrium and 70% habitat preference (each sample contains 70% of all the species randomly). We add some Poisson noise to the data.

```r
df <- read.table('data/demo.n500_p20_eq1.counts.txt', head=F, row.names=1, sep='\t')
## add some poisson noise to the data
aux.poisson <- function(v, depth){
    v <- round(v/sum(v) * depth)
    sapply(v, function(x) rpois(1, x))
}
df.noise <- apply(df, 2, aux.poisson, depth=5000)
```

### Run BEEM

```r
res <- func.EM(df.noise, ncpu=4, scaling=median(biomass.true))
```

### Estimated biomass vs. True biomass

```r
plot(res$trace.m[,30], biomass.true, xlab='BEEM biomass estimation', ylab='True biomass')
```
![](vignettes/biomass_compare.png)


## Citation

A manuscript for BEEM-{Static} is in preparation and please cite BEEM if you use this package:

 - C Li, et al. (2018) An expectation-maximization-like algorithm enables accurate ecological modeling using longitudinal metagenome sequencing data. [BioRxiv](https://www.biorxiv.org/content/early/2018/07/17/288803)
