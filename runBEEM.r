#!/mnt/software/bin/Rscript-3.4.1

args <- commandArgs(trailingOnly=T)

source('/mnt/projects/lich/dream_challenge/BEEM_STATIC/dev/em_functions.r')


folder <- '/mnt/projects/lich/dream_challenge/BEEM_STATIC/simulation/at_eq/'

##for sim in
sim <- args[1] ##'sim11.n500_p20_eq1'
subsample <- as.numeric(args[2])

base <- paste0(folder, sim)

df <- read.table(paste0(base,'.counts.noise.pois_nrep1.txt'), head=F, row.names=1, sep='\t')

## truth to get biomass
df.noNoise <- read.table(paste0(base,'.counts.txt'), head=F, row.names=1, sep='\t')
m.true <- colSums(df.noNoise)

if (subsample == 0){
    subsample <- 500
}

set.seed(0)
idx <- sort(sample(1:ncol(df), subsample))
dat <- df[,idx]


ncpu=10
max.iter <- 40
scaling <- median(m.true)
dev <- Inf

res <- func.EM(dat, ncpu=ncpu, scaling=scaling, dev=dev, max.iter=max.iter)
res$subsample <- idx

saveRDS(res, paste0(base, '.sub', subsample, '.beem.rds'))
