#!/mnt/software/bin/Rscript-3.4.1
library("optparse")
option_list <- list( 
    make_option(c("-i", "--input"), type='character', dest='sim', 
                help="Input prefix"),
    make_option(c("-f", "--folder"), type='character', dest='folder',
                default='/mnt/projects/lich/dream_challenge/BEEM_STATIC/simulation/',
                help="Input folder"),
    make_option(c("-s", "--subsample"), type="integer", dest='subsample',
                default=0,
                help="Subsample from the data"),
    make_option(c("-d", "--dev"), type="double", dest='dev',
                default=Inf,
                help="Deviation for filtering"),  
    make_option("--ncpu", type="integer", default=10,
                help="Number of CPUs",
                dest="ncpu"),
    make_option("--maxIter", type="integer", default=30,
                help="Number of iterations",
                dest="maxIter")    
    )
args <- parse_args(OptionParser(option_list=option_list))

subsample <- args$subsample

if(is.null(args$sim)){
    stop("No input prefix specified!")
}

ncpu <- args$ncpu
max.iter <- args$maxIter
dev <- args$dev

source('/mnt/projects/lich/dream_challenge/BEEM_STATIC/dev/em_functions.r')


base <- paste0(args$folder, args$sim)

df <- read.table(paste0(base,'.counts.noise.pois_nrep1.txt'), head=F, row.names=1, sep='\t')

## truth to get biomass
df.noNoise <- read.table(paste0(base,'.counts.txt'), head=F, row.names=1, sep='\t')
m.true <- colSums(df.noNoise)
scaling <- median(m.true)

set.seed(0)
if(subsample == 0){
    subsample <- ncol(df)
}
idx <- sort(sample(1:ncol(df), subsample))
dat <- df[,idx]

res <- func.EM(dat, ncpu=ncpu, scaling=scaling, dev=dev, max.iter=max.iter)
res$subsample <- idx

saveRDS(res, paste0(base, '.sub', subsample, '.beem.rds'))
