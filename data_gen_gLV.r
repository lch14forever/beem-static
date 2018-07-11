#!/mnt/software/bin/Rscript-3.1.2
load('/mnt/projects/lich/dream_challenge/BEEM_STATIC/dev/dat.gen.model.RData')

suppressMessages(library(deSolve))
suppressMessages(library(reshape2))
args <- commandArgs(trailingOnly=T)
## Number of species
p <- as.numeric(args[1])
## number of individuals
n <- as.numeric(args[2]) ##10
## output folder
out.folder <- args[3]
dir.create(args[3], showWarnings = FALSE)

thres_to_eq <- c(0.2, 0.5) ## > 20% difference, > 10% species
samples_at_equil <- 1*n ##round(0.8 * n)
message(paste0("Number of samples at equilibrium:", samples_at_equil))

#### GLVM specifications #####
## model
lvm <- function(t,x,parms){
    with(as.list(parms,c(x)),{
             if (t >= perturb.s && t <= perturb.e ) {miu=1} else {miu=0}
             x[x < 0.00000001] <- 0
             dx <- alpha.v * x +
                   x * (beta.m %*% x) +
                   gamma * x * miu
             list(dx)
         })
}

n.integrate <- function(time,init.x,model, parms){
    t.out <- seq(time$start,time$end,by=time$steps)
    dat <- as.data.frame(ode(init.x,t.out,lvm,parms))
    return (dat)
}


gr <- gr
## scale interaction such that the scale is not too different
sf <- abs(mean(gr) /mean(self.int))
self.int <- self.int * sf
int <- int[int!=0] * sf


print("Generating model parameters...")
## growth rate
tmp <- abs(rnorm(10*p, mean(gr), sd(gr)))
tmp <- tmp[abs(tmp - mean(gr)) < 5*sd(gr)]
alpha.v <- sample(tmp,p)
gamma.v <- rep(0,p)
perturb <- list(start=1000, end=1000)

## interaction matrix, implemented according to MDSINE's supplement
beta.m <- matrix(0,p,p)
for(j in 1:p){
    temp <- rep(0,p)
    temp[j] <- NA
    temp[!is.na(temp)] <- rbinom(p-1, 1, 0.2)
    idx <- which(temp==1)
    temp[idx] <- rnorm(length(idx), 0, sd(int)/length(idx))
    beta.m[j,] <- temp
}
diag(beta.m) <- -abs(rnorm(p, mean(self.int), sd(self.int)))

## parms
parms <- list(alpha.v=alpha.v, beta.m=beta.m, perturb.s=perturb$start, perturb.e=perturb$end, gamma=gamma.v )

print("Generating simulated dataset...")

counts <- NULL
count <- 0

df <- read.table('/mnt/projects/lich/dream_challenge/BEEM_STATIC/data/curatedMetagenomicData_stool_healthy.v1.7.92.table.s',head=F,row.names=1)

## rank by prevalence
df[df < 0.1] <- 0
df.sorted <- df[order(rowSums(df==0)), ]
## take top p 
df.selected <- df.sorted[1:p,]
## columns with more than 90% 0s are filtered out
df.selected <- df.selected[,colSums(df.selected==0)/p<0.9]
## sparsity 
abs_perc <- as.numeric(colSums(df.selected == 0)/p)
print(summary(abs_perc))

err.to.equil <- NULL

for (i in 1:(1000*n)){
    init.x <- abs(runif(p, 0.001, 30))
    time <- list(start = 0, end = 500, steps = 1)
    ##species to exclude set to 0
    init.x[rbinom(p,1, mean(abs_perc))==1] <- 0
    names(init.x) <- sprintf("sp%03d", 1:p)##c(paste0("sp",1:p))
    while(TRUE){
        out <- n.integrate(time=time,init.x=init.x,model=lvm, parms=parms)
        if( all(abs(out[(nrow(out)),-1]- out[(nrow(out)-1),-1]) < 1e-5) ) break
        time$end <- time$end + 100
    }
    ## if(length(warnings())>0){ # or !is.null(warnings())
    ##     break
    ## }
    end <- out[nrow(out),-1]
    if(count > samples_at_equil){
        ## find the samples satisify this criteria
        tmp <- abs(t(out[,-1]) - as.numeric(end))/as.numeric(end) ## relative error
        idx <- which(apply(tmp>thres_to_eq[1], 2, function(x) sum(x, na.rm = TRUE)/sum(!is.na(x))) > thres_to_eq[2])
        if(length(idx[idx>5])==0){
            message("Skipping")
            next
        }
        ## take a random sample from day 6 to threshold
        res <- out[sample(idx,1),-1]
    }else{
        res <- end
    }    
    end[end < 1e-4] <- 0 
    res[res < 1e-4] <- 0 ## assumed detection limit
    tmp <- as.numeric((res-end)/end)
    tmp[is.na(tmp)] <- 0 ## 0/0 --> no error here
    tmp[is.infinite(tmp)] <- NA ## num/0 --> max error !! Not reaching 0 yet

    err.to.equil <- rbind(err.to.equil, tmp)
    if(any(is.na(res))) {next}
    if( sum(res > 0) < 5 ) {next}
    if(count %%20 ==0) {print(count)}
    count = 1+count
    counts <- rbind(counts, res)
    if (count >= n) break ## generate n replicates
}

print("Writing to output files...")

colnames(counts) <- sprintf("sp_%03d", 1:p)##c(paste0("sp",1:p))
biomass <- data.frame(mass=rowSums(counts))
counts <- rbind(1:nrow(counts), t(counts))
rownames(counts)[1] <- "#OTU ID"
## count table
write.table(counts, paste0(out.folder,  "/counts.txt"), col.names=F, row.names=T, quote=F, sep = "\t")
## error to equilibrium
write.table(t(err.to.equil), paste0(out.folder,  "/err_to_equil.txt"), col.names=F, row.names=F, quote=F, sep = "\t")

## write.table(biomass,  paste0(out.folder, "/biomass.txt"), col.names = T, row.names = F, quote =F, sep = "\t")

## write paramters
parameters.out <- beta.m
rownames(parameters.out) <- colnames(parameters.out) <- sprintf("sp_%03d", 1:p)
parameters.out <- melt(parameters.out)
parameters.out <- rbind(data.frame(Var1=sprintf("sp_%03d", 1:p),Var2=NA,value=alpha.v), parameters.out)[,c(2,1,3)]

colnames(parameters.out)[1:2] <- c('source_taxon', 'target_taxon')

write.table(parameters.out, paste0(out.folder, "/truth.txt"), quote=F, col.names = T, row.names=F , sep ='\t')

