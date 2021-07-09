library(GenomicRanges)
library(rtracklayer)
library(Gviz)
library(GenomicInteractions)
library(biomaRt)
library(Rsamtools)
library(Chicago)

indexBam("/home/master/ra")

genomeAxis <- GenomeAxisTrack(name="MyAxis")
genomeAxis
MALAT1 <- ibed[ibed$bait_name=="EHBP1"]
range <- MALAT1[,c(3,4,5,13)]
#range <- range[range$target.chr=="chr11"]
range <- GRanges(range)

acc <- DataTrack(range,chromosome = "chr2",genome = "hg38")
acc
plotTracks(c(acc,genomeAxis),from = MALAT1$bait.start[1], to = MALAT1$bait.end[1],chromosome = "chr1")

locus <- list(chr = 'chr2', start = 62673851, end = 63046487);
locus.id <- paste0(locus$chr, ':', locus$start, '-', locus$end);

ideogram.track <- IdeogramTrack(genome = 'hg38', chromosome = locus$chr);
genome.axis.track <- GenomeAxisTrack();
intrxn <- ibed[,1:13]
# get counts and interactions involving specific locus
# TO DO: need to switch this to generic other end in case of b2b
# add note about p-value/ q-value
locus.counts <- intrxn[ target.id == locus.id | bait.id == locus.id ];
locus.interactions <- ibed[ p.value < 0.001 & (target.id == locus.id | bait.id == locus.id) ];


# create track for visualising read count
interaction.count.track <- DataTrack(
  range = GRanges(locus.counts$target.id),
  data = locus.counts$count,
  name = 'Reads'
);

# create track for visualising significant interactions
interaction.object <- GenomicInteractions(
  anchor1 = GRanges( locus.interactions$bait.id ),
  anchor2 = GRanges( locus.interactions$target.id ),
  count = -log10(locus.interactions$p.value)
);

interaction.track <- InteractionTrack(
  interaction.object,
  name = 'Chicane'
);

gene.track <- GeneRegionTrack(
  "/home/shashankt693/hg38.gtf",
  chr = locus$chr,
  start = locus$start - 5e5,
  end = locus$end + 7.5e5,
  stacking = 'squish',
  stackHeight = 0.3,
  name = 'Genes'
);

# plot
plotTracks(
  list(
    ideogram.track, 
    genome.axis.track, 
    gene.track,
    interaction.count.track,
    interaction.track
  ),
  sizes = c(0.4, 1, 1, 4, 3),
  type = 'histogram',
  transcriptAnnotation = 'symbol',
  collapseTranscripts = 'longest',
  col = NULL,
  from = locus$start - 6e5,
  to = locus$end + 8e5
)
  
peakReads <- AlignmentsTrack("/home/master/shared_folder/radicl_seq/Data/Processed/D13.bam")
peakReads
plotTracks(peakReads,
           chromosome="chr5",
           from=135312577,
           to=135314146)

##get score functions

ibed$dist <- 0.5*(ibed$bait.end+ibed$bait.start)-0.5*(ibed$target.start+ibed$target.end)
ibed$dist <- ifelse(ibed$target.chr==ibed$bait.chr,ibed$dist,NA)


cis <- ibed[!is.na(ibed$distance)]
q.count <- quantile(cis$count,c(0.30,0.60,0.90,0.95))
cis$q.count <- with(cis,ifelse(count<q.count[1],0,
                                   ifelse(count<q.count[2],1,
                                          ifelse(count<q.count[3],2,
                                                 ifelse(count<q.count[4],3,4)))))

chic_signif_cis <- chic_signif[!is.na(chic_signif$distance)]
q.count <- quantile(chic_signif_cis$count,c(0.30,0.60,0.90,0.95))
chic_signif_cis$q.count <- with(chic_signif_cis,ifelse(count<q.count[1],0,
                               ifelse(count<q.count[2],1,
                                      ifelse(count<q.count[3],2,
                                             ifelse(count<q.count[4],3,4)))))
ggplot(chic_signif_cis %>% arrange(-q.count),aes(x=distance,y=score,size=q.value))+
  geom_point(aes(color=as.factor(q.count)),alpha=0.3)+
  scale_colour_manual(name="percentile of counts",
                      values=c("#69b3a2","black","darkgoldenrod4","coral1","red"),
                      labels=c("0-30","30-60","60-90","90-95","95-100"))+   
  theme_cowplot()


## Correlation
cor_test <- cor.test(ibed[is.na(ibed$distance)]$score,ibed[is.na(ibed$distance)]$q.value,method = "spearman",exact = F)
plt_text <- paste0('r = ', round(cor_test$estimate,3), ', ',
                   paste('p-value < 2.2e-16'))#taken from res of cor_test
setorder(ibed[is.na(ibed$distance)],distance)
ggplot(data=ibed[is.na(ibed$distance)],aes(x = score, 
                                           y = log(q.value,base=10),colour=distance))+ 
  geom_point(alpha=0.2,size=0.1) + 
  annotate('text', x = -Inf, y = Inf, hjust = -.1, vjust = 1  ,label=plt_text)+
  labs(y= "Log 10 P-Value", x = "Log 10 Interaction Counts(Trans)")  +
  theme_cowplot()+
  ggtitle("Trans inetractions")+
  scale_color_viridis()

cor_test <- cor.test(ibed[!is.na(ibed$distance)]$score,ibed[!is.na(ibed$distance)]$q.value,method = "spearman",exact = F)
plt_text <- paste0('r = ', round(cor_test$estimate,3), ', ',
                   paste('p-value < 2.2e-16'))#taken from res of cor_test
setorder(ibed[!is.na(ibed$distance)],distance)
ggplot(ibed[!is.na(ibed$distance)] %>% arrange(distance),aes(x = score, 
                                                             y = log(q.value,base=10),colour=distance))+ 
  geom_point(alpha=0.2,size=0.1) + 
  annotate('text', x = -Inf, y = Inf, hjust = -.1, vjust = 1  ,label=plt_text)+
  labs(y= "Log 10 P-Value", x = "Log 10 Interaction Counts(Cis)")  +
  theme_cowplot()+
  ggtitle("Cis inetractions")+
  scale_color_viridis()

getScores <- function(ibed)
  {
  ## - If method="weightedRelative", we divide by weights (Genovese et al 2006)
  ##Note to self: Algebra is on P96 of my lab notebook.
  
  avgFragLen <- .getAvgFragLength(ibed) ##average fragment length
  
  
    eta.bar <- .getEtaBar(ibed)
    ibed$log.p <- log10(ibed$p.value)
    
    ##Get weights, weight p-values
    message("Calculating p-value weights...")
    ibed[, log.w:= .getWeights(abs(ibed$distance), ibed, eta.bar=eta.bar)]
    ibed[, log.q:= log.p - log.w] ##weighted p-val
    message("Calculating scores...")
    
    ##get score (more interpretable than log.q)
    minval <- .getWeights(0, ibed, eta.bar=eta.bar) ##FIXME could be optimized a *lot*.
    ibed[,score := pmax(- minval - log.q, 0)]
    
  ibed
}

.getAvgFragLength <- function(ibed)
{
  ##Normally, takes rmapfile from cd object.
  ##However, if rmapfile is specified, this overrides cd.
  
  rmap <- unique(ibed[,3:5])
  chrMax <- rmap[,max(target.end),by="target.chr"] ##length of each chr
  sum(as.numeric(chrMax$V1))/nrow(rmap)
}

.getNoOfHypotheses <- function(ibed)
{
  
  ##How many hypotheses are we testing? (algebra on p246 of JMC's lab notebook)
  rmap <- unique(ibed[,3:5])
  chrMax <- rmap[,max(target.end),by="target.chr"] ##length of each chr
  baitmap <- unique(ibed[,6:8])
  nBaits <- table(baitmap$bait.chr)
  avgFragLen <- .getAvgFragLength(ibed) ##average fragment length
  
  chr <- chrMax$target.chr
  ##count # of hypotheses
  Nhyp <- sum(nBaits)*(2*nrow(rmap) - sum(nBaits) - 1)/2L ##number of hypotheses being tested
  Nhyp
}
set <- example@settings
.getEtaBar <- function(cd, includeTrans=TRUE)
{
  
  ##1. Collect parameters
  alpha = set$weightAlpha
  beta = set$weightBeta
  gamma = set$weightGamma
  delta = set$weightDelta
  
  ##2. Get genomic/fragment map information
  rmap <- unique(ibed[,3:5])
  chrMax <- rmap[,max(target.end),by="target.chr"] ##length of each chr
  baitmap <- unique(ibed[,6:8])
  nBaits <- table(baitmap$bait.chr)
  
  
  chr <- as.character(names(nBaits))
  
  avgFragLen <- .getAvgFragLength(ibed)
  
  ##count # of hypotheses
  Nhyp <- .getNoOfHypotheses(ibed)
  
  ##3. Calculate eta.bar
  ##Loop, summing contributions of eta
  
  eta.sigma <- 0 
  for(c in chr)
  {
    ##length of chromosome
    d.c <- as.numeric(chrMax$V1[chrMax$target.chr == c])
    ##no of baits on chromosome
    n.c <- nBaits[c]
    
    for(i in 1:n.c) ##TODO nested for loop can be replaced with lapply & function
    {
      d = d.c*i/n.c
      d.near = min(d, d.c-d) ##dist to nearest chromosome end
      d.other <- seq(from=avgFragLen, to=max(avgFragLen,d.near), by = avgFragLen) ##locations of fragments
      d.other2 <- seq(from=d.near, to=d.c-d.near, by = avgFragLen)
      eta.sigma <- eta.sigma + 2*sum(expit(alpha + beta*log(d.other))) + sum(expit(alpha + beta*log(d.other2)))
      #eta[i] <- 2*sum(expit(alpha + beta*log(d.other))) + sum(expit(alpha + beta*log(d.other2)))
    }
  }
  
  eta.bar <- eta.sigma/Nhyp
  eta.bar
}

.getWeights <- function(dist,ibed,eta.bar)
{
  alpha = set$weightAlpha
  beta = set$weightBeta
  gamma = set$weightGamma
  delta = set$weightDelta
  
  ##4. Calculate weights
  eta <- expit(alpha + beta*log(naToInf(dist)))
  log.w <- log((expit(delta) - expit(gamma))*eta + expit(gamma)) -
    log((expit(delta) - expit(gamma))*eta.bar + expit(gamma))
  
  log.w
}




## Misc functions ---------------------------

wb2b = function(oeID, s, baitmap=NULL){
  # s is the current chicagoData object's settings list
  if (is.null(baitmap)){
    baitmap = .readBaitmap(s)
  }
  which(oeID %in% baitmap[[s$baitmapFragIDcol]])
}

whichbait2bait = function(x, baitmap=NULL){
  stop("whichbait2bait is deprecated. Use wb2b instead")
}

as.dataTableList <- function(cd){ 
  # takes a list of chicagoData objects cd and returns a list of data tables 
  # from the respective cd@x slots
  lapply(cd, function(cdi)cdi@x)
}

geo_mean <- function(data){    
  log_data <- log(data);    
  gm <- exp(mean(log_data[is.finite(log_data)]));    
  return(gm) 
} # http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in

logit <- function(p){log(p/(1-p))}

expit <- function(x){1/(1+exp(-x))}

removeNAs <- function(x) {x[!is.na(x)]}

naToInf <- function(x)
{
  ifelse(is.na(x), Inf, x) ##Convert NAs to infs.
}

ifnotnull = function(var, res){ if(!is.null(var)){res}}

locateFile = function(what, where, pattern){
  message("Locating ", what, " in ", where, "...")
  filename = list.files(where, pattern)
  
  if (length(filename)!=1){
    stop(paste0("Could not unambigously locate a ", what, " in ", where, ". Please specify explicitly in settings\n"))
  }
  
  message("Found ", filename)
  file.path(where, filename)
}


estimateDistFun <- function (cd, method="cubic", plot=TRUE, outfile=NULL) {
  
  # Take the "refBinMean" column of the data x as f(d_b)
  # then interpolate & extrapolate to get f(d).
  # TODO output extra diagnostic information?
  
  if (!method == "cubic"){
    stop ("Unknown method.\n")
  }
  
  # Get f(d_b)
  setkey(cd@x, distbin)
  f.d <- unique(cd@x, by=key(cd@x))[is.na(refBinMean)==FALSE][, c("distbin", "refBinMean"), with=FALSE]

  setDF(f.d) # f.d is tiny, so no need to bother with it being a data.table
  f.d$midpoint <- seq(from=round(cd@settings$binsize/2), by=cd@settings$binsize, length.out=nrow(f.d))
  
  obs.min <- log(min(f.d$midpoint))
  obs.max <- log(max(f.d$midpoint))
  
  #   if(method == "lm") {
  #
  #     ##previously had arguments: n.obs.head=10, n.obs.tail=25
  #
  #     ##On log-scale, do a linear interpolation.
  #     ##Linear models applied to first (n.obs.head) observations, and to last (n.obs.tail) observations.
  #     
  #     ##Interpolation: Estimate f(d) (NB "rule" parameter = 1 forces NAs outside of range)
  #     log.f.obs <- approxfun(log(f.d$midpoint), log(f.d$refBinMean), rule=c(1,1))
  #     
  #     ##Extrapolation: Fit the "head" and "tail" of f using a linear model
  #     head.coef <- coefficients(lm(log(refBinMean)~log(midpoint), data = head(f.d, n.obs.head))) ##Fit for small d
  #     tail.coef <- coefficients(lm(log(refBinMean)~log(midpoint), data = tail(f.d, n.obs.tail))) ##Fit for large d
  #     
  #   }
  
  if(method == "cubic") {
    ##Spline - Cubic fit over observed interval, linear fit elsewhere, assume continuity of f(d) & f'(d).
    distFunParams <- list(method="cubic")
    
    ##cubic fit (quadratic not immensely different TBH)
    f.d.cubic <- lm(log(refBinMean) ~ log(midpoint) + I(log(midpoint)^2) + I(log(midpoint)^3), data = f.d)
    fit <- f.d.cubic$coefficients
    distFunParams[["cubicFit"]] <- fit
    
    ##Extrapolation: Fit the "head" and "tail" of f using continuity
    distFunParams[["obs.min"]] <- log(min(f.d$midpoint))
    distFunParams[["obs.max"]] <- log(max(f.d$midpoint))
    
    beta <- fit[2] + 2*fit[3]*c(obs.min, obs.max) + 3*fit[4]*(c(obs.min, obs.max)^2)
    alpha <- fit[1] + (fit[2] - beta)*c(obs.min, obs.max) + fit[3]*c(obs.min, obs.max)^2 + fit[4]*c(obs.min, obs.max)^3
    
    distFunParams[["head.coef"]] <- c(alpha[1], beta[1])
    distFunParams[["tail.coef"]] <- c(alpha[2], beta[2])
    
  }
  
  cd@params$distFunParams <- distFunParams
  
  if(plot)
  {
    if (!is.null(outfile)){ 
      pdf(outfile)
    }
    
    plotDistFun(cd)
#     my.log.d <- seq(from=obs.min, to=obs.max, length.out = 101)
#     my.d <- exp(my.log.d)
#     plot(my.log.d, log(.distFun(my.d, distFunParams)),
#          type="l",
#          main = "Distance function estimate",
#          xlab = "log distance",
#          ylab = "log f",
#          col = "Red")
     with(f.d, points(log(midpoint), log(refBinMean)))
     legend("topright", legend = c("Data", "Fit"), col = c("Black", "Red"), pch = c(1, NA), lty=c(0,1))
    if (!is.null(outfile)){
      dev.off()
    }
  }
  
  cd
}

.distFun <- function(d, distFunParams)
{
  ##d: distSign column of cd@x
  
  ##distFunParams: list of form
  ##list(method, ...)
  ##so list("cubic", head.coef, tail.coef, obs.min, obs.max, fit)
  ##in future, could support e.g. list("lm", ...)
  
  p <- distFunParams
  
  stopifnot(p[["method"]] == c("cubic")) ##compatibility with future "method"s
  
  ##Transfer parameters
  obs.max <- p[["obs.max"]]
  obs.min <- p[["obs.min"]]  
  head.coef <- p[["head.coef"]]
  tail.coef <- p[["tail.coef"]]
  fit <- p[["cubicFit"]]
  
  ##Put everything together to get the final function
  d <- log(d)
  
  out <- ifelse(d > obs.max,
                tail.coef[1] + d*tail.coef[2], ##Common case evaluated first (large d)
                ifelse(d < obs.min,
                       head.coef[1] + d*head.coef[2],
                       fit[1] + fit[2]*d + fit[3]*(d^2) + fit[4]*(d^3) ##2nd most common case evaluated second (small d)
                )
  )
  
  exp(out)
}




genes <- unique(ibed$bait_name)
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol",
    "entrezgene_id",
    "ensembl_gene_id",
    "gene_biotype"),
  filter = "hgnc_symbol",
  values = genes,
  uniqueRows=TRUE,
  useCache = FALSE)

annotLookup <- as.data.table(annotLookup)

#gene_biotype
setnames(annotLookup,
         c("bait_name","entrezgene_id","ensembl_gene_id","gene_biotype"))