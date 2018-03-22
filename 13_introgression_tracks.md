# HMM

I wrote a R script to do a HMM based on Dstats in windows:

```R
# adapted from http://web.stanford.edu/class/stats366/hmmR2.html
setwd('/projects/2018_HMM')
library("HMM")
# read the data
dat<-read.table("./H3_maura_H1_tog_H2_tonk.concat",header=TRUE)
# reset the max print option so I can see the enture vector
options(max.print=1000000)

# make a histogram of D values
m<-mean(dat$D, na.rm = TRUE)
std<-sqrt(var(dat$D,na.rm = TRUE))
#qplot(dat$D, geom="histogram", binwidth = 0.1) 
hist(dat$D, density=0, breaks=20, prob=TRUE, 
     xlab="D", ylim=c(0, 2), 
     main="normal curve over histogram") 
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")




# set all states to zero
dat$D_state <- "zero"
# update the ones that have a substantial value for D
dat$D_state[which(!is.na(dat$D) & !is.na(dat$dH2H3) & dat$D > 0.8)] <- "pos"
dat$D_state[which(!is.na(dat$D) & !is.na(dat$dH1H3) & dat$D < -0.8)] <- "neg"


# this may be used by similations to simulate the same number of observations
nSim = length(dat$D_state)
# We will have three states
States = c("BABA", "BBAA", "ABBA")
# Symbols are the observations, e.g. the value of D 
# could be D<-0.1, -0.1>D>0.1, or D>0.1
Symbols = matrix(c("neg","zero","pos"))
# set up a 3x3 matrix of transition probabilities
# These are the probabilities that each hidden state changes to another hidden state
transProbs = matrix(c(0.9, 0.05, 0.05, 0.05, 0.9, 0.05, 0.05, 0.05, 0.9), c(length(States),
                                                 length(States)), byrow = TRUE)
# transProbs = matrix(c(0.01, 0.98, 0.01, 0.01, 0.98, 0.01, 0.01, 0.98, 0.01), c(length(States),
#                                                 length(States)), byrow = TRUE)
# not enough mixing
#transProbs = matrix(c(0.005, 0.99, 0.005, 0.005, 0.99, 0.005, 0.005, 0.99, 0.005), c(length(States),
#                                                 length(States)), byrow = TRUE)
# not enough mixing

#transProbs = matrix(c(0.002, 0.996, 0.002, 0.002, 0.996, 0.002, 0.002, 0.996, 0.002), c(length(States),
#                                                 length(States)), byrow = TRUE)
# not enough mixing


# define emmission probabilities
# These are the probabilities that each observed state is emmitted in each hidden state
# This allows, for example, a path within an introgressed region to not look 
# introgressed (e.g. if there is no variation)
#emissionProbs = matrix(c(0.9, 0.05, 0.05, 0.05, 0.9, 0.05, 0.05, 0.05, 0.9),
#                       c(length(States), length(States)), byrow = TRUE)
# too much mixing
emissionProbs = matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8),
                       c(length(States), length(States)), byrow = TRUE)
                       

# define the HMM
hmm = initHMM(States, Symbols, transProbs = transProbs, emissionProbs = emissionProbs)
# This can be used to simulate under the HMM
#sim = simHMM(hmm, nSim)

# I think I can estimate and output ML parameters using this
#bW <- baumWelch(hmm, dat$D_state, maxIterations=100, delta=1E-9, pseudoCount=0)
#bW

# to instead use hmm sims, use this: vit = viterbi(hmm, sim$observation)
# This function gives the ML path for the hmm given the data
vit = viterbi(hmm, dat$D_state)
#vit
f = forward(hmm, dat$D_state)
b = backward(hmm, dat$D_state)
# these are the natural logarithms of the forward probabilities of each state 
# for each position, given an observed sequence 
i <- f[1, nSim] # given an observed sequence, the ln(probabilities of the data) if the last state is BABA, 
j <- f[2, nSim] # given an observed sequence, the ln(probabilities of the data) if the last state is BBAA,
k <- f[3, nSim] # given an observed sequence, the ln(probabilities of the data) if the last state is ABBA,
# below is explained on pg 79 of Durban etal. 2007 for a two state HMM
# This is the demoninator of the Bayes formula, I think, and is the sum of 
# the probabilities of the data over all possible end states
# but instead of calculating this for two states: probObservations = ln (exp(i) + exp(j))
# this is used, which is the same thing but without overflow: probObservations = (i + ln(1 + exp(j - i)))

# From JD:
# The simplest, which is also exact, is:
#   log(a+b) = log(a) + log(1+b/a)
# The natural extension would be
#   log(a+b+c) = log(a) + log(1+b/a) + log(1+c/(b+a))
# Which is the same as this where x = ln(a), y=ln(b), and z=ln(c), and x>y>z
#   (x + log(1 + exp(y - x)) + log(1 + exp(z - a) / (1 + exp(y - x)) ))

# From url
# probObservations = (i + log(1 + exp(j - i)))

# Because I have three states, what I want to calculate is this
# probObservations = log (exp(i) + exp(j) + exp(k))

x <- 0
y <- 0
z <- 0
# first need to rank probabilities so x > y > z
# find highest
if((i>j)&&(i>k)){x <- i;y <- j;z <- k} else if ((j>i)&&(j>k)){x <- j;y <- i;z <- k} else {x <- k;y <- i;z <- j}
x
placeholder<-0
if(z>y){placeholder <- z;z <- y;y <- placeholder} 
y
z

probObservations <- (j + log(1 + exp(i - j)) + log(1 + exp(k-j) / (1+exp(i-j)) ))
probObservations
posterior = exp((f + b) - probObservations)
#posterior[2,]

dat$posteriorBBAA <-posterior[2,] # this is the posterior of BBAA sites
dat$posteriorABBA <-posterior[3,] # this is the posterior of ABBA sites
dat$posteriorBABA <-posterior[1,] # this is the posterior of BABA sites

newdf<-data.frame(dat$chromosome,dat$begin,dat$posteriorBBAA,dat$posteriorABBA,dat$posteriorBABA)

require(reshape2)
reshaped_dat <- melt(newdf, id.vars=c("dat.chromosome", "dat.begin"))

# set up color column
# set all states to zero
reshaped_dat$color <- "red"
# update the ones that have a substantial value for D
reshaped_dat$color[which(reshaped_dat$variable == 'dat.posteriorABBA')] <- "blue"
reshaped_dat$color[which(reshaped_dat$variable == 'dat.posteriorBBAA')] <- "gray90"


library(ggplot2)
pdf("./H3_maura_H1_tog_H2_tonk_HMM.pdf",w=8, h=8, version="1.4", bg="transparent")
c <- ggplot(reshaped_dat, aes(x = dat.begin, y = value, fill = color)) +
  geom_bar(position="stack", stat = "identity") +
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  scale_fill_manual(values=c("blue"="blue", "gray90"="gray90", "red"="red"), 
                      name="Posterior Probability",
                      breaks=c("blue", "gray90", "red"),
                      labels=c("ABBA", "BBAA", "BABA")) +
  facet_grid(dat.chromosome ~ ., scales = "free_x", space = "free_x", 
               margins = FALSE, drop = TRUE, switch = "y")+
  scale_y_continuous(name="Posterior probability", breaks=c(0,1.0))+
  scale_x_continuous(name="Position (Mbp)", 
                      breaks=c(0,50000000,100000000,150000000,200000000),labels=c("0","50","100","150","200"))+
  theme(strip.placement = "outside")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle=180),
        panel.border = element_rect(colour = "black"))
c
dev.off()



# extract windows with high posterior probability for ABBA pattern
ABBA_windowz <- dat[which(dat$posteriorABBA > 0.8),]
#head(ABBA_windowz)
ABBA_windowz[,c("chromosome","begin")]
# extract windows with high posterior probability for BABA pattern
BABA_windowz <- dat[which(dat$posteriorBABA > 0.8),]
#head(BABA_windowz)
BABA_windowz[,c("chromosome","begin")]

# now get annotations of genes in these windows
annotations<-read.table("/projects/2017_SEAsian_macaque_genomz/MacaM_Rhesus_Genome_Annotation_v7.6.8.gtf", sep="\t", header=FALSE)
colnames(annotations) <- c("chr","source","feature","start","end","score","strand","frame","attribute")
annotations_genes_only <- annotations[annotations[,2]=="gene",]

# https://stackoverflow.com/questions/19101849/overlapping-genomic-ranges
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

library(GenomicRanges)
gr0 = with(ABBA_windowz, GRanges(chromosome, IRanges(start=begin, end=end)))
gr1 = with(annotations, GRanges(chr, IRanges(start=start, end=end)))
#  add metadata
mcols(gr1) <- DataFrame(annotations$attribute)

#hits = findOverlaps(gr0, gr1)
hits = findOverlaps(gr1, gr0)
# Then updated the relevant start / end coordinates
#ranges(gr0)[queryHits(hits)] = ranges(gr1)[subjectHits(hits)]
ranges(gr1)[queryHits(hits)] = ranges(gr0)[subjectHits(hits)]
# now gr1 is an object with genomic ranges of annotated genes plus annotations
# in the MacaM genome

# coerce this to a dataframe
df <- data.frame(seqnames=seqnames(gr1),
                 starts=start(gr1)-1,
                 ends=end(gr1),
                 names=c(rep(".", length(gr1))),
                 scores=c(rep(".", length(gr1))),
                 strands=strand(gr1),
                 annotation=list(gr1))
df
# get unique annotations
unique_abba_gene_annotations <- unique(df$annotation.annotations.attribute)


# how many of these ranges are there?
length(gr0)
# what are the data in this object?
seqnames(gr0)
# what are the ranges
ranges(gr0)
# what are the start coordinates
start(gr0)
# what are the end coordinates
end(gr0)
# what are the chrs
names(gr0)
levels(gr0)
```


I looked at three conparisons on sulawesi:
*H3_maura_H1_tog_H2_tonk.concat
*H3_hecki_H1_nigra_H2_nigrescens.concat
*H3_hecki_H1_bru_H2_tonk.concat
and asked whether geneflow (ABBA windows with posterior probability > 0.8) were shared between any pair or all three.  I found two that were shared and both had relatively low divergence. Here they are:

window 1

chr09_62650001
chr09_62700001
chr09_62750001

window 2

chr11_37200001
chr11_37250001
chr11_37300001
chr11_37350001
chr11_37400001
chr11_37450001
chr11_37500001

