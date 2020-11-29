# Twisst

A nice complement to a sliding window analysis is Twisst, which infers phylogenies in genomic windows.

A first step is to phase the data; I did this with Beagle 5.0 like this:
```
#!/bin/sh
#SBATCH --job-name=beagle
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=128gb
#SBATCH --output=beagle.%J.out
#SBATCH --error=beagle.%J.err
#SBATCH --account=def-ben

# sbatch Beagle.sh chr

module load java

java -Xmx12g -jar beagle.18May20.d20.jar gt=../FandM_${1}_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.g
z.recode.vcf.gz out=../FandM_${1}_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz impute=t
rue 
```

Then I made a phased geno file:
```
python ./genomics_general/VCF_processing/parseVCF.py -i FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | gzip > phased_genos/chr01.geno.gz
```
In this directory:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/twisst
```
and only for chrX, I had to convert the male genotypes to diploid like this before phasing:
```
echo "chrX 1 148935249 M 2" > ploidy.txt
bcftools +fixploidy ../all_diploid_haploid_chrX_BSQR_filtered3_noPAR_SNPsonly.vcf.gz.recode.vcf.gz.recode.vcf.gz -Ov -- -p ploidy.txt > chrX_diploid.vcf
```
and then make a geno file from this...
```
python ../genomics_general/VCF_processing/parseVCF.py -i ../all_diploid_haploid_chrX_phased.vcf.gz.vcf.gz | gzip > ../phased_genos/chrX.geno.gz 
```


I had to add a line to the 'phyml_sliding_windows.py' script to let it know where to look for the genomics.py file:

```
sys.path.append("/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_dept
h_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general")
import genomics
```
I also had to install phyml and add it to my path like this:
```
PATH=$PATH:/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/twisst/phyml/src
```
and also load several modules in some weird order...:
```
module load python/3.8
module load StdEnv/2020
module load scipy-stack/2020b
```

and then I made nj trees for each window like this:
```
python3 genomics_general/phylo/phyml_sliding_windows.py -T 10 -g phased_genos/chr02a.geno.gz --prefix phased_genos/chr02a_treez_w50 -w 50 --windType sites --model GTR
```

To get twisst to work I had to update the ete3 package first:
```
pip install --upgrade ete3
```

And now this seems to work (after changing the ../genomics_general/pops_twisst.txt file to not inlcude separate sum and tog populations):
```
python twisst.py -t ../phased_genos/chr01_treez_w50.trees.gz -w ../phased_genos/chr01_output.weights.csv.gz --outputTopos ../phased_genos/topologies_chr01.trees --outgroup papio -g nig -g nge -g hec -g ton -g mau -g bru -g nem -g papio --method complete --groupsFile ../genomics_general/pops_twisst.txt
```

And I'm running twisst with this sbatch script for each chr:
```
#!/bin/sh
#SBATCH --job-name=twisst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=128gb
#SBATCH --output=twisst.%J.out
#SBATCH --error=twisst.%J.err
#SBATCH --account=def-ben

# sbatch Twisst.sh chr
module load StdEnv/2020
module load scipy-stack/2020b
# may have to type this before running: pip install --upgrade ete3

python twisst.py -t ../phased_genos/${1}_treez_w50.trees.gz -w ../phased_genos/${1}_twisstoutput.weights.csv.gz --outputTopos 
../phased_genos/${1}_twisst_topologies.trees --outgroup papio -g nig -g nge -g hec -g ton -g mau -g bru -g nem -g papio --meth
od complete --groupsFile ../genomics_general/pops_twisst.txt --verbose
```


Because Twisst is limited to analysis of 15 topologies or less, I grouped taxa into (i) sumatra, (ii) borneo, (iii) nsulawesi, (iv) ssulawesi, and (v) papio. I combined the functions that Simon Martin included with the Twist package into a modified example file below.  I changes some for the functions and plotting parameters slightly to get the plots I was looking for.

```
# execute these modules if you need to reinstall ape:
# module load nixpkgs/16.09  gcc/7.3.0
# module load gcc/7.3.0 r/4.0.2
# r/3.6.0


rm(list = ls())



simple.loess.predict <- function(x, y, span, new_x=NULL, weights = NULL, max = NULL, min = NULL, family=NULL){
  y.loess <- loess(y ~ x, span = span, weights = weights, family=family)
  if (is.null(new_x)) {y.predict <- predict(y.loess,x)}
  else {y.predict <- predict(y.loess,new_x)}
  if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
  if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
  y.predict
}

smooth.df <- function(x, df, span, new_x = NULL, col.names=NULL, weights=NULL, min=NULL, max=NULL, family=NULL){
  if (is.null(new_x)) {smoothed <- df}
  else smoothed = df[1:length(new_x),]
  if (is.null(col.names)){col.names=colnames(df)}
  for (col.name in col.names){
    print(paste("smoothing",col.name))
    smoothed[,col.name] <- simple.loess.predict(x,df[,col.name],span = span, new_x = new_x, max = max, min = min, weights = weights, family=family)
  }
  smoothed
}

smooth.weights <- function(window_positions, weights_dataframe, span, new_positions=NULL, window_sites=NULL){
  weights_smooth <- smooth.df(x=window_positions,df=weights_dataframe,
                              span=span, new_x=new_positions, min=0, max=1, weights=window_sites)
  
  #return rescaled to sum to 1
  weights_smooth / apply(weights_smooth, 1, sum)
}


stack <- function(mat){
  upper <- t(apply(mat, 1, cumsum))
  lower <- upper - mat
  list(upper=upper,lower=lower)
}

interleave <- function(x1,x2){
  output <- vector(length= length(x1) + length(x2))
  output[seq(1,length(output),2)] <- x1
  output[seq(2,length(output),2)] <- x2
  output
}


sum_df_columns <- function(df, columns_list){
  new_df <- df[,0]
  for (x in 1:length(columns_list)){
    if (length(columns_list[[x]]) > 1) new_df[,x] <- apply(df[,columns_list[[x]]], 1, sum, na.rm=T)
    else new_df[,x] <- df[,columns_list[[x]]]
    if (is.null(names(columns_list)[x]) == FALSE) names(new_df)[x] <- names(columns_list)[x]
  }
  new_df
}


plot.weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,density=NULL,lwd=1,xlim=c(0,230000000),ylim=c(0,1),stacked=FALSE,
                         ylab=NULL, xlab = NULL, main="",xaxt=NULL,yaxt=NULL,bty="n", add=FALSE){
  #get x axis
  x = positions
  #if a two-column matrix is given - plot step-like weights with start and end of each window    
  if (dim(as.matrix(x))[2]==2) {
    x = interleave(positions[,1],positions[,2])
    yreps=2
  }
  else {
    if (is.null(x)==FALSE) x = positions
    else x = 1:nrow(weights_dataframe)
    yreps=1
  }
  
  #set x limits
  if(is.null(xlim)) xlim = c(min(x), max(x))
  
  #if not adding to an old plot, make a new plot
  if (add==FALSE) plot(0, pch = "", xlim = xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main,xaxt='n',yaxt='n',bty=bty, ann=FALSE, frame.plot=TRUE); Axis(side=1, labels=FALSE); Axis(side=2, labels=FALSE)
  
  if (stacked == TRUE){
    y_stacked <- stack(weights_dataframe)
    for (n in 1:ncol(weights_dataframe)){
      y_upper = rep(y_stacked[["upper"]][,n],each=yreps)
      y_lower = rep(y_stacked[["lower"]][,n],each = yreps)
      polygon(c(x,rev(x)),c(y_upper, rev(y_lower)), col = fill_cols[n], density=density[n], border=NA)
    }
  }
  else{
    for (n in 1:ncol(weights_dataframe)){
      y = rep(weights_dataframe[,n],each=yreps)
      polygon(c(x,rev(x)),c(y, rep(0,length(y))), col=fill_cols[n], border=NA,density=density[n])
      lines(x,y, type = "l", col = line_cols[n],lwd=lwd)
    }
  }
}

options(scipen = 7)

#Heres a set of 15 colourful colours from https://en.wikipedia.org/wiki/Help:Distinguishable_colors
topo_cols <- c(
  "#0075DC", #Blue
  "#2BCE48", #Green
  "#FFA405", #Orpiment
  "#5EF1F2", #Sky
  "#FF5005", #Zinnia
  "#005C31", #Forest
  "#00998F", #Turquoise
  "#FF0010", #Red
  "#9DCC00", #Lime
  "#003380", #Navy
  "#F0A3FF", #Amethyst
  "#740AFF", #Violet
  "#426600", #Quagmire
  "#C20088", #Mallow
  "#94FFB5") #Jade




########### Below are some more object-oriented tools for working with standard twisst output files

library(ape)
library(data.table)
library(tools)

#a function that imports weights 
import.twisst <- function(weights_files, window_data_files=NULL, split_by_chrom=TRUE, reorder_by_start=FALSE, na.rm=TRUE, max_window=Inf,
                          lengths=NULL, topos_file=NULL, recalculate_mid=FALSE, names=NULL){
  l = list()
  
  if (length(window_data_files) > 1){
    print("Reading weights and window data")
    l$window_data <- lapply(window_data_files, read.table ,header=TRUE)
    l$weights_raw <- lapply(weights_files, read.table, header=TRUE)
    if (is.null(names) == FALSE) names(l$window_data) <- names(l$weights_raw) <- names
  }
  
  if (length(window_data_files) == 1){
    print("Reading weights and window data")
    l$window_data <- list(read.table(window_data_files, header=TRUE))
    l$weights_raw <- list(read.table(weights_files, header=TRUE))
    if (split_by_chrom == TRUE){
      l$weights_raw <- split(l$weights_raw[[1]], l$window_data[[1]][,1])
      l$window_data <- split(l$window_data[[1]], l$window_data[[1]][,1])
    }
  }
  
  if (is.null(window_data_files) == TRUE) {
    print("Reading weights")
    l$weights_raw <- lapply(weights_files, read.table, header=TRUE)
    n <- nrow(l$weights_raw[[1]])
    l$window_data <- list(data.frame(chrom=rep(0,n), start=1:n, end=1:n))
    if (is.null(names) == FALSE) names(l$window_data) <- names
  }
  
  l$n_regions <- length(l$weights_raw)
  
  if (is.null(names(l$window_data)) == TRUE) {
    names(l$window_data) <- names(l$weights_raw) <- paste0("region", 1:l$n_regions)
  }
  
  print(paste("Number of regions:", l$n_regions))
  
  if (reorder_by_start==TRUE & is.null(window_data_files) == FALSE){
    print("Reordering")
    orders = sapply(l$window_data, function(df) order(df[,2]), simplify=FALSE)
    l$window_data <- sapply(names(orders), function(x) l$window_data[[x]][orders[[x]],], simplify=F)
    l$weights_raw <- sapply(names(orders), function(x) l$weights_raw[[x]][orders[[x]],], simplify=F)
  }
  
  print("Computing summaries")
  
  l$weights <- sapply(l$weights_raw, function(raw) raw/apply(raw, 1, sum), simplify=FALSE)
  
  l$weights_mean <- lapply(l$weights, apply, 2, mean, na.rm=T)
  
  l$weights_overall_mean <- apply(rbindlist(l$weights), 2, mean, na.rm=T)
  
  if (is.null(lengths) == TRUE) l$lengths <- sapply(l$window_data, function(df) tail(df$end,1), simplify=FALSE)
  else l$lengths = lengths
  
  print("Cleaning data")
  
  if (na.rm==TRUE){
    for (i in 1:l$n_regions){
      #remove rows containing NA values
      good_rows = which(is.na(apply(l$weights[[i]],1,sum)) == F &
                          l$window_data[[i]]$end - l$window_data[[i]]$start + 1 <= max_window)
      l$weights[[i]] <- l$weights[[i]][good_rows,]
      l$weights_raw[[i]] <- l$weights[[i]][good_rows,]
      l$window_data[[i]] = l$window_data[[i]][good_rows,]
    }
  }
  
  for (i in 1:length(l$window_data)) {
    if (is.null(l$window_data[[i]]$mid) == TRUE | recalculate_mid == TRUE) {
      l$window_data[[i]]$mid <- (l$window_data[[i]]$start + l$window_data[[i]]$end)/2
    }
  }
  
  l$pos=sapply(l$window_data, function(df) df$mid, simplify=FALSE)
  
  print("Getting topologies")
  
  #attempt to retrieve topologies
  l$topos=NULL
  #first, check if a topologies file is provided
  if (is.null(topos_file) == FALSE) {
    l$topos <- read.tree(file=topos_file)
    if (is.null(names(l$topos)) == TRUE) names(l$topos) <- names(l$weights[[1]])
  }
  else{
    #otherwise we try to retrieve topologies from the (first) weights file
    n_topos = ncol(l$weights[[1]])
    if (file_ext(weights_files[1]) == "gz") cat="gunzip -c" else cat="cat"
    topos_text <- try(system(paste(cat, weights_files[[1]], "2>/dev/null", "| head -n", n_topos), intern = T), silent=TRUE)
    try(l$topos <- read.tree(text = topos_text))
    try(names(l$topos) <- sapply(names(l$topos), substring, 2))
  }
  l
}


smooth.twisst <- function(twisst_object, span=0.05, span_bp=NULL, spacing=NULL) {
  l=list()
  
  l$topos <- twisst_object$topos
  
  l$n_regions <- twisst_object$n_regions
  
  l$weights <- list()
  
  l$lengths = twisst_object$lengths
  
  l$pos <- list()
  
  for (i in 1:l$n_regions){
    if (is.null(span_bp) == FALSE) span <- span_bp/twisst_object$length[[i]]
    
    if (is.null(spacing) == TRUE) spacing <- twisst_object$length[[i]]*span*.1
    
    l$pos[[i]] <- seq(twisst_object$pos[[i]][1], tail(twisst_object$pos[[i]],1), spacing)
    
    l$weights[[i]] <- smooth.weights(twisst_object$pos[[i]], twisst_object$weights[[i]], new_x <- l$pos[[i]], span = span,
                                     window_sites=twisst_object$window_data$sites[[i]])
  }
  l
}


plot.twisst <- function(twisst_object, show_topos=TRUE, ncol_topos=NULL, regions=NULL, ncol_weights=1,
                        cols=topo_cols,xlim=NULL, ylim=NULL, mode=2, rel_height=3, tree_type="clad",
                        concatenate=FALSE, gap=0, include_region_names=FALSE){

  par(oma=c(3,3,3,3))			
  if (mode==3) {
    stacked=TRUE
    fill_cols = cols
    line_cols = NA
    lwd = 0
  }
  
  if (mode==2) {
    stacked=FALSE
    fill_cols = paste0(cols,80)
    line_cols = cols
    lwd=par("lwd")
  }
  
  if (mode==1) {
    stacked=FALSE
    fill_cols = NA
    line_cols = cols
    lwd=par("lwd")
  }
  
  if (is.null(regions)==TRUE) regions <- 1:twisst_object$n_regions
  
  if (concatenate == TRUE) ncol_weights <- 1
  
  if (show_topos==TRUE){
    n_topos <- length(twisst_object$topos)
    
    if (is.null(ncol_topos)) ncol_topos <- n_topos
    
    #if we have too few topologies to fill the spaces in the plot, we can pad in the remainder
    topos_pad <- (n_topos * ncol_weights) %% (ncol_topos*ncol_weights) 
    
    topos_layout_matrix <- matrix(c(rep(1:n_topos, each=ncol_weights), rep(0, topos_pad)),
                                  ncol=ncol_topos*ncol_weights, byrow=T)
  }
  else {
    ncol_topos <- 1
    n_topos <- 0
    topos_layout_matrix <- matrix(NA, nrow= 0, ncol=ncol_topos*ncol_weights)
  }
  
  #if we have too few regions to fill the spaces in the plot, we pad in the remainder
  data_pad <- (length(regions)*ncol_topos) %% (ncol_topos*ncol_weights)
  
  if (concatenate==TRUE) weights_layout_matrix <- matrix(rep(n_topos+1,ncol_topos), nrow=1)
  else {
    weights_layout_matrix <- matrix(c(rep(n_topos+(1:length(regions)), each=ncol_topos),rep(0,data_pad)),
                                    ncol=ncol_topos*ncol_weights, byrow=T)
  }
  
  layout(rbind(topos_layout_matrix, weights_layout_matrix),
         height=c(rep(1, nrow(topos_layout_matrix)), rep(rel_height, nrow(weights_layout_matrix))))
  
  if (show_topos == TRUE){
    if (tree_type=="unrooted"){
      par(mar=c(2,2,2,2), xpd=NA)
      
      for (i in 1:n_topos){
        plot.phylo(twisst_object$topos[[i]], type = tree_type, edge.color=cols[i],
                   edge.width=5, label.offset=0.3, cex=1, rotate.tree = 90)
        mtext(side=3,text=paste0("topo",i), cex=0.75)
      }
    }
    else{
      par(mar=c(1,1,1,1), xpd=NA)
      
      for (i in 1:n_topos){
        plot.phylo(twisst_object$topos[[i]], type = tree_type, edge.color=cols[i],
                   edge.width=5, label.offset=0.3, cex=1, rotate.tree = 0)
        mtext(side=3,text=paste0("topo",i), cex=0.75)
      }
    }
  }
  
  if (is.null(ylim)==TRUE) ylim <- c(0,1)
  
  par(mar=c(0.4,6,0.4,6),xpd=FALSE)
  
  if (concatenate == TRUE) {
    chrom_offsets = cumsum(twisst_object$lengths + gap) - (twisst_object$lengths + gap)
    chrom_ends <- chrom_offsets + twisst_object$lengths
    
    plot(0, pch = "",xlim = c(chrom_offsets[1],tail(chrom_ends,1)), ylim = ylim,
         ylab = "", yaxt = "n", xlab = "", xaxt = "n", bty = "n", main = "")
    
    for (j in regions) {
      if (is.null(twisst_object$window_data[[j]])) positions <- twisst_object$pos[[j]] + chrom_offsets[j]
      else positions <- twisst_object$window_data[[j]][,c("start","end")] + chrom_offsets[j]
      plot.weights(twisst_object$weights[[j]], positions, xlim=xlim,
                   fill_cols = fill_cols, line_cols=line_cols,lwd=lwd,stacked=stacked, add=T)
    }
  }
  else{
    for (j in regions){
      if (is.null(twisst_object$window_data[[j]])) positions <- twisst_object$pos[[j]]
      else positions <- twisst_object$window_data[[j]][,c("start","end")]
      plot.weights(twisst_object$weights[[j]], positions, xlim=xlim, ylim = ylim,
                   fill_cols = fill_cols, line_cols=line_cols,lwd=lwd,stacked=stacked,
                   main = ifelse(include_region_names==TRUE, names(twisst_object$weights)[j], ""))
      }
  }
  mtext("Weights", side=2, line=0, cex=2,outer=TRUE)
  mtext("Position (Mbp)", side=1, line=2, cex=2, outer=TRUE)
}


#function for plotting tree that uses ape to get node positions
draw.tree <- function(phy, x, y, x_scale=1, y_scale=1, method=1, direction="right",
                      col="black", col.label="black", add_labels=TRUE, add_symbols=FALSE,
                      label_offset = 1, symbol_offset=0, col.symbol="black",symbol_bg="NA",
                      pch=19, cex=NULL, lwd=NULL, label_alias=NULL){
  
  n_tips = length(phy$tip.label)
  
  if (direction=="right") {
    node_x = (node.depth(phy, method=method) - 1) * x_scale * -1
    node_y = node.height(phy) * y_scale
    label_x = node_x[1:n_tips] + label_offset
    label_y = node_y[1:n_tips]
    adj_x = 0
    adj_y = .5
    symbol_x = node_x[1:n_tips] + symbol_offset
    symbol_y = node_y[1:n_tips]
  }
  if (direction=="down") {
    node_y = (node.depth(phy, method=method) - 1) * y_scale * 1
    node_x = node.height(phy) * x_scale
    label_x = node_x[1:n_tips]
    label_y = node_y[1:n_tips] - label_offset
    adj_x = .5
    adj_y = 1
    symbol_x = node_x[1:n_tips]
    symbol_y = node_y[1:n_tips] - symbol_offset
  }
  
  #draw edges
  segments(x + node_x[phy$edge[,1]], y + node_y[phy$edge[,1]],
           x + node_x[phy$edge[,2]], y + node_y[phy$edge[,2]], col=col, lwd=lwd)
  
  if (is.null(label_alias) == FALSE) tip_labels <- label_alias[phy$tip.label]
  else tip_labels <- phy$tip.label
  
  if (add_labels=="TRUE") text(x + label_x, y + label_y, col = col.label, labels=tip_labels, adj=c(adj_x,adj_y),cex=cex)
  if (add_symbols=="TRUE") points(x + symbol_x, y + symbol_y, pch = pch, col=col.symbol, bg=symbol_bg)
  
}


#code for plotting a summary barplot
plot.twisst.summary <- function(twisst_object, order_by_weights=TRUE, only_best=NULL, cols=topo_cols,
                                x_scale=0.12, y_scale=0.15, direction="right", col="black", col.label="black",
                                label_offset = 0.05, lwd=NULL, cex=NULL){
  
  # Either order 1-15 or order with highest weigted topology first
  
  if (order_by_weights == TRUE) {
    ord <- order(twisst_object$weights_overall_mean, decreasing=T)
    if (is.null(only_best) == FALSE) ord=ord[1:only_best]
  }
  else ord <- 1:length(twisst_object$topos)
  
  N=length(ord)
  
  #set the plot layout, with the tree panel one third the height of the barplot panel
  layout(matrix(c(2,1)), heights=c(1,3))
  
  par(mar = c(1,4,.5,1))
  
  #make the barplot
  x=barplot(twisst_object$weights_overall_mean[ord], col = cols[ord],
            xaxt="n", las=1, ylab="Average weighting", space = 0.2, xlim = c(0.2, 1.2*N),
	    cex.lab=1.5)
  
  #draw the trees
  #first make an empty plot for the trees. Ensure left and right marhins are the same
  par(mar=c(0,4,0,1))
  plot(0,cex=0,xlim = c(0.2, 1.2*N), xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1), bty="n")
  
  #now run the draw.tree function for each topology. You can set x_scale and y_scale to alter the tree width and height.
  for (i in 1:length(ord)){
    draw.tree(twisst_object$topos[[ord[i]]], x=x[i]+.2, y=0, x_scale=x_scale, y_scale=y_scale,
              col=cols[ord[i]], label_offset=label_offset, cex=cex, lwd=lwd)
  }
  
  #add labels for each topology
  # text(x,.9,names(twisst_object$topos)[ord],col=cols[ord])
}


#code for plotting a summary boxplot
plot.twisst.summary.boxplot <- function(twisst_object, order_by_weights=TRUE, only_best=NULL, cols=topo_cols,
                                        x_scale=0.12, y_scale=0.15, direction="right", col="black", col.label="black",
                                        label_offset = 0.05, lwd=NULL, label_alias=NULL, cex=NULL, outline=FALSE,
                                        cex.outline=NULL, lwd.box=NULL, topo_names=NULL){
  
  # Either order 1-15 or order with highest weigted topology first
  
  if (order_by_weights == TRUE) {
    ord <- order(twisst_object$weights_overall_mean, decreasing=T)
    if (is.null(only_best) == FALSE) ord=ord[1:only_best]
  }
  else ord <- 1:length(twisst_object$topos)
  
  N=length(ord)
  
  #set the plot layout, with the tree panel one third the height of the barplot panel
  layout(matrix(c(2,1)), heights=c(1,3))
  
  par(mar = c(1,4,.5,1))
  
  #make the barplot
  boxplot(as.data.frame(rbindlist(twisst_object$weights))[,ord], col = cols[ord],
          xaxt="n", las=1, xlim = c(.5, N+.5), ylab="Average weighting", outline=outline, cex=cex.outline, lwd=lwd.box)
  
  #draw the trees
  #first make an empty plot for the trees. Ensure left and right marhins are the same
  par(mar=c(0,4,0,1))
  plot(0,cex=0, xlim = c(.5, N+.5), xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1), bty="n")
  
  #now run the draw.tree function for each topology. You can set x_scale and y_scale to alter the tree width and height.
  for (i in 1:N){
    draw.tree(twisst_object$topos[[ord[i]]], i+.1, y=0, x_scale=x_scale, y_scale=y_scale,
              col=cols[ord[i]], label_offset=label_offset, cex=cex, lwd=lwd, label_alias=label_alias)
  }
  
  if (is.null(topo_names)==TRUE) topo_names <- names(twisst_object$topos)
  
  #add labels for each topology
  text(1:N,.9,topo_names[ord],col=cols[ord])
}

#function for subsetting the twisst object by a set of topologies
subset.twisst.by.topos <- function(twisst_object, topos){
  l <- list()
  regions <- names(twisst_object$weights)
  l$window_data <- twisst_object$window_data
  l$n_regions <- twisst_object$n_regions
  l$lengths <- twisst_object$lengths
  l$pos <- twisst_object$pos
  l$weights_raw <- sapply(regions, function(region) twisst_data$weights_raw[[region]][,topos], simplify=F)
  l$weights <- sapply(regions, function(region) twisst_data$weights[[region]][,topos], simplify=F)
  l$weights_mean <- sapply(regions, function(region) twisst_data$weights_mean[[region]][topos], simplify=F)
  l$weights_overall_mean <- twisst_data$weights_overall_mean[topos]
  l$topos <- twisst_data$topos[topos]
  l
}

#function for subsetting the twisst object by specific regions
subset.twisst.by.regions <- function(twisst_object, regions){
  l <- list()
  regions <- names(twisst_object$weights[regions])
  l$window_data <- twisst_object$window_data[regions]
  l$n_regions <- length(regions)
  l$lengths <- twisst_object$lengths[regions]
  l$pos <- twisst_object$pos[regions]
  l$weights_raw <- twisst_data$weights_raw[regions]
  l$weights <- twisst_data$weights[regions]
  l$weights_mean <- twisst_data$weights_mean[regions]
  l$weights_overall_mean <- apply(rbindlist(l$weights), 2, mean, na.rm=T)
  l$topos <- twisst_data$topos
  l
}

################################# overview #####################################

# The main data produced by Twisst is a weights file which has columns for each
# topology and their number of observations of that topology within each
# genealogy. Weights files produced by Twisst also contain initial comment lines
# speficying the topologies.

# The other data file that may be of interest is the window data. That is, the
# chromosome/scaffold and start and end positions for each of the regions or
# windows represented in the weights file.

# Both of the above files can be read into R, manipulated and plotted however
# you like, but I have written some functions to make these tasks easier.
# These functions are provided in the script plot_twisst.R

################### load helpful plotting functions #############################
setwd(".")
#source("plot_twisst.R")

############################## input files ######################################

# It is possible to import one or more weights files at a time.
# Here we just import 1

#weights file with a column for each topology
#weights_file <- "chr01_nemsulatwisstoutput.weights.csv.gz"
weights_file <- c(
"chr01_nemsulatwisstoutput.weights.csv.gz",
"chr02a_nemsulatwisstoutput.weights.csv.gz",
"chr02b_nemsulatwisstoutput.weights.csv.gz",
"chr03_nemsulatwisstoutput.weights.csv.gz",
"chr04_nemsulatwisstoutput.weights.csv.gz",
"chr05_nemsulatwisstoutput.weights.csv.gz",
"chr06_nemsulatwisstoutput.weights.csv.gz",
"chr07_nemsulatwisstoutput.weights.csv.gz",
"chr08_nemsulatwisstoutput.weights.csv.gz",
"chr09_nemsulatwisstoutput.weights.csv.gz",
"chr10_nemsulatwisstoutput.weights.csv.gz",
"chr11_nemsulatwisstoutput.weights.csv.gz",
"chr12_nemsulatwisstoutput.weights.csv.gz",
"chr13_nemsulatwisstoutput.weights.csv.gz",
"chr14_nemsulatwisstoutput.weights.csv.gz",
"chr15_nemsulatwisstoutput.weights.csv.gz",
"chr16_nemsulatwisstoutput.weights.csv.gz",
"chr17_nemsulatwisstoutput.weights.csv.gz",
"chr18_nemsulatwisstoutput.weights.csv.gz",
"chr19_nemsulatwisstoutput.weights.csv.gz",
"chrX_nemsulatwisstoutput.weights.csv.gz")
#weights_file <- "chr18_twisstoutput.weights.csv.gz"


# It is not necessary to import window data files, but if you do there should be one for
# each weights file

#coordinates file for each window
#window_data_file <- "chr01_treez_w50.data.tsv"
window_data_file <- c(
"chr01_treez_w50.data.tsv",
"chr02a_treez_w50.data.tsv",
"chr02b_treez_w50.data.tsv",
"chr03_treez_w50.data.tsv",
"chr04_treez_w50.data.tsv",
"chr05_treez_w50.data.tsv",
"chr06_treez_w50.data.tsv",
"chr07_treez_w50.data.tsv",
"chr08_treez_w50.data.tsv",
"chr09_treez_w50.data.tsv",
"chr10_treez_w50.data.tsv",
"chr11_treez_w50.data.tsv",
"chr12_treez_w50.data.tsv",
"chr13_treez_w50.data.tsv",
"chr14_treez_w50.data.tsv",
"chr15_treez_w50.data.tsv",
"chr16_treez_w50.data.tsv",
"chr17_treez_w50.data.tsv",
"chr18_treez_w50.data.tsv",
"chr19_treez_w50.data.tsv",
"chrX_treez_w50.data.tsv")
#window_data_file <- "chr18_treez_w50.data.tsv"


################################# import data ##################################

# The function import.twisst reads the weights, window data  files into a list object
# If there are multiple weights files, or a single file with different chromosomes/scaffolds/contigs
# in the window data file, these will be separated when importing.

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)


#get list of the most abundant topologies (top 5 in this example)
top5_topos <- order(twisst_data$weights_overall_mean, decreasing=T)[1:5]
ordered_twisst <- order(twisst_data$weights_overall_mean, decreasing=T)
#interesting_topos <- order(twisst_data$weights_overall_mean, decreasing=T)[c(1:3,6:7)]

#subset twisst object for these
top5_topos_twisst <- subset.twisst.by.topos(twisst_data, top5_topos)
ordered_topos_twisst <- subset.twisst.by.topos(twisst_data, ordered_twisst)
#this can then be used in all the same plotting functions below


############################## combined plots ##################################
# there are a functions available to plot both the weightings and the topologies
#a summary plot shows all the topologies and a bar plot of their relative weightings
png("chr_all_summary.png",width = 3000, height = 2000, res = 200, bg = "transparent")
plot.twisst.summary(ordered_topos_twisst, order_by_weights=TRUE, lwd=3, cex=1.5,
                    x_scale=0.05, y_scale=0.10,)
dev.off()


#png("chr_all_top5summary.png",width = 2000, height = 2000, res = 200, bg = "transparent")
#plot.twisst.summary(top5_topos_twisst, order_by_weights=TRUE, lwd=3, cex=0.7,
#                    x_scale=0.05, y_scale=0.10,)
#dev.off()

# this is an attempt to make a stacked bar plot - difficult to interpret
#wgts <- twisst_data$weights
#wgtss <- as.data.frame(wgts)
#start <- twisst_data$window_data[["chrX"]]["start"]
#weightsdf <- cbind(start,wgtss)
#barplot(t(weightsdf[2:16]), col=rainbow(16), border=NA,space=0.05)



#or plot ALL the data across the chromosome(s)
#Note, this is not recommended if there are large numbers of windows.
# instead, it is recommended to first smooth the weghtings and plot the smoothed values
# plot.twisst(twisst_data, mode=1, show_topos=TRUE)


# make smooth weightings and plot those across chromosomes
smooth_top5_topos_twisst <- smooth.twisst(top5_topos_twisst, span_bp = 200000, spacing = 1000)
smooth_twisst_data <- smooth.twisst(ordered_topos_twisst, span_bp = 200000, spacing = 1000)

#png("Summary_chrX_smooth5interesting.png",width = 900, height = 600,res = 200, bg = "transparent")

#plot.twisst(interesting_topos_twisst, mode=3) # mode 2 overlays polygons, 
                                          # mode 3 would stack them
#dev.off()

png("chr_all_smooth_twisst_all.png",width = 3000, height = 4000, res = 200, bg = "transparent")

plot.twisst(smooth_twisst_data, mode=3) # mode 2 overlays polygons,
                                          # mode 3 would stack them
dev.off()

#png("chr_all_smooth_twisst_top5.png",width = 3000, height = 4000, res = 200, bg = "transparent")

#plot.twisst(smooth_top5_topos_twisst, mode=3) # mode 2 overlays polygons,
                                          # mode 3 would stack them
#dev.off()


png("chr_all_smooth_twisst_concat.png",width = 3500, height =3500, res = 200, bg = "transparent")
plot.twisst(smooth_twisst_data, mode=3, xlim = c(0,230000000), show_topos=FALSE)
dev.off()
```
