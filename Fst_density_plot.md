# Fst density plot.

I modified my permutation scripts so that they output Fst for each genomic window and add whether or not the window has a gene, and whether it has an N_interact gene.  I made density plots for each comparison using the R script below:
```R
setwd('/Users/Shared/Previously Relocated Items/Security/projects/2017_SEAsian_macaque_genomz/Fst_windows/FST_densities')

# I made the permutation perl scripts print out TajD and FW_H values of genomic windows in one 
# column and whether or not the window contains an N_interact gene in the other column

# this script will make faceted density plots for each of the 5 populations for all N_interact

# read in the TajD data 
bor_mau <- read.table("./nem_mau_windowstats.concat_FST__density.txt", header = T, sep="\t")
bor_ton <- read.table("./nem_ton_windowstats.concat_FST__density.txt", header = T, sep="\t")
bor_hec <- read.table("./nem_hec_windowstats.concat_FST__density.txt", header = T, sep="\t")
bor_nig <- read.table("./nem_nig_windowstats.concat_FST__density.txt", header = T, sep="\t")
mau_ton <- read.table("./mau_ton_windowstats.concat_FST__density.txt", header = T, sep="\t")
mau_hec <- read.table("./mau_hec_windowstats.concat_FST__density.txt", header = T, sep="\t")
mau_nig <- read.table("./mau_nig_windowstats.concat_FST__density.txt", header = T, sep="\t")
ton_hec <- read.table("./ton_hec_windowstats.concat_FST__density.txt", header = T, sep="\t")
ton_nig <- read.table("./ton_nig_windowstats.concat_FST__density.txt", header = T, sep="\t")
hec_nig <- read.table("./hec_nig_windowstats.concat_FST__density.txt", header = T, sep="\t")

# define a pairwise variable
bor_mau$pairwise <- "bor_mau"
bor_ton$pairwise <- "bor_ton"
bor_hec$pairwise <- "bor_hec"
bor_nig$pairwise <- "bor_nig"
mau_ton$pairwise <- "mau_ton"
mau_hec$pairwise <- "mau_hec"
mau_nig$pairwise <- "mau_nig"
ton_hec$pairwise <- "ton_hec"
ton_nig$pairwise <- "ton_nig"
hec_nig$pairwise <- "hec_nig"

# some quick statistics 
# independent 2-group Mann-Whitney U Test
wilcox.test(bor_mau$Fst~bor_mau$containsgenes)
wilcox.test(bor_ton$Fst~bor_ton$containsgenes)
wilcox.test(bor_hec$Fst~bor_hec$containsgenes)
wilcox.test(bor_nig$Fst~bor_nig$containsgenes)
wilcox.test(mau_ton$Fst~mau_ton$containsgenes)
wilcox.test(mau_hec$Fst~mau_hec$containsgenes)
wilcox.test(mau_nig$Fst~mau_nig$containsgenes)
wilcox.test(ton_hec$Fst~ton_hec$containsgenes)
wilcox.test(ton_nig$Fst~ton_nig$containsgenes)
wilcox.test(hec_nig$Fst~hec_nig$containsgenes)
# where y is numeric and A is A binary factor

# all are significant.
# sanity check to make sure that Fst is higher for genic compared to non-genic windows
only_nongenez <- hec_nig[which(hec_nig$containsgenes == 0),] ; mean(only_nongenez$Fst, na.rm = T)
only_genez <- hec_nig[which(hec_nig$containsgenes == 1),]; mean(only_genez$Fst, na.rm = T) 

only_nongenez <- bor_nig[which(bor_nig$containsgenes == 0),] ; mean(only_nongenez$Fst, na.rm = T)
only_genez <- bor_nig[which(bor_nig$containsgenes == 1),]; mean(only_genez$Fst, na.rm = T) 

# select only windows with genes
only_genez <- my_data[which(my_data$containsgenes != 0),] 
# independent 2-group Mann-Whitney U Test
wilcox.test(only_genez$TajD~only_genez$containsNinteractgenez)
wilcox.test(only_genez$FW_H~only_genez$containsNinteractgenez)
# where y is numeric and A is A binary factor



# now rbind them together
all_my_data <- rbind(bor_mau,bor_ton,bor_hec,bor_nig,mau_ton,mau_hec,mau_nig,ton_hec,ton_nig,hec_nig)

# define a grouping variable
all_my_data$group <- "NA"
# numbering the categories like this makes the one sided p-value sensible
all_my_data$group[which((all_my_data$containsgenes == 0)&(all_my_data$containsNinteractgenez == 0))] <- "No genes"
all_my_data$group[which((all_my_data$containsgenes == 1)&(all_my_data$containsNinteractgenez == 0))] <- "Other genes"
all_my_data$group[which((all_my_data$containsgenes == 1)&(all_my_data$containsNinteractgenez == 1))] <- "Ninteract genes"
table(all_my_data$group)

# re-order the species
all_my_data$pairwise_f = factor(all_my_data$pairwise, levels=c('bor_mau','bor_ton','bor_hec','bor_nig',
                                                             'mau_ton','mau_hec','mau_nig','ton_hec',
                                                             'ton_nig','hec_nig'), ordered = T)
# reorder the group
all_my_data$group_f = factor(all_my_data$group, levels=c('No genes','Other genes','Ninteract genes'), ordered = T)


#Plot.
library(ggplot2)
png(paste("Fst_density_allNinteract.png",sep=""),
    width = 300, height = 200, units='mm', res = 100)
      
      ggplot(all_my_data) + 
      geom_density(aes(x=Fst, colour=group_f),show_guide=FALSE)+
      stat_density(aes(x=Fst, colour=group_f),
                   geom="line",position="identity")+
      scale_color_manual(values = c("black", "blue", "red")) +
      facet_wrap(vars(pairwise_f), nrow=2, ncol = 5, scales = "free_x") +
      theme_bw() +
      # remove legend key border color & background
      theme(legend.key=element_blank()) +
      theme(legend.box.background = element_blank())+
      #theme(legend.key = element_rect(colour = NA, fill = NA)) +
      theme(legend.title = element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(strip.background = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 20),
          strip.text = element_text(size = 20),
          legend.text = element_text(size = 20)) +
      scale_x_continuous(limits = c(0.05, .99),
        labels = scales::number_format(accuracy = 0.1)) 
dev.off()

```
