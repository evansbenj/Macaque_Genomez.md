# Here is a plot that makes the pi fig
```R
library(ggplot2)
library(reshape)
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2017_SEAsian_macaque_genomz/Fst_windows')

# this script will make a faceted plot of pi vs position for all of the interesting genes
window <- 10000000
options(scipen = 999) 

# initilize a dataframe for plotting
big_df <- data.frame(scaffold=character(),
                 start=numeric(), 
                 end=numeric(), 
                 Species=character(),
                 color=character(),
                 gener=character(),
                 value=numeric()
                 ) 
# all 11
#gene_vector <- c("NDUFAF3","NDUFA2","TACO1","HARS2","COX16","C11orf83","MRPS21","UQCRC1","ATP5B","MRPL48","NDUFB11")
#chr_vector <- c("chr03","chr05","chr17","chr05","chr14","chr11","chr01","chr03","chr12","chr11","chrX")
#beginning_vector <- c(104846075,138178553,56920802,138224999,132330617,11484179,124258900,104359835,54910245,65083808,47194137)
#end_vector <- c(104847927,138180781,56928897,138233197,132354632,11486231,124272013,104371673,54920096,65165565,47196604)

# just 6 good ones
#gene_vector <- c("NDUFAF3","NDUFA2","HARS2","MRPS21","UQCRC1","ATP5B")
#chr_vector <- c("chr03","chr05","chr05","chr01","chr03","chr12")
#beginning_vector <- c(104846075,138178553,138224999,124258900,104359835,54910245)
#end_vector <- c(104847927,138180781,138233197,124272013,104371673,54920096)

### NDUFAF3 and UQCRC1 are in the same place ### 
### NDUFA2 and HARS2 are in the same place ### 

# just 4 good regions with 6 genes
gene_vector <- c("NDUFAF3, UQCRC1","NDUFA2,HARS2","MRPS21","ATP5B")
chr_vector <- c("chr03","chr05","chr01","chr12")
beginning_vector <- c(104846075,138178553,124258900,54910245)
end_vector <- c(104847927,138180781,124272013,54920096)

genez_df <- cbind(gene_vector,chr_vector,beginning_vector,end_vector)
v_line_info <- cbind(paste(gene_vector,chr_vector,sep="; "),beginning_vector,end_vector)
colnames(v_line_info)[1] <- "gene"
colnames(v_line_info)[2] <- "start"
v_line_info<-as.data.frame(v_line_info)
v_line_info$start <- as.numeric(as.character(v_line_info$start))


# additional lines (same beginning for singletons, different for pairs)
gene_vector3 <- c("NDUFAF3, UQCRC1","NDUFA2,HARS2")
chr_vector3 <- c("chr03","chr05")
beginning_vector3 <- c(104359835,138224999)
end_vector3 <- c(104371673,138233197)
genez_df3 <- cbind(gene_vector3,chr_vector3,beginning_vector3,end_vector3)
v_line_info3 <- cbind(paste(gene_vector3,chr_vector3,sep="; "),beginning_vector3,end_vector3)
colnames(v_line_info3)[1] <- "gene"
colnames(v_line_info3)[2] <- "start"
v_line_info3<-as.data.frame(v_line_info3)
v_line_info3$start <- as.numeric(as.character(v_line_info3$start))

for(i in 1:nrow(genez_df)) {
    row <- as.vector(genez_df[i,])
    #  print(row[4]) # 1 is the gene name, 2 is the chr, 3 is the beginning, 4 is the end
    
    # open the data for this gene
    mydat1<-read.table(paste("./fst_bor/popgenWindows_",eval(row[2]), "_nem_hec_windowstats.csv", sep=""),sep=",", header=T)
    mydat2<-read.table(paste("./fst_bor/popgenWindows_",eval(row[2]), "_nem_mau_windowstats.csv", sep=""),sep=",", header=T)
    mydat3<-read.table(paste("./fst_bor/popgenWindows_",eval(row[2]), "_nem_ton_windowstats.csv", sep=""),sep=",", header=T)
    mydat4<-read.table(paste("./fst_bor/popgenWindows_",eval(row[2]), "_nem_nig_windowstats.csv", sep=""),sep=",", header=T)
    
    # merge the dfs for this gene
    temp_dat1 <- merge(mydat1, mydat2, by = c("scaffold","start","end"), all=T)
    temp_dat2 <- merge(mydat3, mydat4, by = c("scaffold","start","end"), all=T)
    all_dat <- merge(temp_dat1, temp_dat2, by = c("scaffold","start","end"), all=T)
    
    # reassign columnnames
    colnames(all_dat)[6] <- "nem"
    colnames(all_dat)[7] <- "hec"
    colnames(all_dat)[13] <- "mau"
    colnames(all_dat)[25] <- "nig"
    colnames(all_dat)[19] <- "ton"
    
    # subset the data
    # First and Third Column with All rows
    my_little_bit <- all_dat[,c("scaffold","start","end","nem","hec","ton","mau","nig")]
    my_plot_data <- subset(my_little_bit, (start >= as.numeric(row[3])-window)&(end <= as.numeric(row[4])+window))
    my_ready_to_plot <- melt(my_plot_data, id=c("scaffold","start","end"))
    colnames(my_ready_to_plot)[4] <- "Species"
    
    # add a colum for colors to match the depth plot
    my_ready_to_plot$color <- NA
    my_ready_to_plot$color[my_ready_to_plot$Species == "nem"] <- "blue"
    my_ready_to_plot$color[my_ready_to_plot$Species == "mau"] <- "orange"
    my_ready_to_plot$color[my_ready_to_plot$Species == "ton"] <- "black"
    my_ready_to_plot$color[my_ready_to_plot$Species == "hec"] <- "maroon"
    my_ready_to_plot$color[my_ready_to_plot$Species == "nig"] <- "purple"

    # add a gene column for use in faceting
    my_ready_to_plot$gene <- paste(row[1],row[2],sep="; ")
    
    # add each gene to a big_df dataframe
    big_df <- rbind(big_df,my_ready_to_plot) 
}

View(big_df)

p <- ggplot(big_df, aes(x=start/1000000, y=value, by=Species, color=Species)) +
  geom_smooth(span=.1,  size=0.5) + #se=F,
  labs( x = "Position (Mb)", y = expression(pi)) +
  ylim(0,0.1) +
  scale_y_continuous(breaks = c(0,0.05,0.10)) +
  facet_wrap( ~ gene, scales = "free_x") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.text = element_text (size = 20),
        axis.text.x = element_text(size = 20),
        axis.title = element_text (size = 20),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        legend.justification = "top")
  

png(filename="pi_closeup.png",
    width = 500, height = 300, units='mm', res = 100)  
    p + geom_vline(data = v_line_info, mapping = aes(xintercept = start/1000000),  color = "black", size=1) +
   geom_vline(data = v_line_info3, mapping = aes(xintercept = start/1000000),  color = "black", size=1)
dev.off()
```
