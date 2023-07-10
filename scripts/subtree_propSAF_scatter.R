#### PACKAGES ####
library(ggplot2)

#### SET UP CLADES AND THEME ####
clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

theme_violin_bar <- theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background= element_blank(),
                          axis.line = element_line(colour = "black"),
                          axis.text.x = element_blank(),
                          axis.text.y = element_text(size = 10),
                          axis.title = element_text(face = "bold",
                                                    size = 10),
                          legend.title = element_blank(),
                          plot.title = element_text(face = "bold",
                                                    size = 17,
                                                    hjust=0.5))
#### SET UP DATA ####
clades_observed <- data.frame()
clades_null <- as.data.frame(matrix(ncol=2,nrow=5))

#loop through csv files for clades
for(i in 1:5){
  #get raw data
  raw.dat <- read.csv(paste0("../outputs/subtrees/proportions_",clades[i],"_raw.csv"))[,-1]
  
  clades_observed <- rbind(clades_observed,
                           cbind(raw.dat[1:100,1],
                                 clades[i]))
  
  clades_null[i,] <- c(mean(raw.dat[101:200,1]),clades[i])
}

colnames(clades_null) <- c("mean","clade")
colnames(clades_observed) <- c("prop","clade")
clades_observed$prop <- as.numeric(clades_observed$prop)
clades_null$mean <- as.numeric(clades_null$mean)

clades_observed$clade <- factor(clades_observed$clade,
                                   levels=unique(clades_observed$clade))
clades_null$clade <- factor(clades_null$clade,
                               levels=unique(clades_null$clade))
levels(clades_observed$clade) <- levels(clades_null$clade) <- c("Carnivora",
                                                                "Artiodctlya",
                                                                "Yangochiroptera",
                                                                "Rodentia",
                                                                "Primatomorpha")

#### PLOT ####
subtree_scatter <- ggplot()+
  geom_jitter(data=clades_observed,
              mapping = aes(x=clades_observed$clade,
                  y=clades_observed$prop,
                  col=clades_observed$clade),
              alpha=0.5,size=0.75,position=position_jitterdodge())+
  geom_point(mapping=aes(x=factor(clades_null$clade),y=clades_null$mean,
                        color=clades_null$clade),
            show.legend = F,shape=95,size=15)+
  scale_y_continuous(name="Proportion SAF")+
  xlab("Clade")+
  scale_color_viridis_d()+
  theme_violin_bar

ggsave(subtree_scatter,
       filename = paste0("../figures/clades_scatter.pdf"),
       width = 7,
       height = 7,
       units = "in")














