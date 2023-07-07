#### PACKAGES ####
library(ggplot2)

#### SET UP CLADES AND THEME ####
clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

#basic density theme
theme_density <- theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background= element_blank(),
                       axis.line = element_line(colour = "black"),
                       axis.text.x = element_text(size = 15),
                       axis.text.y = element_text(size = 15),
                       axis.title = element_text(face = "bold",
                                                 size = 15),
                       legend.title = element_blank(),
                       plot.title = element_text(face = "bold",
                                                 size = 17,
                                                 hjust=0.5))

for(i in 1:5){
  
  #### LOAD DATA ####
  hpd.intervals <- read.csv(paste0("../outputs/subtrees/HPD_",clades[i],"_intervals.csv"))[,-1]
  raw.dat <- read.csv(paste0("../outputs/subtrees/proportions_",clades[i],"_raw.csv"))[,-1]
  
  #### PLOT ####
  SAF.overlap <- ggplot()+
    geom_density(aes(x=raw.dat$proportion,fill=raw.dat$category),
                 alpha=0.5)+
    geom_line(mapping=aes(x=hpd.intervals$x,y=hpd.intervals$y,
                          color=hpd.intervals$category),
              show.legend = F)+
    scale_y_continuous("Density")+
    scale_fill_viridis_d()+
    scale_color_viridis_d()+
    xlab("Proportion SAF")+
    labs(subtitle = paste0(clades[i]))+
    theme_density
  
  plot(SAF.overlap)
  
  #### SAVE PLOT ####
  ggsave(SAF.overlap,
         filename = paste0("../figures/observed_null_overlap_",clades[i],".pdf"),
         width = 7,
         height = 7,
         units = "in")
  
}



