# Environment prep --------------------------------------------------------
EnvPrep<-function (){
  library(tidyverse)
  library(ggpubr)
  library(Hmisc)
  library(ggplot2)
  library(ggprism)
  library(patchwork)
  library(scales)
  library(MASS)
  select<-dplyr::select
}

#first, load in the datasets and merge them, adding necessary metadata
EnvPrep()
#assign the treatments
puncta<-read_csv("Puncta.csv")


#puncta analysis
puncta<-puncta%>%mutate(norm_Intensity_MeanIntensity_OrigRed = rescale(Intensity_MeanIntensity_OrigRed, to = c(0, 1)))
puncta<-puncta%>%mutate(norm_Intensity_MeanIntensity_OrigGreen = rescale(Intensity_MeanIntensity_OrigGreen, to = c(0, 1)))

puncta_fl<-puncta%>%dplyr::select(AreaShape_Area,AreaShape_MaxFeretDiameter,norm_Intensity_MeanIntensity_OrigGreen,norm_Intensity_MeanIntensity_OrigRed)%>%
                                pivot_longer(cols=c(norm_Intensity_MeanIntensity_OrigGreen,
                                   norm_Intensity_MeanIntensity_OrigRed),
                                 names_prefix="norm_Intensity_MeanIntensity_Orig",
                                 names_to="Color")%>%
                                  mutate(Component=ifelse((Color=="Red"),"H-Dot","Dex"))
puncta_fl$Component=as.factor(puncta_fl$Component)
#Calculate the range
xlim <- range(puncta_fl$AreaShape_Area)
ylim <-range(puncta_fl$value)

#Genrate the kernel density for each group
newplot_data <- puncta_fl %>% group_by(Component) %>% do(Dens=kde2d(.$AreaShape_Area, .$value, n=100, lims=c(xlim,ylim)))
#Transform the density in  data.frame
newplot_data  %<>%  do(Component=.$Component, V=expand.grid(.$Dens$x,.$Dens$y), 
                    Value=c(.$Dens$z)) %>% do(data.frame(Component=.$Component,x=.$V$Var1, y=.$V$Var2, Value=.$Value))

ggplot(newplot_data, aes(x,y,z=Value, fill=Component, alpha=..level..))  + 
  stat_contour(geom="polygon", colour="grey")+
  scale_fill_manual(values=c("#05f62b", "#f605d0", "#66CC99"))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme_minimal()+
  xlab("Puncta area")+
  ylab("Normalized Signal")
  scale_alpha_continuous(range = c(0, 0.5))
  
ggplot(puncta_fl, aes(x=AreaShape_MaxFeretDiameter, y=value) ) +
  stat_density_2d(aes(alpha = ..level.., fill=Component), geom = "polygon", colour="white")+
  scale_fill_manual(values=c("#05f62b", "#f62b05", "#66CC99"))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+

  guides(alpha="none")+
 # scale_alpha_continuous(range = c(0, 1))+
  #xlim(0,90)+
  theme_minimal()+
  theme(legend.position = c(0.15, 0.85),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), )+
  xlab("Puncta Max Feret Diameter (px)")+
  
  ylab("Normalized Signal")


dodge<-position_dodge2(width=.2)
ggplot(data=filter(puncta_fl,(AreaShape_MaxFeretDiameter>13)&(AreaShape_MaxFeretDiameter<30)), aes(x=Component, y=value))+
  geom_point(position=dodge, aes(colour=Component), show.legend=FALSE)+
  scale_colour_manual(values=c("#05f62b", "#f62b05", "#66CC99"))+
  labs(x="Component", y="Puncta Intensity")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        panel.border = element_rect(colour = "grey", fill=NA, size=1))
  

#now plot fluorescence in cytoplasm, nucleus, and puncta.
cytoplasms<-read_csv("Cytoplasms.csv")
nuclei<-read_csv("Nuclei.csv")

cytoplasms<-cytoplasms%>%mutate(norm_Intensity_MeanIntensity_OrigRed = rescale(Intensity_MeanIntensity_OrigRed, to = c(0, 1)))
cytoplasms<-cytoplasms%>%mutate(norm_Intensity_MeanIntensity_OrigGreen = rescale(Intensity_MeanIntensity_OrigGreen, to = c(0, 1)))
nuclei<-nuclei%>%mutate(norm_Intensity_MeanIntensity_OrigRed = rescale(Intensity_MeanIntensity_OrigRed, to = c(0, 1)))
nuclei<-nuclei%>%mutate(norm_Intensity_MeanIntensity_OrigGreen = rescale(Intensity_MeanIntensity_OrigGreen, to = c(0, 1)))

cytoplasms<-cytoplasms%>%mutate(RedPctTotal=Intensity_MeanIntensity_OrigRed/sum(Intensity_MeanIntensity_OrigRed))
cytoplasms<-cytoplasms%>%mutate(GreenPctTotal=Intensity_MeanIntensity_OrigGreen/sum(Intensity_MeanIntensity_OrigGreen))
nuclei<-nuclei%>%mutate(RedPctTotal=Intensity_MeanIntensity_OrigRed/sum(Intensity_MeanIntensity_OrigRed))
nuclei<-nuclei%>%mutate(GreenPctTotal=Intensity_MeanIntensity_OrigGreen/sum(Intensity_MeanIntensity_OrigGreen))
puncta<-puncta%>%mutate(RedPctTotal=Intensity_MeanIntensity_OrigRed/sum(Intensity_MeanIntensity_OrigRed))
puncta<-puncta%>%mutate(GreenPctTotal=Intensity_MeanIntensity_OrigGreen/sum(Intensity_MeanIntensity_OrigGreen))


puncta_in<-puncta%>%select(Intensity_MeanIntensity_OrigRed,Intensity_MeanIntensity_OrigGreen,RedPctTotal,GreenPctTotal)
cytoplasm_in<-cytoplasms%>%select(Intensity_MeanIntensity_OrigRed,Intensity_MeanIntensity_OrigGreen,RedPctTotal,GreenPctTotal)
nuclei_in<-nuclei%>%select(Intensity_MeanIntensity_OrigRed,Intensity_MeanIntensity_OrigGreen,RedPctTotal,GreenPctTotal)
nuclei_in$compartment<-"n"
puncta_in$compartment<-"p"
cytoplasm_in$compartment<-"c"

compartment.intensity<-bind_rows(puncta_in,nuclei_in,cytoplasm_in)

ggplot(compartment.intensity, aes(x=compartment, y=RedPctTotal))+
  geom_boxplot()+
  theme_minimal()

ggplot(compartment.intensity, aes(x=compartment, y=GreenPctTotal))+
  geom_boxplot()+
  theme_minimal()
