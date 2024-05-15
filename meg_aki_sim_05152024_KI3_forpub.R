library(tidyverse)
library(units)
library(rxode2)
library(cowplot)
library(ggforce)
library(ggsci)


MWCilDexHDot<-16.367
ng_dose<-50*MWCilDexHDot
DOSE=100#ng_dose

Cpl0<-14.39 #nmol/mL estimated plasma concentration at time zero
Vlumen<-8.3 #nL
Vpt<-91 #nl
Vplasma<-112.5*0.025*0.4 #mL per mouse 

gfr=40.9 #nL/min
meg=105.5#gfr/1.5
k10=6
k12=meg
k21=meg
zero<-0
casenum<-10
theta_multi<-data.frame(
  Kgfr=c(gfr,gfr*0.5, gfr*0.14,gfr*0.5,gfr*0.5,gfr*0.14,gfr*0.14,gfr, gfr,gfr),
  K12=c(meg,meg,meg,meg*0.8,meg*0.5,meg*0.8,meg*0.5,meg*0.8,meg*0.5,zero),  
  K21=c(meg,meg,meg,meg*0.8,meg*0.5,meg*0.8,meg*0.5,meg*0.8,meg*0.5,meg),
  K10=rep(k10,casenum),
  Vlumen=rep(8.3,casenum),
  Vpt=rep(91,casenum),
  Vplasma=rep(Vplasma,casenum),
  Cpl0=rep(Cpl0,casenum)
)

case_labels<-data.frame(
  case=factor(c("100% GFR, 100% megalin",
               "50% GFR, 100% megalin",
                "14% GFR, 100% megalin",
               "50% GFR, 80% megalin",
               "50% GFR, 50% megalin",
              "14% GFR, 80% megalin",
               "14% GFR, 50% megalin",
              "100% GFR, 80% megalin",
              "100% GFR, 50% megalin",
              "100% GFR, ZERO megalin"
              ),ordered=TRUE),
  GFR_Lrp2=factor(c("100:100",
                    "50:100",
                    "14:100",
                    "50:80",
                    "50:50",
                    "14:80",
                    "14:50",
                    "100:80",
                    "100:50",
                    "100:0"),ordered=TRUE),
  id=seq(1:casenum))




mod <- rxode2({
  # initial values for all "compartments"
  tubular_lumenA1(0) = 0;
  PT_cellA2(0) = 0;
  aucPTC_A2(0) = 0;
  plasma(0)=Cpl0/Vplasma;
 
  #
  Clumen= tubular_lumenA1/Vlumen;
  CPTcell= PT_cellA2/Vpt;
  Cplasma=plasma/Vplasma;
  
  # differential equations
  d/dt(plasma)          =(-Kgfr*Cplasma);
  d/dt(tubular_lumenA1) = (Kgfr*Cplasma) + (K21*CPTcell) -((K12 + K10)*Clumen);
  #d/dt(tubular_lumenA1) = (Kgfr) + (K21*PT_cellA2) -((K12 + K10)*tubular_lumenA1);
  d/dt(PT_cellA2)       = K12*Clumen  -K21*CPTcell;
  d/dt(aucPTC_A2)       =  CPTcell;
})




# state variables and their amounts at time 0 (the use of names is
# encouraged, but not required)
#inits <- c(tubular_lumenA1=0, PT_cellA2=0, aucPTC_A2=0)
events <- et(amountUnits = "nL", timeUnits = "hour")
events <- events |> 
  et(time = 0, amt = DOSE, addl = 0, ii = 0, cmt = "plasma")
events <- events |> et(time = c(seq(0,4,0.001)))#,seq(7,168,0.25)))
events 
events <- events |> et(id=1:nrow(case_labels))
events


#run it
out <- rxSolve(mod, theta_multi, events)

pct<-function(x,y)(x/y)


pct_cmax<-out%>%
  mutate_at(c("plasma", "tubular_lumenA1","PT_cellA2","aucPTC_A2"), ~pct(.,max(.[id==1])))%>%
  left_join(case_labels,by="id")%>%dplyr::arrange("GFR:LRP2")


#remove the units
pct_cmax_nu<-drop_units(pct_cmax)



# PLOTS FOR PAPER ---------------------------------------------------------
pct_cmax_nu<-pct_cmax_nu%>%mutate(case=factor(case,ordered=TRUE,levels=
                                                c("100% GFR, 100% megalin",
                                                  "100% GFR, 80% megalin",
                                                  "100% GFR, 50% megalin",
                                                  "50% GFR, 100% megalin",
                                                  "50% GFR, 80% megalin",
                                                  "50% GFR, 50% megalin",
                                                  "14% GFR, 100% megalin",
                                                  "14% GFR, 80% megalin",
                                                  "14% GFR, 50% megalin")))
pct_cmax_nu<-pct_cmax_nu%>%mutate(case=fct_relevel(case,"100% GFR, 100% megalin",
                                                   "100% GFR, 80% megalin",
                                                   "100% GFR, 50% megalin",
                                                   "50% GFR, 100% megalin",
                                                   "50% GFR, 80% megalin",
                                                   "50% GFR, 50% megalin",
                                                   "14% GFR, 100% megalin",
                                                   "14% GFR, 80% megalin",
                                                   "14% GFR, 50% megalin",
                                                   
                                                   ))

meg_varies<-pct_cmax_nu%>%filter(id%in%(c(1,8,9,10)))
gfr_varies<-pct_cmax_nu%>%filter(id%in%(c(1,2,3)))
all_plots<-pct_cmax_nu%>%filter(id%in%(c(1,2,3,8,9,10)))


addSmallLegend <- function(myPlot, pointSize = 1, textSize = 1, spaceLegend = 1, titleSize=1) { ##after https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"),
          plot.title = element_text(size = titleSize, face = "bold"))
}



pt<-ggplot(all_plots,aes(x=time,y=PT_cellA2,color=factor(GFR_Lrp2)))+
  #geom_line(linewidth=1.5)+
  geom_point(size=0.7)+
  #geom_smooth(method = "lm", formula = y ~ poly(x, 18), se = FALSE)+
  scale_color_aaas()+
  ylab("PT cell amount (% Control Cmax)") +
  theme_bw()+
  xlim(0,4)+ 
  ylim(0,1)+
  xlab("time (h)")+
  labs(colour = "GFR:Megalin (%)") +
  theme(legend.position=c(0.6,0.3))+
  ggtitle("All models")

pt


pt_meg_varies<-ggplot(meg_varies,aes(x=time,y=PT_cellA2,color=factor(GFR_Lrp2)))+
  #geom_line(linewidth=1.5)+
  geom_point(size=1)+
  scale_color_aaas()+
  ylab("PT cell amount (% Control Cmax)") +
  theme_bw()+
  xlim(0,4)+ 
  ylim(0,1)+
  xlab("time (h)")+
  labs(colour = "GFR:Megalin (%)") +
  theme(legend.position=c(0.6,0.3))+
  ggtitle("Megalin varies, GFR constant")

pt_meg_varies

pt_gfr_varies<-ggplot(gfr_varies,aes(x=time,y=PT_cellA2,color=factor(GFR_Lrp2)))+
  #geom_line(linewidth=1.5)+
  geom_point(size=1)+
  scale_color_aaas()+
  ylab("PT cell amount (% Control Cmax)") +
  theme_bw()+
  xlim(0,4)+ 
  ylim(0,1)+
  xlab("time (h)")+
  labs(colour = "GFR:Megalin (%)") +
  theme(legend.position=c(0.6,0.3))+
  ggtitle("GFR varies, megalin constant")

pt_gfr_varies

pt<-addSmallLegend(pt,1.5,10,0.5)
ggsave("pt.tiff",pt, units="in", width=2.2, height=3.2)

pt_meg_varies<-addSmallLegend(pt_meg_varies,1.5,10,0.5,8)
ggsave("pt_meg_varies.tiff",pt_meg_varies, units="in", width=2.2, height=3.2)

pt_gfr_varies<-addSmallLegend(pt_gfr_varies,1.5,10,0.5,8)
ggsave("pt_gfr_varies.tiff",pt_gfr_varies, units="in", width=2.2, height=3.2)
