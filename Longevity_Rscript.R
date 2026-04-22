#Load packages

library(lme4)
library(lmerTest)
library(carData)
library(car)
library(emmeans)
library(multcomp)
library(effects)
library(Hmisc)
library(ggbeeswarm)
library(ggplot2)

#Set contrasts to sum to zero
options(contrasts=c('contr.sum','contr.poly'))


#Import data
lifetime <- read.csv("Longevity_Data.csv", header = TRUE, strip.white = TRUE)


#Convert variables to factors
lifetime$hapgroup<- as.factor(lifetime$hapgroup)
lifetime$hapdup<- as.factor(lifetime$hapdup)
lifetime$location<- as.factor(lifetime$location)
lifetime$block<- as.factor(lifetime$block)
lifetime$bio<- as.factor(lifetime$bio)
lifetime$vial<- as.factor(lifetime$vial)
lifetime$temp<- as.factor(lifetime$temp)
lifetime$sex<- as.factor(lifetime$sex)

#Fitting a linear mixed effects model with fixed and random effects using lmer with random slopes to avoid pseudoreplication.
#The model includes explicity nested structure for haplotypes present within haplogroups
library(lme4)
lifetime_model1<-lmer (age ~ (hapgroup/hapdup) + sex + temp + block              # Main effects          
                       + sex:temp                                                # 2-way interactions
                       + (hapgroup/hapdup):sex                                   # 2-way interactions
                       + (hapgroup/hapdup):temp                                  # 2-way interactions                       
                       + (hapgroup/hapdup):sex:temp                              # 3-way interactions
                       + (1|vial)                                                # Random intercept for vial (Cohort ID)                                               
                       + (1 +sex + temp + sex:temp|bio),                         # Random intercept and slope for bio (Strain ID)
                       data=lifetime,                                            # Data frame
                       REML=T)

summary(lifetime_model1)
Anova(lifetime_model1, type="III", test="F")


## Plotting the haplotype and temperature interaction graph
## Emmeans were calculated using the emmeans package and were used later for plotting
emm_options(pbkrtest.limit=5000)

hap_temp_emm<- emmeans(lifetime_model1, specs = c("hapgroup","hapdup", "sex","temp","block"))
hap_temp_emm_dataframe<- as.data.frame(hap_temp_emm)


hap_temp_emm_dataframe$haplotype<-paste(hap_temp_emm_dataframe$hapgroup, hap_temp_emm_dataframe$hapdup, sep = " ")
hap_temp_emm_dataframe$haplotype<- as.factor(hap_temp_emm_dataframe$haplotype)

hap_temp_emm_dat<-emmip(lifetime_model1, temp ~ hapgroup/hapdup, CIs = TRUE, plotit = FALSE)
str(hap_temp_emm_dat)

colnames(hap_temp_emm_dataframe)[which(names(hap_temp_emm_dataframe)=="emmean")]<-"yvar"
colnames(hap_temp_emm_dataframe)[which(names(hap_temp_emm_dataframe)=="haplotype")]<-"xvar"

hap_temp_plot<-ggplot(data=hap_temp_emm_dat, aes(x=xvar, y=yvar, color=xvar)) +
  facet_grid(.~temp)+
  geom_errorbar(aes(ymin=yvar-SE, ymax=yvar+SE), width=.3, color="black")+
  geom_point(size=7, color="black")+
  geom_quasirandom(data = hap_temp_emm_dataframe, shape=18, size = 6, dodge.width = 0.1, alpha=0.5)+
  xlab("Mitochondrial haplotypes")+
  ylab("Longevity (Emmeans)")+
  theme_bw(base_size = 14, base_family = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
hap_temp_plot + guides(color=FALSE) + theme(axis.title.x = element_text(face="bold", colour="black", size=20),
                                            axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
                                            axis.title.y = element_text(face="bold", colour="black", size=20),
                                            axis.text.y = element_text(vjust=0.5, size=16),
                                            strip.text.x = element_text(size=20, face = "bold"))



##Plotting the haplogroup and temperature interaction
## Emmeans were calculated using the emmeans package and were used later for plotting
hapgroup_temp_emm<- emmeans(lifetime_model1, specs = c("hapgroup","hapdup", "sex","temp","block"))
hapgroup_temp_emm_dataframe<- as.data.frame(hap_temp_emm)

hapgroup_temp_emm_dat<-emmip(lifetime_model1, temp ~ hapgroup, CIs = TRUE, plotit = FALSE)
str(hapgroup_temp_emm_dat)

colnames(hapgroup_temp_emm_dataframe)[which(names(hapgroup_temp_emm_dataframe)=="emmean")]<-"yvar"

hapgroup_temp_plot<-ggplot(data=hapgroup_temp_emm_dat, aes(x=hapgroup, y=yvar, color=hapgroup)) +
  facet_grid(.~temp)+
  geom_errorbar(aes(ymin=yvar-SE, ymax=yvar+SE), width=.1, color="black")+
  geom_point(size=6, color="black")+
  geom_quasirandom(data = hapgroup_temp_emm_dataframe, shape=18, size = 6, dodge.width = 0.3, alpha=0.5)+
  xlab("Mitochondrial haplogroups")+
  ylab("Longevity (Emmeans)")+
  theme_bw(base_size = 14, base_family = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
hapgroup_temp_plot + guides(color=FALSE) + theme(axis.title.x = element_text(face="bold", colour="black", size=20),
                                                 axis.text.x  = element_text(size=20),
                                                 axis.title.y = element_text(face="bold", colour="black", size=20),
                                                 axis.text.y = element_text(vjust=0.5, size=16),
                                                 strip.text.x = element_text(size=20, face = "bold"))



##Calculating effect sizes using Omega squared values
library(effectsize)
omega_squared(lifetime_model1, alternative="two.sided")



##Supplementary material

#Creating a composite term with bio.temp to resolve singularity error
lifetime$bio.temp<- paste(lifetime$bio, lifetime$temp, sep=".")
lifetime$bio.temp<-as.factor(lifetime$bio.temp)


##The random intercept model to counter singularity error
lifetime_model2<-lmer (age ~ (hapgroup/hapdup) + sex + temp + block              # Main effects          
                       + sex:temp                                                # 2-way interactions
                       + (hapgroup/hapdup):sex                                   # 2-way interactions
                       + (hapgroup/hapdup):temp                                  # 2-way interactions                       
                       + (hapgroup/hapdup):sex:temp                              # 3-way interactions
                       + (1|vial)                                                # Random intercept for vial (Cohort ID)              
                       + (1|bio.temp:sex) + (1|bio.temp:block),                  # Random intercept for composite term::bio.temp(StrainID.temp)
                       data=lifetime,                                            # Data frame
                       REML=T)

summary(lifetime_model2)
Anova(lifetime_model2, type="III",test="F")
