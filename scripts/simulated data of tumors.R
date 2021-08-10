## Data is from:
## A Near-Infrared Phosphorescent
## Nanoprobe Enables Quantitative, Longitudinal Imaging of
## Tumor Hypoxia Dynamics during Radiotherapy
## --
## Zheng et. al.
## doi: 10.1158/0008-5472.CAN-19-0530
library(here)
library(ggplot2)
library(tidyverse)
library(mgcv)
library(patchwork)
#####
data<-read.csv(here("data","tumor_data.csv"))
data$Group<-as.factor(data$Group)
#plot data
ggplot(data,aes(x=Day,y=Volume))+geom_line()+facet_wrap(~Group)

#This function simulates data for the tumor data using default parameters of 10 observations per timepoint (as in the
# paper)and Standard deviation (sd) of ~15%.
#Because tumor volume cannot go below 0%, data is  generated with a cutoff value of 0.0001 (the "StO2_sim")

simulate_data <- function(dat, n = 10, sd = 15) {
    dat_sim <- dat %>%
        slice(rep(1:n(), each = n)) %>%
        group_by(Group, Day) %>%
        mutate(
            Vol_sim = pmax(rnorm(n, Volume, sd), 0.0001),
            subject=rep(1:10),
            subject=factor(paste(subject, Group, sep = "-"))
        ) %>%
        ungroup()

    return(dat_sim)
}

n <- 10 #number of observations (from paper)
sd <- 40 #mm3 approximate sd from paper
dat_sim <- simulate_data(data, n, sd)
breaks_s<-c(0,2,4,6,8,10,12,14)

#plotting simulated data
f1<-ggplot(dat_sim, aes(x = Day, y = Vol_sim, color = Group)) +
    geom_point(show.legend=FALSE,size=1.5,alpha=0.5)+
    stat_summary(aes(y = Vol_sim,
                     group=Group),
                 fun=mean, geom="line",size=1,show.legend = FALSE)+
    labs(y=expression(atop(Volume (mm^3),'(simulated)')))+
    theme_classic()+
    theme(
        axis.text=element_text(size=22)
    )+
    scale_x_continuous(breaks=breaks_s)

f1

##### fit GAM and rm-ANOVA

#GAM for StO2

m1 <- gam(Vol_sim ~ Group+s(Day, by = Group, k = 10),
          method='REML',
          data  = dat_sim)



#linear model
lm1<-lm(Vol_sim ~ Day + Group + Day * Group, data = dat_sim)


#creates a dataframe using the length of the covariates for the GAM
gam_predict <- expand_grid(Group = factor(c("T1", "T2")),
                           Day = seq(0, 15, by = 0.1),
                           subject=factor(rep(1:10)))

#creates a dataframe using the length of the covariates for rm-ANOVA
lm_predict<-expand_grid(Group = factor(c("T1", "T2")),
                        Day = c(0:15),
                        subject=factor(rep(1:10)),
)
lm_predict$subject<-factor(paste(lm_predict$subject, lm_predict$Group, sep = "-"))

#adds the predictions to the grid and creates a confidence interval for GAM
gam_predict<-gam_predict%>%
    mutate(fit = predict(m1,gam_predict,se.fit = TRUE,type='response')$fit,
           se.fit = predict(m1, gam_predict,se.fit = TRUE,type='response')$se.fit)

#using lm
lm_predict<-lm_predict%>%
    mutate(fit = predict(lm1,lm_predict,se.fit = TRUE,type='response')$fit,
           se.fit = predict(lm1, lm_predict,se.fit = TRUE,type='response')$se.fit)

#plot smooths and confidence interval for GAM
f3<-ggplot(data=dat_sim, aes(x=Day, y=Vol_sim, group=Group)) +
    geom_point(aes(color=Group),size=1.5,alpha=0.5,show.legend = FALSE)+
    geom_ribbon(aes( x=Day,ymin=(fit - 2*se.fit),
                     ymax=(fit + 2*se.fit),
                     fill=Group
    ),
    alpha=0.3,
    data=gam_predict,
    show.legend=FALSE,
    inherit.aes=FALSE) +
    geom_line(aes(y=fit,
                  color=Group),
              size=1,data=gam_predict,
              show.legend = FALSE)+
    #facet_wrap(~Group)+
    labs(y=expression(atop(StO[2],'complete')))+
    scale_x_continuous(breaks=breaks_s)+
    theme_classic()+
    theme(
        axis.text=element_text(size=22)
    )

f3

#plot linear fit for rm-ANOVA
f4<-ggplot(data=dat_sim, aes(x=Day, y=Vol_sim, group=Group)) +
    geom_point(aes(color=Group),size=1.5,alpha=0.5,show.legend = FALSE)+
    geom_ribbon(aes( x=Day,ymin=(fit - 2*se.fit),
                     ymax=(fit + 2*se.fit),fill=Group),
                alpha=0.3,
                data=lm_predict,
                show.legend = FALSE,
                inherit.aes=FALSE) +
    geom_line(aes(y=fit,
                  color=Group),
              size=1,data=lm_predict,
              show.legend = FALSE)+
    #facet_wrap(~Group)+
    labs(y=expression(paste('Volume'(mm^3))))+
    scale_x_continuous(breaks=breaks_s)+
    theme_classic()+
    theme(
        axis.text=element_text(size=22)
    )
f4

