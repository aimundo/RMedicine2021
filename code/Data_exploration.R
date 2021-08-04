library(here)
library(tidyverse)
library(broom)
library(mgcv)
library(ggsci)
library(ggpubr)
data<-read.csv(here('data','mice_weights.csv')) #load dataframe
data<-data[,-c(40:83)] #remove blank columns
data$ID<-as.factor(data$ID) #make ID a factor
data_long<-pivot_longer(data=data,!ID,values_to="weight") #pivot to long format for analysis

data_long<-data_long %>% mutate(Day=as.numeric(str_replace(name,"X",""))) #remove the Xs from "name" and name the column "Day"

data_long<-data_long %>%
    select(-("name")) #remove the "name" column

#make a plot
ggplot(data_long,aes(x=Day,y=weight))+geom_point()

test<-gam(weight~s(Day,k=20),data=data_long)
plot(test)
