# Plan

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = F,eval = T,warning = F,message = F)

```

```{r library}
library(xtable)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
library(gridExtra)
library(Hmisc)
library(readstata13)
library(survival)
library(ggplot2)
library(survminer)
library(incidence)
library(rms)
library(skimr)
library(finalfit)
library(mice)
```

```{r options}
options( tinytex.verbose = TRUE)
```

```{r }
data<-read.dta13("updata3.dta",nonint.factors = T,missing.type = T)
data[is.na(data)] <- NA  #575
```

```{r data}

attach(data)

data<-data %>%
  mutate(t=dateinstitution)%>%           
  mutate(t=coalesce(t,datelost))%>%
  mutate(t=coalesce(t,datedeath))%>%
  mutate(t=coalesce(t,datelastknownalive))

# sum(is.na(data$t))

data$cens<-ifelse(is.na(dateinstitution),0,1)


data$tt<-as.Date(as.character(data$t), format="%Y-%m-%d")-
                as.Date(as.character(data$datevisitbl), format="%Y-%m-%d")



attach(data)

```

```{r data2.1}

data2<- subset(data,!(diagnosiscat %in% c("Control")))
```

```{r agebl}

my_skim <- skim_with(base = sfl(n = length))

data%>%
  group_by(diagnosiscat)%>%
  select(agebl)%>%
  my_skim()%>%
  transmute(Variable=skim_variable,Group=diagnosiscat,n=n,Mean=numeric.mean,SD=numeric.sd,Min=numeric.p0,Median=numeric.p50,Max=numeric.p100,IQR=numeric.p75-numeric.p50)%>%
  kable(caption='\\label{tab:agebl} Summary statistics of the age at baseline visit.',booktabs=T,format = "latex",linesep="",digits=2) %>%
  kable_styling(font_size = 11,latex_options = "HOLD_position")

```


```{r ageboxplot, fig.height=5, fig.width=10, fig.align = "center",fig.cap = "\\label{ } Boxplots of patients' age at diagnosis and at onset of parkinsonian symptoms.", fig.pos ='H' }

plot1<-ggplot(data=data2,aes(y=agedx,x=diagnosiscat))+
  geom_boxplot()+
  ylab("Age at diagnosis of parkinsonian syndrome")+
  xlab("Diagnoses group")


plot2<-ggplot(data=data2,aes(agesymptomonset,x=diagnosiscat))+
  geom_boxplot()+
  ylab("Age at onset of parkinsonian symptoms")+
  xlab("Diagnoses group")

grid.arrange(plot1, plot2,ncol=2, nrow=1)

```

```{r sex}

data$sex.mf<-ifelse(sex==0,"male","female")


data%>%
  group_by(diagnosiscat,sex.mf)%>%
  summarise(n=n())%>%
  mutate(percentage=paste0(round(100*n/sum(n),2),'%'))%>%
  kable(caption='\\label{tab:sex} Summary statistics of the sex.',booktabs=T,format = "latex",linesep="",digits=2) %>%
  kable_styling(font_size = 11,latex_options = "HOLD_position")

```

```{r packyears}

my_skim <- skim_with(base = sfl(n = length))

data%>%
  group_by(diagnosiscat)%>%
  select(packyears)%>%
  my_skim()%>%
  transmute(Variable=skim_variable,Group=diagnosiscat,Mean=numeric.mean,SD=numeric.sd,Median=numeric.p50,Max=numeric.p100,IQR=numeric.p75-numeric.p50)%>%
  kable(caption='\\label{tab:packyears} Summary statistics of smoking history in pack-years.',booktabs=T,format = "latex",linesep="",digits=2) %>%
  kable_styling(font_size = 10,latex_options = "HOLD_position")

```

```{r  eversmok}

data$eversmok.yn<-ifelse(eversmok==0,"no","yes")

data%>%
  group_by(diagnosiscat,eversmok.yn)%>%
  summarise(n=n())%>%
  mutate(percentage=paste0(round(100*n/sum(n),2),'%'))%>%
  transmute(eversmok= eversmok.yn, n=n,percentage=percentage)%>%
  kable(caption='\\label{tab:eversmok} Summary statistics of ever smoking or not.',booktabs=T,format = "latex",linesep="",digits=2) %>%
  kable_styling(font_size = 11,latex_options = "HOLD_position")

```


```{r livesalone}

data$livesalone.yn<-ifelse(livesalone==1,"yes","no")

data%>%
  group_by(diagnosiscat,livesalone.yn)%>%
  summarise(n=n())%>%
  mutate(percentage=paste0(round(100*n/sum(n),2),'%'))%>%
  
  kable(caption='\\label{tab:livesalone} Summary statistics of whether participants live alone or not.',booktabs=T,format = "latex",linesep="",digits=2) %>%
  kable_styling(font_size = 11,latex_options = "HOLD_position")

```

```{r KM}

KM2<-survfit(Surv(tt/365,cens)~diagnosiscat,data=data,type="fh")

KM2

```

```{r KMplot,out.width="60%", fig.align = "center",fig.cap = "\\label { } Kaplan-Maier survival curve.", fig.pos ='H'}



ggsurvplot(survfit(Surv(tt/365,cens)~diagnosiscat,data=data),xlab="Survival days",data,surv.median.line = "hv",legend="top")



```


```{r missplot,out.width="60%", fig.align = "center",fig.cap = "\\label { } Cluster analysis showing which predictors tend to be missing on the same patients.", fig.pos ='H'}

plot(naclus(data2[,c(6:18)]))

```

 
```{r mistable1}

explanatory<-c("packyears","mmsebltotal","mmseyr1total")
dependent<- "cens"

fit1<-data2%>%
  ff_glimpse(dependent,explanatory)

data.frame(Variable=c("packyears","mmsebltotal","mmseyr1total"),n=c(308,289,257),missing_n=c(7,26,58),missing_percent=c("2.2%","8.3%","18.4%")) %>%
  kable(caption='\\label{tab:mistable1}  Summary of missing data.',booktabs=T,format = "latex",linesep="",digits=2) %>%
  kable_styling(font_size = 10,latex_options = "HOLD_position")
  

```

```{r misspattern,out.width="65%", fig.align = "center",fig.cap = "\\label { } Missing data pattern.", fig.pos ='H'}

data2 %>% 
  missing_pattern(explanatory)
```

```{r logistic}

data2$mispackyears<-ifelse(is.na(data2$packyears),1,0)  # if has missing then is 1, otherwise is 0

data2$mismmsebltotal<-ifelse(is.na(data2$mmsebltotal),1,0)

# data2$mismmseyr1total<-ifelse(is.na(data2$mmseyr1total),1,0)

data2<- data2 %>%
  mutate(mispmm = ifelse((data2$mispackyears|data2$mismmsebltotal)==1,1,0)) # missing either one will be 1


# sum(data2$mispmm) 

data2$NelsonAalen<-nelsonaalen(data2,tt,cens)  # Nelson-Aalen estimator

res.glm<-glm(mispmm~agebl+sex+livesalone+charlsonbl+updrsblpart3+hybl+sebl+NelsonAalen,data=data2,family =  binomial)

summary(res.glm)

```

```{r fcs}

fcs<-mice(data2[,6:18],seed=81420,m=10.5,print=F)
```

```{r densityplot1}

# xyplot(fcs, mmsebltotal ~ mmseyr1total | as.factor(.imp))

densityplot(fcs,~ mmsebltotal)

```

```{r densityplot2}

imp.packyears<-fcs[["imp"]][["packyears"]]

imp.packyears<-imp.packyears%>%
  rename(
    ten=10,
    nine=9
  )

packyearsnozero<-subset(data2,packyears>0)

plot(density(packyearsnozero$packyears),col="blue",main=" ", xlab = "packyears")
lines(density(imp.packyears$ten),col="red")


```


```{r merge1}

merge1<-merge(data2,imp.packyears,by="row.names")%>%
  mutate(packyears=ten)

merge1<-merge1[,2:32]

data.imp<-bind_rows(data2,merge1)%>%
  filter(is.na(packyears) == F)

```


```{r merge2}

imp.mmsebltotal<-fcs[["imp"]][["mmsebltotal"]]%>%
   rename(
    ten=10,
    nine=9
  )


merge2<-merge(data2,imp.mmsebltotal,by="row.names")%>%
  mutate(mmsebltotal=ten)

merge2<-merge2[,2:32]

data.imp2<-bind_rows(data.imp,merge2)%>%
  filter(is.na(mmsebltotal) == F)

```

```{r cox}
res.cox<-coxph(Surv(tt,cens)~agebl+sex+packyears+livesalone+charlsonbl+updrsblpart3+hybl+sebl+mmsebltotal,data=data.imp2)

summary(res.cox)
```

```{r testcox}
test.ph<-cox.zph(res.cox,transform = "km")
test.ph
```
