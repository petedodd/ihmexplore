library(here)
library(data.table)
library(ggplot2)
library(ggrepel)
library(scales)
library(ggpubr)

## utilities
rh <- function(x) fread(here(x))
absspace <- function(x,...) {             #works
  format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}
`%ni%` <- Negate(`%in%`)

## age key
hagz <- c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65plus")
agekey <- c("<5 years","5-9 years","10-14 years","15-19 years",
            "20-24 years","25-29 years","30-34 years","35-39 years",
            "40-44 years","45-49 years","50-54 years","55-59 years",
            "60-64 years","65-69 years","70-74 years","75-79 years",
            "80+ years")
agekey <- data.table(age_group_name=agekey)
agekey[,age:=hagz[c(1,rep(2:7,each=2),rep(8,4))]]


## colots
who.col <- "#E69F00"
ihme.col <- "#56B4E9"
clz <- c('who'=who.col,'ihme'=ihme.col)
CLZ <- c('WHO'=who.col,'IHME'=ihme.col)

## load data
load(here('data/prev.rda'))
## load(here('data/agekey.rda'))
load(here('data/hbc30key.rda'))
load(here('data/N.rda'))

## IHME CSV data
HT <- rh('indata/IHME-GBD_2019_DATA-8e052381-1.csv') #tot TB
HY <- rh('indata/IHME-GBD_2019_DATA-56b8db41-1.csv') #by t for PI ratio
HD <- rh('indata/IHME-GBD_2019_DATA-a0de4e33-1.csv') #disaggregated 19

## WHO CSV data
TBA <- rh('indata/TB_burden_age_sex_2020-10-15.csv')
TBC <- rh('indata/TB_burden_countries_2020-10-15.csv')
TBN <- rh('indata/TB_notifications_2020-10-15.csv')


## ------------- data work on incidence/notifications by age

## --- WHO
## incidence estimates
AS <- TBA[iso3 %in% hbc30key$iso3 &
          age_group %in% agekey$age &
          risk_factor=='all' &
          sex!='a',
          .(iso3,sex=toupper(sex),age=age_group,best)]
AS$age <- factor(AS$age,levels=hagz,ordered=TRUE)
AS$sex <- factor(AS$sex,levels=c('M','F'),ordered=FALSE)

## notifications
nmz <- c('iso3','year',grep('newrel_m|newrel_f',names(TBN),value=TRUE))
NAS <- TBN[,..nmz]
NAS <- NAS[iso3 %in% hbc30key$iso3 & year==2019]
NAS <- melt(NAS,id=c('iso3','year'))
NAS[,sex:=ifelse(grepl('m',variable),'M','F')]
NAS[,agegp:=gsub('newrel_m|newrel_f','',variable)]
NAS <- NAS[!is.na(as.integer(agegp))]
## key to change ages
akey <- data.table(age=hagz)
akey[,agegp:=gsub('-','',age)]
akey[agegp=='65plus',agegp:='65']
NAS <- merge(NAS,akey,by='agegp',all.x=TRUE,all.y=FALSE)
NAS <- NAS[!is.na(age)]
NAS$age <- factor(NAS$age,levels=hagz,ordered=TRUE)
NAS$sex <- factor(NAS$sex,levels=c('M','F'),ordered=FALSE)

## join
AS <- merge(AS,NAS[,.(iso3,age,sex,newrel=value)],
            all=TRUE,
            by=c('iso3','sex','age'))
AS[,table(iso3)]

## IHME
ASI <- HD[measure_name=='Incidence'&
                         metric_name=='Number',
          .(location_name,val,sex_name,age_group_name=age_name,age_id,
            cause_name)]
ASI <- merge(ASI,hbc30key,by='location_name',all.x=TRUE)
ASI <- merge(ASI,agekey,by='age_group_name',all.x=TRUE)
unique(ASI[is.na(age),.(age_group_name,age)]) #check drop
ASI <- ASI[!is.na(age)]           #drop
ASI[,table(cause_name)] #check

## aggregate
ASI <- ASI[,.(ihme=sum(val)),
           by=.(iso3,
                name=location_name,
                sex=ifelse(grepl('F',sex_name),'F','M'),
                age)]

## join all
ALL <- merge(AS,ASI,by=c('iso3','sex','age'))
ALL$age <- factor(ALL$age,levels=hagz,ordered=TRUE)
ALL$sex <- factor(ALL$sex,levels=c('M','F'),ordered=FALSE)


## -------------- incidence & notifications plot -------------
## incidence
F1 <- ggplot(data=ALL,aes(x=age,y=newrel,fill=sex)) +
  coord_flip() +
  facet_wrap(~name,scales='free')+
  geom_bar(stat='identity',aes(y=ifelse(sex=='M',newrel,-newrel)))+
  xlab('Age group')+ylab('Incidence or notifications')+
  scale_y_continuous(labels = absspace)+
  geom_bar(stat='identity',
           aes(x=age,y=ifelse(sex=='M',best,-best)),
           fill='transparent',col=1)+
  geom_point(aes(x=age,y=ifelse(sex=='M',ihme,-ihme)),
             size=2,col=1,
             shape=ifelse(ALL$ihme<ALL$newrel,16,1))+ #using solid points for CDR>1
  theme_light()+theme(legend.position = 'none')

ww <- 18
ggsave(F1,file=here('plots/aF1.pdf'),h=ww*0.8,w=ww)



## same as above but singling out 1 country (BGD)
tmp <- ALL[iso3=='BGD']

F1 <- ggplot(data=tmp,aes(x=age,y=newrel,fill=sex)) +
  coord_flip() +
  facet_wrap(~name,scales='free')+
  geom_bar(stat='identity',aes(y=ifelse(sex=='M',newrel,-newrel)))+
  xlab('Age group')+ylab('Incidence or notifications')+
  scale_y_continuous(labels = absspace)+
  geom_bar(stat='identity',
           aes(x=age,y=ifelse(sex=='M',best,-best)),
           fill='transparent',col=1)+
  geom_point(aes(x=age,y=ifelse(sex=='M',ihme,-ihme)),
             size=2,col=1,
             shape=ifelse(tmp$ihme<tmp$newrel,16,1))+
  theme_light()+theme(legend.position = 'none')

ggsave(F1,file=here('plots/F1.pdf'),h=6,w=6)


## ------------- uncertainty data
## IHME
UI <- HT[location_name %in% hbc30key$location_name
         & metric_name=='Number'
         & cause_name=='Tuberculosis'
         & measure_name=='Incidence'
         & sex_name=='Both'
         & year==2019
         & age_name=="All ages"
        ,.(location_name,ihme.frac.unc=(upper-lower)/val)]
UI <- merge(UI,hbc30key,by='location_name')

## WHO
UW <- TBC[iso3 %in% hbc30key$iso3
          & year==2019
         ,.(iso3,who.frac.unc=(e_inc_num_hi-e_inc_num_lo)/e_inc_num)]

## merge
UB <- merge(UI,UW,by='iso3')

## --- figure 2B
m <- 1.3
F2b <- ggplot(UB,aes(who.frac.unc,ihme.frac.unc,label=iso3))+
  geom_point()+
  geom_text_repel(max.overlaps = 20)+
  coord_fixed()+
  scale_x_continuous(limits=c(0,m),label=percent)+
  scale_y_continuous(limits=c(0,m),label=percent)+
  geom_abline(slope=1,intercept = 0,lty=2,col=2)+
  theme_classic()+grids()+
  xlab('WHO fractional incidence uncertainty')+
  ylab('IHME fractional incidence uncertainty')
sf <- 1.5
ggsave(F2b,file=here('plots/F2b.pdf'),w=7*sf,h=7*sf)



## ----------------- P:N ratios
## 
## 2 A) Male to female ratios of empirical TB prevalence:notification ratios vs of GBD
## estimates prevalence:incidence ratio

## --- empirical 15plus prev rate: 15plus N rate by sex
## notifications
EP <- TBN[,.(iso3,year,newrel_m15plus,newrel_f15plus)]
EP <- EP[iso3 %in% hbc30key$iso3]
EP <- merge(EP,
            N[age=='15plus',.(iso3,year,PopMale,PopFemale)],
            by=c('iso3','year')) #merge pop denom for CNR
EP[,c('nm','nf'):=.(newrel_m15plus/PopMale,
                    newrel_f15plus/PopFemale)]

## prevalence
P <- dcast(prev,iso3 + year ~ sex,value.var='prev') #prev data
EP <- merge(EP,P,by=c('iso3','year'),all.x = FALSE,all.y=TRUE)
EP[,MFpn:=(m/nm) / (f/nf)]

## --- IHME prev:inc by sex, 15plus
PI <- HY[location_name %in% hbc30key$location_name
         & metric_name=='Number'
         & cause_name=='Tuberculosis'
         & measure_name %in% c('Incidence','Prevalence')
         & sex_name %in% c('Male','Female')
        ,.(location_name,year,age_name,sex_name,measure_name,val)]
PI <- merge(PI,hbc30key,by='location_name')
## right country-years
PI <- PI[paste0(iso3,year) %in% prev[,unique(paste0(iso3,year))]]

## reshape & calculate
PI <- dcast(PI,
            location_name+iso3+measure_name+sex_name~age_name,
            value.var = 'val')
PI[,o15:=`All ages`-`0-14 years`]
PI <- dcast(PI,location_name+iso3+sex_name~measure_name,value.var ='o15')
PI[,PtoI:=Prevalence/Incidence]
PI <- dcast(PI,location_name+iso3~sex_name,value.var = 'PtoI')
PI[,MFPtoI:=Male/Female]

## --- merge
PNB <- merge(PI[,.(iso3,location_name,MFPtoI)],
             EP[!is.na(MFpn),.(iso3,year,MFpn)],
             by='iso3')


## --- plot
m2 <- 2.25
F2a <- ggplot(PNB,aes(MFpn,MFPtoI,label=paste0(iso3,', ',year)))+
  geom_point()+
  geom_text_repel(max.overlaps = 20)+
  coord_fixed()+
  scale_x_continuous(limits=c(0,m2))+
  scale_y_continuous(limits=c(0,m2))+
  geom_abline(slope=1,intercept = 0,lty=2,col=2)+
  theme_classic()+grids()+
  xlab('Empirical M:F ratio of P:N ratio')+
  ylab('IHME M:F ratio of P:I ratio')

sf <- 1.5
ggsave(F2a,file=here('plots/F2a.pdf'),w=7*sf,h=7*sf)

## combined figure 2
sf2 <- 1
F2 <- ggarrange(F2a,F2b,labels = c('A','B'))
ggsave(F2,file=here('plots/F2.pdf'),w=14*sf,h=7*sf)
