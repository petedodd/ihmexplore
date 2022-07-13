library(here)
library(data.table)
library(ggplot2)

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
load(here('data/agekey.rda'))
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
ggplot(data=ALL,aes(x=age,y=newrel,fill=sex)) +
  coord_flip() +
  facet_wrap(~name,scales='free')+
  geom_bar(stat='identity',aes(y=ifelse(sex=='M',newrel,-newrel)))+
  xlab('Age group')+ylab('Incidence or notifications')+
  scale_y_continuous(labels = absspace)+
  geom_bar(stat='identity',
           aes(x=age,y=ifelse(sex=='M',best,-best)),
           fill='transparent',col=1)+
  geom_point(aes(x=age,y=ifelse(sex=='M',ihme,-ihme)),
             shape=1,col=1,size=2)+
  theme_light()+theme(legend.position = 'none')+
  ggtitle('INCIDENCE, 2019: new & relapse notifications (red/right=M, green/left=F); WHO incidence estimates (bars); IHME incidence estimates (circles)')

ww <- 18
ggsave(here('plots/notinc.pdf'),h=ww*0.8,w=ww)
