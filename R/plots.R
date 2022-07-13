library(here)
library(data.table)
library(ggplot2)

## utilities
rh <- function(x) fread(here(x))
absspace <- function(x,...) {             #works
  format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}
`%ni%` <- Negate(`%in%`)

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
## WHO
AS <- TBA[iso3 %in% hbc30key$iso3 & age_group %in% agekey$age & risk_factor=='all' & sex!='a',
          .(iso3,sex=toupper(sex),age=age_group,best)]
hagz <- c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65plus")
AS$age <- factor(AS$age,levels=hagz,ordered=TRUE)
AS$sex <- factor(AS$sex,levels=c('M','F'),ordered=FALSE)



AS <- merge(AS,
            T1[,.(iso3=`IHME Location ID`,country=`Location Name`)],
            all.x=TRUE,all.y=FALSE,by='iso3')
AS[iso3=='MOZ' & age %ni% c('0-4','5-14'),newrel:=NA]

## IHME
ASI <- T5[year_id==2019]
ASI[,age_group_id:=as.integer(age_group_id)]
unique(ASI[,.(age_group_id,age_group_name)])
ASI[age_group_id <= 5,age:='0-4']
ASI[age_group_id %in% (6:7),age:='5-14']
ASI[age_group_id %in% (6:7+2),age:='15-24']
ASI[age_group_id %in% (6:7+2*2),age:='25-34']
ASI[age_group_id %in% (6:7+2*3),age:='35-44']
ASI[age_group_id %in% (6:7+2*4),age:='45-54']
ASI[age_group_id %in% (6:7+2*5),age:='55-64']
ASI[age_group_id >=18,age:='65plus']
unique(ASI[,.(age_group_id,age_group_name,age)])
ASI[sex=='Male',sex:='M']
ASI[sex=='Female',sex:='F']
ASI$age <- factor(ASI$age,levels=hagz,ordered=TRUE)
ASI$sex <- factor(ASI$sex,levels=c('M','F'),ordered=FALSE)
ASI <- ASI[,.(inc.ihme=sum(val)),by=.(iso3,sex,age)]

## join
ASA <- merge(AS,ASI,by=c('iso3','sex','age'))



## -------------- incidence & notifications plot -------------
## incidence
ggplot(data=ASA,aes(x=age,y=newrel,fill=sex)) +
  coord_flip() +
  facet_wrap(~country,scales='free')+
  geom_bar(stat='identity',aes(y=ifelse(sex=='M',newrel,-newrel)))+
  xlab('Age group')+ylab('Incidence or notifications')+
  scale_y_continuous(labels = absspace)+
  geom_bar(stat='identity',
           aes(x=age,y=ifelse(sex=='M',inc,-inc)),
           fill='transparent',col=1)+
  geom_point(aes(x=age,y=ifelse(sex=='M',inc.ihme,-inc.ihme)),
             shape=1,col=1,size=2)+
  theme_light()+theme(legend.position = 'none')+
  ggtitle('INCIDENCE, 2019: new & relapse notifications (red/right=M, green/left=F); WHO incidence estimates (bars); IHME incidence estimates (circles)')

ww <- 18
ggsave(here('plots/02asINCcompare.pdf'),h=ww*0.8,w=ww)
