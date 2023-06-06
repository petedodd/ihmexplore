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
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

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
load(here('data/hbc30key.rda'))
load(here('data/N.rda'))
load(here('data/aprev.Rdata'))

## IHME CSV data
HT <- rh('indata/IHME-GBD_2019_DATA-8e052381-1.csv') #tot TB
HD <- rh('indata/IHME-GBD_2019_DATA-a0de4e33-1.csv') #disaggregated 19

## new
HS <- fread(here('indata/IHME-GBD_2019_DATA-b4d5febd-1.csv'))
HM <- fread(here('indata/IHME-GBD_2019_DATA-561f3057-1.csv'))
HX <- fread(here('indata/IHME-GBD_2019_DATA-3feca1e5-1.csv'))
HY <- Reduce(merge,list(HS[,.(location_name,sex_name,age_name,measure_name,year,ds=val)],
                        HM[,.(location_name,sex_name,age_name,measure_name,year,dm=val)],
                        HX[,.(location_name,sex_name,age_name,measure_name,year,dx=val)]))
HY[,val:=ds+dm+dx] #sum over DST


## prevalence snapshots
flz <- dir(path=here('indata/ihmeprev'),full.names = TRUE)
IP <- lapply(flz,fread)
IP <- rbindlist(IP)
IP <- merge(IP,agekey,by.x='age_name',by.y = 'age_group_name')
IP <- IP[,.(val=sum(val)),by=.(location_name,year,age)]
IP[,agey:=age]
IP[age=='65plus',agey:='65+']

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
ALL$agey <- ALL$age
ALL[age=='65plus',agey:='65+']

## -------------- incidence & notifications plot -------------
## incidence
aF1 <- ggplot(data=ALL,aes(x=agey,y=newrel,fill=sex)) +
  coord_flip() +
  facet_wrap(~name,scales='free')+
  geom_bar(stat='identity',aes(y=ifelse(sex=='M',newrel,-newrel)))+
  xlab('Age group')+ylab('Incidence or notifications')+
  scale_y_continuous(labels = absspace)+
  geom_bar(stat='identity',
           aes(x=agey,y=ifelse(sex=='M',best,-best)),
           fill='transparent',col=1)+
  geom_point(aes(x=agey,y=ifelse(sex=='M',ihme,-ihme)),
             size=2,col=1,
             shape=ifelse(is.na(ALL$newrel) | ALL$ihme>ALL$newrel,1,16))+ #using solid points for CDR>1
  theme_light()+theme(legend.position = 'none')

ww <- 16
ggsave(aF1,file=here('plots/aF1.pdf'),h=ww*0.8,w=ww)


## same as above but singling out 1 country (BGD)
tmp <- ALL[iso3=='BGD']

F1 <- ggplot(data=tmp,aes(x=agey,y=newrel,fill=sex)) +
  coord_flip() +
  facet_wrap(~name,scales='free')+
  geom_bar(stat='identity',aes(y=ifelse(sex=='M',newrel,-newrel)))+
  xlab('Age group')+ylab('Incidence or notifications')+
  scale_y_continuous(labels = absspace)+
  geom_bar(stat='identity',
           aes(x=agey,y=ifelse(sex=='M',best,-best)),
           fill='transparent',col=1)+
  geom_point(aes(x=agey,y=ifelse(sex=='M',ihme,-ihme)),
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

## --- figure 2A
m <- 1.3
F2a <- ggplot(UB,aes(who.frac.unc,ihme.frac.unc,label=iso3))+
  geom_point()+
  geom_text_repel(max.overlaps = 50)+
  coord_fixed()+
  scale_x_continuous(limits=c(0,m),label=percent)+
  scale_y_continuous(limits=c(0,m),label=percent)+
  geom_abline(slope=1,intercept = 0,lty=2,col=2)+
  theme_classic()+grids()+
  xlab('WHO fractional incidence uncertainty')+
  ylab('IHME fractional incidence uncertainty')

sf <- 1.0
ggsave(F2a,file=here('plots/F2a.pdf'),w=7*sf,h=7*sf)



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
         ## & metric_name=='Number'
         ## & cause_name=='Tuberculosis'
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


## IHME prevalence in 15 plus for relevant years
UP <- HY[location_name %in% hbc30key$location_name
         ## & metric_name=='Number'
         ## & cause_name=='Tuberculosis'
         & measure_name=='Prevalence'
        ,.(location_name,year,age_name,sex_name,val)]
UP <- dcast(data=UP,formula = location_name + year + sex_name ~ age_name,
            value.var = 'val')
UP[,prev.num.15plus:=`All ages` - `0-14 years`]
UP <- merge(UP,hbc30key,by='location_name')
UP <- dcast(data=UP,formula = iso3+location_name+year ~ sex_name,
            value.var = 'prev.num.15plus')

## join in notification data and compute ratios
UP <- merge(EP,UP,by=c('iso3','year'),all.x=TRUE,all.y=FALSE)
UP[,MFpn.ihme:=(Male/newrel_m15plus)/(Female/newrel_f15plus)]

## --- merge
PNB <- merge(UP[!is.na(MFpn.ihme),.(iso3,year,location_name,MFpn.ihme)],
             EP[!is.na(MFpn),.(iso3,year,MFpn)],
             by= c('iso3','year'))

## --- figure 2B
m2 <- 2.25
F2b <- ggplot(PNB,aes(MFpn,MFpn.ihme,label=paste0(iso3,', ',year)))+
  geom_point()+
  geom_text_repel(max.overlaps = 20)+
  coord_fixed()+
  scale_x_continuous(limits=c(0,m2))+
  scale_y_continuous(limits=c(0,m2))+
  geom_abline(slope=1,intercept = 0,lty=2,col=2)+
  geom_vline(xintercept = 1,lty=2,col='darkgrey')+
  geom_hline(yintercept = 1,lty=2,col='darkgrey')+
  theme_classic()+grids()+
  xlab('Empirical M:F ratio of prevalence/notifications')+
  ylab('IHME M:F ratio of prevalence/incidence')

sf <- 1.0
ggsave(F2b,file=here('plots/F2b.pdf'),w=7*sf,h=7*sf)


## combined figure 2
sf2 <- 0.77
F2 <- ggarrange(F2a,F2b,labels = c('A','B'))
ggsave(F2,file=here('plots/F2.pdf'),w=14*sf2,h=7*sf2)



## ======= other plots for appendix ===

## --- WHO vs IHME over time: incidence and mortality

## IHME - sum over both sexes
IC <- HY[location_name %in% hbc30key$location_name
         & measure_name=='Incidence'
         & age_name=="All ages"
        ,.(ihme.inc=sum(val)),
         by=.(location_name,year)]

IC <- merge(IC,hbc30key,by='location_name')

## WHO
IW <- TBC[iso3 %in% hbc30key$iso3
         ,.(iso3,year,who.inc=e_inc_num)]

## merge
IC <- merge(IC,IW,by=c('iso3','year'))
IC[year==2019,iso:=iso3]
IC[,pcd:=round(abs(ihme.inc/who.inc-1)*1e2)] #% difference
IC[!is.na(iso),isop:=paste0(iso,'\n',pcd,'%')]

## incidence
m <- max(IC$ihme.inc,IC$who.inc)
aF2a <- ggplot(IC,
               aes(who.inc,ihme.inc,label=isop,col=year,group=iso3))+
  geom_point()+ geom_line()+
  scale_x_sqrt(label=comma,limits = c(0,m))+
  scale_y_sqrt(label=comma,limits = c(0,m))+
  scale_color_viridis_c(direction = -1)+
  coord_fixed()+geom_abline(intercept = 0,slope=1,col=2)+
  geom_text_repel(max.overlaps = 30)+
  theme_light()+
  xlab('WHO incidence estimates') +
  ylab('IHME incidence estimates')+
  ggtitle('INCIDENCE estimates (sqrt scale)')


ggsave(aF2a,file=here('plots/aF2a.pdf'),h=10,w=10)


## IHME - sum over both sexes
MC <- HY[location_name %in% hbc30key$location_name
         & measure_name=='Deaths'
         & age_name=="All ages"
        ,.(ihme.inc=sum(val)),
         by=.(location_name,year)]

MC <- merge(MC,hbc30key,by='location_name')

## WHO
MW <- TBC[iso3 %in% hbc30key$iso3
         ,.(iso3,year,who.inc=e_mort_num)]

## merge
MC <- merge(MC,MW,by=c('iso3','year'))
MC[year==2019,iso:=iso3]
MC[,pcd:=round(abs(ihme.inc/who.inc-1)*1e2)] #% difference
MC[!is.na(iso),isop:=paste0(iso,'\n',pcd,'%')]

## mortality
m <- max(MC$ihme.inc,MC$who.inc,na.rm=TRUE)
aF2b <- ggplot(MC,
               aes(who.inc,ihme.inc,label=isop,col=year,group=iso3))+
  geom_point()+ geom_line()+
  scale_x_sqrt(label=comma,limits = c(0,m))+
  scale_y_sqrt(label=comma,limits = c(0,m))+
  scale_color_viridis_c(direction = -1)+
  coord_fixed()+geom_abline(intercept = 0,slope=1,col=2)+
  geom_text_repel(max.overlaps = 30)+
  theme_light()+
  xlab('WHO mortality estimates') +
  ylab('IHME mortality estimates')+
  ggtitle('MORTALITY estimates (sqrt scale)')

ggsave(aF2b,file=here('plots/aF2b.pdf'),h=10,w=10)


## --- duration by age/sex
HD[,unique(location_name)]
HD[,unique(age_name)]

ASIP <- HD[measure_name %in% c('Incidence','Prevalence')&
          metric_name=='Number',
          .(location_name,val,sex_name,age_group_name=age_name,age_id,
            measure_name)]

ASIP <- merge(ASIP,hbc30key,by='location_name',all.x=TRUE)
ASIP <- merge(ASIP,agekey,by='age_group_name',all.x=TRUE)
unique(ASIP[is.na(age),.(age_group_name,age)]) #check drop
ASIP <- ASIP[!is.na(age)]           #drop
ASIP[,table(measure_name)] #check

## aggregate
ASIP <- ASIP[,.(ihme=sum(val)),
             by=.(iso3,
                  name=location_name,
                  sex=ifelse(grepl('F',sex_name),'F','M'),
                  age,measure_name)]

ASIP <- dcast(ASIP,iso3+name+sex+age~measure_name,value.var = 'ihme')
ASIP$age <- factor(ASIP$age,levels=hagz,ordered=TRUE)
ASIP$sex <- factor(ASIP$sex,levels=c('M','F'),ordered=FALSE)
ASIP$agey <- ASIP$age
ASIP[age=='65plus',agey:='65+']


aF5 <-ggplot(ASIP,
             aes(agey,Prevalence/Incidence,col=sex,lty=sex,group=paste0(name,sex)))+
  geom_line()+
  facet_wrap(~name,scales='fixed')+
  scale_y_sqrt()+
  theme_light()+
  rot45+
  xlab('Age group') +
  ylab('IHME Prevalence/Incidence (square root scale)')

ggsave(aF5,file=here('plots/aF5.pdf'),h=15,w=15)


## --- duration over time
DY <- HY[location_name %in% hbc30key$location_name
         & age_name=="All ages"
        ,.(val=sum(val)),
         by=.(location_name,year,measure_name,sex=sex_name)]

DY <- dcast(DY,
            location_name+year+sex~measure_name,
            value.var='val')

aF6 <- ggplot(DY,
              aes(year,Prevalence/Incidence,col=sex,lty=sex))+
  geom_line()+
  facet_wrap(~location_name,scales='fixed')+
  scale_y_sqrt()+
  theme_light()+
  xlab('Year') +
  ylab('IHME Prevalence/Incidence (square root scale)')

ggsave(aF6,file=here('plots/aF6.pdf'),h=15,w=15)

## --- CFR vs age/sex
CFM <- HD[measure_name %in% c('Incidence','Deaths')&
          metric_name=='Number',
          .(location_name,val,sex_name,age_group_name=age_name,age_id,
            measure_name)]

CFM <- merge(CFM,hbc30key,by='location_name',all.x=TRUE)
CFM <- merge(CFM,agekey,by='age_group_name',all.x=TRUE)
unique(CFM[is.na(age),.(age_group_name,age)]) #check drop
CFM <- CFM[!is.na(age)]           #drop
CFM[,table(measure_name)] #check

## aggregate
CFM <- CFM[,.(ihme=sum(val)),
             by=.(iso3,
                  name=location_name,
                  sex=ifelse(grepl('F',sex_name),'F','M'),
                  age,measure_name)]

CFM <- dcast(CFM,iso3+name+sex+age~measure_name,value.var = 'ihme')
CFM$age <- factor(CFM$age,levels=hagz,ordered=TRUE)
CFM$sex <- factor(CFM$sex,levels=c('M','F'),ordered=FALSE)
CFM$agey <- CFM$age
CFM[age=='65plus',agey:='65+']


## TODO WHO include
aF4 <- ggplot(CFM,aes(agey,Deaths/Incidence,col=sex,group=paste(name,sex)))+
  geom_line() +
  facet_wrap(~name)+
  scale_y_continuous(label=percent)+
  theme_light()+rot45+
  xlab('Age group')+ylab('Mortality/Incidence (Case Fatality Ratio)')

ggsave(aF4,file=here('plots/aF4.pdf'),h=15,w=15)

## --- CFR vs CDR
CFY <- HY[age_name=='All ages' &
          measure_name!='Prevalence'
         ,.(val=sum(val)),
          by=.(location_name,measure_name,year)]
CFY <- dcast(CFY,location_name+year~measure_name,value.var = 'val' )
CFY[,CFR:=Deaths/Incidence]
CFY <- merge(CFY,hbc30key,by='location_name')

CFY <- merge(CFY,TBN[,.(iso3,year,c_newinc)],
             by=c('iso3','year'),all.x=TRUE,all.y=FALSE)
CFY[,CDR:=c_newinc/Incidence]
CFY[year==2019,iso:=iso3]

m <- max(CFY$CDR)

aF8 <- ggplot(CFY,
              aes(CFR,CDR,label=iso,col=year,group=iso3))+
  geom_hline(yintercept=1,col='darkgrey')+
  geom_point()+ geom_line()+
  scale_x_continuous(label=percent,limits = c(0,1))+
  scale_y_continuous(label=percent,limits = c(0,m))+
  scale_color_viridis_c(direction = -1)+
  coord_fixed()+
  geom_abline(intercept = 0,slope=1,col=2)+
  geom_text_repel(max.overlaps =  20 )+
  theme_light()+
  ylab('Case Detection Ratio (New & relapse notifications/Incidence)') +
  xlab('Case Fatality Ratio (Mortality/Incidence)')

ggsave(aF8,file=here('plots/aF8.pdf'),w=10,h=15)



##  --- HIV specific P:N ratio
## check data:
HD[metric_name=='Number' &
   sex_name=='Both' &
   age_name=="All ages" & 
   measure_name=='Incidence',sum(val)/1e6,by=cause_name][order(V1)]
HD[metric_name=='Number' &
   sex_name=='Both' &
   age_name=="All ages" & 
   measure_name=='Deaths',sum(val)/1e5,by=cause_name][order(V1)]



HD[,HIV:='-ve']
HD[grepl('HIV',cause_name),HIV:='+ve']


HH <- HD[metric_name=='Number' &
         sex_name=='Both' &
         age_name=="All ages" &
         measure_name!='Deaths',
         .(val=sum(val)),
         by=.(location_name,HIV,measure_name)]

HH <- dcast(HH,location_name+HIV~measure_name,value.var = 'val')
HH[,PI:=Prevalence/Incidence]

aF7 <- ggplot(HH,aes(location_name,PI,col=HIV))+
  geom_point(shape=1,size=2,stroke=1)+
  scale_x_discrete(limits=rev)+
  scale_color_manual(values=c('blue','red'))+
  coord_flip()+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = 'top')+
  xlab('Country')+
  ylab('Implied duration in years (Prevalence/Incidence)')+
  ylim(c(0,2.5))

ggsave(aF7,file=here('plots/aF7.pdf'),h=8,w=10)


## prevalence survey by age/sex
aprev[,agey:=gsub('plus','+',age)]
aprev[agey!='65+',
      agey:=paste(substr(agey,start=1,stop=2),
                  substr(agey,start=3,stop=4),sep='-')]

IP <- merge(IP,hbc30key,by='location_name',all.y=FALSE,all.x=TRUE)
IP[,age:=NULL]
IP <- merge(IP,aprev,by=c('iso3','year','agey'),all.x=FALSE,all.y=FALSE)
IP[,isoy:=paste0(iso3,', ',year)]

aF3 <- ggplot(IP,aes(agey,y=prev,ymin=prev.lo,ymax=prev.hi,group=isoy))+
  geom_point()+geom_errorbar(width=0)+geom_line()+
  geom_point(aes(y=1e5*val/pop),shape=4,col=2,stroke=1)+
  geom_line(aes(agey,1e5*val/pop,group=isoy),col=2)+
  facet_wrap(~isoy,scales='free')+
  theme_light()+rot45+
  xlab('Age group (years)')+
  ylab('TB prevalence per 100,000')

ggsave(aF3,file=here('plots/aF3.pdf'),h=10,w=10)



## WHO vs IHME by HIV
HHI <- TBC[iso3 %in% hbc30key$iso3
           & year==2019
          ,.(iso3,e_inc_num,e_inc_tbhiv_num)]
HHI[is.na(e_inc_tbhiv_num),e_inc_tbhiv_num:=0]
HHI <- HHI[,.(iso3,`-ve`=e_inc_num-e_inc_tbhiv_num,`+ve`=e_inc_tbhiv_num)]
HHI <- melt(HHI,id='iso3')
names(HHI)[3] <- 'WHO'
names(HHI)[2] <- 'HIV'
HHI <- merge(HHI,hbc30key,by='iso3')

HHI <- merge(HHI,HH[,.(location_name,HIV,IHME=Incidence)],by=c('location_name','HIV'))

aF9 <- ggplot(HHI,aes(HIV,abs(WHO/IHME-1),col=HIV,label=iso3))+
  geom_line(aes(group=iso3),col='grey',alpha=0.7)+
  geom_point(shape=1)+
  geom_boxplot(outlier.shape=NA,alpha=0.0)+
  geom_text_repel(show.legend=FALSE)+
  scale_colour_manual(values=c('blue','red'))+
  scale_y_sqrt()+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = 'top')+
  xlab('Country')+
  ylab('Incidence estimates in 2019: |WHO/IHME-1| (square root scale) ')

ggsave(aF9,file=here('plots/aF9.pdf'),h=7,w=7)
