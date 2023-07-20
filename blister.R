# > version
# _                           
# platform       x86_64-w64-mingw32          
# arch           x86_64                      
# os             mingw32                     
# system         x86_64, mingw32             
# status                                     
# major          3                           
# minor          3.1                         
# year           2016                        
# month          06                          
# day            21                          
# svn rev        70800                       
# language       R                           
# version.string R version 3.3.1 (2016-06-21)
# nickname       Bug in Your Hair  
#Script for extracting dryness indexes associated with fire from water balance
#Fire polygons are from MTBS data base
#Water balance is daily using Daymet climate input to water balance model "Daymet_Penman_Hamon_Batch.xlsm"
#Steps
#Identify ignition date in MTBS database for polygon
#Find that date in the daily water balance file associated with same polygon
#Compute backward sums or other stats for each water balance metric of interest

#set default R package library path
#.libPaths(c("C:\\Users\\dthoma\\R packages", .libPaths()))

library(xts)
library(TTR)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(devtools)
library(easyGgplot2)
library(extrafont)
library(data.table)
library(lubridate)
library(lme4)
require(randomForest)
require(MASS)#Package which contains the Boston housing dataset
library(corrplot)
library(MuMIn)
library(lmtest)
#install.packages("devtools")
# devtools::install_github("cardiomoon/ggiraphExtra")
# library(ggiraphExtra)
#loadfonts()#only have to do this once
# If you want to output to .ps files instead of .pdf, use:
# loadfonts(device="postscript")
# After the fonts are registered with R's PDF device, you can create figures with them. 
fonts()
windowsFonts(Times=windowsFont("TT Times New Roman"))

#next line is very slow, so only do once
#font_import(paths = NULL, recursive = TRUE, prompt = TRUE,pattern = NULL)

####################master################################
##working with "master" data base depracated b/c attempts to unambiguously sort infected from not infected trees
##using boolian strings in Excel was not successful###
##bounce down to sectoin called transition###########
# setwd("C:\\David\\Water balance\\blister rust\\")
# master<-read.csv("wbp_sites.csv",header=TRUE, stringsAsFactors=FALSE,na.strings=c(""," ","na","no trees"))#remove all the garbage that can't be plotted
# head(master)
# #replace "no trees" with na
# #master[master == "no trees"] <- NA
# #drop unused columns
# master<-master[,c(1:6,9:17)];head(master)#drop unused columns
# #rename column headers
# names(master)[]<-c("site","lat","lon","slope","aspect","whc","elev","T1","T1count","T2","T2count","T3","T3count","T4","T4count")
# head(master);str(master)
# unique(master$T2)
# #master[,8:15] <- lapply(master, function(x) as.numeric(as.character(x)))
# meltmaster<-melt(master,id=c("site","lat","lon","slope","aspect","whc","elev","T1count","T2count","T3count","T4count"))
# head(meltmaster)
# 
# 
# #elevation effect?
# eplot<-ggplot(meltmaster,na.rm=TRUE, aes(x=elev, y=value)) + geom_point(size = 5) + # geom_bar(stat="identity", position=position_dodge())+#
#   facet_wrap(~variable, scales="free")+geom_smooth(method=lm)+ 
#   theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   theme(axis.text=element_text(size=24,face="bold",family="Courier New"),
#         axis.title=element_text(size=28,face="bold",family="Courier New"),
#         title=element_text(size=20,face="bold",family="Courier New"))+ labs(title = "", x = NULL, y = "Proportion infected", color = NULL)
# eplot
# 
# fit1<-lm(data=master, T1~elev);summary(fit1)
# plot(master$elev,master$T1)
# abline(fit1)


########################################################
###############transition###############################
setwd("C:\\David\\Water balance\\blister rust\\")
transition<-read.csv("transition.csv",header=TRUE, stringsAsFactors=TRUE,na.strings=c(""," ","na","no trees"))#remove all the garbage that can't be plotted
head(transition)
names(transition)
#replace "no trees" with na
#trans[trans == "no trees"] <- NA
#drop unused columns
transition<-transition[,c(1,5,6,7,8,9,12:14,15,18,21,22,25,26)];head(transition)#drop unused columns
#rename column headers
names(transition)[]<-c("panel","site","treeid","year","infected","trans","height_class","recentdbh","dbh_class","status","elev","long","lat","slope","aspect")
head(transition)
str(transition)
unique(transition$trans)

#the transitions data set has multiple observations for each tree.  We need to pull out just one observation for each tree
#that can help determine whether it was infected or not, so each tree appears in the resulting data set just once
length(unique(transition$treeid))
# infatap<-subset(transition, trans=="Y"|trans=="first");head(infatap)#infected anytime anyplace
# deduped<-unique(infatap[,c(3,6)]);head(deduped)#remove duplicates eg. a site could be "Y" for several years
# a<-!(infatap %in% transition)#index of not infected at any time or place
# notinfatap<-transition[a];head(notinfatap)#data set of not infected anytime or place
# #combine the infected trees with non infected trees into one data set
# inf_status<-rbind(infatap,notinfatap);head(inf_status)
# nrow(inf_status)


#because climate data end with 2017 ignore trees observed in 2018
transub1<-subset(transition, year<2019);head(transub1);length(unique(transub1$treeid))

#pull out unique and recent records for any tree ever infected
i="10596.1.59"
inf<-NULL
uninf<-NULL
infamb<-NULL
#infatap<-subset(transition, trans=="Y"|trans=="first");head(infatap)#inf anytime anyplace
#create a dataframe of infected and uninfected trees based on all observations
#rule 1: tree must be visited at least twice = 4 years between observations.  eliminates new trees since 2014
#rule 2: any visit marked as infected based on 3 indicators, ie. any visit "Y" for infected
#rule 3: any visit marked "N".  These are what's left after pulling "Y"'s in previous if then step
#the "any" Y or N aspect is important in deciding both Y and N status ultimately and it depends on the evaluation order here with Y's first
#NA almost always means tree died, so N's with NA's mean a tree was infected before it died
for(i in unique(transub1$treeid)){
  sub<-subset(transub1, treeid==i);sub
  subinf<-sub$infected;subinf
  if(length(subinf)<2){#must have 2 or more visits across 4 years to determine BR status
    mxyr<-max(sub$year);mxyr#this is the most recent observation
  sub1<-subset(sub,year==mxyr);sub1
  infamb<-rbind(sub1,infamb)#build set of trees  with ambiguous status   
  #test<-unique(subinf);test
    }else if(any(subinf=="Y")){
    mxyr<-max(sub$year);mxyr#this is the most recent observation
    sub1<-subset(sub,year==mxyr);sub1
    inf<-rbind(sub1,inf)#build set of infected trees    
    }else if(any(subinf=="N")){#if the only unique infection status is "N", then uninfected through time so find max year and capture that record
    mxyr<-max(sub$year);mxyr#this is the most recent observation used to pull just one observation record
    sub1<-subset(sub,year==mxyr);sub1
    uninf<-rbind(sub1,uninf)#build set of unabmiguously uninfected trees 
    }
    }
    
    
head(inf);nrow(inf)#;write.csv(inf,"inf.csv")#tree was recorded as infected at any one of 2 or more site visits
head(uninf);nrow(uninf)#;write.csv(uninf,"uninf.csv")#tree was uninfected while it was alive 
head(infamb);nrow(infamb)#;write.csv(infamb,"infamb.csv")#tree visited only once so not used in analysis of blister rust distribution

tree_cnt<-nrow(inf)+nrow(uninf)+nrow(infamb) 
tree_cnt
inf_pct<-nrow(inf)/tree_cnt;inf_pct
uninf_pct<-nrow(uninf)/tree_cnt;uninf_pct
amb_pct<-nrow(infamb)/tree_cnt;amb_pct

   
head(uninf);unique(uninf$infected)
head(infamb)

#check all trees accounted for in one of two classes infected or uninfected
length(unique(transub1$treeid))#this and next row should be equal
sum(nrow(uninf)+nrow(inf)+nrow(infamb))

test<-subset(uninf, infected=="Y");test#this should be empty


#create a reclassified infection category and append it to uninf data frame
#this helps make calls of infected vs. uninfected unambiguous since some trees are NA after they die and some infected become uninfected over time
infrc<-rep("N",nrow(uninf));infrc
uninf2<-cbind(uninf,infrc)
head(uninf2);nrow(uninf2)

head(inf);unique(inf$infected)
test<-subset(inf, infected=="N");test#this is the list of trees where infection occurred on some visit but may not be infected last visit
#create a reclassified infection category for the infected trees and append it to uninf data frame
infrc<-rep("Y",nrow(inf));infrc
inf2<-cbind(inf,infrc)
head(inf2);nrow(inf2)

head(infamb);unique(infamb$infected)
test<-subset(infamb, infected=="N");test#this is the list of trees where infection occurred on some visit but may not be infected last visit
#create a reclassified infection category for the infected trees and append it to uninf data frame
infrc<-rep("NA",nrow(infamb));infrc
infamb2<-cbind(infamb,infrc)
head(infamb2);nrow(infamb2)

#logic check should be true
length(unique(transub1$treeid))==nrow(uninf2)+nrow(inf2)+nrow(infamb2)

#now merge infected trees with uninfected trees to form a single data frame for analysis
br<-rbind(inf2,uninf2)#,infamb2)
nrow(br)


#plot spatial patterns
#where did infenction occur in each year?  year here is last observed year, not first year infection detected
ggplot(br,aes(x=long, y=lat)) + geom_point() + facet_wrap(~infrc)
ggplot(br,aes(x=long, y=lat)) + geom_point(aes(colour = (aspect)))
ggplot(br,aes(x=long, y=lat)) + geom_point(aes(colour = (elev)))#lower elev up north, higher elevation down south

#ggplot(br,aes(x=long, y=lat, color=infrc)) + geom_point() + facet_wrap(~year)#here color doesn't mean anything except last year visited

# #some unnecessary subsetting of just live or just infected going on here
# subyn<-subset(br, status=="L" & infrc=="Y" | infrc=="N");head(subyn)#selecting only live trees not necessary
# suby<-subset(br, infrc=="Y" );head(suby)#status=="L" & #selecting only y infected not necessary
# 
# #create elevation bins by 200m intervals
# subyn$elbin<-.bincode(subyn$elev, breaks=c(2400, 2600,2800,3000,3200,3400), TRUE);head(subyn)
# suby$elbin<-.bincode(suby$elev, breaks=c(2400, 2600,2800,3000,3200,3400), TRUE);head(suby)

#create elevation bins by 200m intervals
br$elbin<-.bincode(br$elev, breaks=c(2400, 2600,2800,3000,3200,3400), TRUE);head(br)

#calculate standard errors of infected vs. non-infected in the br$infrc2 column which is numberic form of infrc column
br$infrc2<-as.numeric(br$infrc);head(br);tail(br)
unique(br$infrc2)

max(br$elev)
min(br$elev)
# convert to character vector
elevc <- as.character(br$elev)
# fill any refs with "999"
#temp[grep("REF", toupper(temp)) ] <- "999"
# use cut to get desired categories
br$elev_bin<-cut(as.numeric(elevc), breaks=c(2400,2600,2800,3000,3200), labels=c("2400-2600 m", "2600-2800 m", "2800-3000 m",">3000 m"),include.lowest=T)
head(br)
write.csv(br,"br.csv")
write.csv(br,"br.csv")

names(br)
elev_stats<-tapply(br$dbh_class,br$elev_bin,count);elev_stats

#br %>% group_by(elbin,dbh_class) %>% summarize(inf=count(infrc))#, sum=sum(dt))
as.data.frame(tapply(br$infrc, c(br$dbh_class,br$elev_bin),summary))

#summarize using data.table #https://s3.amazonaws.com/assets.datacamp.com/img/blog/data+table+cheat+sheet.pdf
dt<-data.table(br);dt
dt[ ,count(infrc)]

dt[ ,count(infrc), by=.(elev_bin,dbh_class)]

#in order to plot bars in geom_bar side-by-side we have to melt the data frame
#aql <- melt(airquality, id.vars = c("month", "day"),variable.name = "climate_variable", value.name = "climate_value")

br_melt<-melt(br, id.vars = c("panel","site","treeid","year" ,"infected","trans" ,"height_class","recentdbh","dbh_class" ,"status",
                              "elev", "long","lat", "slope","aspect","elbin", "infrc2","elev_bin"))#,"infrc"
head(br_melt)
#checks out in Excel br_chart_check.xlxs
ggplot(br_melt, aes(x = elev_bin, fill = dbh_class)) +  geom_bar(width=.5, position = "dodge")  +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Infected")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

#checks out in Excel br_chart_check.xlxs
ggplot(br, aes(x = elev_bin, fill = infrc)) +  geom_bar(width=.5, position = "dodge") + facet_wrap(~dbh_class)  +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Infected")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

br_dbh<-subset(br, dbh_class!="no DBH data")
unique(br_dbh$dbh_class)

#shift to a smaller bin width using continuous elevation data
ggplot(br_dbh, aes(x = elev, fill = infrc)) +  geom_histogram(binwidth=100, position = "dodge") + 
  facet_wrap(~dbh_class)  +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Infected")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

color_table <- tibble(
  infrc = c("Y", "N","NA"),
  Color = c("dark red", "dark green","blue"))

ggplot(br_dbh, aes(x=elev, fill=infrc)) +geom_histogram(aes(y=..count../sum(..count..)),binwidth=50, position="dodge")+ #aes(y = (..count..)/sum(..count..)),
  scale_y_continuous(labels=scales::percent) +#
  ylab("Percent")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Infected")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

ggplot(br_dbh, aes(x=elev, fill=dbh_class)) +geom_histogram(aes(y=..count../sum(..count..)),binwidth=50, position="dodge")+ #aes(y = (..count..)/sum(..count..)),
  scale_y_continuous(labels=scales::percent) +#
  ylab("Percent")+facet_wrap(~infrc)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Size Class")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))
#following chart is percentage infected/uninfected by size class and elevation.
#percentage is percentage across the entire sampled population of trees classed as infected or not
#a few trees ~25 were dropped because the did not have a dbh measurement and thus no dbh_class
#converting this to proportion infected within size class would require some more calculations
#not sure how to fix truncated density plot at high elevation.
ggplot(br_dbh, aes(x = elev, fill = infrc)) +  
  geom_histogram(binwidth= 50, position = "dodge",aes(y=..count../sum(..count..)*100)) + 
  #scale_y_continuous(labels=scales::percent)+ 
  xlab("Elevation (m)") + ylab("Percent")+
  scale_x_continuous(limits=c(2400,3200))+
  #stat_bin(aes(y=..density..), breaks = seq(min(br$elev), max(br$elev), by = 100), color="black") +
  geom_density(aes(y = ..count..), alpha=0.2)+ 
  #geom_line(stat="density", size = 1) +
  #facet_wrap(~dbh_class)  +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Infected")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 25))

#code for testing other x-axis variables as plots
ggplot(br_dbh, aes(x = lat, fill = infrc)) +  
  geom_histogram(binwidth= 0.1, position = "dodge",aes(y=..count../sum(..count..)*100)) + 
  #scale_y_continuous(labels=scales::percent)+ 
  xlab("Elevation (m)") + ylab("Percent")+
  scale_x_continuous()+
  #stat_bin(aes(y=..density..), breaks = seq(min(br$elev), max(br$elev), by = 100), color="black") +
  #geom_density(aes(y=..count../sum(..count..)*3000), alpha=0.2)+ 
  #geom_line(stat="density", size = 1) +
  #facet_wrap(~dbh_class)  +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Infected")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 25))

# library(scales)
# 
# ggplot(br_dbh, aes(x = elev, fill = infrc)) +  
#   geom_histogram(binwidth= 50, position = "dodge",aes(y = ..density..)) + 
#   scale_y_continuous(labels = percent_format())+ xlab("Elevation (m)") + ylab("Percent")+
#     scale_x_continuous(limits=c(2400,3200))+
#   #stat_bin(aes(y=..density..), breaks = seq(min(br$elev), max(br$elev), by = 100), color="black") +
#   geom_density(alpha = 0.2)+
#   #geom_line(stat="density", size = 1) +
#   facet_wrap(~dbh_class)  +
#    theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   scale_fill_discrete(name = "Infected")+
#   #theme(legend.position="none")+#removes legend+
#   theme(text = element_text(size = 25))



br_inf<-subset(br_dbh, infrc=="Y");head(br_inf)

# ggplot(br_inf, aes(x=elev, fill=dbh_class)) + geom_density(alpha =0.2)+ #
#   scale_y_continuous() +#
#   ylab("Percent")+
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   scale_fill_discrete(name = "Infected")+
#   #theme(legend.position="none")+#removes legend+
#   theme(text = element_text(size = 30))


# ggplot(br, aes(x=elev)) +geom_histogram(aes(y=..count../sum(..count..)),binwidth=50)+ #aes(y = (..count..)/sum(..count..)),
#   scale_y_continuous(labels=scales::percent) + facet_wrap(~infrc)+#
#   ylab("Percent")+
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   scale_fill_discrete(name = "Infected")+
#   #theme(legend.position="none")+#removes legend+
#   theme(text = element_text(size = 30))
# 
# 
# ggplot(br, aes(x=elev, fill=infrc)) +geom_density(aes(y = (..count..)/sum(..count..)),alpha=0.2)+ 
#   scale_y_continuous(labels=scales::percent) +
#   ylab("Relative frequencies")+
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   scale_fill_discrete(name = "Infected")+
#   #theme(legend.position="none")+#removes legend+
#   theme(text = element_text(size = 30))
# 
# 
# br_inf<-subset(br, infrc=="Y");head(br_inf)#subset only infected class
# ggplot(br_inf, aes(x=elev, fill=dbh_class)) +geom_histogram(binwidth=50)+ #aes(y = (..count..)/sum(..count..)),
#   scale_y_continuous() +#labels=scales::percent
#   ylab("Count")+
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   scale_fill_discrete(name = "Infected")+
#   #theme(legend.position="none")+#removes legend+
#   theme(text = element_text(size = 30))
# 
# ggplot(br_inf, aes(x=elev, fill=dbh_class)) +geom_density(aes(y = (..count..)/sum(..count..)),alpha=0.4)+ #
#   scale_y_continuous(labels=scales::percent) +scale_fill_manual( values = c("dark red","blue","orange","purple","blue"))+
#   ylab("Relative frequencies")+
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   scale_fill_discrete(name = "Infected")+
#   #theme(legend.position="none")+#removes legend+
#   theme(text = element_text(size = 30))

ggplot(br_inf, aes(x=lat, fill=dbh_class)) +geom_density(aes(y = (..count..)/sum(..count..)),alpha=0.4)+ 
  scale_y_continuous(labels=scales::percent) +
  ylab("Relative frequencies")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Infected")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

# ggplot(br_inf, aes(x=lat, fill=infrc)) +geom_histogram(aes(y = (..count..)/sum(..count..)),binwidth=1)+ 
#   scale_y_continuous(labels=scales::percent) +
#   ylab("Relative frequencies")+
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   scale_fill_discrete(name = "Infected")+
#   #theme(legend.position="none")+#removes legend+
#   theme(text = element_text(size = 30))


#chart proportion across all elevation bands
# ggplot(br, aes(infrc, fill=infrc)) +geom_bar(aes(y = (..count..)/sum(..count..)))+ 
#   scale_y_continuous(labels=scales::percent) +
#   ylab("relative frequencies")+ facet_grid(~elev_bin)
  
#chart proportion of infected within each elevation band
ggplot(transform(br), aes(x= infrc,  group= elev_bin)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat = "count",position="dodge") +
  scale_fill_manual(values = color_table$Color) +
  labs(y = "Percent", fill="infected")+
  facet_grid(~elev_bin) + xlab("Infection Status")+
  scale_y_continuous(labels = scales::percent)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position="none")+#removes legend+
theme(text = element_text(size = 30))


# ggplot(transform(br), aes(x= lat,  fill= infrc)) + 
#   geom_density(alpha = 0.1) +
#   scale_fill_manual(values = color_table$Color) +
#   labs(y = "Percent", fill="infected")+
#   #facet_grid(~elev_bin) + xlab("Infection Status")+
#   scale_y_continuous(labels = scales::percent)+
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   #theme(legend.position="none")+#removes legend+
#   theme(text = element_text(size = 30))


# #panel by year doesn't make sense if year = last time visited
# ggplot(transform(br, panel = c("2004-2007", "2008-2011", "2012-2015","2016-2018")[as.numeric(panel)]),
#        aes(x= infrc,  group= elbin)) + 
#   geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat = "count",position="dodge") +
#   labs(y = "Percent", fill="infected") +
#   facet_grid(~panel) + xlab("Infection Status")+
#   scale_y_continuous(labels = scales::percent)
# 
# ggplot(transform(br, panel = c("2004-2007", "2008-2011", "2012-2015","2016-2018")[as.numeric(panel)]),
#        aes(x= elbin,  group= infrc)) + 
#   geom_histogram(position="dodge")+ facet_grid(~panel) 
# 
# ggplot(transform(br, panel = c("2004-2007", "2008-2011", "2012-2015","2016-2018")[as.numeric(panel)]),
#        aes(x= elbin,  group= infrc)) + 
#   geom_density(position="dodge")+ facet_grid(~panel) 
# 
# 
# ggplot(br,aes(x= elbin, group= panel)) + 
#   geom_density(position="dodge")+ facet_grid(~panel) 
# 
# ggplot(transform(br, panel = c("2004-2007", "2008-2011", "2012-2015","2016-2018")[as.numeric(panel)]),
#        aes(x= elbin, group= panel)) + labs(y = "Density Live Infected", fill="infrc")+
#   geom_density(position="dodge")+ facet_grid(~panel) 
# 
# ggplot(transform(subyn, panel = c("2004-2007", "2008-2011", "2012-2015","2016-2018")[as.numeric(panel)]),
#        aes(x= elbin, group= panel)) + labs(y = "Density Infected", fill="infrc")+
#   geom_density(position="dodge")+ facet_grid(~panel)
# 
# ggplot(transform(suby, panel = c("2004-2007", "2008-2011", "2012-2015","2016-2018")[as.numeric(panel)]),
#        aes(x= elbin, group= panel)) + 
#   geom_bar(aes(y = ..prop..), stat = "count",position="dodge") +
#   labs(y = "Percent Live Infected", fill="infrc")+ facet_grid(~panel) 
# 
# ggplot(transform(subyn, panel = c("2004-2007", "2008-2011", "2012-2015","2016-2018")[as.numeric(panel)]),
#        aes(x= elbin, group= panel)) + 
#   geom_bar(aes(y = ..prop..), stat = "count",position="dodge") +
#   labs(y = "Percent Infected", fill="infrc")+ facet_grid(~panel)

# sumtab<-NULL
# for(yr in 2004:2018){#unique(transition$year)
#   subn<-subset(transition, year==yr & trans=="N" & status=="L");head(subn)
#   countn<-nrow(subn)
#   suby<-subset(transition, year==yr &trans=="Y"& status=="L");head(suby)
#   county<-nrow(suby)
#   subp<-subset(transition, year==yr &trans=="prune"& status=="L");head(subp)
#   countp<-nrow(subp)
#   subf<-subset(transition, year==yr &trans=="first"& status=="L");head(subf)
#   countf<-nrow(subf)
#   total<-countn+county+countp+countf
#   countn_pct<-countn/total
#   county_pct<-county/total
#   countp_pct<-countp/total
#   countf_pct<-countf/total
#   rowsum<-cbind(yr,countn,countn_pct,county,county_pct,countp,countp_pct,countf,countf_pct);rowsum
#   sumtab<-rbind(sumtab,rowsum)
# }
# head(sumtab)
# str(sumtab)
# sumtab<-as.data.frame(sumtab)#have to convert numeric matrix to df then convert yr to factor for ggplotting
# str(sumtab)
# sumtab[,1]<-as.factor(sumtab[,1]);str(sumtab)
# #mquant<-melt(nonforest_ecdfquants[,1:6],id=c("maj_class","dec"));head(mquant)
# msumtab<-melt(sumtab,id=c("yr","countn","county","countp","countf"));head(msumtab)
# ggplot(msumtab,aes(x=yr, y=value))+geom_point()+ facet_wrap(~variable)

#plot spatial patterns
#where did infenction occur in each year?
# head(transition)
# suby<-subset(transition, status=="L" & infected=="Y" );head(suby)
# mtrans<-melt(suby, id=c("site", "treeid","year","infected", "elev","aspect", "height_class", "recentdbh",  "dbh_class", "status","long","lat"));head(mtrans)
# ggplot(mtrans,aes(x=long, y=lat)) + geom_point(aes(colour = factor(year), size=aspect))
# ggplot(mtrans,aes(x=long, y=lat)) + geom_point() + facet_wrap(~year)
# 
# subyn<-subset(transition, status=="L" & infected=="Y" | infected=="N");head(subyn)
# #create elevation bins by 200m intervals
# subyn$elbin<-.bincode(subyn$elev, breaks=c(2400, 2600,2800,3000,3200,3400), TRUE);head(subyn)
# suby$elbin<-.bincode(suby$elev, breaks=c(2400, 2600,2800,3000,3200,3400), TRUE);head(suby)
# 
# ggplot(transform(subyn, panel = c("2004-2007", "2008-2011", "2012-2015","2016-2018")[as.numeric(panel)]),
#        aes(x= infected,  group= elbin)) + 
#   geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat = "count",position="dodge") +
#   labs(y = "Percent", fill="infected") +
#   facet_grid(~panel) +
#   scale_y_continuous(labels = scales::percent)
# 
# ggplot(transform(suby, panel = c("2004-2007", "2008-2011", "2012-2015","2016-2018")[as.numeric(panel)]),
#        aes(x= elbin,  group= infected)) + 
#   geom_histogram(position="dodge")+ facet_grid(~panel) 
#     
# aes(y = ..prop.., fill = factor(..x..)), stat = "count",position="dodge") +
#   labs(y = "Percent", fill="infected") +
#  
# 
#   scale_y_continuous(labels = scales::percent)


#######################################################################################
#################add the water balance data to the transition data file################
#this is where we add water balance information to the transition data frame###########
#first identify the site name and year of interest, then open the site's water balance
#file and extract the desired water balance data 
#attach that data for each site/year row combination to the transitions data frame
#new columns can be added as needed for different water balance variables and periods
#of interest, including lags###########################################################
setwd("C:\\David\\Water balance\\blister rust\\daily")#find the water balance data
head(br)
n=10596.1
sumstack<-NULL #empty variable for holding summarized water balance data
#count of unique sites
#nsites<-unique(transition$site);nsites;length(nsites)
nsites<-unique(br$site);nsites;length(nsites)
#nsites<- nsites[nsites!= 8386.2];length(nsites)#this site missing b/x long missing in site file
for (n in nsites){#this loop takes a long time - 30 min to an hour
  site<-paste(n,".csv",sep="");site
  # tree<-transition[n,2];tree#build a date value for year of observation linked to last day of observation year Dec 31
  # dt<-transition[n,3];dt#year of interest for a site.  use as an index in the water balance file to start compiling wb data
  # dtx<-as.Date(paste(dt,"-12-31",sep=""));dtx#format date as POSIX so you and add and subtract number of days that bracket time period of interest
  # #define periods of interest for extracting wb 
  start_date <- '1980-01-01';start_date#beginning of daymet record
  #end_date<-paste(dt,"-31-12",sep="");end_date#year of observation.  collect water balance data for this year and earlier
  #pull the water balance data and do some formatting
  sitewb<-read.csv(site, stringsAsFactors=FALSE);head(sitewb)
  nrow(sitewb); names(sitewb)
  colnames(sitewb)[1:26]<-c("date","year","month","day","SVP","RH","VPD","DRO","P","TEMP","F","RAIN","SNOW","PACK","MELT","W","PET","W_PET","SOIL","SOIL_Delta",
                            "AET","RUNOFF","D","GDD","site","Dew")#rename the last 3 columns
  names(sitewb)
  corewb<-sitewb[,c(-1,-25)];head(corewb)#drop dates and strings prior to converting to xts
  #keep any strings out of the xts matrix or it will convert everything to strings with obnoxious quotes
  sitecol<-sitewb[,25]#in case you need to add the site information back into a data frame
  f2<-as.Date(sitewb[,1], format = "%m/%d/%Y");head(f2);length(f2)
  wbxts<-xts(coredata(corewb),order.by=f2);head(wbxts)#build xts data structure 
  doy<-yday(wbxts);head(doy)
  mth<-month(wbxts);head(mth)
  yr<-year(wbxts);head(yr)
  wbxts<-cbind(wbxts,yr,mth,doy);head(wbxts);ncol(wbxts)
  names(wbxts)[25:27]<-c("yr","mth","doy");head(wbxts)#rename the doy column in the xts
  
  ########################################
  #begin summarizing water balance data variables by periods of interest
  #the strategy is to roll daily up to monthly, then monthly up to annual
  #if interest is in a chunk of a year, still have to roll up to a single annual value
  #representing the summary of each chunk as an annual value for that chunk
  #after summarizing across the entire daymet record then pull out years of interest based on when trees were measured
  #for instance if a tree was measured in 2010 pull annual summaries for 2010 and earlier
  ############Annual summary###############
  historic<-wbxts[paste("2000-01-01","2019-01-01", sep="/"), c("yr","mth","doy","SVP","RH","VPD","P","TEMP","RAIN","PACK","W","PET","W_PET","SOIL",
                                                               "SOIL_Delta","AET","RUNOFF","D","GDD","Dew")]
  head(historic);tail(historic)
  #annual summary from daily values
  a_means<-apply.yearly(historic[,c("SVP","RH","VPD","TEMP","SOIL")], colMeans);a_means
  a_sums<-apply.yearly(historic[,c("P","RAIN","W","PET","W_PET","AET","RUNOFF","D","GDD","Dew")], colSums);a_sums
  a_max<- apply.yearly(historic[,c("PACK")], max);a_max 
  a_summary<-cbind(a_means,a_sums,a_max);head(a_summary)
  colnames(a_summary)[]<-c("aSVP","aRH","aVPD","aTEMP","aSOIL","aP","aRAIN","aW","aPET","aW_PET","aAET","aRUNOFF","aD","aGDD","aDew","aPACK");head(a_summary)
  #monthly summary from daily values
  m_means<-apply.monthly(historic[,c("SVP","RH","VPD","TEMP","SOIL")], colMeans);m_means
  m_sums<-apply.monthly(historic[,c("P","RAIN","W","PET","W_PET","AET","RUNOFF","D","GDD","Dew")], colSums);m_sums
  m_max<- apply.monthly(historic[,c("PACK")], max);m_max 
  m_summary<-cbind(m_means,m_sums,m_max);head(m_summary)
  colnames(m_summary)[]<-c("mSVP","mRH","mVPD","mTEMP","mSOIL","mP","mRAIN","mW","mPET","mW_PET","mAET","mRUNOFF","mD","mGDD","mDew","mPACK");head(m_summary)  
  
  #specified periods according to hypotheses of infection, in this case humidity in august and september
  #eg Aug, Sept relative humidity. this can be modified to select periods between days of year in doy column
  fall_mean<-subset(historic, mth==8 | mth==9, select=c(SVP,RH,VPD,TEMP,SOIL));head(fall_mean)
  fall_sum<-subset(historic, mth==8 | mth==9, select=c(RAIN,AET,D,GDD,Dew,PET));head(fall_sum)
  #summarize monthly chunks to yearly by setting end points for summary function to "years"
  ep<-endpoints(fall_mean, on = "years", k = 1)#summarize 3 month blocks by year
  # Now calculate the weekly mean and display the results
  fall_mean1<-period.apply(fall_mean[], INDEX = ep, FUN = colMeans);head(fall_mean1)
  fall_sum1<-period.apply(fall_sum[], INDEX = ep, FUN = colSums);head(fall_sum1)
  as_summary<-cbind(fall_mean1,fall_sum1);head(as_summary)#'ao' stands for august - october period
  colnames(as_summary)[]<-c("asSVP","asRH","asVPD","asTEMP","asSOIL","asRAIN","asAET","asD","asGDD","asDew","asPET")
  
  #specified periods according to hypotheses of infection, in this case humidity in May
  #eg May relative humidity. this can be modified to select periods between days of year in doy column
  spring_mean<-subset(historic, mth==5 | mth==5, select=c(SVP,RH,VPD,TEMP,SOIL));head(spring_mean)
  spring_sum<-subset(historic, mth==5 | mth==5, select=c(RAIN,AET,D,GDD,Dew,PET));head(spring_sum)
  #summarize monthly chunks to yearly by setting end points for summary function to "years"
  ep<-endpoints(spring_mean, on = "years", k = 1)#summarize 3 month blocks by year
  # Now calculate the weekly mean and display the results
  spring_mean1<-period.apply(spring_mean[], INDEX = ep, FUN = colMeans);head(spring_mean1)
  spring_sum1<-period.apply(spring_sum[], INDEX = ep, FUN = colSums);head(spring_sum1)
  s_summary<-cbind(spring_mean1,spring_sum1);head(s_summary)#'s' stands for spring period
  colnames(s_summary)[]<-c("mSVP","mRH","mVPD","mTEMP","mSOIL","mRAIN","mAET","mD","mGDD","mDew","mPET")#m for may
  
  
  #daily rolling sum of 2 day relative humidity
  period_sums<-rollapply(fall_mean[,c("RH")], 2, sum, by = 1, by.column = TRUE);head(period_sums);tail(period_sums)#fall_mean holds RH column
  ep<-endpoints(period_sums, on = "years", k = 1);ep#summarize 3 month blocks by year
  # a function to count 2 day rolling sums greater than 0.97*2 on an annual basis:
  f <- function(x) {
    result<-sum(x > 0.7*2, na.rm=TRUE)#na.rm needed because first result in 2 day rolling sum is NA
  }
  fall_rhcount<-period.apply(period_sums, INDEX=ep,FUN=f);fall_rhcount
  #use for checking result sum(period_sums[,1] >(0.97*2), na.rm=T)#0.97 RH two days in a row meet conditoins for basidiospore transport to pine
  colnames(fall_rhcount)[]<-"dRHc"
  as_summary<-cbind(fall_rhcount,as_summary);head(as_summary)
  
  #make the summary xts files into a data frame that strips the time stamps.  this is necessary because annual summaries end Dec 31, and monthly summaries end Octber 31
  #joining them results in two separate rows for each year of annual or monthly chunked data
  yr<-year(a_summary);yr#make a column of year for each row
  #now create a massive vector of all the different periods of summarized data
  #then rbind them as the loop progresses.
  #each pass through this loop creates values you want to append to the historic data file
  #NOTE: can only append data summarized on annual time step to field visit file.  fall_summary is an annual data set of only fall time period
  #first strip off time by converting to data frame, then identify the year that goes with field visit year extracted earlier as dt 
  vec<-cbind.data.frame(site,yr,coredata(a_summary),coredata(s_summary), coredata(as_summary));head(vec)
  sumstack<-rbind(sumstack,vec)
}#compile water balance summary by period for next row 

setwd("C:\\David\\Water balance\\blister rust\\")#find the water balance data
save.image("C:/David/Water balance/blister rust/.RData")
#some formatting steps
sumstack_backup<-sumstack# incase you need to retrieve it without having to run through the entire loop above
head(sumstack)
site2<-as.data.frame(unlist(strsplit(as.character(sumstack$site), ".csv")))
colnames(site2)[]<-"site2";head(site2)
sumstack<-cbind(site2,sumstack)
sumstack<-sumstack[,-2]#drop the site column with .csv extension
colnames(sumstack)[1:2]<-c("site","year")


head(sumstack);nrow(sumstack);length(unique(sumstack$year))*length(unique(sumstack$site))#these should be equal to nrow sumstack
head(br);nrow(br)
#all a column for each of 4 lag years prior to site visits. These columns are used to match corresponding water balance data for a given lag
# transition$yearL1<-transition$year-1;head(transition)
# transition$yearL2<-transition$year-2;head(transition)
# transition$yearL3<-transition$year-3;head(transition)
# transition$yearL4<-transition$year-4;head(transition)

#pull out the most recent observation for each tree
#this is the data set of infected / not-infected trees we evaluate as a response to climate 
#transub<-subset(transub1, infrc=="Y" | infrc=="N");head(transub);unique(transub$infrc)

i="10596.1.55"
holder<-NULL
for(i in unique(br$treeid)){
  test<-subset(br,treeid==i)
  test2<-max(test$year)
  test3<-subset(test, year==test2)
  holder<-rbind(holder,test3) 
  
}
head(holder);nrow(holder)
#check for duplicates
n_occur <- data.frame(table(holder$treeid))
holder[holder$treeid %in% n_occur$Var1[n_occur$Freq > 1],]#result should be empty or zero indicating no duplicates

#with the subset of trees identified, match to them the long-term climate data for their location summarized over
#the period of time from year 2000 through the time when they were last observed
head(sumstack)#annualized climate data for all transects in one data frame
head(holder)#the trees of interest

#NOTE: holder and br are equal based on the way year in br was recorded from the last observation
#join tree data with long-term average climate data starting in 2000.
#Variables have previously been summed or averaged for each year according to the appropriate method for each
# here annual average climate for each stand location is calculated across years by averaging all values except for counts of dew and days RH> threshold. 
#as is average or sum of august to september climate data then averaged across years
unique(br$infrc);count(br$infrc=="Y")#False is not infected, true is infected tree count
#subset(br, infrc=="N")#make sure there are Y's and N's in the resulting dataframe
head(sumstack)#has all years of climate data for each stand starting at 2000
i = "958.1.64"  
treeclim<-NULL
for(i in unique(br$treeid)){#unique(nsites)
  test<-subset(br,treeid==i);test
  yr<-test$year
  st<-test$site
  substack<-subset(sumstack,site==st );head(substack);names(substack)#if want to limit to yrs before last sampling & year<yr+1
  #next line results in a vector so must convert to data frame and transpose to get cbind to work
  subavg<-t(as.data.frame(apply(substack[,c(3:16,18:28,31:39,41)],2,mean)));subavg# average SVP, RH, VPD, T, Soil
  subsum<-t(as.data.frame(apply(substack[,c(17,30,40)],2,sum)));subsum# sum DEW, dRHc
  #subsum<-sum(substack[,c(17,19,29)]);subsum#count of days relative humidity > threshold and days with dew
  #cbind(x, t = t[names(x)])
  subs<-cbind(test,subavg, row.names = NULL,subsum)
  treeclim<-rbind(treeclim,subs)
}
names(treeclim)
#colnames(treeclim)[46] <-"dRHc"
head(treeclim);tail(treeclim);max(treeclim$dRHc)
unique(treeclim$infrc)
unique(treeclim$dbh_class)
# #match these data by site and year to the transitions data frame  
# transwb<-merge(transition, sumstack, by.x=c("site","year"), by.y=c("site","year"))# all.x=TRUE
# head(transwb) # Okay if not all Y's match because some trees not visited all yrs
# nrow(transwb)#after merge should equal nrow(transition) minus 2018 trees that don't have assoicated daymet data b/c it stops at Dec 2017
# nrow(transition)
# 
# #now create the lagged data sets
# transwbL1<-merge(transition, sumstack, by.x=c("site","yearL1"), by.y=c("site","year"));head(transwbL1)# all.x=TRUE
# transwbL2<-merge(transition, sumstack, by.x=c("site","yearL2"), by.y=c("site","year"));head(transwbL2)# all.x=TRUE
# transwbL3<-merge(transition, sumstack, by.x=c("site","yearL3"), by.y=c("site","year"));head(transwbL1)# all.x=TRUE
# transwbL4<-merge(transition, sumstack, by.x=c("site","yearL4"), by.y=c("site","year"));head(transwbL1)# all.x=TRUE

save.image("C:/David/Water balance/blister rust/.RData")

# #subsample just the transition yes or transition n rows of data
# transwb_tmp<-subset(transwb, trans=="N" | trans=="Y")
# transwb_yn<-subset(transwb_tmp, status=="L");unique(transwb_yn$trans);unique(transwb_yn$status)
# transwb_tmp<-subset(transwbL1, trans=="N" | trans=="Y")
# transwb_ynL1<-subset(transwb_tmp, status=="L")
# transwb_tmp<-subset(transwbL3, trans=="N" | trans=="Y")
# transwb_ynL3<-subset(transwb_tmp, status=="L")
# transwb_tmp<-subset(transwbL4, trans=="N" | trans=="Y")
# transwb_ynL4<-subset(transwb_tmp, status=="L")
# 
# head(transwb_yn);nrow(transwb_yn)
# treeid<-(paste(transwb_yn[,1],".",transwb_yn[,3],sep=""));length(treeid);str(treeid)
# transwb_yn<-cbind(as.data.frame(treeid),transwb_yn);head(transwb_yn)



#########################################################################################
#########################exploratory plotting#############################################
##########################################################################################
#just do this once because variables correlated in one year will be correlated in all others
#determine if variables are correlated.  Drop one from each pair if correlation > 0.75
#write.csv(cor(transwb_yn[,c(15:39)]),"varcorel.csv")
names(treeclim)
write.csv(cor(treeclim[,c(20:57)], use = "complete.obs"),"varcorel.csv")#complete obs necessary if any of the site wb missing b/c if failed to run (happens if bad lat long submitted to Daymet server)
varcor<-cor(treeclim[,c(20:57)], use = "complete.obs")
corrplot(varcor,method="circle", type="lower")
corrplot(varcor,method="number", type="lower")


#corrplot requires a correlation matrix as input
# varcor<-cor(transwb_yn[,c(15:39)])
# corrplot(varcor,method="circle", type="lower")
# corrplot(varcor,method="number", type="lower")
#based on r^2>0.75 these variables can be dropped
#Drop aSVP","aRH","aVPD","aGDD","aD","aSOIL","aW_PET","aRUNOFF","aoSVP","aoGDD","aoSOIL", and a few categorical variables for plotting box plots
names(treeclim)
#set up a data frame for exploratory plotting that has only numbers in vlaues column so it works with ggplot
# drops <- c("aSVP","aRH","aVPD","aT","aSOIL","aP","aRAIN","aW","aPET","aW_PET","aRUNOFF","aD","aGDD","aoGDD","aoSOIL","aoAET","aoSVP",
#            "panel","infected","trans","height_class","recentdbh", "dbh_class","status")
# #subtrans<-transwb_yn[ , !(names(transwb_yn) %in% drops)]
# subtreeclim<-treeclim[ , !(names(treeclim) %in% drops)]
# 
# treemelt<-melt(subtreeclim, id=c("site","treeid","infrc","year","elev","long","lat"));head(treemelt)
# ggplot(data=treemelt, aes(x=elev, y=value,group = round_any(elev, 100, floor)))+
#   geom_boxplot()+facet_wrap(~variable, scales="free")
# ggplot(data=treemelt, aes(x=long, y=value,group = round_any(long, 0.5, floor)))+
#   geom_boxplot()+facet_wrap(~variable, scales="free")
# ggplot(data=treemelt, aes(x=lat, y=value,group = round_any(lat, 0.5, floor)))+
#   geom_boxplot()+facet_wrap(~variable, scales="free")
#########################################################################################
#########################################################################################
#automate removal of correlated variables by setting diagonal and upper triangle to zero
tmp<-cor(treeclim[,c(20:57)], use = "complete.obs")
write.csv(tmp,"tmp.csv")
tmp[!lower.tri(tmp)] <- 0
tmp.new <- tmp[,!apply(tmp,2,function(x) any(abs(x) > 0.7))];tmp.new
#############################glmer modeling##############################################
#########################################################################################
#drop the more esoteric variable in pairs with r^2>0.49, r =0.7 , see varcorel.xlsx for colored matrix
drops <- c("aSVP","aRH","aVPD","aTEMP","aSOIL","aP","aRAIN","aW","aPET","aW_PET","aAET","aRUNOFF",
           "aD","aGDD","aDew","aoGDD","aoSVP","aoD","aoSOIL","aoAET","aoVPD","aoDew","aoPET","dRHc")

##################################################################
#########################Lag0#####################################
#swap out the lagged data sets in the next line;transwb_ynL1,transwb_ynL2,transwb_ynL3,transwb_ynL4
#subtrans<-transwb_yn[ , !(names(transwb_yn) %in% drops)]
treeclim2<-treeclim[ , !(names(treeclim) %in% drops)]
unique(treeclim$infrc)
unique(treeclim2$dbh_class)
a<-table(treeclim2$dbh_class)
a#count of trees in each dbh class
head(treeclim2)
names(treeclim2)
#summarize using data.table #https://s3.amazonaws.com/assets.datacamp.com/img/blog/data+table+cheat+sheet.pdf
dt<-data.table(treeclim2);dt
dt[ ,count(infrc), by=.(elev_bin,dbh_class)]


#treeclim4<-subset(treeclim2, dbh_class!="no DBH data")
#exploratory plotting

color_table <- tibble(
  infrc = c("Y", "N"),
  Color = c("dark red", "dark green"))

#size class distribution counts
ggplot(treeclim2, aes(x= elev,  fill= dbh_class)) + 
  geom_histogram(binwidth= 25,position="dodge") +facet_grid(~dbh_class)+
  #scale_fill_manual(values = color_table$Color) + 
  labs(y = "Count")+theme(legend.position="right")+
  xlab("Infection Status")+
  #scale_y_continuous(labels = scales::percent)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

#size class distribution percent
ggplot(treeclim2, aes(x= elev,  fill= dbh_class)) + 
  geom_histogram(aes(y=..count../sum(..count..)),binwidth=25, position="dodge")+facet_grid(~dbh_class)+
  #scale_fill_manual(values = color_table$Color) + #
  labs(y = "Percent")+theme(legend.position="right")+# facet_grid(~dbh_class)+
  #xlab("Infection Status")+
  #scale_y_continuous(labels = scales::percent)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

#chart proportion of infected within each elevation band
ggplot(treeclim2, aes(x= elev,  fill= infrc)) + 
  geom_histogram(binwidth= 25,position="dodge") +
  scale_fill_manual(values = color_table$Color) + facet_grid(~dbh_class)+
  labs(y = "Count")+theme(legend.position="right")+
   xlab("Infection Status")+
  #scale_y_continuous(labels = scales::percent)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

ggplot(treeclim2, aes(x= elev,  fill= infrc)) + 
  geom_histogram(aes(y=..count../sum(..count..)),binwidth=25, position="dodge") +
  scale_fill_manual(values = color_table$Color) + facet_grid(~dbh_class)+
  labs(y = "Percent")+theme(legend.position="right")+
  xlab("Infection Status")+
  #scale_y_continuous(labels = scales::percent)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

# 
ggplot(treeclim2, aes(x=elev, fill=infrc)) + geom_histogram(binwidth= 25,position="dodge") +  #facet_grid(~dbh_class)+
  geom_density(aes(y=..count../sum(..count..)*100000),alpha=0.1)+ 
  scale_color_brewer(palette = "Set3")

#  library(scales)
# 
#  ggplot(treeclim2, aes(x = elev, fill = infrc)) +  
#    geom_histogram(binwidth= 25, position = "dodge",aes(y = ..density..)) + 
#    scale_y_continuous(labels = percent_format())+ xlab("Elevation (m)") + ylab("Percent")+
#      scale_x_continuous(limits=c(2400,3200))+
#    theme_bw() + theme(panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#    #theme(legend.position="none")+#removes legend+
#    theme(text = element_text(size = 30))

ggplot(treeclim2, aes(x= elev,  fill= infrc)) + 
  geom_histogram(aes(y=..count../sum(..count..)),binwidth=25, position="dodge")+
  scale_fill_manual(values = color_table$Color) + facet_grid(~elev_bin)+
  labs(y = "Percent")+theme(legend.position="right")+# facet_grid(~dbh_class)+
  #xlab("Infection Status")+
  #scale_y_continuous(labels = scales::percent)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

ggplot(treeclim2, aes(x= elev,  fill= infrc)) + 
  geom_density(aes(y=..count../sum(..count..)*10), alpha= 0.5)+
  scale_fill_manual(values = color_table$Color) + #facet_grid(~elev_bin)+
  labs(y = "Percent")+theme(legend.position="right")+# facet_grid(~dbh_class)+
  #xlab("Infection Status")+
  scale_y_continuous(labels = scales::percent)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))






# Correlation panel

# Create the plots
#pairs(treeclim2[,c(20:41)], lower.panel = NULL)#slow with lots of variables
cor(treeclim2[,c(20:41)],use = "complete.obs")
#corrplot(treeclim2[,c(20:24)],use = "complete.obs",method="number")

names(treeclim2)
treeclim2_melt<-melt(treeclim2, id=c("panel","site","treeid","year","infected","trans","height_class","recentdbh","dbh_class","status",
                                     "long","lat","slope","aspect","infrc","elbin","infrc2","elev_bin","elev" ))

treeclim2_melt<-melt(treeclim2, id=c("panel","site","treeid","year","infected","trans","height_class","recentdbh","dbh_class","status",
                                     "long","lat","slope","aspect","infrc","elbin","infrc2","elev_bin","elev" ))
head(treeclim2_melt)

treeclim22_melt<-subset(treeclim2_melt, variable=="aPACK"|variable=="asRH"|variable=="asTEMP"|variable=="asRAIN")
head(treeclim22_melt)

library(ggpubr)# may have to adjust equation text in adobe
ggplot(treeclim22_melt, aes(x=elev, y= value))+ geom_point()+geom_smooth(method=lm, se=FALSE) + facet_wrap(~variable, scales="free")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #theme(legend.position="none")+#removes legend+ 
  theme(text = element_text(size = 30)) + xlab("Elevation (m)") +ylab("")+
  stat_cor(label.y = 7) +
  stat_regline_equation(label.y = 4)


treeclim2_infected<-subset(treeclim2, infrc=="Y");unique(treeclim2_infected$infrc)#subset of infected trees only
head(treeclim2_infected)

#order levels used in legend so they are increasing
unique(treeclim2$dbh_class)
treeclim2$dbh_class <- factor(treeclim3$dbh_class, levels = c("<=2.5cm", ">2.5<=10cm", ">10<=30cm",">30cm"))
names(treeclim2)
#exploratory plotting swap out x= with weather variables summarized in treeclim2
treeclim2<-subset(treeclim2,dbh_class!="no DBH data")
ggplot(treeclim2, aes(x=elev, fill=dbh_class)) +geom_density(aes(y = (..count..)/sum(..count..)),alpha=0.6)+ 
  scale_y_continuous(labels=scales::percent) + facet_grid(~infrc)+
  scale_x_discrete(limits=c(2600,2800,3000))+#comment out this line if plotting lat, long
  ylab("Relative frequencies")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual( name = "Size Class",values = c('#d7191c','#fdae61','#2c7bb6','#abd9e9'))+
  xlab("Elevation (m)")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))+  theme(panel.spacing.x=unit(2, "lines"),panel.spacing.y=unit(1, "lines"))


#rescale RH from 0 to 1 to 0 to 100:  this avoids warning about some variables being on very different scales
treeclim2$asRH<-(treeclim2$asRH)*100;head(treeclim2$aoRH)
#scale all numeric variables before statistical analysis
names(treeclim2)

#replace 1,2 in infrc2 with 0,1 so logistic models work 
treeclim2$infrc2[treeclim2$infrc=="Y"]<-1
treeclim2$infrc2[treeclim2$infrc=="N"]<-0

treeclim3<-subset(treeclim2,dbh_class!="no DBH data");head(treeclim3)
write.csv(treeclim3,"treeclim3.csv")

#infection stats
dt<-data.table(treeclim3);dt
dt[ ,count(infrc), by=.(elev_bin,dbh_class)]
dt[ ,count(infrc)]

#exploratory plotting swap out x= with weather variables summarized in treeclim2

#order levels used in legend so they are increasing
unique(treeclim3$dbh_class)
treeclim3$dbh_class <- factor(treeclim3$dbh_class, levels = c("<=2.5cm", ">2.5<=10cm", ">10<=30cm",">30cm"))

ggplot(treeclim3, aes(x=elev, fill=dbh_class)) +geom_density(aes(y = (..count..)/sum(..count..)),alpha=0.6)+ 
  scale_y_continuous(labels=scales::percent) + scale_fill_manual(values = c("#E495A5","#ABB065","#39BEB1","#ACA4E2"))+
  facet_grid(~infrc)+
  ylab("Relative frequencies")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Size Class")+ 
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))



#scale numeric variables by subtracting mean values
scaled<- scale(treeclim3[c(11:15,20:41)], center = TRUE, scale = TRUE);head(scaled)
treeclim4<-cbind(treeclim3[,c(1:10,16:19)],scaled);head(treeclim4)

#can add models of spring conditions here too, but Ribes is ubiquitous so not likely relations with spring matter even if they are good.

m1<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m2<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+elev,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m3<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+elev+long,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m4<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+elev+long+lat,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m5<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+elev+long+lat+slope,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m6<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+elev+long+lat+slope+aspect,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m8<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+aPACK+asRH+asTEMP+asRAIN,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m9<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+aPACK+asRH+asTEMP+asRAIN,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m10<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+aPACK+asRH+asTEMP,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m11<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+aPACK+asRH,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m12<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+aPACK,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m13<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+long,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m14<-glmer(data=treeclim4, infrc2~ (1|site)+elev,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m15<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+asRH+asTEMP,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m16<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+elev+long+asTEMP+asRH,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m17<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+asRH+asTEMP+asRH*asTEMP,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m19<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+asRH,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m20<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+asTEMP,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m21<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+asRAIN,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m22<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+aPACK,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))

mtest<-glmer(data=treeclim4, infrc2~ (1|site)+recentdbh+asRAIN+elev,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))

setwd("C:\\David\\Water balance\\blister rust\\")
model_rank<-model.sel(m1,m2,m3,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m17,m19,m20,m21,mtest)#keep 18 and 19 out b/c spatial pattern in climate double counted if have both position and climate in same model
model_rank
write.csv(model_rank,"model_rank.csv")
model_anova<-anova(m1,m2,m3,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m17);model_anova
write.csv(model_anova,"model_anova.csv")

# Nakagawa & Schielzeth's (2013) conditional and marginal r^2 values
r.squaredGLMM(m17)
r.squaredGLMM(m15)
r.squaredGLMM(m10)
r.squaredGLMM(m8)
r.squaredGLMM(m9)
r.squaredGLMM(m3)
r.squaredGLMM(m4)
r.squaredGLMM(m6)
r.squaredGLMM(m2)
r.squaredGLMM(m13)
r.squaredGLMM(m19)
r.squaredGLMM(m11)
r.squaredGLMM(m21)
r.squaredGLMM(m20)
r.squaredGLMM(m1)
r.squaredGLMM(m12)
r.squaredGLMM(m14)



summary(treeclim4)#use to 4generate simple summary stats including counts of infected/uninfected

#summarize using data.table #https://s3.amazonaws.com/assets.datacamp.com/img/blog/data+table+cheat+sheet.pdf
dt<-data.table(treeclim2);dt
dt[ ,count(infrc), by=.(elev_bin,dbh_class)]

head(treeclim3)
#chart proportion of infected within each elevation band
ggplot(transform(treeclim3), aes(x= infrc,  group= elev_bin)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat = "count",position="dodge") +
  scale_fill_manual(values = color_table$Color) +
  labs(y = "Percent", fill="infected")+
  facet_grid(~elev_bin) + xlab("Infection Status")+
  scale_y_continuous(labels = scales::percent)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

#The website http://colorbrewer2.org/ provides a large variety of customizable palettes that are colorblind- and photocopy-safe. 
ggplot(treeclim3, aes(x=asRH, fill=dbh_class)) +geom_density(aes(y = (..count..)/sum(..count..)),alpha=0.6)+ 
  scale_y_continuous(labels=scales::percent) + facet_wrap(~infrc)+ 
    ylab("Relative Frequency")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual( name = "Size Class",values = c('#d7191c','#fdae61','#2c7bb6','#abd9e9'))+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))
ggplot(treeclim3, aes(x=asTEMP, fill=dbh_class)) +geom_density(aes(y = (..count..)/sum(..count..)),alpha=0.6)+ 
  scale_y_continuous(labels=scales::percent) + facet_wrap(~infrc)+
  ylab("Relative Frequency")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual( name = "Size Class",values = c('#d7191c','#fdae61','#2c7bb6','#abd9e9'))+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

ggplot(treeclim3, aes(x=long, fill=dbh_class)) +geom_density(aes(y = (..count..)/sum(..count..)),alpha=0.4)+ 
  scale_y_continuous(labels=scales::percent) + facet_wrap(~infrc)+
  ylab("Relative frequencies")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Size Class")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

ggplot(treeclim3, aes(x=lat, fill=dbh_class)) +geom_density(aes(y = (..count..)/sum(..count..)),alpha=0.4)+ 
  scale_y_continuous(labels=scales::percent) + facet_wrap(~infrc)+
  ylab("Relative frequencies")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Size Class")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

ggplot(treeclim3, aes(x=elev, fill=dbh_class)) +geom_density(aes(y = (..count..)/sum(..count..)),alpha=0.4)+ 
  scale_y_continuous(labels=scales::percent) + facet_wrap(~infrc)+
  ylab("Relative frequencies")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name = "Size Class")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))

plot(treeclim3$lat,treeclim4$asRH)
plot(treeclim4$long,treeclim4$asRH)

ggplot(treeclim3, aes(x = recentdbh, fill = infrc)) +geom_density(aes(y = (..count..)/sum(..count..)),alpha=0.6)+ 
  scale_y_continuous(labels=scales::percent)+#geom_bar(width=.5, position = "dodge")  +
  ylab("Relative Frequency")+ xlab("DBH (cm)")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual( name = "Infected",values = c("red","green"))+
  #scale_fill_discrete(name = "Infected")+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 30))



#Using predict with a fitted model to chart effects using the unscaled dataframe
m17.5<-glmer(data=treeclim3, infrc2~ (1|site)+recentdbh+asRH+asTEMP+asRH*asTEMP,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))

min(treeclim3$asRH);max(treeclim3$asRH)
min(treeclim3$asTEMP);max(treeclim3$asTEMP)

#change variable ranges here to identify how Temperature, dbh or site affect prob of infection
newdata2<-data.frame(asRH=rep(seq(from=36,to =71),length.out=50),asTEMP=rep(seq(from=7,to =12),length.out=50),#mean(treeclim3$aoTEMP),
                     recentdbh=rep(seq(from=20,to =20),length.out=50),site=rep(c(958.1),length.out=50));newdata2

newdata1 <- cbind(newdata2, predict(m17.5, newdata2, type ="response"));head(newdata1)
colnames(newdata1)[5]<-"fitinfrc";newdata1

ggplot(newdata1, aes(x=asRH, y=fitinfrc)) + #geom_point() + #facet_wrap(~dbh_class)+
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)+
  xlab("asRH")+ylab("Probability of Infection")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #theme(legend.position="none")+#removes legend+
  theme(text = element_text(size = 20))
#warning message thrown here can be ignored
#https://stackoverflow.com/questions/12953045/warning-non-integer-successes-in-a-binomial-glm-survey-packages

#map of relative humidity and Temperature in August to October period
#note over plotting obscures patterns better to plot in ARCmap and size yes or no's differently
ggplot(treeclim3,aes(x=long, y=lat)) + geom_point(aes(colour = (asRH)))#lower elev up north, higher elevation down south
ggplot(treeclim3,aes(x=long, y=lat)) + geom_point(aes(colour = (asTEMP)))#lower elev up north, higher elevation down south
#note in next figure cant see overplotted symbols export to Arcmap and size symbols differently
ggplot(treeclim3,aes(x=long, y=lat)) + geom_point(aes(colour = (infrc)))#lower elev up north, higher elevation down south
infected<-(subset(treeclim3,infrc2==1))#count of infected trees
notinfected<-(subset(treeclim3,infrc2==0))#count of uninfected trees
unique(infected$treeid);nrow(infected)
unique(notinfected$treeid);nrow(notinfected)
length(unique(infected$lat))


#how does infection vary by dbh when Temp and location fixed
#build data frame that varies by dbh and temperature to show how prob of infection varies with dbh and temp according to RH
#Here 3 is a fixed temperature and number following "_" is dbh
unique(treeclim3$site)

min(treeclim3$asRH);max(treeclim3$asRH)
min(treeclim3$asTEMP);max(treeclim3$asTEMP)

newdata3_2.5<-data.frame(asRH=rep(seq(from=36,to =71),length.out=50),asTEMP=rep(seq(from=7,to =12),length.out=50),#mean(treeclim3$asTEMP),
                     recentdbh=rep(seq(from=2.5,to =2.5),length.out=50),site=rep(c(-1.4593),length.out=50));newdata3_2.5
newdata3_5<-data.frame(asRH=rep(seq(from=36,to =71),length.out=50),asTEMP=rep(seq(from=7,to =12),length.out=50),#mean(treeclim3$asTEMP),
                         recentdbh=rep(seq(from=5,to =5),length.out=50),site=rep(c(-1.4593),length.out=50));newdata3_5
newdata3_15<-data.frame(asRH=rep(seq(from=36,to =71),length.out=50),asTEMP=rep(seq(from=7,to =12),length.out=50),#mean(treeclim3$asTEMP),
                         recentdbh=rep(seq(from=15,to =15),length.out=50),site=rep(c(-1.4593),length.out=50));newdata3_15
newdata3_35<-data.frame(asRH=rep(seq(from=36,to =71),length.out=50),asTEMP=rep(seq(from=7,to =12),length.out=50),#mean(treeclim3$asTEMP),
                         recentdbh=rep(seq(from=35,to =35),length.out=50),site=rep(c(-1.4593),length.out=50));newdata3_35
newsize<-rbind(newdata3_2.5,newdata3_5,newdata3_15,newdata3_35);head(newsize)
#fit newsize
newdata1 <- cbind(newsize, predict(m17, newsize, type ="response"))
colnames(newdata1)[5]<-"fitinfrc";head(newdata1);tail(newdata1)

#newdatamelt<-melt(newdata1);head(newdatamelt)
newdatacast<-melt(newdata1, id.vars=c("asRH","asTEMP","site","recentdbh"), variable.name="fitinfrc",value.name="value");head(newdatacast);tail(newdatacast)
#size<-melt(newdata1, id.vars = c("site","aoTEMP","recentdbh", "aoRH"));head(size);tail(size)

facet_names <- c(`7` = "asTEMP = 7",
  `8` = "asTEMP = 8",
  `9` = "asTEMP = 9",
  `10` = "asTEMP = 10",
  `11` = "asTEMP = 11",
  `12` = "asTEMP = 12")

ggplot(newdatacast, aes(x=asRH, y=value, color=recentdbh)) + facet_wrap(asTEMP~., labeller= as_labeller(facet_names))+#facet_wrap(~asTEMP)+ # geom_point() 
  stat_smooth(aes(group=recentdbh),method="glm", method.args=list(family="binomial"), se=FALSE, size=2)+#facet_wrap(~recentdbh)+
  xlab("Aug-Sep Relative Humidity (%)")+ylab("Probability of Infection")+ 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(color = "DBH (cm)")+#no idea why this works to change legend title, but does
  theme(text = element_text(size = 20))
  # scale_colour_gradientn(name = "Size Class",colors = c('#fdae61','#abd9e9','#2c7bb6','#d7191c'))


plot(treeclim4$long,treeclim4$asRH)
fit2<-lm(treeclim4$asRH~treeclim4$long)
summary(fit2)
abline(fit2)

plot(treeclim4$long,treeclim4$asTEMP)
fit2<-lm(treeclim4$asTEMP~treeclim4$long)
summary(fit2)
abline(fit2)

plot(treeclim4$long,treeclim4$recentdbh)
fit2<-lm(treeclim4$recentdbh~treeclim4$long)
summary(fit2)
abline(fit2)

plot(treeclim4$elev,treeclim4$recentdbh)
fit2<-lm(treeclim4$recentdbh~treeclim4$elev)
summary(fit2)
abline(fit2)

save.image()
############################################################
############################################################
#map of blister rust probability############################
############################################################

#https://stackoverflow.com/questions/30829178/predict-with-glmer-where-new-data-is-a-raster-stack-of-fixed-efefcts
library(raster)


setwd("C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\test")

####################test plots#########################
test <-raster("C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\test\\tmin\\DAYMET.003_tmax_doy2000214_aid0001.tif")
plot(test)
# Plot distribution of raster values 
testhist<-hist(test,breaks=5, main="Histogram Temperature Model GYE", col="wheat3", xlab= "Deg C")  # label the x-axis
test <-raster("C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\test\\tmax\\DAYMET.003_tmax_doy2000233_aid0001.tif")
plot(test)
dir <- "C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\test\\tmax"
test <-raster("C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\tmin\\DAYMET.003_tmin_doy2000240_aid0001.tif")
plot(test)
dir <- "C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\test\\tmin"
test <-raster("C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\tmin\\DAYMET.003_tmin_doy2000240_aid0001.tif")
plot(test)
dir <- "C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\test\\vp"


######################stack raster data for analysis##################################
#read in the daymet climate data files as raster stacks of daily values
setwd("C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\tmax")
files <- list.files(pattern = ".tif")
tmax<-stack(files)

setwd("C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\tmin")
files <- list.files(pattern = ".tif")
tmin<-stack(files)

setwd("C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\vp")
files <- list.files(pattern = ".tif")# or files<-list.files
vp<-stack(files)

#https://gis.stackexchange.com/questions/152252/creating-a-new-raster-layer-based-on-mean-of-other-raster-layers-ignoring-na-va
#s <- stack(r2000, r2001, r2002, r2003, r2004, r2005)
#x <- reclassify(s, cbind(0, NA))

#compute daily mean values for tmean, svp, rh as daily raster stacks
tmean<-(tmax+tmin)/2#this is a raster stack
svp<-610.8*10^((17.27 * tmean)/(237.3+tmean))
rh <- (vp/svp)*100

#plot(rh)#plot of panels one for each day
head(treeclim3)
head(newdata2)#data frame made earlier using observed ranges of climate variables that works with predict in model m17.5

#compute long-term mean relative humidity as single layer raster by averaging through time on the raster stacks
asRH1<- mean(rh)
asTEMP1<-mean(tmean)

#identify min and max values used at the stand levels for building model of probability.  don't predict outside this range
min(treeclim3$asRH);max(treeclim3$asRH)
min(treeclim3$asTEMP);max(treeclim3$asTEMP)

asRH <- clamp(asRH1, lower=0, upper=70.87, useValues=FALSE)#set to range of values in data set used to build model
asTEMP <- clamp(asTEMP1, lower=7.67, upper=12.03, useValues=FALSE)#
#vp_avg<-mean(vp)
plot(asRH)
plot(asTEMP)
#plot(vp_avg)
#build a constant raster of recentdbh by taking a random temperature layer, multiplying by zero then adding dbh of your choice
test <-raster("C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\tmax\\DAYMET.003_tmax_doy2000215_aid0001.tif")
plot(test)
recentdbh<-(test*0) + 30#DBh of tree size you want to model probability of infection for
plot(recentdbh)
site<-(test*0) -1.4593#+ 38.1# this is the intercept from m17.5
plot(site)

new_rast<-stack(recentdbh, asRH,asTEMP,site)#
names(new_rast)#rename the raster layers because the retain old garbage names
names(new_rast)<-c("recentdbh","asRH","asTEMP","site")#rename raster layers to exactly match climate variables used for making prediction that match glmer model
names(new_rast)
 
  #top model from glmer blister script
#m17.5<-glmer(data=treeclim3, infrc2~ (1|site)+recentdbh+asRH+asTEMP+asRH*asTEMP,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))

x <- predict(new_rast, m17.5, allow.new.levels=TRUE, type ="response",progress="text",re.form=~0)#,re.form=~0
plot(x)
y <- predict(new_rast, m17.5, allow.new.levels=FALSE, type ="response",progress="text")#const=(data.frame(site=1707.1)
plot(y)

library(rgdal)
setwd("C:/David/Water balance/blister rust/")
writeRaster(x,filename="prob_map30cm.tif",format="GTiff",overwrite=TRUE)
writeRaster(asRH,filename="asRH_map.tif",format="GTiff",overwrite=TRUE)
writeRaster(asTEMP,filename="asTEMP_map.tif",format="GTiff",overwrite=TRUE)

save.image()

# tmax_avg <- mean(tmax, na.rm=TRUE)
# tmin_avg <- mean(tmin, na.rm=TRUE)
# vp_avg <- mean(vp,na.rm=TRUE)
# gdd <- ((tmax-tmin)/2)-10
# 
# 
# #build a constant raster of recentdbh by taking a random temperature layer, multiplying by zero then adding dbh of your choice
# test <-raster("C:\\David\\Water balance\\blister rust\\daymet\\Blister rust daymet grids\\test\\tmax\\DAYMET.003_tmax_doy2000215_aid0001.tif")
# plot(test)
# test10<-(test*0) + 10
# plot(test10)
# 
# 
# #rename raster stacks to match variable names in top model
# asRH<-rh_mean
# asTEMP<-tmean
# 
# interact<-asRH*asTEMP
# 
# head(newdata2)#data frame that works with predict in model m17.5
# 
# x <- predict(asRH,asTEMP, m17.5, const=(data.frame(Year=10)))
# 
# 
# 
# m <- lmer(pa ~ red + blue + (1 | Year), data=v)
# 
# # here adding Year as a constant, as it is not a variable (RasterLayer) in the RasterStack object
# x <- predict(logo, m, const=(data.frame(Year=2000)))
# 


# Logistics Regression
#https://www.datacamp.com/community/tutorials/logistic-regression-R
#http://www.shizukalab.com/toolkits/plotting-logistic-regression-in-r
#http://www.cookbook-r.com/Statistical_analysis/Logistic_regression/



#https://cran.r-project.org/web/packages/ggiraphExtra/vignettes/ggPredict.html
fit7=glm(infrc2~elev+lat,data=treeclim4,family=binomial)
summary(fit7)
#https://stackoverflow.com/questions/36942443/plotting-a-multiple-logistic-regression-for-binary-and-continuous-values-in-r
newdata = expand.grid(elev=seq(2400,3400, length.out=100), lat=42:46, dbh_class=c(A,B,C,D))
newdata$prob = predict(fit7, newdata, type="response");head(newdata)
ggplot(newdata, aes(lat, prob, color=factor(dbh_class), group=dbh_class)) +  geom_line()



new.data <- expand.grid(elev=seq(2400,3400, length.out=100), lat=42:46)
preds <- predict(fit7, newdata = new.data, type = 'response',se = TRUE)
new.data$pred.full <- preds$fit

new.data$ymin <- new.data$pred.full - 2*preds$se.fit
new.data$ymax <- new.data$pred.full + 2*preds$se.fit  

ggplot(treeclim4,aes(x = elev, y = infrc2)) + 
  facet_wrap(~dbh_class) + 
  geom_point() + 
  geom_ribbon(data = new.data,aes(y = pred.full, ymin = ymin, ymax = ymax),alpha = 0.25) +
  geom_line(data = new.data,aes(y = pred.full),colour = "blue")






plot(treeclim4$dbh_class)
treeclim5<-subset(treeclim4, recentdbh<30 & dbh_class!="no DBH data");unique(treeclim5$dbh_class)




glm.fit <- glm(infrc ~ elev, data = treeclim2, family = binomial)
summary(glm.fit)
plot(treeclim2$elev,treeclim2$infrc,xlab="elev",ylab="Inf status")
curve(predict(glm.fit,data.frame(elev=x),type="resp"),add=TRUE)

predict(glm.fit)


#m2<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+dbh_class,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m5<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+long,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m6<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+lat,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m7<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+slope,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m8<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+aspect,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m9<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+aTEMP,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# 
# m10<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+aAET,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m11<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+ aPACK,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m12<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+aoRH,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# #m13<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+aoVPD,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m14<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+aoTEMP,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m15<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+aoRAIN,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# #m16<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+aoD,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m17<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+dRHc,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m25<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+aAET+aPACK+aoRH+aoTEMP+aoRAIN+dRHc,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
# m26<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))

m18<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+elev+long+aTEMP+aoRH+aoRAIN,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m19<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+elev+long+aTEMP+aoRH,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m20<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+elev+long+aTEMP,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m21<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+elev+long,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m22<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+aoTEMP*aoRH,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m23<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+aoTEMP*aoRAIN,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m24<-glmer(data=treeclim2scale, infrc~ (1|site)+dbh_class+elev+long+aoTEMP*aoRH,family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))


model_rank<-model.sel(m0,m1,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24);model_rank
anova(m0,m1,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18)

names(treeclim2scale)
treeclim2scale_melt<-melt(treeclim2scale,id.vars=c("panel","site","treeid","year","infected","trans","height_class","recentdbh","dbh_class","status",
                                       "infrc","elbin","infrc2", "elev_bin","elev") );head(treeclim2scale_melt);tail(treeclim2scale_melt)

ggplot(treeclim2scale_melt, aes(x=elev,y=value))+ geom_density() + facet_wrap(~variable, scales="free")

#stopped here
m1<-glmer(data=subtrans, trans~ height_class +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m2<-glmer(data=subtrans, trans~ dbh_class +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m3<-glmer(data=subtrans, trans~ elev +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m4<-glmer(data=subtrans, trans~ long +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m5<-glmer(data=subtrans, trans~ lat +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m6<-glmer(data=subtrans, trans~ slope +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m7<-glmer(data=subtrans, trans~ aspect +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m8<-glmer(data=subtrans, trans~ aTEMP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m9<-glmer(data=subtrans, trans~ aP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m10<-glmer(data=subtrans, trans~ aRAIN +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m11<-glmer(data=subtrans, trans~ aPET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m12<-glmer(data=subtrans, trans~ aAET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m13<-glmer(data=subtrans, trans~ aoD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m14<-glmer(data=subtrans, trans~ aPACK+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m15<-glmer(data=subtrans, trans~ dRHc +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m16<-glmer(data=subtrans, trans~ aoSVP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m17<-glmer(data=subtrans, trans~ aoRH +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m18<-glmer(data=subtrans, trans~ aoVPD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m19<-glmer(data=subtrans, trans~ aoTEMP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m20<-glmer(data=subtrans, trans~ aoSOIL +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m21<-glmer(data=subtrans, trans~ aoRAIN +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m22<-glmer(data=subtrans, trans~ aoAET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m23<-glmer(data=subtrans, trans~ aoD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m24<-glmer(data=subtrans, trans~ aoGDD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#single variable rankings
model_rank<-model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m17,m19,m21,m22);model_rank
anova(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m17,m19,m21,m22)
#select significant variables from single variable rankings and combine into a kitchen sink model
#Lag0 kitchen sink model
m25<-glmer(data=subtrans, trans~ dbh_class+elev+long +aspect + aTEMP+aP+aAET+aPACK+aoRH+aoAET+aoRAIN+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#the model with the significant terms is now the top model.  you can test with 
model_rank<-model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m17,m19,m21,m22,m25);model_rank
#but which of the variables are most important.  start dropping them to figure it out.
#start dropping terms from kitchen sink model one at a time
m26<-glmer(data=subtrans, trans~ dbh_class+elev+long+aspect +aTEMP+aPET+aoD+aPACK+aoRH+aoTEMP+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m27<-glmer(data=subtrans, trans~ dbh_class+elev+long+aspect +aTEMP+aPET+aoD+aPACK+aoRH+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m28<-glmer(data=subtrans, trans~ dbh_class+elev+long+aspect +aTEMP+aPET+aoD+aPACK+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m29<-glmer(data=subtrans, trans~ dbh_class+elev+long+aspect +aTEMP+aPET+aoD+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m30<-glmer(data=subtrans, trans~ dbh_class+elev+long+aspect +aTEMP+aPET+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m31<-glmer(data=subtrans, trans~ dbh_class+elev+long+aspect +aTEMP+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m32<-glmer(data=subtrans, trans~ dbh_class+elev+long+aspect +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m33<-glmer(data=subtrans, trans~ dbh_class+elev+long +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m34<-glmer(data=subtrans, trans~ dbh_class+elev+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m35<-glmer(data=subtrans, trans~ dbh_class+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))

#then recalculate model ranks
model_rank<-model.sel(m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35);model_rank
anova(m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35)

##################################################################
#########################Lag1#####################################
subtrans<-transwb_ynL1[ , !(names(transwb_ynL1) %in% drops)]
names(subtrans)
#rescale RH from 0 to 1 to 0 to 100:  this avoids warning about some variables being on very different scales
subtrans$aoRH<-(subtrans$aoRH)*100;head(subtrans$aoRH)
m41<-glmer(data=subtrans, trans~ aTEMP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m42<-glmer(data=subtrans, trans~ aP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m43<-glmer(data=subtrans, trans~ aRAIN +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m44<-glmer(data=subtrans, trans~ aPET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m45<-glmer(data=subtrans, trans~ aAET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m46<-glmer(data=subtrans, trans~ aoD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m47<-glmer(data=subtrans, trans~ aPACK+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m48<-glmer(data=subtrans, trans~ dRHc +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m16<-glmer(data=subtrans, trans~ aoSVP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m49<-glmer(data=subtrans, trans~ aoRH +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m18<-glmer(data=subtrans, trans~ aoVPD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m50<-glmer(data=subtrans, trans~ aoTEMP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m20<-glmer(data=subtrans, trans~ aoSOIL +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m51<-glmer(data=subtrans, trans~ aoRAIN +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m52<-glmer(data=subtrans, trans~ aoAET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m23<-glmer(data=subtrans, trans~ aoD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m24<-glmer(data=subtrans, trans~ aoGDD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#single variable rankings
model_rank<-model.sel(m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52);model_rank
anova(m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52)
#select significant variables from single variable rankings and combine into a kitchen sink model
#Lag0 kitchen sink model
m53<-glmer(data=subtrans, trans~ dbh_class+elev+long +aspect + aP+aAET+dRHc+aoRH+aoAET+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#the model with the significant terms is now the top model.  you can test with 
model_rank<-model.sel(m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53);model_rank
#but which of the variables are most important.  start dropping them to figure it out.
#start dropping terms from kitchen sink model one at a time
m54<-glmer(data=subtrans, trans~ dbh_class+elev+long +aspect + aP+aAET+dRHc+aoRH+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m55<-glmer(data=subtrans, trans~ dbh_class+elev+long +aspect + aP+aAET+dRHc+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m56<-glmer(data=subtrans, trans~ dbh_class+elev+long +aspect + aP+aAET+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m57<-glmer(data=subtrans, trans~ dbh_class+elev+long +aspect + aP+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#stop reducing models here b/c m57 without aTEMP is same as m31
#then recalculate model ranks
model_rank<-model.sel(m53,m54,m55,m56,m57);model_rank
anova(m53,m54,m55,m56,m57)


##################################################################
#########################Lag2#####################################
subtrans<-transwb_ynL2[ , !(names(transwb_ynL2) %in% drops)]
names(subtrans)
#rescale RH from 0 to 1 to 0 to 100:  this avoids warning about some variables being on very different scales
subtrans$aoRH<-(subtrans$aoRH)*100;head(subtrans$aoRH)
m70<-glmer(data=subtrans, trans~ aTEMP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m71<-glmer(data=subtrans, trans~ aP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m72<-glmer(data=subtrans, trans~ aRAIN +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m73<-glmer(data=subtrans, trans~ aPET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m74<-glmer(data=subtrans, trans~ aAET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m75<-glmer(data=subtrans, trans~ aoD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m76<-glmer(data=subtrans, trans~ aPACK+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m77<-glmer(data=subtrans, trans~ dRHc +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m16<-glmer(data=subtrans, trans~ aoSVP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m78<-glmer(data=subtrans, trans~ aoRH +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m18<-glmer(data=subtrans, trans~ aoVPD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m79<-glmer(data=subtrans, trans~ aoTEMP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m20<-glmer(data=subtrans, trans~ aoSOIL +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m80<-glmer(data=subtrans, trans~ aoRAIN +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m81<-glmer(data=subtrans, trans~ aoAET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m23<-glmer(data=subtrans, trans~ aoD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m24<-glmer(data=subtrans, trans~ aoGDD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#single variable rankings
model_rank<-model.sel(m70,m71,m72,m73,m74,m75,m76,m77,m78,m79,m80,m81);model_rank
anova(m70,m71,m72,m73,m74,m75,m76,m77,m78,m79,m80,m81)
#select significant variables from single variable rankings and combine into a kitchen sink model
#Lag0 kitchen sink model
m82<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aRAIN+aAET+aoD+aoRH+aoRAIN+aoAET+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#the model with the significant terms is now the top model.  you can test with 
#model_rank<-model.sel(m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53);model_rank
#but which of the variables are most important.  start dropping them to figure it out.
#start dropping terms from kitchen sink model one at a time
m82<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aRAIN+aAET+aoD+aoRH+aoRAIN++(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m83<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aRAIN+aAET+aoD+aoRH+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m84<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aRAIN+aAET+aoD+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m85<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aRAIN+aAET+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m86<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aRAIN+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m87<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + (1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))


##################################################################
#########################Lag3#####################################
subtrans<-transwb_ynL3[ , !(names(transwb_ynL3) %in% drops)]
names(subtrans)
#rescale RH from 0 to 1 to 0 to 100:  this avoids warning about some variables being on very different scales
subtrans$aoRH<-(subtrans$aoRH)*100;head(subtrans$aoRH)
m90<-glmer(data=subtrans, trans~ aTEMP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m91<-glmer(data=subtrans, trans~ aP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m92<-glmer(data=subtrans, trans~ aRAIN +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m93<-glmer(data=subtrans, trans~ aPET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m94<-glmer(data=subtrans, trans~ aAET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m95<-glmer(data=subtrans, trans~ aoD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m96<-glmer(data=subtrans, trans~ aPACK+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m97<-glmer(data=subtrans, trans~ dRHc +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m16<-glmer(data=subtrans, trans~ aoSVP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m98<-glmer(data=subtrans, trans~ aoRH +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m18<-glmer(data=subtrans, trans~ aoVPD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m99<-glmer(data=subtrans, trans~ aoTEMP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m20<-glmer(data=subtrans, trans~ aoSOIL +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m100<-glmer(data=subtrans, trans~ aoRAIN +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m101<-glmer(data=subtrans, trans~ aoAET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m23<-glmer(data=subtrans, trans~ aoD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m24<-glmer(data=subtrans, trans~ aoGDD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#single variable rankings
model_rank<-model.sel(m90,m91,m92,m93,m94,m95,m96,m97,m98,m99,m100,m101);model_rank
anova(m90,m91,m92,m93,m94,m95,m96,m97,m98,m99,m100,m101)
#select significant variables from single variable rankings and combine into a kitchen sink model
#Lag0 kitchen sink model
m102<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+aAET+aPACK+aoRH+aoTEMP+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#the model with the significant terms is now the top model.  you can test with 
#model_rank<-model.sel(m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53);model_rank
#but which of the variables are most important.  start dropping them to figure it out.
#start dropping terms from kitchen sink model one at a time
m103<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+aAET+aPACK+aoRH+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m104<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+aAET+aPACK+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m105<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+aAET+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m106<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m107<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + (1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#stop reducing models here b/c m58 without aTEMP is same as m31
#then recalculate model ranks
model_rank<-model.sel(m102,m103,m104,m105,m106,m107);model_rank
anova(m102,m103,m104,m105,m106,m107)

##################################################################
#########################Lag4#####################################
subtrans<-transwb_ynL4[ , !(names(transwb_ynL4) %in% drops)]
names(subtrans)
#rescale RH from 0 to 1 to 0 to 100:  this avoids warning about some variables being on very different scales
subtrans$aoRH<-(subtrans$aoRH)*100;head(subtrans$aoRH)
m120<-glmer(data=subtrans, trans~ aTEMP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m121<-glmer(data=subtrans, trans~ aP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m122<-glmer(data=subtrans, trans~ aRAIN +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m123<-glmer(data=subtrans, trans~ aPET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m124<-glmer(data=subtrans, trans~ aAET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m125<-glmer(data=subtrans, trans~ aoD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m126<-glmer(data=subtrans, trans~ aPACK+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m127<-glmer(data=subtrans, trans~ dRHc +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m16<-glmer(data=subtrans, trans~ aoSVP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m128<-glmer(data=subtrans, trans~ aoRH +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m18<-glmer(data=subtrans, trans~ aoVPD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m129<-glmer(data=subtrans, trans~ aoTEMP +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m20<-glmer(data=subtrans, trans~ aoSOIL +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m130<-glmer(data=subtrans, trans~ aoRAIN +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m131<-glmer(data=subtrans, trans~ aoAET +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m23<-glmer(data=subtrans, trans~ aoD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#m24<-glmer(data=subtrans, trans~ aoGDD +(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#single variable rankings
model_rank<-model.sel(m120,m121,m122,m123,m124,m125,m126,m127,m128,m129,m130,m131);model_rank
anova(m120,m121,m122,m123,m124,m125,m126,m127,m128,m129,m130,m131)
#select significant variables from single variable rankings and combine into a kitchen sink model
#Lag0 kitchen sink model
m132<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+aAET+aRAIN+aoD+dRHc+aoRH+aoRAIN+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
#the model with the significant terms is now the top model.  you can test with 
#model_rank<-model.sel(m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53);model_rank
#but which of the variables are most important.  start dropping them to figure it out.
#start dropping terms from kitchen sink model one at a time
m133<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+aAET+aRAIN+aoD+dRHc+aoRH+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m134<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+aAET+aRAIN+aoD+dRHc+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m135<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+aAET+aRAIN+aoD+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m136<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+aAET+aRAIN+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m137<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+aAET+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))
m138<-glmer(data=subtrans, trans~ dbh_class+elev+long + aspect + aP+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))

#stop reducing models here b/c m58 without aTEMP is same as m31
#then recalculate model ranks
model_rank<-model.sel(m132,m133,m134,m135,m136,m137,m138);model_rank
anova(m132,m133,m134,m135,m136,m137,m138)

#giant lag comparison
model_rank_lag<-model.sel(m25,m26,m27,m53,m54,m55,m56,m57,m72,m75,m78,m79,m80,m81,m102,m132)


##########end here#######################################
#below is an attempt to run all the coefficients through glmer in one call
names(subtrans)
df<-transwb_yn[,c(1,3,4:23)];head(df)
models <- lapply(paste(names(df), '~ .'), function(f){
  eval(substitute(data = df, glmer(trans~(1|site)+(1|tree)+ ., family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap")), 
                  list(frm = as.formula(f))))
})
models


models <- vector(mode = "list", length = ncol(transwb_yn)-15)
test_cols<-seq(from = 15, to = ncol(transwb_yn), by = 1);test_cols
var<-NULL
coeffs_table_nonforest<-NULL
for (j in test_cols){#note model number reflects panel order in ggplot figures following the model comparison where 
  var<-colnames(transwb_yn[j]);var
  models[[j-2]] <-glmer(data=transwb_yn, trans~ var+colnames(transwb_yn[j])+(1|site)+(1|tree),family=binomial(link = "logit"),nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))#linear model as exponential fuction
  var[j-2]<-colnames(transwb_yn)[j]
  coeffs<-as.data.frame(summary(models[[j-2]])[10], stringsAsFactors = FALSE)
  coeffs2<-cbind("nonforest",var[j-2],coeffs[1,1],coeffs[1,2])
  coeffs_table_nonforest<-rbind(coeffs_table_nonforest,coeffs2)
}
colnames(coeffs_table_nonforest)[1:4]<-c("maj_class","variable","a","b");coeffs_table_nonforest
model_rank<-model.sel(models);nonforest_model_rank;var




### Added by shuysman <2023-07-19 Wed> investigating issue reproducing models from paper

### Blister rust probability models
### Reflect m17 and m17.5 from paper.
### m17.5 appears to be the model actually used to create the
### probabiltiy x asRH/asTEMP/recentDBH plots

predict_br_17 <- function (rh, tavg, dbh = 1) {
    ## m17 from paper
    intercept <- -1.459
    coef_dbh <- 0.062
    coef_rh <- 0.738
    coef_tavg <- 0.541
    coef_interaction <- -0.435

    logit <- intercept + (coef_dbh * dbh) + (coef_rh * rh) + (coef_tavg * tavg) + (coef_interaction * rh * tavg)
    odds <- exp(logit)
    p <- odds / (1 + odds)
    return(p)
}

curve(predict_br_17(rh = x, tavg = 8, dbh = 10), from = 0, to = 1)

predict_br_17.5 <- function (rh, tavg, dbh = 1) {
    ## m17.5 from blister.R
    intercept <- -53.51028
    coef_dbh <- 0.06177
    coef_rh <- 0.78306
    coef_tavg <- 4.71199
    coef_interaction <- -0.06966

    logit <- intercept + (coef_dbh * dbh) + (coef_rh * rh) + (coef_tavg * tavg) + (coef_interaction * rh * tavg)
    odds <- exp(logit)
    p <- odds / (1 + odds)
    return(p)
}

png("img/m17_5_curve.png")
m17.5_curve <- curve(predict_br_17.5(rh = x, tavg = 8, dbh = 35), from = 40, to = 70)
dev.off()

## Can also use the model objects straight from the paper
## attach("./blister.RData")


library(reshape2) ## for melt()

ggpredict(m17.5)
mydf <- ggpredict(m17.5, terms = c("asRH"), back.transform = TRUE)
ggplot(mydf, aes(x, predicted, colour = group)) + geom_line()



### Reproduce asRH x infection probability plot from paper
### Need to use m17.5 to match the figure from the paper
newdata3_2.5 <- data.frame(
    asRH = rep(seq(from = 36, to = 71), length.out = 50),
    asTEMP = rep(seq(from = 7, to = 12), length.out = 50),
    recentdbh = rep(seq(from = 2.5, to = 2.5), length.out = 50)
)
newdata3_2.5

newdata3_5 <- data.frame(
    asRH = rep(seq(from = 36, to = 71), length.out = 50),
    asTEMP = rep(seq(from = 7, to = 12), length.out = 50),
    recentdbh = rep(seq(from = 5, to = 5), length.out = 50)
)
newdata3_5


newdata3_15 <- data.frame(
    asRH = rep(seq(from = 36, to = 71), length.out = 50),
    asTEMP = rep(seq(from = 7, to = 12), length.out = 50),
    recentdbh = rep(seq(from = 15, to = 15), length.out = 50)
)
newdata3_15

newdata3_35 <- data.frame(
    asRH = rep(seq(from = 36, to = 71), length.out = 50),
    asTEMP = rep(seq(from = 7, to = 12), length.out = 50),
    recentdbh = rep(seq(from = 35, to = 35), length.out = 50)
)
newdata3_35

newsize <- rbind(newdata3_2.5, newdata3_5, newdata3_15, newdata3_35)
head(newsize)
## fit newsize
predictions <- predict(m17.5, newdata = newsize, type = "response", re.form = ~0)
newdata1 <- cbind(newsize, predictions)
colnames(newdata1)[4] <- "fitinfrc"
head(newdata1)
tail(newdata1)

# newdatamelt<-melt(newdata1);head(newdatamelt)
newdatacast <- melt(newdata1, id.vars = c("asRH", "asTEMP", "recentdbh"), variable.name = "fitinfrc", value.name = "value")
head(newdatacast)
tail(newdatacast)
# size<-melt(newdata1, id.vars = c("site","aoTEMP","recentdbh", "aoRH"));head(size);tail(size)

facet_names <- c(
    `7` = "asTEMP = 7",
    `8` = "asTEMP = 8",
    `9` = "asTEMP = 9",
    `10` = "asTEMP = 10",
    `11` = "asTEMP = 11",
    `12` = "asTEMP = 12"
)

ggplot(newdatacast, aes(x = asRH, y = value, color = recentdbh)) +
    facet_wrap(asTEMP ~ ., labeller = as_labeller(facet_names), scales = "fixed") + # facet_wrap(~asTEMP)+ # geom_point()
    stat_smooth(aes(group = recentdbh), method = "glm", method.args = list(family = "binomial"), se = FALSE, size = 2) + # facet_wrap(~recentdbh)+
    xlab("Aug-Sep Relative Humidity (%)") +
    ylab("Probability of Infection") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
    ) +
    labs(color = "DBH (cm)") + # no idea why this works to change legend title, but does
    theme(text = element_text(size = 20))
# scale_colour_gradientn(name = "Size Class",colors = c('#fdae61','#abd9e9','#2c7bb6','#d7191c'))
ggsave("img/m17_5_plot.png", width = 12)

