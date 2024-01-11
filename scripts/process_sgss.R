library(lubridate)
library(here)
library(RColorBrewer)

sgss_ages<- read.csv(here("data","sgss_all_age.csv"), header=TRUE)#, sep=,)
head(sgss_ages)

dates2<-lubridate::ymd("2014-01-01") + weeks(sgss_ages$Reporting_week - 1) + years(sgss_ages$Year - 2014)

data<-t(as.matrix(sgss_ages[,3:6]))
rownames(data) <- c("a0_4","a15_64","a5_14","a65p")


coul <- brewer.pal(4,'Spectral') 
# Transform this data in %
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})

agg_age<-rowSums(data, na.rm = T)/sum(rowSums(data, na.rm = T))
barplot(agg_age)
pie(agg_age,labels = rownames(data), main="Age distribution of Noro cases reported, England & Wales (2014-2023)")

total_cases<-colSums(data, na.rm = T)


# Make a stacked barplot--> it will be in %!
windows()
sgss_vec<-seq(dates2[2], dates2[length(dates2)], by='6 months')
barplot(data_percentage, col=coul , border="white",space = 0,
        legend.text = TRUE, 
        xaxt = "n",
        args.legend = list(x = "topright",
                           inset = c(- 0.05, -0.15), cex=0.7))
axis(side= 1, at=seq(1,length(dates2), by =24), 
     labels = dates2[c(seq(1,length(dates2), by =24))],
     las=2,cex.axis=0.8, lwd=0, lwd.tick=1)


windows()
barplot(data, col=coul , border="white", space = 0,
        legend.text = TRUE, 
        xaxt = "n",
        args.legend = list(x = "topright",
                           inset = c(- 0.05, -0.15), cex=0.7))
axis(side= 1, at=seq(1,length(dates2), by =24), 
     labels = dates2[c(seq(1,length(dates2), by =24))],
     las=2,cex.axis=0.8, lwd=0, lwd.tick=1)









####### Strains

sgss_str<- read.csv(here("data","sgss_strain_individual.csv"), header=TRUE)#, sep=,)
head(sgss_str)
table(sgss_str$Characterisation.results)
id1<-which(sgss_str$Characterisation.results=="GI indeterminate")
sgss_str<-sgss_str[-id1,]
id2<-which(sgss_str$Characterisation.results=="GII indeterminate")
sgss_str<-sgss_str[-id2,]
id3<-which(sgss_str$Characterisation.results=="GII.4 and GI.3")
sgss_str<-sgss_str[-id3,]
id4<-which(sgss_str$Characterisation.results=="GII.4 and other GI")
sgss_str<-sgss_str[-id4,]
id5<-which(sgss_str$Characterisation.results=="Mixed indeterminate")
sgss_str<-sgss_str[-id5,]
id6<-which(sgss_str$Characterisation.results=="Other GII and other GI")
sgss_str<-sgss_str[-id6,]
table(sgss_str$Characterisation.results)
head(sgss_str)

sgss_str<-sgss_str %>%
  group_by(Year,ISOweek,Characterisation.results) %>%
  summarize(n_observations = n())
head(sgss_str)

sgss_str


sgss_str_wide <- spread(sgss_str, Characterisation.results, n_observations)
head(sgss_str_wide)
dates<-lubridate::ymd("2018-01-01") + weeks(sgss_str_wide$ISOweek - 1) + 
  years(sgss_str_wide$Year - sgss_str_wide$Year[1])

data<-t(as.matrix(sgss_str_wide[,3:6]))
rownames(data) <- c("GI.3","GII.4","O-GI","O-GII")
id<-which(is.na(data[1,]))
data[1,id]<-0
id<-which(is.na(data[2,]))
data[2,id]<-0
id<-which(is.na(data[3,]))
data[3,id]<-0
id<-which(is.na(data[4,]))
data[4,id]<-0

coul <- brewer.pal(4,'Spectral') 
# Transform this data in %
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})

# Aggreagted
agg_strain<-rowSums(data, na.rm = T)/sum(rowSums(data, na.rm = T))
barplot(agg_strain)
pie(agg_strain,labels = rownames(data), main="Genotype distribution of Noro cases reported, England & Wales (2014-2023)")



# Make a stacked barplot--> it will be in %!
windows()
barplot(data_percentage, col=coul , border="white",space = 0,
        legend.text = TRUE, 
        xaxt = "n",
        args.legend = list(x = "topright",
                           inset = c(- 0.05, -0.15), cex=0.7))
axis(side= 1, at=seq(1,length(dates), by =10), 
     labels = dates[c(seq(1,length(dates), by =10))],
     las=2,cex.axis=0.8, lwd=0, lwd.tick=1)



windows()
barplot(data, col=coul , border="white",space = 0,
        legend.text = TRUE, 
        xaxt = "n",
        args.legend = list(x = "topright",
                           inset = c(- 0.05, -0.15), cex=0.7))
axis(side= 1, at=seq(1,length(dates), by =10), 
     labels = dates[c(seq(1,length(dates), by =10))],
     las=2,cex.axis=0.8, lwd=0, lwd.tick=1)





















sgss<- read.csv(here("data","sgss_weekly_cases.csv"), header=TRUE)#, sep=,)

sgss<-sgss %>% 
  arrange(week)
sgss_start<-as.Date('2010-07-01')
sgss_end<-as.Date('2011-06-20')

sgss_vec<-seq(ymd(sgss_start), ymd(sgss_end), by='1 week')
wknum<-lubridate::isoweek(sgss_vec)
span_sgss<-which(days_vec%in%sgss_vec)
sgss$cases<-round(sgss$cases)
sgss$day<-span_sgss


# Create data
set.seed(1124)
data <- matrix(sample(1:30,15) , nrow=3)
colnames(data) <- c("A","B","C","D","E")
rownames(data) <- c("var1","var2","var3")

# create color palette:
library(RColorBrewer)
coul <- brewer.pal(3, "Pastel2") 

# Transform this data in %
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})

# Make a stacked barplot--> it will be in %!
barplot(data_percentage, col=coul , border="white", xlab="group")
