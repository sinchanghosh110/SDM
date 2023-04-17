setwd("D:/SDM_blue_tailed_bee_eater/Prey")
getwd()

library(dismo)
library(rgbif)
library(sp)
library(raster)
library(spatialEco)
#789

sname=c(789,797,1457,1458)
years=as.character(seq(from=2010,to=2018,length.out=9))
mn=as.character(seq(from=1,to=12,length.out=12))
years
r_w=array(data=NA,dim = c(12,9,4))
sn_w=array(data=NA,dim = c(12,9,4))
sp_w=array(data=NA,dim = c(12,9,4))
r_w
for (i in 2:2) {
 for (j in 7:7) {
    for (k in 12:12) {
  
  prey<-occ_data(orderKey=sname[i], year = years[j], month=mn[k],limit=10000)
  birds<-occ_data(scientificName="Merops philippinus", year = years[j], month=mn[k],limit=10000)
  
  data1=prey$data
  data1$lon=data1$decimalLongitude
  data1$lat=data1$decimalLatitude
  data1$year=data1$year
  data1$month=data1$month
  prey1=data1[!is.na("lon") & !is.na("lat") & !is.na("year")& !is.na("month")]
  prey1=prey1[,c("lon","lat")]
  prey1=data.frame(prey1)
  prey1=na.omit(prey1)
  prey1$x=rep(1,nrow(prey1))
  
  data2=birds$data
  data2$lon=data2$decimalLongitude
  data2$lat=data2$decimalLatitude
  data2$year=data2$year
  data2$month=data2$month
  birds1=data2[!is.na("lon") & !is.na("lat") & !is.na("year")& !is.na("month")]
  birds1=birds1[,c("lon","lat")]
  birds1=na.omit(birds1)
  birds1$x=rep(1,nrow(birds1))
  
  
  coordinates(prey1)=~lon+lat
  proj4string(prey1)<-CRS("+proj=longlat +datum=WGS84")
  coordinates(birds1)=~lon+lat
  proj4string(birds1)<-CRS("+proj=longlat +datum=WGS84")
  
  rp=raster(prey1);res(rp)=0.00001
  p_1=gridSample(prey1,rp,n=1)
  b_1=gridSample(birds1,rp,n=1)
 
  p_1=data.frame(p_1);b_1=data.frame(b_1)
  p_1$x=rep(1,nrow(p_1));b_1$x=rep(1,nrow(b_1))
  colnames(p_1)=colnames(b_1)=c("lon","lat","x")
  
  
  r_obj <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, resolution=c(10,10))
  r_p<- rasterize(x=p_1[, 1:2], # lon-lat data
                             y=r_obj, 
                  field=p_1[, 3],
                             fun=mean) # aggregate
  
  r_b<- rasterize(x=b_1[, 1:2], # lon-lat data
                  y=r_obj, 
                  field=b_1[, 3],
                  fun=mean) # aggregate
  
  r_c=raster(r_p)
  
  r_p[]=ifelse(is.na(r_p[])==T,0,1)
  r_b[]=ifelse(is.na(r_b[])==T,0,1)
  r_c[]=ifelse(r_b[]==r_p[],r_b[],0)
  #r_w[k,j,i]=sum(r_c[]==1)/sum(r_b[]==1)
  sn_w[k,j,i]=sum(r_c[]==1)/(sum(r_p[]==1))
  sp_w[k,j,i]=sum(r_c[]==0)/(sum(r_c[]==0)+sum((r_c[]-r_b[])==-1))
 #
  print(sp_w[k,j,i])
  
  #r_c[]=ifelse(r_b[]==r_p[]| r)
  #r.cor <- rasterCorrelation(r_b, r_p, s = 11, type = "pearson")
  #plot(r.cor)
                
  #writeRaster(r_c,paste(sname[i], year = years[j], month=mn[k], sep =""),format="GTiff",overwrite=TRUE)

  }}}
r_w
sn_w
sp_w
#######
psn=read.csv("Prey_sensitivity.csv", header = T)
psp=read.csv("Prey_specificity.csv", header = T)
names(psp)
psn=psn[,-15]
psp=psp[,-15]
View(psn1)
#cr=numeric(12)
psn1=psn[1:9,-c(1,2)]
psn2=psn[10:18,-c(1,2)]
psn3=psn[19:27,-c(1,2)]
psp1=psp[1:9,-c(1,2)]
psp2=psp[10:18,-c(1,2)]
psp3=psp[19:27,-c(1,2)]
dim(psn1)
csn1=csn2=csn3=csp1=csp2=csp3=NULL
#psn1=cbind.data.frame(psn1,years)
#as.numeric(psp3[1,])
for (i in 1:12) {
  csn1[i]=list(as.numeric(psn1[,i]),csn1[i-1])
  csn2[i]=list(as.numeric(psn2[,i]),csn2[i-1])
  csn3[i]=list(as.numeric(psn3[,i]),csn3[i-1])
  csp1[i]=list(as.numeric(psp1[,i]),csp1[i-1])
  csp2[i]=list(as.numeric(psp2[,i]),csp2[i-1])
  csp3[i]=list(as.numeric(psp3[,i]),csp3[i-1])
  
}
psp3
#csn1
print(bartlett.test(csn1))
print(bartlett.test(csn2))
print(bartlett.test(csn3))
print(bartlett.test(csp1))
print(bartlett.test(csp2))
print(bartlett.test(csp3))






###################
for (i in 1:9) {
  csn1[i]=chisq.test(psn1[i,-c(1,2)])
  print(chisq.test(psn1[i,-c(1,2)]))
  csn2[i]=chisq.test(psn2[i,-c(1,2)])
  print(chisq.test(psn2[i,-c(1,2)]))
  csn3[i]=chisq.test(psn3[i,-c(1,2)])
  print(chisq.test(psn3[i,-c(1,2)]))
  
  csp1[i]=chisq.test(psp1[i,-c(1,2)])
  print(chisq.test(psp1[i,-c(1,2)]))
  csp2[i]=chisq.test(psp2[i,-c(1,2)])
  print(chisq.test(psp2[i,-c(1,2)]))
  csp3[i]=chisq.test(psp3[i,-c(1,2)])
  print(chisq.test(psp3[i,-c(1,2)]))
}

psn1=psn1[,-c(1,2)];psn2=psn2[,-c(1,2)];psn3=psn3[,-c(1,2)]
psp1=psp1[,-c(1,2)];psp2=psp2[,-c(1,2)];psp3=psp3[,-c(1,2)]
tpr1=1-psp1;tpr2=1-psp2;tpr3=1-psp3
rocc1=NULL


for (i in 1:9) {
  rocc1_1=lm(as.numeric(psn1[1,])~as.numeric(tpr1[1,]))
  
}
rocc1
plot(rocc1_1)
tpsn1=c(psn1[1,], psn1[2,], psn1[3,],psn1[4,], psn1[5,],psn1[6,], psn1[7,,],psn1[8,],psn1[9,])
tpsn2=c(psn2[1,], psn2[2,], psn2[3,],psn2[4,], psn2[5,],psn2[6,], psn2[7,,],psn2[8,],psn2[9,])
tpsn3=c(psn3[1,], psn3[2,], psn3[3,],psn3[4,], psn3[5,],psn3[6,], psn3[7,,],psn3[8,],psn3[9,])

tpsn1=as.numeric(tpsn1)
tpsn2=as.numeric(tpsn2)
tpsn3=as.numeric(tpsn3)
pcan=prcomp(x=cbind.data.frame(tpsn1,tpsn2,tpsn3))
biplot(pcan)
psn$Prey
plot(tpsn1, pch=".", cex=3, ylab = "Sensitivity for Odonata")
lines(tpsn1)

plot(tpsn2, pch=".", cex=3, ylab = "Sensitivity for Lepidoptera")
lines(tpsn2)

plot(tpsn3, pch=".", cex=3, ylab = "Sensitivity for Hymenoptera")
lines(tpsn3)


tpsp1=c(psp1[1,], psp1[2,], psp1[3,],psp1[4,], psp1[5,],psp1[6,], psp1[7,,],psp1[8,],psp1[9,])
tpsp2=c(psp2[1,], psp2[2,], psp2[3,],psp2[4,], psp2[5,],psp2[6,], psp2[7,,],psp2[8,],psp2[9,])
tpsp3=c(psp3[1,], psp3[2,], psp3[3,],psp3[4,], psp3[5,],psp3[6,], psp3[7,,],psp3[8,],psp3[9,])

tpsp1=as.numeric(tpsp1)
tpsp2=as.numeric(tpsp2)
tpsp3=as.numeric(tpsp3)

pcap=prcomp(x=cbind.data.frame(tpsp1,tpsp2,tpsp3))
biplot(pcap)
psp$Prey
plot(tpsp1, pch=".", cex=3, ylab = "specificity for Odonata")
lines(tpsp1)

plot(tpsp2, pch=".", cex=3, ylab = "specificity for Lepidoptera")
lines(tpsp2)

plot(tpsp3, pch=".", cex=3, ylab = "specificity for Hymenoptera")
lines(tpsp3)

#plot(sort(c(1-tpsp1)[1:12],tpsn1[1:12]),type="l", lty=2, lwd=1, ylab = "ROC curve")

#tpsn1=numeric(108)
rocc1=lm(tpsn1~1-tpsp1)
rocc1
lines(predict(rocc1))
plot(tpsn1~c(1-tpsp1))
abline(fit)

library(ggplot2)
length(tpsn1)
#create scatterplot with fitted regression line
ggplot(data.frame(tpsn1,seq(1,108)), aes(x = seq(1,108), y = tpsn1)) + 
  geom_point() +
  xlab("Time point") +
  ylab("Sensitivity to Odonata") +
  stat_smooth(method = "gam")

ggplot(data.frame(tpsn2,seq(1,108)), aes(x = seq(1,108), y = tpsn2)) + 
  geom_point() +
  xlab("Time point") +
  ylab("Sensitivity to Lepidoptera") +
  stat_smooth(method = "gam")

ggplot(data.frame(tpsn3,seq(1,108)), aes(x = seq(1,108), y = tpsn3)) + 
  geom_point() +
  xlab("Time point") +
  ylab("Sensitivity to Hymenoptera") +
  stat_smooth(method = "gam")


ggplot(data.frame(tpsp1,seq(1,108)), aes(x = seq(1,108), y = tpsp1)) + 
  geom_point() +
  xlab("Time point") +
  ylab("Specificity to Odonata") +
  stat_smooth(method = "gam")

ggplot(data.frame(tpsp2,seq(1,108)), aes(x = seq(1,108), y = tpsp2)) + 
  geom_point() +
  xlab("Time point") +
  ylab("Specificity to Lepidoptera") +
  stat_smooth(method = "gam")

ggplot(data.frame(tpsp3,seq(1,108)), aes(x = seq(1,108), y = tpsp3)) + 
  geom_point() +
  xlab("Time point") +
  ylab("Specificity to Hymenoptera") +
  stat_smooth(method = "gam")
#Moving average
tpsn1_ma=arima(tpsn1, order = c(0, 0, 0),seasonal = list(order = c(2L, 2L, 2L)))
tpsn1_ma
residualsn1=residuals(tpsn1_ma)
tpsn1_fitted <- tpsn1 - residualsn1
ts.plot(tpsn1, ylab="Sensitivity to Odonata")
points(tpsn1_fitted, type = "l", col = 2, lty = 2, lwd=3)
legend(x = "topright",          # Position
       legend = c("Original", "Fitted moving average"),  # Legend texts
       lty = c(1, 2),           # Line types
       col = c("black", "red"),           # Line colors
       lwd = c(2,3)) 

tpsn2_ma=arima(tpsn2, order = c(0, 0, 0),seasonal = list(order = c(2L, 2L, 2L)))
tpsn2_ma
residualsn2=residuals(tpsn2_ma)
tpsn2_fitted <- tpsn2 - residualsn2
ts.plot(tpsn2, ylab="Sensitivity to Lepidoptera")
points(tpsn2_fitted, type = "l", col = 2, lty = 2, lwd=3)
legend(x = "topright",          # Position
       legend = c("Original", "Fitted moving average"),  # Legend texts
       lty = c(1, 2),           # Line types
       col = c("black", "red"),           # Line colors
       lwd = c(2,3)) 

tpsn3_ma=arima(tpsn3, order = c(0, 0, 0),seasonal = list(order = c(2L, 2L, 2L)))
tpsn3_ma
residualsn3=residuals(tpsn3_ma)
tpsn3_fitted <- tpsn3 - residualsn3
ts.plot(tpsn3, ylab="Sensitivity to Hymenoptera")
points(tpsn3_fitted, type = "l", col = 2, lty = 2, lwd=3)
legend(x = "topright",          # Position
       legend = c("Original", "Fitted moving average"),  # Legend texts
       lty = c(1, 2),           # Line types
       col = c("black", "red"),           # Line colors
       lwd = c(2,3)) 

tpsp1_ma=arima(tpsp1, order = c(0, 0, 0),seasonal = list(order = c(2L, 2L, 2L)))
tpsp1_ma
residualsp1=residuals(tpsp1_ma)
tpsp1_fitted <- tpsp1 - residualsp1
ts.plot(tpsp1, ylab="Specificity to Odonata")
points(tpsp1_fitted, type = "l", col = 2, lty = 2, lwd=3)
legend(x = "topright",          # Position
       legend = c("Original", "Fitted moving average"),  # Legend texts
       lty = c(1, 2),           # Line types
       col = c("black", "red"),           # Line colors
       lwd = c(2,3)) 

tpsp2_ma=arima(tpsp2, order = c(0, 0, 0),seasonal = list(order = c(2L, 2L, 2L)))
tpsp2_ma
residualsp2=residuals(tpsp2_ma)
tpsp2_fitted <- tpsp2 - residualsp2
ts.plot(tpsp2, ylab="Specificity to Lepidoptera")
points(tpsp2_fitted, type = "l", col = 2, lty = 2, lwd=3)
legend(x = "topright",          # Position
       legend = c("Original", "Fitted moving average"),  # Legend texts
       lty = c(1, 2),           # Line types
       col = c("black", "red"),           # Line colors
       lwd = c(2,3)) 

tpsp3_ma=arima(tpsp3, order = c(0, 0, 0),seasonal = list(order = c(2L, 2L, 2L)))
tpsp3_ma
residualsp3=residuals(tpsp3_ma)
tpsp3_fitted <- tpsp3 - residualsp3
ts.plot(tpsp3, ylab="Specificity to Hymenoptera")
points(tpsp3_fitted, type = "l", col = 2, lty = 2, lwd=3)
legend(x = "topright",          # Position
       legend = c("Original", "Fitted moving average"),  # Legend texts
       lty = c(1, 2),           # Line types
       col = c("black", "red"),           # Line colors
       lwd = c(2,3)) 
