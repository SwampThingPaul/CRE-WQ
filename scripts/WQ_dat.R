## Interactive water quality tool Demo
## Base code
## 
##
## Code was compiled by Paul Julian
## contact infor: pjulian@ufl.edu

#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(AnalystHelper);#devtools::install_github("SwampThingPaul/AnalystHelper")
library(plyr)
library(reshape)
library(Hmisc)
library(zoo)
library(RcppRoll)


#Spatial
library(rgeos)
library(rgdal)
library(leaflet)

library(tmap)
#reminder https://rstudio.github.io/leaflet/shapes.html

#Paths
setwd("D:/_GitHub/WQ_Demo")

paths=paste0(getwd(),c("/Exports/","/Plots/","/Data/","/GIS"))
#Folder.Maker(paths);#One and done. Creates folders in working directory.
export.path=paths[1]
plot.path=paths[2]
data.path=paths[3]
GIS.path=paths[4]


# Spatial Data ------------------------------------------------------------
utm17=CRS("+proj=utm +zone=17 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

counties=spTransform(readOGR(GIS.path,"Florida_Counties"),utm17)
wbids=spTransform(readOGR(GIS.path,"FDEP_WBID_LeeCo"),utm17)
#monitoring=spTransform(readOGR(GIS.path,"Environmental_Monitoring_Stations"),utm17)

lee.co=subset(counties,COUNTYNAME=="LEE")
wbids.cal=subset(wbids,WBID%in%c("3240C","3240B","3240A"))
#monitoring.cal=subset(monitoring,STATION%in%c("S79",paste0("CES0",2:8))&ACTIVITY_T=="Chemistry"&ACTIVITY_S=="Surface Water Grab")
#wq.mon=subset(monitoring,ACTIVITY_T=="Chemistry"&ACTIVITY_S=="Surface Water Grab")
#writeOGR(monitoring.cal, GIS.path, "Caloosa_wq", driver="ESRI Shapefile")
#writeOGR(wq.mon, GIS.path, "WMD_wq", driver="ESRI Shapefile")

monitoring.cal=spTransform(readOGR(GIS.path,"Caloosa_wq"),utm17)

tmap_mode("view")
tm_basemap(leaflet::providers$Esri.WorldImagery,alpha=0.9)+
  tm_shape(wbids)+tm_borders("grey")+
  tm_shape(wbids.cal)+tm_polygons(id="WBID",alpha=0.5,border.col="red")+
  tm_shape(monitoring.cal)+tm_dots(id="STATION",col="yellow",size=0.25)+
  tm_shape(lee.co)+tm_borders("black",lty=2,lwd=2)
  
  

# Estuary Water Quality ---------------------------------------------------
#https://www.leegov.com/naturalresources/Documents/Water%20Quality%20Status/TMDL%20update%20to%20BoCC%20Dec52017%20.pdf
#https://floridadep.gov/sites/default/files/tidal-caloosa-nutr-tmdl_0.pdf
estuary.site.list=data.frame(Station.ID=c("CES02","CES03","CAL 01","CES04","CAL 03","CES05","CAL 05","CES06","CAL 07","CES07","CAL 09","CES08","CAL 17","S79"),
                     ALIAS=c("CES02",rep("CES03",2),rep("CES04",2),rep("CES05",2),rep("CES06",2),rep("CES07",2),rep("CES08",2),"S79"),
                     WBID=c(rep("3240C",3),rep("3240B",2),rep("3240C",8),NA))
estuary.site.list$CES=ifelse(substring(estuary.site.list$Station.ID,1,3)=="CES"|estuary.site.list$Station.ID=="S79",1,0)
estuary.site.list=subset(estuary.site.list,CES==1)
wq.params=data.frame(Test.Number=c(61,179,21,18,197,25,98,9,8,7,80),param=c("Chla","Chla","TKN","NOx","K_PAR","TP","Sal","SPC","DO","Temp","TN"))
wq.params=subset(wq.params,param%in%c("TP","TN","TKN","NOx","Sal","SPC","DO","Temp","Chla"))

dates=as.Date(c("2010-05-01","2018-04-30"))
est.dat=data.frame()
for(i in 1:nrow(estuary.site.list)){
tmp=DBHYDRO_WQ(dates[1],dates[2],estuary.site.list$Station.ID[i],wq.params$Test.Number)
est.dat=rbind(est.dat,tmp)
print(i)
}
est.dat=subset(est.dat,Collection.Method=="G")



S79.flow=SFWMD.DBHYDRO.Data.daily(dates[1],dates[2],"00865")
S79.flow$Date.EST=date.fun(S79.flow$Date)
S79.flow$Q.cfs.S79=S79.flow$Data.Value

est.dat$WY=WY(est.dat$Date.EST)
est.dat=merge(est.dat,wq.params,"Test.Number")
est.dat=merge(est.dat,estuary.site.list,"Station.ID")

est.dat.xtab=cast(est.dat,WBID+Station.ID+DateTime.EST+WY~param,value="HalfMDL",mean)
est.dat.xtab$Date.EST=date.fun(est.dat.xtab$DateTime.EST)
est.dat.xtab$TP=est.dat.xtab$TP*1000
est.dat.xtab$TN=with(est.dat.xtab,TN_Combine(NOx,TKN,TN))
est.dat.xtab$Sal.calc=with(est.dat.xtab,round(SalinityCalc(SPC,Temp),2))
est.dat.xtab$DO.persat=with(est.dat.xtab,DO_PerSat(Temp,DO,Sal.calc))
#Predominately feshwater is SPC<4580 (see SFER)
est.dat.xtab$DO.WQS.TOD=with(est.dat.xtab,ifelse(SPC<4580&Station.ID%in%c("S79","CES02"),DO.TOD.WQS.stream(est.dat.xtab$DateTime.EST),42.0));
est.dat.xtab$month=as.numeric(format(est.dat.xtab$Date.EST,"%m"))
est.dat.xtab$season=FL.Hydroseason(est.dat.xtab$Date.EST)
est.dat.xtab=merge(est.dat.xtab,S79.flow[,c("Date.EST","Q.cfs.S79")],all.x=T)

N.scn=cast(est.dat.xtab,Station.ID+WY~season,value="TP",fun.aggregate = function(x) N(x))
N.scn$NTotal=with(N.scn,A_Wet+B_Dry);
N.scn$SeasonScreen=with(N.scn,ifelse(B_Dry>=1&A_Wet>=1,1,0));
N.scn$NScreen=with(N.scn,ifelse(NTotal>=4,1,0));
N.scn$UseData=with(N.scn,ifelse((SeasonScreen+NScreen)==2,"Yes","No"));

est.dat.xtab=merge(est.dat.xtab,N.scn[,c("Station.ID","WY","UseData")],c("Station.ID","WY"))

WY.GM.dat=ddply(est.dat.xtab,c("Station.ID","WY","UseData"),summarise,TP.GM=exp(mean(log(TP),na.rm=T)),TN.GM=exp(mean(log(TN),na.rm=T)),Chla.GM=exp(mean(log(Chla),na.rm=T)),TP.FWM=wtd.mean(TP,Q.cfs.S79,na.rm=T),TN.FWM=wtd.mean(TN,Q.cfs.S79,na.rm=T),mean.DO=mean(DO.persat,na.rm=T),mean.Sal=mean(Sal.calc,na.rm=T),DO.ExceedPer=(sum(DO.persat<DO.WQS.TOD,na.rm=T)/N(DO.persat))*100)
#WY.GM.dat$TP.Cusum=with(WY.GM.dat,ave(TP.GM,Station.ID,FUN=function(test) round(cumsum(base::scale(test)),1)))
WY.GM.dat$LT.AVG.TP=with(WY.GM.dat,ave(TP.GM,Station.ID,FUN=function(x) roll_meanr(x,n=5)))
WY.GM.dat$LT.AVG.TN=with(WY.GM.dat,ave(TN.GM,Station.ID,FUN=function(x) roll_meanr(x,n=5)))
WY.GM.dat$LT.AVG.Chla=with(WY.GM.dat,ave(Chla.GM,Station.ID,FUN=function(x) roll_meanr(x,n=5)))
WY.GM.dat

subset(WY.GM.dat,Station.ID==site.val)

est.dat.xtab=est.dat.xtab[order(est.dat.xtab$Station.ID,est.dat.xtab$Date.EST),]

#https://ca.dep.state.fl.us/mapdirect/?focus=nutrientcriteria
data.frame(WBID=c("3240C","3240B","3240A"),TP=c(86,55,40),TN=c(NA,NA,NA),Chla=c(4.2,6.5,5.6),Name=c("Upper Caloosa","Middle Caloosa","Lower Caloosa"))

par(family="serif",mar=c(2,2,1,2),oma=c(1,2.5,1,0.1));
layout(matrix(1:10,2,5,byrow=F));
i=1
for(i in 1:8){
  site.val=as.character(estuary.site.list$Station.ID[i])
  ylim.val=c(10,1000);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
  xlim.val=date.fun(c("2010-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"3 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
  plot(TP~Date.EST,est.dat.xtab,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA,log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(est.dat.xtab,Station.ID==site.val),lines(Date.EST,TP,lty=2,col="grey",lwd=2))
  with(subset(est.dat.xtab,Station.ID==site.val&Q.cfs.S79>0),points(Date.EST,TP,pch=21,bg="dodgerblue1"))
  with(subset(est.dat.xtab,Station.ID==site.val&Q.cfs.S79<=0),points(Date.EST,TP,pch=22,bg="indianred1"))
  axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),cex=0.9);  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=2,line=2.5,"TP (\u03BCg L\u207B\u00B9)")
  mtext(side=1,line=2,"Date (Month-Year)")
  legend("topleft",legend=c("Q > 0 at S-79","Q = 0 at S-79"),pch=c(21,22),pt.bg=c("dodgerblue1","indianred1"),
       pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)

  ylim.val=c(40,200);ymaj=c(50,log.scale.fun(ylim.val,"major"));ymin=log.scale.fun(ylim.val,"minor")
  xlim.val=c(2011,2018);by.x=3;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
  plot(TP.GM~WY,WY.GM.dat,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA,log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(WY.GM.dat,Station.ID==site.val),lines(WY,LT.AVG.TP,col="red",lwd=2))
  with(subset(WY.GM.dat,Station.ID==site.val),pt_line(WY,TP.GM,2,"grey",2,21,"olivedrab1",cex=1.25))
  #with(subset(WY.GM.dat,Station.ID==site.val),pt_line(WY,TP.FWM,2,"grey",2,23,"skyblue",cex=1.25))
    axis_fun(1,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=2,line=2.5,"GM TP (\u03BCg L\u207B\u00B9)")
  mtext(side=1,line=2,"Water Year")
  legend("topleft",legend=c("Annual GM","5-WY Mean"),pch=c(21,NA),lty=c(2,1),col=c("grey","red"),lwd=c(1,2),pt.bg=c("olivedrab1",NA),pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)

  ylim.val=c(0,5);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);#;ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
  xlim.val=date.fun(c("2010-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"3 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
  plot(TN~Date.EST,est.dat.xtab,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(est.dat.xtab,Station.ID==site.val),lines(Date.EST,TN,lty=2,col="grey",lwd=2))
  with(subset(est.dat.xtab,Station.ID==site.val&Q.cfs.S79>0),points(Date.EST,TN,pch=21,bg="dodgerblue1"))
  with(subset(est.dat.xtab,Station.ID==site.val&Q.cfs.S79<=0),points(Date.EST,TN,pch=22,bg="indianred1"))
  axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),cex=0.9);  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=2,line=2,"TN (mg L\u207B\u00B9)")
  mtext(side=1,line=2,"Date (Month-Year)")
  #legend("topleft",legend=c("Q > 0 at S-79","Q = 0 at S-79"),pch=c(21,22),pt.bg=c("dodgerblue1","indianred1"),
  #       pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
  
  ylim.val=c(0,5);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);#;ymaj=c(50,log.scale.fun(ylim.val,"major"));ymin=log.scale.fun(ylim.val,"minor")
  xlim.val=c(2011,2018);by.x=3;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
  plot(TP.GM~WY,WY.GM.dat,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(WY.GM.dat,Station.ID==site.val),lines(WY,LT.AVG.TN,col="red",lwd=2))
  with(subset(WY.GM.dat,Station.ID==site.val),pt_line(WY,TN.GM,2,"grey",2,21,"olivedrab1",cex=1.25))
  #with(subset(WY.GM.dat,Station.ID==site.val),pt_line(WY,TN.FWM,2,"grey",2,23,"skyblue",cex=1.25))
  axis_fun(1,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=2,line=2,"GM TN (mg L\u207B\u00B9)")
  mtext(side=1,line=2,"Water Year")

  
  ylim.val=c(0,100);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);#;ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
  xlim.val=date.fun(c("2010-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"3 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
  plot(Chla~Date.EST,est.dat.xtab,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(est.dat.xtab,Station.ID==site.val),lines(Date.EST,Chla,lty=2,col="grey",lwd=2))
  with(subset(est.dat.xtab,Station.ID==site.val&Q.cfs.S79>0),points(Date.EST,Chla,pch=21,bg="dodgerblue1"))
  with(subset(est.dat.xtab,Station.ID==site.val&Q.cfs.S79<=0),points(Date.EST,Chla,pch=22,bg="indianred1"))
  axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),cex=0.9);  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=2,line=2.5,"Chlorophyll-A (\u03BCg L\u207B\u00B9)")
  mtext(side=1,line=2,"Date (Month-Year)")
  
  ylim.val=c(0,20);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);#;ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
  xlim.val=c(2011,2018);by.x=3;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
  plot(Chla.GM~WY,WY.GM.dat,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(WY.GM.dat,Station.ID==site.val),lines(WY,LT.AVG.Chla,col="red",lwd=2))
  with(subset(WY.GM.dat,Station.ID==site.val),pt_line(WY,Chla.GM,2,"grey",2,21,"olivedrab1",cex=1.25))
  #with(subset(WY.GM.dat,Station.ID==site.val),pt_line(WY,TP.FWM,2,"grey",2,23,"skyblue",cex=1.25))
  axis_fun(1,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=2,line=2.5,"GM Chl-A (\u03BCg L\u207B\u00B9)")
  mtext(side=1,line=2,"Water Year")
  
  
  ylim.val=c(0,200);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);#;ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
  xlim.val=date.fun(c("2010-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"3 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
  plot(DO.persat~Date.EST,est.dat.xtab,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(est.dat.xtab,Station.ID==site.val),lines(Date.EST,DO.persat,lty=2,col="grey",lwd=2))
  with(subset(est.dat.xtab,Station.ID==site.val&Q.cfs.S79>0),points(Date.EST,DO.persat,pch=21,bg="dodgerblue1"))
  with(subset(est.dat.xtab,Station.ID==site.val&Q.cfs.S79<=0),points(Date.EST,DO.persat,pch=22,bg="indianred1"))
  with(subset(est.dat.xtab,Station.ID==site.val),lines(Date.EST,DO.WQS.TOD,col="red",lty=1,lwd=2))
  axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),cex=0.9);  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=2,line=2.5,"DO (% Sat.)")
  mtext(side=1,line=2,"Date (Month-Year)")
  legend("topleft",legend=c("WQS"),lty=1,lwd=2,col="red",
         ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
  
  ylim.val=c(0,200);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);#;ymaj=c(50,log.scale.fun(ylim.val,"major"));ymin=log.scale.fun(ylim.val,"minor")
  xlim.val=c(2011,2018);by.x=3;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
  plot(mean.DO~WY,WY.GM.dat,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(WY.GM.dat,Station.ID==site.val),pt_line(WY,mean.DO,2,"grey",2,21,ifelse(DO.ExceedPer>10,"red","olivedrab1"),cex=ifelse(DO.ExceedPer>10,2,1.25)))
  #with(subset(WY.GM.dat,Station.ID==site.val),pt_line(WY,TN.FWM,2,"grey",2,23,"skyblue",cex=1.25))
  axis_fun(1,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=2,line=2.5,"Mean DO (% Sat.)")
  mtext(side=1,line=2,"Water Year")
  legend("topleft",legend=c("Exceed WQS", "Achieve WQS"),pch=c(21,21),pt.bg=c("red","olivedrab1"),
         pt.cex=c(2,1.25),ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)

  
  ylim.val=c(0,50);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);#;ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
  xlim.val=date.fun(c("2010-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"3 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
  plot(Sal.calc~Date.EST,est.dat.xtab,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(est.dat.xtab,Station.ID==site.val),lines(Date.EST,Sal.calc,lty=2,col="grey",lwd=2))
  with(subset(est.dat.xtab,Station.ID==site.val&Q.cfs.S79>0),points(Date.EST,Sal.calc,pch=21,bg="dodgerblue1"))
  with(subset(est.dat.xtab,Station.ID==site.val&Q.cfs.S79<=0),points(Date.EST,Sal.calc,pch=22,bg="indianred1"))
  axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),cex=0.9);  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=2,line=2,"Salinity")
  mtext(side=1,line=2,"Date (Month-Year)")
  #legend("topleft",legend=c("Q > 0 at S-79","Q = 0 at S-79"),pch=c(21,22),pt.bg=c("dodgerblue1","indianred1"),
  #       pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
  
  ylim.val=c(0,40);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);#;ymaj=c(50,log.scale.fun(ylim.val,"major"));ymin=log.scale.fun(ylim.val,"minor")
  xlim.val=c(2011,2018);by.x=3;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
  plot(mean.Sal~WY,WY.GM.dat,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(WY.GM.dat,Station.ID==site.val),pt_line(WY,mean.Sal,2,"grey",2,21,"olivedrab1",cex=1.25))
  axis_fun(1,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=2,line=2,"Mean Salinity")
  mtext(side=1,line=2,"Water Year")
  mtext(side=3,outer=T,site.val[i])
}