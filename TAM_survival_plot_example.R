#Function for round off p value
p.text<-function(p) ifelse(p>0.99,'>0.99',
                           ifelse(p<0.001,'p<0.001',paste0('p=',round(p,3))))

#Function for generate survival plot
KMplot<-function(datin,time,event,group,titlein='',
                 palette.value=c( "#1B1919FF","#925E9FFF"),
                 break.time.value=24,legend.pos='right',legend.title='',legend.labs=NULL,
                 x.lab="Months",y.lab,ci.TF=FALSE,pval.xy=c(5,0.1),
                 HR.x=5, HR.y=0.3,p.value=TRUE,add.HR=TRUE){
  max.os<-max(datin[,time],na.rm=TRUE)
  #KM plot
  fit<-eval(parse(text=paste0("survfit(Surv(",time,",", event,")~",group,", data=datin,type='kaplan-meier',conf.type='log')")))
  fit.sum<-summary(fit,seq(0,max.os,by=break.time.value))
  #rho = 0 this is the log-rank or Mantel-Haenszel test, and with rho = 1 it is equivalent to the Peto & Peto modification of the Gehan-Wilcoxon test.
  lr.test<-eval(parse(text=paste0("survdiff(Surv(",time,",", event,")~",group,", data=datin,rho=0)")))
  
  lr.p<-pchisq(lr.test$chisq,df=length(lr.test[[1]])-1,lower.tail=FALSE)
  # lr.pt<-paste0('Log-rank test: ',p.text(lr.p))
  lr.pt<-p.text(lr.p)
  if (p.value==FALSE) {
    lr.pt=FALSE 
  } else {
    p.logrank.text<-paste0(lr.pt,' (Log-rank)')
  }
  cat('----Summary statistics for survival--\n')
  print(fit)
  #print(fit.sum)
  cat('\n\n')
  cat('----Log-rank test--\n')
  print(lr.test) 
  
  if (add.HR==TRUE){
    #Add HR------------------
    ##coxph function
    f.coxph<-eval(parse(text=paste0('coxph(Surv(',time,',',event,')~',group,',data=datin)')))
    #print(f.coxph)
    sum.coxph<-summary(f.coxph)
    HR.CI<-sum.coxph$conf.int
   
    HR.CI<-HR.CI[,c('exp(coef)','lower .95', 'upper .95')]
    HR.text<-paste0('HR (95% CI)\n',round(HR.CI[1],3),' (',round(HR.CI[2],3),', ',round(HR.CI[3],3),')')
    add.HR.text<-''
  
  }
  
  #survival plot
  ggsurv <- ggsurvplot(
    fit,
    data = datin, 
    size = 0.6,  # change line size
    palette = palette.value,
    conf.int = FALSE,                        # Add confidence interval
    pval=FALSE,
    surv.plot.height = 0.60,                 # the height of the survival plot on the grid
    censor.shape = "+",
    censor.size=2,
    legend =legend.pos,  # "top", "bottom", "left", "right", "none"
    font.legend=10,
    legend.title = legend.title,
    legend.labs = legend.labs,  # Change legend labels
    ylab=y.lab,
    xlab =x.lab,                   # customize X axis label
    break.time.by=break.time.value,
    risk.table = TRUE,                       # Add risk table
    risk.table.col = "strata",               # Risk table color by groups
    risk.table.height = 0.25,                # Useful to change when you have multiple groups
    risk.table.y.text.col = TRUE,            # colour risk table text annotations
    risk.table.y.text = TRUE,               # show bars instead of names in text annotations
    risk.table.fontsize = 3,
    tables.theme = theme_cleantable(plot.title = element_text(hjust = -0.5,size=5))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            #legend.text=element_text(size=8),
            plot.title=element_text(size=10)),  #text size for number at risk
    
    title=titlein
    
  )
  
  if (p.value==TRUE) {
    ggsurv$plot <- ggsurv$plot +
      annotate("text", x = pval.xy[1], y = pval.xy[2], label =p.logrank.text ,  hjust = 0,vjust=1,size=3.5)
  }
  if (add.HR==TRUE){
    ggsurv$plot <- ggsurv$plot +
      theme(axis.text.x=element_text(size=9),  #to control axis text on plot
            axis.text.y=element_text(size=9))+
      annotate("text", x = HR.x, y = HR.y, label = HR.text,  hjust = 0,vjust=1,size=3.5)+
      scale_y_continuous(breaks = seq(0,1,by=0.25), labels = seq(0,100,by=25))
  } else {
    ggsurv$plot <- ggsurv$plot +
      theme(axis.text.x=element_text(size=9),  #to control axis text on plot
            axis.text.y=element_text(size=9))+
      scale_y_continuous(breaks = seq(0,1,by=0.25), labels = seq(0,100,by=25))
  }
  
  print(ggsurv)
}

#############
#Example

library(survival)
library(ggplot2)
library(ggsurvfit)
library(survminer)
library(rms)
data(df_lung)
#data set name
#datin: your data include time, event, group variables
datin<-df_lung
#variable name for a binary group
group<-'sex'
#variable name for time to event or lost follow-up
time<-'time'  
# variable name for event
event<-'status'
#
y.lab<-"Survival (%)"

ddist<-datadist(datin)
options(datadist="ddist")
KMplot(datin=datin,time=time,event=event,group=group,
       break.time.value=3,
       legend.pos='right',  #range from 0 to 1
       legend.labs=levels(datin[,group]),
       x.lab= "Months",y.lab=y.lab,
       HR.x=1, HR.y=0.2, #position for HR
       ci.TF=FALSE,
       pval.xy=c(1,0.1), #position for p value
       p.value=TRUE 
      # add.HR=add.HR
       )
