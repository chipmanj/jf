# 06/03/2016   v0.03 - Allow to change y-axis label
#                    - Thicken lines and circle
#                    - set las=1
#                    - no need to input data, just use y
# 06/06/2016   v0.05 - unfinished, genearlizing to any amount of quantiles
# Plot quartiles
#_______________



plotoe <- function(p,y,nquants=5,title,ylab.outcome,rd.digits,add.details=TRUE) {
 par(mar=c(7,4,4,2)+.1, las=1)

#  if(sum(is.na(data[,prob]))>0) {
#  warning(paste('removing',sum(is.na(data[,prob])),'records with missing observations.'),call.=FALSE)
#  data = data[!is.na(data[,prob]),]  
#  }

 rd <- function(obj,d=3) formatC(obj,format="f",digits=rd.digits)
 
 # Break data into quintiles
 # quint <- quantile(data[,p],probs=seq(.2,1,by=.2))
 quants <- quantile(p,probs=seq(0,1,length.out=(nquants + 1)))
 qs     <- cut2(p,cuts=quantile(p,seq(0,1,length.out=(nquants + 1))))
 
 # Proportion brca+ within each quintile and 95% margin of error
 pj  <- tapply(y, qs, mean)
 mej <- 1.96 * sqrt(pj*(1-pj)/table(qs))
                
 # o / e
 oe         <- tapply(y, qs, sum) / tapply(p, qs, sum)
 oe.overall <- sum(y) / sum(p)
 
 # Median of quintile
 xmedian <- tapply(p,qs, median)

 # O/E Curve
 plot(x=xmedian,y=pj,xlim=c(0,1),ylim=c(0,1),xlab='risk (at median of quintile)',ylab="",
      main=paste(title,"\n",'Proportion',ylab.outcome,'with 95% CI'),lwd=2)
 apply(as.matrix(cbind(x=xmedian,p=pj,me=mej)),1,
       function(CI) lines(x=rep(CI["x"],2),y=CI["p"] + c(-1,1)*CI["me"], 
                          col="blue", lwd=2))
 lines(x=c(0,1),y=c(0,1),col='grey')

 if(add.details==TRUE){
   mtext(at =  .05, side=1, line=4, text = 'median of quintile exp:',cex=.75,adj=1)
   mtext(at = .1, side=1, line=4, text = paste('q1: ', rd(xmedian[1],rd.digits)),cex=.75,adj=0)
   mtext(at = .3, side=1, line=4, text = paste('q2: ', rd(xmedian[2],rd.digits)),cex=.75,adj=0)
   mtext(at = .5, side=1, line=4, text = paste('q3: ', rd(xmedian[3],rd.digits)),cex=.75,adj=0)
   mtext(at = .7, side=1, line=4, text = paste('q4: ', rd(xmedian[4],rd.digits)),cex=.75,adj=0)
   mtext(at = .9, side=1, line=4, text = paste('q5: ', rd(xmedian[5],rd.digits)),cex=.75,adj=0)
   
   mtext(at =  .05, side=1, line=5, text = 'o/e per quintile:',cex=.75,adj=1)
   mtext(at = .1, side=1, line=5, text = paste('q1: ', rd(oe[1],rd.digits)),cex=.75,adj=0)
   mtext(at = .3, side=1, line=5, text = paste('q2: ', rd(oe[2],rd.digits)),cex=.75,adj=0)
   mtext(at = .5, side=1, line=5, text = paste('q3: ', rd(oe[3],rd.digits)),cex=.75,adj=0)
   mtext(at = .7, side=1, line=5, text = paste('q4: ', rd(oe[4],rd.digits)),cex=.75,adj=0)
   mtext(at = .9, side=1, line=5, text = paste('q5: ', rd(oe[5],rd.digits)),cex=.75,adj=0)
   
   mtext(at =  .05, side=1, line=6, text = 'o/e overall:',cex=.75,adj=1)
   mtext(at =  .1,  side=1, line=6, text = rd(oe.overall,rd.digits),cex=.75,adj=0)
 }

 par(mar=c(5,4,4,2)+.1)
}


#check=VVc[!is.na(VVc$bp1.01) & !is.na(VVc$brca.res),c('bp1.01','brca.res'),]
#plotfun(VVc[!is.na(VVc$bp1.01),],prob='bp1.01',y='brca.res',title='test')
#plotfun(VVc[!is.na(VVc$bpAIM),],prob='bpAIM',y='brca.res',title='test')
#plotfun(VVc[!is.na(VVc$bp2.07),],prob='bp2.07',y='brca.res',title='test')
