
#### Functions ####

# single Marker functions 

singmarktab= function(markers, status, outdir=NULL){
	aucs=NULL
	for (k in 1:ncol(markers)){
		tmp=singmark(markers[,k], status, name=colnames(markers)[k],plot.save=outdir)
		aucs=c(aucs,tmp$auc)
	}
	cbind("Marker"=colnames(markers)[order(aucs,decreasing=T)],"AUC"=round(aucs[order(aucs,decreasing=T)],3))
}

singmark = function(marker, status, name="Marker", plots=TRUE, plot.save=NULL){
	probs=NULL
	for (j in 1:length(marker)){
		fit=glm(status[-j]~marker[-j],family=binomial)
		logOdds=fit$coef[1]+fit$coef[2]*marker[j]
  		prob=exp(logOdds)/(1+exp(logOdds))
  		probs=c(probs,prob)
	}
	pred = prediction(probs,status)
	perf = performance(pred,"tpr","fpr")
	auc = performance(pred,"auc")@'y.values'[[1]]
	
	if(plots){
		if(!is.null(plot.save)){pdf(paste(plot.save,name,'.pdf',sep=''),height=5,width=10)}
		par(mfrow=c(1,2))
		plot(perf)
		title(paste("ROC Curve:",name))
		text(.8,.2, paste("AUC: ",round(auc,3),sep=''))

		plot(rep(0,length(probs)),probs, type='n',xaxt='n',ylab="Predicted Probability",xlab="Disease Status",ylim=c(0,1),xlim=c(-.5,1.5))
		axis(1,0:1,levels(status))
		title(name)
		points(jitter(rep(0,sum(status==levels(status)[1])),amount=.2),1-probs[status==levels(status)[1]],col=2,pch=19)
		points(jitter(rep(1,sum(status==levels(status)[2])),amount=.2),1-probs[status==levels(status)[2]],col=4,pch=19)
		if(!is.null(plot.save)){dev.off()}
	}
	return(list(auc=auc,probs=probs))
}


# Multiple Marker functions 

markerfreq=function(modsaucs,marks=NULL){ 
	s= length(modsaucs$models)
	tab=table(unlist(modsaucs$models))
	tab/s
	tab[order(tab,decreasing=T)]/s

	for (i in 1:s){
	cat(modsaucs$AUCs[i],'\t',modsaucs$models[[i]],"\n")
	}
	tab1=matrix(cbind(tab[order(tab,decreasing=T)]/s,order(tab,decreasing=T)),ncol=2)
	if(!is.null(marks)){rownames(tab1)=marks[order(tab,decreasing=T)]}
	tab1
}

markerselect = function(markers, status, threshold=0.03){
	bestAUC=0
	AUCs=NULL
	models=list()
	for (s in 1:ncol(markers)){
 		var=s
 		probs=NULL
 		for (z in 1:nrow(markers)){
		    fit=glm(status[-z]~markers[-z,s],family=binomial)
		    logOdds=fit$coef[1]+fit$coef[2]%*%markers[z,s]
    		prob=exp(logOdds)/(1+exp(logOdds))
    		probs=c(probs,prob)
  		}
 		pred = prediction(probs, as.character(status))
 		oldAUC=newAUC = performance(pred,"auc")@'y.values'[[1]]

		while(newAUC>=oldAUC){
  			oldAUC=newAUC
  			add=NULL
  			for (v in (1:ncol(markers))[-var]){
	    		probs=NULL
    			for (i in 1:nrow(markers)){
	      			fit=glm(status[-i]~as.matrix(markers[-i,c(var,v)]),family=binomial)
	      			logOdds=fit$coef[1]+as.numeric(markers[i,c(var,v)])%*%fit$coef[-1]
    	  			prob=exp(logOdds)/(1+exp(logOdds))
      	  			probs=c(probs,prob)
    			}
    		pred = prediction(probs, status)
    		tmpAUC = performance(pred,"auc")@'y.values'[[1]]
    		if(tmpAUC>newAUC+threshold){add=v;newAUC=tmpAUC}
  			}
   			if(is.null(add)){
    			break
  			}else{
  				var=c(var,add)
  			}
 		}

	cat("Done adding variables for iteration",s,"; Current AUC:",newAUC,"; Variables:",var,"\n")
	if (newAUC>=bestAUC){
   	bestAUC=newAUC
   	bestvars=var
 	}
 	models[[s]]=var
 	AUCs=c(AUCs,newAUC)
	}
	return(list(models=models,AUCs=AUCs))
}
	
modelplots = function(models,mark_dat,status,outdir=NULL){	
 for (k in models){
  markers=unlist(strsplit(colnames(mark_dat)[k],"marker."))
  markers=markers[markers!=""]
  probs=NULL
  for (i in 1:nrow(mark_dat)){
    fit=glm(status[-i]~as.matrix(mark_dat[-i,k]),family=binomial)
    logOdds=fit$coef[1]+as.matrix(mark_dat[i,k])%*%fit$coef[-1]
    prob=exp(logOdds)/(1+exp(logOdds))
    probs=c(probs,prob)
  }
   pred = prediction(probs, as.character(status))
   perf = performance(pred,"tpr","fpr")
   plot(perf)

   auc = performance(pred,"auc")@'y.values'[[1]]
   cat("Variables:",k,"; peptides:",colnames(mark_dat)[k],"; AUC:",auc,'\n',sep=' ')

  pdf(paste(outdir,"markers_",paste(markers,collapse="_"),'.pdf',sep=''), height=5, width=10)
  par(mfrow=c(1,2))
  pred = prediction(probs, as.character(status))
  perf = performance(pred,"tpr","fpr")
  plot(perf)
  title(paste("ROC plot: Markers ",paste(markers,collapse=", "),sep='')) 
  text(.6,.2, paste("AUC:",round(auc,3)))

  plot(status, probs, type='n',xaxt='n',ylab="Predicted Probability",xlab="Disease Status",ylim=c(0,1),xlim=c(-1.5,.5))
  axis(1,-1:0,c('Control',"Alzheimer's"))
  title(paste("Markers ",paste(markers,collapse=", "),sep=''))
  points(jitter(rep(-1,nrow(mark_dat))[status=="Control"],amount=.2),1-probs[status=="Control"],col=2,pch=19)
  points(jitter(rep(0,nrow(mark_dat))[status=="Case"],amount=.2),1-probs[status=="Case"],col=4,pch=19)
  dev.off()

}
}

modelplot = function(model,mark_dat,status){	
 k=model
  markers=unlist(strsplit(colnames(mark_dat)[k],"marker."))
  markers=markers[markers!=""]
  probs=NULL
  for (i in 1:nrow(mark_dat)){
    fit=glm(status[-i]~as.matrix(mark_dat[-i,k]),family=binomial)
    logOdds=fit$coef[1]+as.matrix(mark_dat[i,k])%*%fit$coef[-1]
    prob=exp(logOdds)/(1+exp(logOdds))
    probs=c(probs,prob)
  }
   pred = prediction(probs, as.character(status))
   perf = performance(pred,"tpr","fpr")
   plot(perf)

   auc = performance(pred,"auc")@'y.values'[[1]]
   cat("Variables:",k,"; peptides:",colnames(mark_dat)[k],"; AUC:",auc,'\n',sep=' ')

  par(mfrow=c(1,2))
  pred = prediction(probs, as.character(status))
  perf = performance(pred,"tpr","fpr")
  plot(perf)
  title(paste("ROC plot: Markers ",paste(markers,collapse=", "),sep='')) 
  text(.6,.2, paste("AUC:",round(auc,3)))

  plot(status, probs, type='n',xaxt='n',ylab="Predicted Probability",xlab="Disease Status",ylim=c(0,1),xlim=c(-1.5,.5))
  axis(1,-1:0,c('Control',"Alzheimer's"))
  title(paste("Markers ",paste(markers,collapse=", "),sep=''))
  points(jitter(rep(-1,nrow(mark_dat))[status=="Control"],amount=.2),1-probs[status=="Control"],col=2,pch=19)
  points(jitter(rep(0,nrow(mark_dat))[status=="Case"],amount=.2),1-probs[status=="Case"],col=4,pch=19)
}


	