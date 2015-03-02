meanpower<-as.matrix(read.table(file="meanpower.txt"))[1,] 
erp_signal<-as.matrix(read.table(file="signal.txt"))[1,] 
rec_time<-as.matrix(read.table(file="time.txt"))[1,] 
library(zoo)
#----Functions To Be Used For The Test-----#
#Adopted for R from http://www.cs.bris.ac.uk/~rafal/phasereset/
#Generates noise at the spectrum of human EEG
noiseR <- function(frames,srate,meanpower){    
  sumsig<-20
  freq=0;     
  signal<-rep(0,frames)
  range = 1:frames 
  for (i in 1:sumsig) {
    freq <- freq + (4*runif(1))
    freqamp = meanpower[min (ceiling(freq), 125)] / meanpower[1]
    phase = runif(1)*2*pi
    signal<- signal+ sin (c(1:frames)/srate*2*pi*freq + phase) * freqamp      
  }
  return(signal)
} 

#Creates an EEG data set with only noise
#As a matrix, where each row represents a trial
#numoftrials: number of trials for the EEG
#size: length of the eeg trials (timepoints)
#k: k factor of which to scale the noise 
createnoiseeeg<-function(numoftrials,size,k){
  EEG<-matrix(NA,numoftrials,size)
  for(i in 1:numoftrials){
    EEG[i,]<-k*noiseR(size,2500,meanpower)
  }
  return (EEG)
}

#Creates an EEG data set : Signal + noise
#As a matrix, where each row represents a trial
#numoftrials: number of trials for the EEG
#signal: signal to be used at each trial
#size: length of the eeg trials (timepoints)
#k: k factor of which to scale the noise 
createnoiseeeg_withsignal<-function(numoftrials,signal,k){
  size<-length(signal)
  EEG<-matrix(NA,numoftrials,size)
  for(i in 1:numoftrials){
    EEG[i,]<-k*noiseR(size,2500,meanpower)+signal
  }
  return (EEG)
}

#Create an ERP - average across trials
#eeg: the EEG data, as a matrix 
#where the rows represent trials and the columns time points
erp<-function(eeg){
  return(colSums(eeg)/dim(eeg)[1])
}

#Peak to Peak measurment
#signal : the signal to apply p2p 
#start: where to start looking for the  maximum window
#max_end: where to stop looking for the  maximum
#end: where to stop looking for the minimum
#size: the size of the window for which the minimum/maximum is applied
peak2peak_stat<-function(signal,start,max_end,end,size){    
  win_avg_max<-rollapply(signal[start:max_end],size,mean)
  max_roll<-which.max(win_avg_max)
  max_latency<-(start+max_roll-1)+round(size/2) 
  max_avg_peak<- max(win_avg_max)
  win_avg_min<-rollapply(signal[max_latency:end],size,mean)
  min_avg_peak<- win_avg_min[which.min(win_avg_min)]
  return (max_avg_peak - min_avg_peak)
}
#----------Test ----------------------------------------------------#
#Goal: to compare the p-value geneate from bootstrap resampling and permutation under the null hypothesis (EEG data sets consisting only of noise)
#---------Set up test variables ------------------------------------#
trials<-50 # number of trials per condition
len<-2200  #number of time points per trial (based on the time file)
ns<-10 #scaling factor of the noise -- not important for this test
p2p_start<-801 #(which time==801) Selecting realistic time periods to search for p2p based on the time availabe here, in order to keep consistent for tests with signal 
p2p_max_end<-2001 #(which time==801)
p2p_end<-2200 #end of trial 
p2p_win<-200 #size of p2p internal window 200 tp=100 ms 
numberoftests<-10000 #How many test to perform
numberofreps<-1000 #How many resampling per test - should not set to less than 1000
#---------End of Variables Configuration----------------------------#
#-------------------------------------------------------------------#
#---------Initialise matrices to store results----------------------#
p_values_perm<-rep(NA,numberoftests)
p_values_z<-rep(NA,numberoftests)
trueobserved_p2p<-rep(NA,numberoftests)
#---------Start the Test ------------------------------------------#
pb <- txtProgressBar(min=0, max=numberoftests, initial=0, style=3)
for(k in 1:numberoftests){    
    fc<-createnoiseeeg(trials,len,ns) #first condition EEG set
    sc<-createnoiseeeg(trials,len,ns) #second condition EEG set
    
    fc_erp<-erp(fc) #first condition ERP
    sc_erp<-erp(sc) #second condition ERP
    alltogether<-rbind(fc,sc) #set of all trials

    fc_erp_p2p<-peak2peak_stat(fc_erp,p2p_start,p2p_max_end,p2p_end,p2p_win) #p2p for first condition on ERP
    sc_erp_p2p<-peak2peak_stat(sc_erp,p2p_start,p2p_max_end,p2p_end,p2p_win) #p2p for second condition on ERP
    #*****#
    st<-fc_erp_p2p-sc_erp_p2p
    trueobserved_p2p[k]<-st #true observed difference of p2p -- to be used for permutation
    #*****#
    boot_st<-array(NA,numberofreps) #array to store resampling results from bootstrap
    perm_st<-array(NA,numberofreps) #array to store resampling results from permutation

    for(i in 1:numberofreps){ # start resamplings
      
      bootsample1<-sample(trials,replace=T) 
      bootsample2<-sample(trials,replace=T)      
      fc_boot<-fc[bootsample1,] #Select randomly from available trials from the first condition eeg set to create bootstraped first condition
      sc_boot<-sc[bootsample2,] #Create second condition bootstraped EEG set      
      fc_boot_erp<-erp(fc_boot) #Bootstraped first condition ERP
      sc_boot_erp<-erp(sc_boot) #Bootstraped second condition ERP
 
      fc_boot_erp_p2p<-peak2peak_stat(fc_boot_erp,p2p_start,p2p_max_end,p2p_end,p2p_win) #
      sc_boot_erp_p2p<-peak2peak_stat(sc_boot_erp,p2p_start,p2p_max_end,p2p_end,p2p_win)
      boot_st[i]<-fc_boot_erp_p2p-sc_boot_erp_p2p #store bootstraped p2p difference
            
      #Randomly divide the available trials to first and second condition
      permsample<-sample(trials*2)      
      fc_perm_erp<-erp(alltogether[permsample[1:trials],]) # ERP from random trials to create first condition
      sc_perm_erp<-erp(alltogether[permsample[(trials+1):(2*trials)],]) #ERP from rest trials to create second condition
      fc_perm_erp_p2p<-peak2peak_stat(fc_perm_erp,p2p_start,p2p_max_end,p2p_end,p2p_win)
      sc_perm_erp_p2p<-peak2peak_stat(sc_perm_erp,p2p_start,p2p_max_end,p2p_end,p2p_win)
      perm_st[i]<-fc_perm_erp_p2p-sc_perm_erp_p2p #store p2p differece from permutation resampling
      
    }
        
    z<-1-length(boot_st[boot_st>0])/numberofreps #calculate boostrap p-value
    p_values_z[k]<-z #store bootstrap pvalue
    
    p_value_perm<-length(perm_st[perm_st>=st])/numberofreps
    p_values_perm[k]<-p_value_perm #store permutation p-value
    setTxtProgressBar(pb, k)
    
}
  
#create data.frame with the result 
results<-data.frame(p_values_perm=p_values_perm,p_values_z=p_values_z)
#save results as R data
save(results, file = paste('permvboot',numberoftests,'_nr',numberofreps,'_trials',trials ,'_p2p',p2p_win,'noise_',ns,'.Rda', sep=''))
#save results as csv file
write.csv(results,file = paste('permvboot',numberoftests,'_nr',numberofreps,'_trials',trials ,'p2p',p2p_win,'noise_',ns,'.csv', sep=''))

#--------Visualise Results
pdf("results.pdf") 
par(xpd=TRUE)
hist(p_values_z,main=paste0("Histogram of p-values \n",numberoftests,", tests ",numberofreps," resamplings\n pure noise"),col=rgb(1,0,0,1/4),xlab='p-values')
hist(p_values_perm,col=rgb(0,0,1,1/4),add=TRUE)
legend(-0.01,-15,c("boot","perm"), lty = c(1,1),lwd=2,col=c("red","blue"),bty ="n")
dev.off()