setwd("F:/R.bao")
rm(list=ls())
######Package loading######
library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)
library(KEGGgraph)
library(openxlsx)
######Data Preparation#######
pedges<-read.csv("./data/pedges.csv",row.names = 1) #47020/2 Edge data start to end
#GSE30550
dsubjects <- c("Subject 01","Subject 05","Subject 06","Subject 07","Subject 08","Subject 10","Subject 12","Subject 13","Subject 15")#º≤≤° 8 13»± ß
nsubjects <- c("Subject 02","Subject 03","Subject 04","Subject 09","Subject 11","Subject 14","Subject 16","Subject 17")#’˝≥£ 17»± ß
subjects <- c(dsubjects,nsubjects)
hours<-c("Baseline","Hour 00","Hour 005","Hour 012",
         "Hour 021","Hour 029","Hour 036","Hour 045",
         "Hour 053","Hour 060","Hour 069","Hour 077",
         "Hour 084","Hour 093","Hour 101","Hour 108")
allgvs<-read.csv("./data/allgvs.csv",row.names = 1) #3867/268 exp
tts<-read.csv("./data/tts_30550.csv",row.names = 1) #268/3  sample_id
#GSE52428
dsubjects_H1N1<-c("H1N1_002","H1N1_003","H1N1_006","H1N1_007","H1N1_008","H1N1_009","H1N1_010","H1N1_012","H1N1_013","H1N1_017","H1N1_020","H1N1_021") #≤°Ã¨12
nsubjects_H1N1<-c("H1N1_001","H1N1_004","H1N1_005","H1N1_011","H1N1_014","H1N1_015","H1N1_016","H1N1_018","H1N1_019","H1N1_022","H1N1_023","H1N1_024") #Œﬁ÷¢◊¥12
subjects_H1N1<-c(dsubjects_H1N1,nsubjects_H1N1)
dsubjects_H3N2<-c("H3N2_001","H3N2_005","H3N2_006","H3N2_007","H3N2_008","H3N2_010","H3N2_012","H3N2_013","H3N2_015") #≤°Ã¨9
nsubjects_H3N2<-c("H3N2_002","H3N2_003","H3N2_004","H3N2_009","H3N2_011","H3N2_014","H3N2_016","H3N2_017") #Œﬁ÷¢◊¥ 8
subjects_H3N2<-c(dsubjects_H3N2,nsubjects_H3N2)
hours_52428<-c("Baseline","Hour 00","Hour 005","Hour 012",
         "Hour 021","Hour 029","Hour 036","Hour 045",
         "Hour 053","Hour 060","Hour 069","Hour 077",
         "Hour 084","Hour 093","Hour 101","Hour 108")
allgvs_52428<-read.csv("./data/allgvs_52428.csv",row.names = 1) #3905/649 exp
allgvs_52428<-log2(allgvs_52428+1) #normalization
tts_52428<-read.csv("./data/tts_52428.csv",row.names = 1) #649/3 sample_id
####function####
regulation_strength<-function(ii,subjects,tts,allgvs,pedges){
  tt1 <- tts[which(tts[,1] == subjects[[ii]]),]
  samsone <- tt1[,3]
  same <- allgvs[,samsone]
  evalues <- NULL
  for(j in 1:length(samsone)){
    tgvs <- same[,j] #sample tt1's exp at time j 
    names(tgvs) <- rownames(same)
    fv <- tgvs[pedges[,1]] #upstream gene
    tv <- tgvs[pedges[,2]] #downstream gene
    fv[is.na(fv)] <- 0    #nan
    tv[is.na(tv)] <- 0
    evalue <- as.numeric(tv)-as.numeric(fv) #regulation strength
    evalues <- cbind(evalues,evalue) 
  }
  return(evalues)
}
diff_regulation_strength<-function(evalues){
  refev <- apply(evalues[,1:4],1,mean) #47020 
  #refev<-evalues[,1]
  edgesv <- abs(evalues-refev) #the difference of regulation strength between T and ref
  return(edgesv)
}
diff_regulation_strength_order<-function(edgesv){
  diff_score_order<-NULL #order
  for(i in 1:dim(edgesv)[2]){
    diff_score_order<-cbind(diff_score_order,order(-edgesv[,i]))
  }
  return(diff_score_order)
}
diff_regulation_strength_topn<-function(topn,edgesv,diff_score_order){
  edgesv_topn<-NULL
  for(i in 1:dim(edgesv)[2]){
    d<-diff_score_order[1:topn,i]
    edgesv_topn<-cbind(edgesv_topn,edgesv[d,i])
  }
  return(edgesv_topn)
}
get_score<-function(ii,subjects,tts,edgesv_topn,hours){
  tt1 <- tts[which(tts[,1] == subjects[[ii]]),]
  score <- apply(edgesv_topn,2,mean)  #50 16 missed
  score<-as.data.frame(t(score))
  colnames(score)<-tt1[,2]
  score_t <- as.data.frame(matrix(nrow = 0,ncol = length(hours)))
  colnames(score_t)<-hours
  score<-merge(score_t,score,all=T) 
  score<-score[match(hours,colnames(score))]
  w<-which(is.na(match(hours,tt1[,2]))) #Fill in missing values
  if(length(w)==0){}
  else{
    for(h in 1:length(w)){
      if(w[h]==1){score[w[h]]<-0}
      else{score[w[h]]<-score[w[h]-1]}
    }
  }
  return(score)
}

dscore<-function(ii,subjects,tts,allgvs,pedges,topn,hours){
  
  evalues<-regulation_strength(ii,subjects,tts,allgvs,pedges)
  edgesv<-diff_regulation_strength(evalues) #the difference
  diff_score_order<-diff_regulation_strength_order(edgesv) #gene's
  edgesv_topn<-diff_regulation_strength_topn(topn,edgesv,diff_score_order)
  score<-get_score(ii,subjects,tts,edgesv_topn,hours)
  
  return(score)
}
dscore_add<-function(score){
  score_num<-as.numeric(score)
  score_add<-NULL
  for(hh in 1:length(score_num)){
    score_add<-c(score_add,abs(score_num[hh]-mean(score_num[1:hh])) )
  }
  return(score_add)
}

####K=50#####test GSE30550####
scores <- NULL
scores_add<-NULL
for(ii in 1:length(subjects)){
  score<-dscore(ii,subjects,tts,allgvs,pedges,topn = 50,hours)
  score_add<-dscore_add(score) #add  
  scores <- rbind(scores,score)
  scores_add<-rbind(scores_add,score_add)
}
rownames(scores) <- subjects
colnames(scores_add)<-hours
rownames(scores_add)<-subjects
#write.csv(scores,"allscroes_50.csv")
write.csv(scores_add,"allscroes_add_50.csv")

#### K ######Process of experiment####
topn_list<-c(seq(1,100,by=1))
topn_list<-c(seq(1,9,by=1),seq(10,99,by=10),seq(100,999,by=100),seq(1000,9999,by=1000),seq(10000,dim(pedges)[1],by=10000))#dim(pedges)[1]

threshold_list<-NULL
for(m in 1:length(subjects)){ #sample m
  tt1 <- tts[which(tts[,1] == subjects[[m]]),]
  evalues<-regulation_strength(m,subjects,tts,allgvs,pedges)
  edgesv<-diff_regulation_strength(evalues) 
  diff_score_order<-diff_regulation_strength_order(edgesv) 
threshold<-NULL
  for(t in 1:length(topn_list)){
    topn<-topn_list[t]
    edgesv_topn<-diff_regulation_strength_topn(topn,edgesv,diff_score_order)
    score<-get_score(m,subjects,tts,edgesv_topn,hours)
    score_add<-dscore_add(score) 
    threshold<-c(threshold,score_add[which.max(score_add)])
  }
threshold_list<-rbind(threshold_list,threshold)
}
rownames(threshold_list)<-subjects
colnames(threshold_list)<-topn_list
write.csv(t(threshold_list),"K_list.csv")

####Dysregulated DNBs(K=10)####
#This process was repeated for cross-validation with the replacement dataset
#H1N1
#subjects<-subjects_H1N1
#tts<-tts_52428
#allgvs<-allgvs_52428
#hours<-hours_52428
#dsubjects<-dsubjects_H1N1
#nsubjects<-nsubjects_H1N1
topn<-10
scores <- NULL
scores_add<-NULL
KDNs<-NULL
for(ii in 1:length(dsubjects)){
  evalues<-regulation_strength(ii,subjects,tts,allgvs,pedges)
  edgesv<-diff_regulation_strength(evalues) 
  diff_score_order<-diff_regulation_strength_order(edgesv) 
  edgesv_topn<-diff_regulation_strength_topn(topn,edgesv,diff_score_order)
  score<-get_score(ii,subjects,tts,edgesv_topn,hours)
  score_add<-dscore_add(score) 
  timei<-which.max(score_add) 
  KDN<-diff_score_order[1:topn,timei]
  KDNs<-cbind(KDNs,KDN) 
}
colnames(KDNs)<-dsubjects
KDNs_cut<-as.numeric(names(which(table(KDNs)>1))) 
union(pedges[KDNs_cut,]$fromgene,pedges[KDNs_cut,]$togene)
library("clusterProfiler")
geneIDs<-union(pedges[KDNs_cut,]$fromgene,pedges[KDNs_cut,]$togene)
gene.df <- bitr(geneIDs, fromType = "ENTREZID", 
                toType = c("ENTREZID", "SYMBOL"), 
                OrgDb = org.Hs.eg.db)

print(gene.df)#Dysregulated_DNBs
write.csv(KDNs_cut,"Dysregulated_DNBs_10.csv")
####test####
KDNs_cut<-read.csv("Dysregulated_DNBs_10.csv",row.names = 1)
KDNs_edge<-KDNs_cut$x 

scores <- NULL
scores_add<-NULL
for(ii in 1:length(subjects)){
  tt1 <- tts[which(tts[,1] == subjects[[ii]]),]
  evalues<-regulation_strength(ii,subjects,tts,allgvs,pedges)
  edgesv<-diff_regulation_strength(evalues) #µ˜øÿ≤Ó“Ï
  edgesv_topn<-edgesv[KDNs_edge,]
  score<-get_score(ii,subjects,tts,edgesv_topn,hours)
  score_add<-dscore_add(score) #‘ˆ¡ø ’“µ„ #print(which.max(score_add)) >„–÷µ
  scores <- rbind(scores,score)
  scores_add<-rbind(scores_add,score_add)
 
}
rownames(scores) <- subjects
colnames(scores_add)<-hours
rownames(scores_add)<-subjects
#write.csv(scores,"scroes_10_cut_reverse.csv")
write.csv(scores_add,"scroes.csv")

#threshold
t_test_datamax<-apply(scores_add[nsubjects,],1,max)
threshold_value<-t.test(t_test_datamax)$conf.int[2]
print(threshold_value) #0.8131584 cut
for(ii in 1:length(subjects)){
  w<-which(scores_add[ii,]>threshold_value)[1]
  if(!is.na(w)){print(c(subjects[ii],as.numeric(w),hours[w]))} else{print(c(subjects[ii],"F"))}
}
########## end ################################
