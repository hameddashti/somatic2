{
  library(tidyverse)
}
#=======Read files========
path<-'/home/san/halinejad/Desktop/Dashti/somatic2'
pcnt<-paste0(path,'/AGRE_cnv_control.csv')
pcs<-paste0(path,'/AGRE_cnv_case.csv')
input_case<- read.csv(pcs,header = T)
input_control<- read.csv(pcnt,header = T)
t1<-subset(input_case,input_case$CNV.Type=='Dup')
t1<-as.data.frame(t1$Patient.Id)
t1<-unique(t1)
no_case_Dup<- nrow(t1)
t1<-subset(input_control,input_control$CNV.Type=='Dup')
t1<-as.data.frame(t1$Patient.Id)
t1<-unique(t1)
no_control_Dup<- nrow(t1)
input_case<- subset(input_case,input_case$chr.1=='Y')
input_control<-subset(input_control,input_control$chr.1=='Y')
st_case<-min(input_case$start)
en_case<-max(input_case$end)
st_control<-min(input_control$start)
en_control<-max(input_control$end)
st<-min(st_case,st_control)
en<-max(en_control,en_case)
ln<-en-st+1
rm(en_case,st_case,en_control,st_control,t1)
#==========Choose Duplications================
case_Dup<-subset(input_case,input_case$CNV.Type=='Dup')
control_Dup<-subset(input_control,input_control$CNV.Type=='Dup')
rm(input_case,input_control)
gc()
#==========Create Matrix=======================
chrY_dup_cnv<-matrix(0, nrow = ln, ncol = 6)
#==========Fill matrix=========================
for (i in 1:ln) {
  chrY_dup_cnv[i,1]<-i+st-1
}
case_Dup <- as.matrix(case_Dup)
if (nrow(case_Dup) != 0) {
  for (i in 1:nrow(case_Dup)) {
    k1 <- as.integer(case_Dup[i, 3])
    k2 <- as.integer(case_Dup[i, 4])
    k1 <- k1 - st + 1
    k2 <- k2 - st + 1
    chrY_dup_cnv[k1:k2, 2] <- chrY_dup_cnv[k1:k2, 2] + 1
  }
}
control_Dup<-as.matrix(control_Dup)
if(nrow(control_Dup)!=0){
  for (i in 1:nrow(control_Dup)){
    k1<-as.integer(control_Dup[i,3])
    k2<-as.integer(control_Dup[i,4])
    k1<-k1-st+1
    k2<-k2-st+1
    chrY_dup_cnv[k1:k2,3]<-chrY_dup_cnv[k1:k2,3]+1
  }
}
chrY_dup_cnv<-subset(chrY_dup_cnv,chrY_dup_cnv[,2]>=5)
if(nrow(chrY_dup_cnv)!=0){
  for (i in 1:nrow(chrY_dup_cnv)){
    m<-matrix(c(chrY_dup_cnv[i,2],chrY_dup_cnv[i,3],no_case_Dup-chrY_dup_cnv[i,2],no_control_Dup -chrY_dup_cnv[i,3]),nrow = 2)
    chrY_dup_cnv[i,6] <-fisher.test(m,alternative = "two.sided",conf.level =0.9)$p.value
    chrY_dup_cnv[i,4] <-fisher.test(m,alternative = "greater",conf.level =0.9)$p.value
    chrY_dup_cnv[i,5] <-fisher.test(m,alternative = "less",conf.level =0.9)$p.value
  }
}
path<-'/home/san/halinejad/Desktop/Dashti/somatic2'
path<-paste0(path,'/Result/chrY_dup.csv')
write.csv(chrY_dup_cnv,path)
#======================Significant regions=======================
significant_pval<- 0.0001609981
significant_file_Dup<- as.data.frame(subset(chrY_dup_cnv,chrY_dup_cnv[,4]<=significant_pval))
significant_regions_chrY_dup_cnv<- significant_file_Dup%>%
  as.data.frame() %>%
  mutate(
    jump = V1 - c(-1, V1[-length(V1)]),
    region = cumsum(jump != 1)
  ) %>%
  group_by(region) %>%
  summarize(
    start = min(V1),
    end = max(V1),
    min_pval = min(V4),
    max_pval = max(V4),
    mean_pval = mean(V4),
    min_case = min(V2),
    max_case = max(V2),
    mean_case = mean(V2),
    min_control = min(V3),
    max_control = max(V3),
    mean_control = mean(V3),
  )

path<-'/home/san/halinejad/Desktop/Dashti/somatic2'
path<-paste0(path,'/Result/regions_chrY_dup.csv')
write.csv(significant_regions_chrY_dup_cnv,path)