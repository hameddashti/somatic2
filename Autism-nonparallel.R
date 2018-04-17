{
  library(parallel)
  library(tidyverse)
}
#=======Read files========
n1 <- readline(prompt="Enter file address: ")
no_cores<-readline(prompt="Enter number of CPUs: ")
genome<-readline(prompt="Enter the source of genome(hg19/hg18): ")
no_cores<-as.double(no_cores)
input_case<- read.csv(n1,header = T)
input_control<- read.csv('/Users/soudbeh/Desktop/ all/Autism/COOPER_cnv_control.csv',header = T)
#=======find total number of Duplication and Deletion========
no_case_dup<- nrow(subset(input_case,input_case$cnv_type=='Dup'))
no_case_del<-nrow(subset(input_case,input_case$cnv_type=='Del'))
no_control_dup<- nrow(subset(input_control,input_control$cnv_type=='Dup'))
no_control_del<- nrow(subset(input_control,input_control$cnv_type=='Del'))
#=======hg18 and hg19==========
hg18<- data.frame(chr22=247249719,
                  chr2=242951149,
                  chr3=199501827,
                  chr4=191273063,
                  chr5=180857866,
                  chr6=170899992,
                  chr7=158821424,
                  chr8=146274826,
                  chr9=140273252,
                  chr10=135374737,
                  chr11=134452384,
                  chr12=132349534,
                  chr13=114142980,
                  chr14=106368585,
                  chr15=100338915,
                  chr16=88827254,
                  chr17=78774742,
                  chr18=76117153,
                  chr19=63811651,
                  chr20=62435964,
                  chr21=46944323,
                  chr22=49691432,
                  chrX=154913754,
                  chrY=57772954,
                  chrM=16571
)
hg19<- data.frame(chr1=249250621,
                  chr2=243199373,
                  chr3=198022430,
                  chr4=191154276,
                  chr5=180915260,
                  chr6=171115067,
                  chr7=159138663,
                  chr8=146364022,
                  chr9=141213431,
                  chr10=135534747,
                  chr11=135006516,
                  chr12=133851895,
                  chr13=115169878,
                  chr14=107349540,
                  chr15=102531392,
                  chr16=90354753,
                  chr17=81195210,
                  chr18=78077248,
                  chr19=59128983,
                  chr20=63025520,
                  chr21=48129895,
                  chr22=51304566,
                  chrX=155270560,
                  chrY=59373566,
                  chrM=16569
)
#==========seprate Dup and Del================
case_del<-subset(input_case,input_case$cnv_type=='Del')
case_dup<-subset(input_case,input_case$cnv_type=='Dup')
control_del<-subset(input_control,input_control$cnv_type=='Del')
control_dup<-subset(input_control,input_control$cnv_type=='Dup')
rm(input_case,input_control)
gc()

#==========Create Matrix=======================
if (genome == 'hg18'){
  chr22_del<-matrix(0, nrow = hg18$chr22, ncol = 6)
}else{
  chr22_del<-matrix(0, nrow = hg19$chr22, ncol = 6)
}
if (genome == 'hg18'){
  chr22_dup<-matrix(0, nrow = hg18$chr22, ncol = 6)
}else{
  chr22_dup<-matrix(0, nrow = hg19$chr22, ncol = 6)
}
#chr22_case_no_del<-nrow(subset(case_del,case_del$chr=='chr22'))
#chr22_control_no_del<-nrow(subset(control_del,control_del$chr=='chr22'))
#chr22_case_no_dup<-nrow(subset(case_dup,case_dup$chr=='chr22'))
#chr22_control_no_dup<-nrow(subset(control_dup,control_dup$chr=='chr22'))
#==========Define function count===============
{
  count_del1<-function(x){
    chr22_del[x,1]<<-x
  }
  count_dup1<-function(x){
    chr22_dup[x,1]<<-x
  }
  count_del2<-function(x){
    chr22_del[x,2]<<-nrow(subset(case_del,case_del$start<=x & case_del$end >= x))
  }
  count_dup2<-function(x){
    chr22_dup[x,2]<<-nrow(subset(case_dup,case_dup$start<=x & case_dup$end >= x))
  }
  count_del3<-function(x){
    chr22_del[x,3]<<-nrow(subset(control_del,control_del$start<=x & control_del$end >= x))
  }
  count_dup3<-function(x){
    chr22_dup[x,3]<<-nrow(subset(control_dup,control_dup$start<=x & control_dup$end >= x))
  }
  count_del4<-function(x){
    fisher_matrix <<- matrix(c(chr22_del[x,2],chr22_del[x,3],no_case_del- chr22_del[x,2],no_control_del -chr22_del[x,3]),nrow = 2,ncol = 2)
    test <<- fisher.test(fisher_matrix,alternative = "two.sided",conf.level =0.9)
    chr22_del[x,4]<<- test$p.value
  }
  count_dup4<-function(x){
    fisher_matrix <<- matrix(c(chr22_dup[x,2],chr22_dup[x,3],no_case_dup- chr22_dup[x,2],no_control_dup -chr22_dup[x,3]),nrow = 2,ncol = 2)
    test <<-fisher.test(fisher_matrix,alternative = "two.sided",conf.level =0.9)
    chr22_dup[x,4]<<- test$p.value
  }
  count_del5<-function(x){
    fisher_matrix <<- matrix(c(chr22_del[x,2],chr22_del[x,3],no_case_del- chr22_del[x,2],no_control_del -chr22_del[x,3]),nrow = 2,ncol = 2)
    test<<-fisher.test(fisher_matrix,alternative = "greater",conf.level =0.9)
    chr22_del[x,5]<<- test$p.value
  }
  count_dup5<-function(x){
    fisher_matrix <<- matrix(c(chr22_dup[x,2],chr22_dup[x,3],no_case_dup- chr22_dup[x,2],no_control_dup -chr22_dup[x,3]),nrow = 2,ncol = 2)
    test<<-fisher.test(fisher_matrix,alternative = "greater",conf.level =0.9)
    chr22_dup[x,5]<<-test$p.value
  }
  count_del6<-function(x){
    fisher_matrix <<- matrix(c(chr22_del[x,2],chr22_del[x,3],no_case_del - chr22_del[x,2],no_control_del -chr22_del[x,3]),nrow = 2,ncol = 2)
    test<-fisher.test(fisher_matrix,alternative = "less",conf.level =0.9)
    chr22_del[x,6]<<- test$p.value
  }
  count_dup6<-function(x){
    fisher_matrix <<- matrix(c(chr22_dup[x,2],chr22_dup[x,3],no_case_dup - chr22_dup[x,2],no_control_dup -chr22_dup[x,3]),nrow = 2,ncol = 2)
    test<<-fisher.test(fisher_matrix,alternative = "less",conf.level =0.9)
    chr22_dup[x,6]<<- test$p.value
  }
}

#============chromosome count==================

if (genome == 'hg18'){
  start_time<-Sys.time()
  t<-hg18$chr22
  lapply(chr22_del[,1],count_del1)
  lapply(chr22_dup[,1],count_dup1)
  lapply(chr22_del[,2],count_del2)
  lapply(chr22_dup[,2],count_dup2)
  lapply(chr22_del[,3],count_del3)
  lapply(chr22_dup[,3],count_dup3)
  lapply(chr22_del[,4],count_del4)
  lapply(chr22_dup[,4],count_dup4)
  lapply(chr22_del[,5],count_del5)
  lapply(chr22_dup[,5],count_dup5)
  lapply(chr22_del[,6],count_del6)
  lapply(chr22_dup[,6],count_dup6)
  end_time<-Sys.time()
  print(start_time - end_time )
}else{
  start_time<-Sys.time()
  t<-hg19$chr22
  lapply(chr22_del[,1],count_del1)
  lapply(chr22_dup[,1],count_dup1)
  lapply(chr22_del[,2],count_del2)
  lapply(chr22_dup[,2],count_dup2)
  lapply(chr22_del[,3],count_del3)
  lapply(chr22_dup[,3],count_dup3)
  lapply(chr22_del[,4],count_del4)
  lapply(chr22_dup[,4],count_dup4)
  lapply(chr22_del[,5],count_del5)
  lapply(chr22_dup[,5],count_dup5)
  lapply(chr22_del[,6],count_del6)
  lapply(chr22_dup[,6],count_dup6)
  end_time<-Sys.time()
  print(start_time - end_time )
}
write.csv(chr22_del,file = '../chr22_del.csv',col.names = T)
write.csv(chr22_dup,file = '../chr22_dup.csv',col.names = T)

#=============Permutation================
number_of_permutation<- readline(prompt="Enter number of permutation: ")
number_of_permutation<-as.integer(number_of_permutation)
max_cnv_chr22<- max(chr22_del[,2],chr22_dup[,2],chr22_del[,3],chr22_dup[,3])
#----Del----
ASD_random_del_chr22<- matrix(, nrow =no_case_del+no_control_del , ncol = 2)
ASD_random_del_chr22[1:no_case_del,1]<-'case'
tt<-no_case_del+1
ttt<-nrow(ASD_random_del_chr22)
ASD_random_del_chr22[tt:ttt,1]<-'control'
CNVarray_report <- matrix(,nrow=number_of_permutation, ncol = max_cnv_chr22)
for (numberCNV in 1:max_cnv_chr22) {
  ASD_random_del_chr22[,2]<- runif(ttt,0,1)
  ASD_random_del_chr22<-ASD_random_del_chr22[order(ASD_random_del_chr22[,2]),]
  for (i in 1:number_of_permutation) {
    r_case_control<- sample(1:ttt, 1)
    ASD<-ASD_random_del_chr22[1:r_case_control,]
    CNV_case_positive<- nrow(subset(ASD,ASD[,1]=='case'))
    CNV_control_positive <-nrow(subset(ASD,ASD[,1]=='control'))
    CNV_case_negative <- no_case_del -CNV_case_positive
    CNV_control_negative <- no_control_del - CNV_control_positive
    contingency_table<-matrix(c(CNV_case_positive,CNV_control_positive,CNV_case_negative,CNV_control_negative),nrow = 2)
    test<-fisher.test(contingency_table,alternative = "greater",conf.level =0.9)
    pval<-test$p.value
    CNVarray_report[i,numberCNV]<-pval
  }
}
#----Dup----
ASD_random_dup_chr22<- matrix(, nrow =no_case_dup+no_control_dup , ncol = 2)
ASD_random_dup_chr22[1:no_case_dup,1]<-'case'
tt<-no_case_dup+1
ttt<-nrow(ASD_random_dup_chr22)
ASD_random_dup_chr22[tt:ttt,1]<-'control'
CNVarray_report <- matrix(,nrow=number_of_permutation, ncol = max_cnv_chr22)
for (numberCNV in 1:max_cnv_chr22) {
  ASD_random_dup_chr22[,2]<- runif(ttt,0,1)
  ASD_random_dup_chr22<-ASD_random_dup_chr22[order(ASD_random_dup_chr22[,2]),]
  for (i in 1:number_of_permutation) {
    r_case_control<- sample(1:ttt, 1)
    ASD<-ASD_random_dup_chr22[1:r_case_control,]
    CNV_case_positive<- nrow(subset(ASD,ASD[,1]=='case'))
    CNV_control_positive <-nrow(subset(ASD,ASD[,1]=='control'))
    CNV_case_negative <- no_case_dup -CNV_case_positive
    CNV_control_negative <- no_control_dup - CNV_control_positive
    contingency_table<-matrix(c(CNV_case_positive,CNV_control_positive,CNV_case_negative,CNV_control_negative),nrow = 2)
    test<-fisher.test(contingency_table,alternative = "greater",conf.level =0.9)
    pval<-test$p.value
    CNVarray_report[i,numberCNV]<-pval
  }
}
#================Significant Region========================================
ci <- readline(prompt="Enter C.I.(99.95,99.99,100): ")
ci<-as.double(ci)
significant_pval<- as.integer(number_of_permutation*(1.0 - ci/100)+1)
CNVarray_report<- CNVarray_report[order(CNVarray_report[,1]),]
significant_pval1<- CNVarray_report[significant_pval,1]
#---del---
significant_regions_chr22_del<-data.frame(region=integer(),start=integer(),end=integer(),cnv_type=character(),min_pvalue=double(),max_pvalue=double(),mean_pvalue=double(),min_case=integer(),max_cas=integer(),min_control=integer(),max_control=integer())
significant_file_del<- as.data.frame(subset(chr22_del,chr22_del[,4]<=significant_pval1))
significant_regions_chr22_del<- significant_file_del%>%
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
#---dup---
significant_regions_chr22_dup<-data.frame(region=integer(),start=integer(),end=integer(),cnv_type=character(),min_pvalue=double(),max_pvalue=double(),mean_pvalue=double(),min_case=integer(),max_cas=integer(),min_control=integer(),max_control=integer())
significant_file_dup<- as.data.frame(subset(chr22_dup,chr22_dup[,4]<=significant_pval1))
significant_regions_chr22_dup<- significant_file_dup%>%
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





