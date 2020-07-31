BiocInstaller::biocLite("RTCGA.clinical")
rm(list = ls())
library(RTCGA.clinical) 

cl_df <- LUAD.clinical
tmp=as.data.frame(colnames(cl_df))
as.character(tmp$`colnames(cl_df)`[grep("days",tmp$`colnames(cl_df)`)])

cl_df_LUAD=as.data.frame(cl_df[c( 'patient.bcr_patient_barcode',
                                 "patient.vital_status",
                                 "patient.follow_ups.follow_up.days_to_death",
                                 "patient.follow_ups.follow_up.days_to_last_followup",
                                 "patient.gender",
                                 "patient.age_at_initial_pathologic_diagnosis")])
names(cl_df_LUAD) <- c("TCGA_id","vital_status","days_to_death","days_to_last_followup","gender","age")
cl_df_LUAD$stage <- "LUAD"


if(T){
  cl_df <- LUSC.clinical
  
  cl_df_LUSC=as.data.frame(cl_df[c( 'patient.bcr_patient_barcode',
                                   "patient.vital_status",
                                   "patient.follow_ups.follow_up.days_to_death",
                                   "patient.follow_ups.follow_up.days_to_last_followup",
                                   "patient.gender",
                                   "patient.age_at_initial_pathologic_diagnosis")])
  names(cl_df_LUSC) <- c("TCGA_id","vital_status","days_to_death","days_to_last_followup","gender","age")
  cl_df_LUSC$stage <- "LUSC"
}


cl_df <- rbind(cl_df_LUSC,cl_df_LUAD)

cl_df$fustat <-NA
cl_df$futime <-NA

cl_df[is.na(cl_df)] =0

for (i in 1:nrow(cl_df)) {
  if(cl_df$vital_status[i]=="alive"){
    cl_df$fustat[i]= 0
    cl_df$futime[i]= cl_df$days_to_last_followup[i]
  }else{
    cl_df$fustat[i]=1
    cl_df$futime[i]=cl_df$days_to_death[i]
  }
}

cl_df <- cl_df[,c("TCGA_id","fustat", "futime","gender", "age", "stage")]
cl_df$TCGA_id  <- toupper(cl_df$TCGA_id)
save(cl_df,file = "cl_df.Rdata")


rm(list = ls())
load("NSCLC_exprSet_plot.Rda")
load("cl_df.Rdata")
test <- NSCLC_exprSet[1:10,1:10]

index <- grep("SMAD4",colnames(NSCLC_exprSet))
rt <- NSCLC_exprSet[,c(1,2,index)]

library(dplyr)
rt$TCGA_id <- substring(rownames(rt),1,12)
rt_plot <- rt %>% 
  arrange(desc(SMAD4)) %>% 
  distinct(TCGA_id,.keep_all = T) %>% 
  select(-stage) %>% 
  inner_join(cl_df,by="TCGA_id")
rt_plot <- rt_plot[,c(3,4,5,6,7,8,1,2)]
rt_plot$futime <- as.numeric(rt_plot$futime)
#save(rt_plot,file = "rt_plot.Rdata")
load(file = "rt_plot.Rdata")


library(ggpubr)
my_comparisons <- list(
  c("female", "male")
)
ggboxplot(
  rt_plot, x = "gender", y = "SMAD4",
  color = "gender", palette = "lacent",
  add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")


data$age <- ifelse(rt_plot$age >= 50,"age >= 50","age <50")
clinical_plot("age")

clinical_plot("stage")

clinical_plot("sample")


#### ºÏ²¢
if(T){
  ## gender
  p1 <- clinical_plot("gender")
  ## age
  p2 <- clinical_plot("age")
  ## stage
  p3 <- clinical_plot("stage")
  ## sample
  p4 <- clinical_plot("sample")
  
  library(cowplot)
  plot_grid(p1,p2,p3,p4,labels = LETTERS[1:4])
}