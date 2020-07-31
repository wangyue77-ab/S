

load(file = "rt_plot.Rdata")
 
library(dplyr)

library(tidyr)

library(tidyverse)

library(survival)

library(survminer)

if(T){
  rt <- rt_plot[,c(3,2,8)]
  res.cut <- surv_cutpoint(rt, 
                           time =names(rt)[1], 
                           event = names(rt)[2], 
                           variables = names(rt)[3], 
  
  
  sur.cut <- surv_cutpoint(data,   time= 'futime',
                           event = 'fustat' ,
                           variables = 'SMAD4' )
  
  surv <- surv_categorize(sur.cut)
  
  table(sur.cat$SMAD4)
  
  #or
  risk <- surv_categorize(res.cut)
  rt$risk <- risk$SMAD4
  rt$riskScore <- rt$SMAD4
  cutoff <- sort(rt$riskScore)[sum(rt$risk=="low")]
  
  my.surv <- Surv(rt$futime, rt$fustat)
  group <- rt$risk
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
}
  



ggsurvplot(fit, data = survival_dat)

if(T){
  
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  x = summary(coxph(Surv(futime, fustat)~riskScore, data = rt))
  HR = signif(x$coef[2], digits=2)
  up95 = signif(x$conf.int[,"upper .95"],2)
  low95 = signif(x$conf.int[,"lower .95"], 2)
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  
  rt <- rt[order(rt[,"riskScore"],decreasing = T),]
  ggsurvplot(fit, data = survival_dat ,
             #ggtheme = theme_bw(), 
             conf.int = F, 
             #conf.int.style = "step",
             censor = F, 
             palette = c("#D95F02","#1B9E77"), 
             font.legend = 11,
            
             legend.labs=c(paste0(">",round(cutoff,2),"(",sum(rt$risk=="high"),")"),
                           paste0("<",round(cutoff,2),"(",sum(rt$risk=="low"),")")),
             
             pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                        paste("p = ",round(p.val,3), sep = "")),
                          HR, CI, sep = "\n"))
}
