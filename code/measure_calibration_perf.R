##Script to measure the performance of the SDM calibration.
library(dismo)
library(readr)
library(dplyr)
library(foreach)
library(doParallel)


proj_dir <- "~/code/WesternNA_SDM/"
model_dir <- "scratch/models/"
out_dir <- "scratch/models_calib/"
setwd(proj_dir)


pres <- read_csv(paste(proj_dir,"./data/occurences_final_3_1_2018.csv",sep=""))
spp_samp <- unique(pres$species)

# spp_samp <- c("Primula parryi", "Sidalcea oregana","Phacelia corymbosa", "Rubus parviflorus", "Mertensia paniculata",
#              "Vicia americana","Mimulus cardinalis","Boechera holboellii", "Dasiphora fruticosa", "Geum triflorum",
#              "Phlox diffusa", "Ivesia santolinoides", "Heliomeris multiflora", "Dodecatheon redolens", "Hulsea algida",
#              "Ligusticum grayi","Collinsia torreyi","Cardamine bellidifolia","Rhodiola integrifolia","Ranunculus inamoenus",
#              "Minuartia obtusiloba","Viola bakeri","Angelica breweri","Penstemon davidsonii","Solidago multiradiata")
overwrite=FALSE
#spp_samp <- spp_samp[1:7]

calib_stats <- foreach(i=1:length(spp_samp),.combine='rbind',.packages=c("sdm","dismo")) %do% {
  
  print(paste("Evaluating model for",spp_samp[i]))
  
  model_file <- paste("sdm_",gsub(" ","_",spp_samp[i]),"_calib.Rdata", sep="")
  local_path <- paste(model_dir,model_file,sep="")
  
  if(!file.exists(paste(model_dir,model_file,sep="")) | overwrite==TRUE){
    s3_path <- paste("s3://sdmdata/models_calibrated/",model_file, sep="")
    s3_dl_string <- paste("~/miniconda2/bin/aws s3 cp",s3_path,local_path)
    system(s3_dl_string,wait=TRUE)
  }
  
  ##Loads model
  spp_model <- readRDS(paste(proj_dir,local_path,sep=""))
  
  ##Extracts test stats.
  eval_data <- as.tbl(spp_model$test_data)
  if(sum(eval_data$PRES)==0){
    out_stats <- data.frame(spp=spp_samp[i],
                            orig_np=sum(spp_model$corr_data$TR_PRES==1),
                            orig_na=sum(spp_model$corr_data$TR_PRES==0),
                            orig_nmodels=sum(!is.na(as.numeric(spp_model$orig_stats$AUC))),
                            orig_AUC_te=weighted.mean(as.numeric(spp_model$orig_stats$AUC),
                                                      w=as.numeric(spp_model$orig_stats$ensemble_weight)),
                            eval_np=0,
                            eval_na=length(eval_data$PRES),
                            eval_AUC=NA,
                            eval_cal=NA,
                            corr_AUC=NA,
                            corr_cal=NA,
                            calib_AUC=NA,
                            calib_cal=NA)
  }else{
    eval_orig <- evaluate(p=eval_data$Ens_pred[eval_data$PRES==1],
                          a=eval_data$Ens_pred[eval_data$PRES==0])
    eval_orig
    
    eval_corr <- evaluate(p=eval_data$Corr_pred_prob[eval_data$PRES==1],
                          a=eval_data$Corr_pred_prob[eval_data$PRES==0])
    eval_corr
    
    eval_calib <- evaluate(p=eval_data$Calib_pred[eval_data$PRES==1],
                           a=eval_data$Calib_pred[eval_data$PRES==0])
    eval_calib
    
    ##Measures calibration
    calib_orig <- calibration(x=eval_data$PRES,
                              p=eval_data$Ens_pred)
    
    calib_corr <- calibration(x=eval_data$PRES,
                              p=eval_data$Corr_pred_prob)
    
    calib_calib <- calibration(x=eval_data$PRES,
                               p=eval_data$Calib_pred)
    
    out_stats <- data.frame(spp=spp_samp[i],
                            orig_np=sum(spp_model$corr_data$TR_PRES==1),
                            orig_na=sum(spp_model$corr_data$TR_PRES==0),
                            orig_nmodels=sum(!is.na(as.numeric(spp_model$orig_stats$AUC))),
                            orig_AUC_te=weighted.mean(as.numeric(spp_model$orig_stats$AUC),
                                                      w=as.numeric(spp_model$orig_stats$ensemble_weight)),
                            eval_np=eval_orig@np,
                            eval_na=eval_orig@na,
                            eval_AUC=eval_orig@auc,
                            eval_cal=calib_orig@statistic,
                            corr_AUC=eval_corr@auc,
                            corr_cal=calib_corr@statistic,
                            calib_AUC=eval_calib@auc,
                            calib_cal=calib_calib@statistic)
  }
  (out_stats)
}
write.csv(out_stats,"./output/sdm_calib_stats.csv",row.names=FALSE)


calib_stats$AUC_imp <- calib_stats$calib_AUC - calib_stats$eval_AUC
View(calib_stats)
summary(calib_stats)
plot(density(calib_stats$AUC_imp[calib_stats$eval_np>100],na.rm=TRUE))
abline(v=0,lty=2)

library(ggplot2)
library(gridExtra)
ggplot(subset(calib_stats,eval_np>20))+
  geom_point(aes(x=eval_AUC,y=AUC_imp,size=log(eval_np)))+
  #scale_y_continuous(limits=c(-0.05,0.15))+
  scale_size_continuous(range=c(0.0001,3))+
  theme_bw()

p1 <- ggplot(subset(calib_stats,eval_np>20))+
  geom_point(aes(x=eval_AUC,y=calib_AUC,size=(eval_np/(eval_np+eval_na))),
             fill="white",shape=21)+
  geom_abline(aes(intercept=0,slope=1),linetype="dotted")+
  scale_y_continuous("AUC (Calibrated Model)",limits=c(0.5,1))+
  scale_x_continuous("AUC (Original Model)",limits=c(0.5,1))+
  scale_size_continuous("Prevalence",range=c(0.0001,5),breaks=c(0.001,0.01,0.1))+
  theme_bw()+
  theme(panel.grid=element_blank())

p2 <- ggplot(subset(calib_stats,eval_np>20))+
  geom_point(aes(x=eval_cal,y=calib_cal,size=(eval_np/(eval_np+eval_na))),
             fill="white",shape=21)+
  geom_abline(aes(intercept=0,slope=1),linetype="dotted")+
  scale_y_continuous("Calibration (Calibrated Model)",limits=c(0.4,1))+
  scale_x_continuous("Calibration (Original Model)",limits=c(0.4,1))+
  scale_size_continuous("Prevalence",range=c(0.0001,5),breaks=c(0.001,0.01,0.1))+
  theme_bw()+
  theme(panel.grid=element_blank())

pdf("./figs/calibration_AUC.pdf",width=10,height=4)
grid.arrange(p1,p2,ncol=2)
dev.off()

p2 <- ggplot(subset(calib_stats,eval_np>20))+
  geom_point(aes(x=(orig_np/(orig_np+orig_na)),y=orig_AUC_te),
             fill="white",shape=21)+
  #geom_abline(aes(intercept=0,slope=1),linetype="dotted")+
  scale_y_continuous("AUC (Calibrated Model)",limits=c(0.4,1))+
  #scale_x_continuous("Calibration (Original Model)",limits=c(0.4,1))+
  scale_size_continuous("Prevalence",range=c(0.0001,5),breaks=c(0.001,0.01,0.1))+
  theme_bw()+
  theme(panel.grid=element_blank())
