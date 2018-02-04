##Script to calibrate and geographically correct SDMs
library(readr)
library(dplyr)
library(mgcv)
library(foreach)
library(doParallel)
library(doRNG)
library(raster)
library(sdm)
##library(openblasctl)

proj_dir <- "~/code/WesternNA_SDM/"
model_dir <- "scratch/models/"
out_dir <- "scratch/models_calib/"
setwd(proj_dir)

##Downloads presence point data from Amazon S3 if it doesn't already exist.
if(!file.exists("./data/occurences_final_1_9_2018.csv")){
  aws_dl <- "~/.local/bin/aws s3 cp s3://sdmdata/occurences/occurences_final_1_9_2018.tar.gz ./data/occurences_final_1_9_2018.tar.gz"
  system(paste("cd",proj_dir,"&&",aws_dl),wait=TRUE)
  tar_dl <- "tar -xf ./data/occurences_final_1_9_2018.tar.gz -C ./data/"
  system(paste("cd",proj_dir,"&&",tar_dl),wait=TRUE)
}

##Downloads presence-absence point data from Amazon S3 if it doesn't already exist.
if(!file.exists("./data/calibration_plotdata_1_26_2018.csv")){
  aws_dl2 <- "~/.local/bin/aws s3 cp s3://sdmdata/occurences/calibration_plotdata_1_26_2018.tar.gz ./data/calibration_plotdata_1_26_2018.tar.gz"
  system(paste("cd",proj_dir,"&&",aws_dl2),wait=TRUE)
  tar_dl2 <- "tar -xf ./data/calibration_plotdata_1_26_2018.tar.gz -C ./data/"
  system(paste("cd",proj_dir,"&&",tar_dl2),wait=TRUE)
}

##Reads in presence data for species list.
pres <- read_csv(paste(proj_dir,"data/occurences_final_1_9_2018.csv",sep=""))
spp <- unique(pres$species)
gen <- unique(pres$genus)

##Reads in presence-absence data.
pres_abs <- read_csv(paste(proj_dir,"data/calibration_plotdata_1_26_2018.csv",sep=""))
pres_abs <- filter(pres_abs,AcceptedGe %in% gen)

pres_abs$PlotAreaHa[is.na(pres_abs$PlotAreaHa)] <- 0.04
pres_abs <- pres_abs[which(complete.cases(pres_abs[,19:50])),]
pres_abs$Loc <- paste("W",round(pres_abs$X/200),"N",round(pres_abs$Y/200),sep="")

##Computes SDM corrections and calibrations for each focal species.
set.seed(37)

#cl <- makeCluster(3,outfile="sdm_messages.log")
#registerDoParallel(cl)
overwrite <- TRUE
# 
spp <- c("Xerophyllum tenax","Primula parryi", "Sidalcea oregana","Phacelia corymbosa", "Rubus parviflorus", "Mertensia paniculata",
        "Vicia americana","Mimulus cardinalis","Boechera holboellii", "Dasiphora fruticosa", "Geum triflorum",
        "Phlox diffusa", "Ivesia santolinoides", "Heliomeris multiflora", "Dodecatheon redolens", "Hulsea algida",
        "Ligusticum grayi","Collinsia torreyi","Cardamine bellidifolia","Rhodiola integrifolia","Ranunculus inamoenus",
        "Minuartia obtusiloba","Viola bakeri","Angelica breweri","Penstemon davidsonii","Solidago multiradiata")
#spp <- spp[1:7]

##PROBLEM: sdm predict() function fails when run in parallel with %dopar% or %dorng%

cals <- foreach(i=1:length(spp),.packages=c("dplyr")) %do% {
  setwd(proj_dir)
  library(sdm)
  ##Defines local utility functions
  logit <- function(x){log(replace(x,x < 1e-4, 1e-4)/(1-replace(x,x < 1e-4, 1e-4)))} 
  inv_logit <- function(x){exp(x)/(1+exp(x))}
  wfun <- function(x){weighted.mean(x,w=spp_weights)}
  
  ##Prevents multithreaded linear algebra library from implicitly parallelizing.
  #openblas_set_num_threads(1)
  
  ##Checks to see if file already exists on S3
  fit_file <- paste("sdm_",gsub(" ","_",spp[i]),"_calib.Rdata",sep="")
  file_exists_aws <- system(paste("~/miniconda2/bin/aws s3 ls",paste("s3://sdmdata/models_calibrated/",
                                  fit_file,sep="")),wait=TRUE)
  
  if(file_exists_aws == 0 & overwrite==FALSE){
    print(paste("Output exists on S3, skipping..."))
    cat(paste("Output file",paste(fit_file),"exists, skipping...\n"),
        file="./scratch/sdm_progress.log",append=TRUE)
    
  }else{
    print(paste("Calibrating model for ",spp[i],"(",i,"of",length(spp),")"))
    cat(paste("Calibrating model for ",spp[i],"(",i,"of",length(spp),") on",
              Sys.time(),"\n"),file="./scratch/sdm_progress.log",append=TRUE)
    
    ##Downloads model from S3
    model_file <- paste("sdm_",gsub(" ","_",spp[i]),".Rdata", sep="")
    s3_path <- paste("s3://sdmdata/models/",model_file, sep="")
    local_path <- paste(model_dir,model_file,sep="")
    s3_dl_string <- paste("~/miniconda2/bin/aws s3 cp",s3_path,local_path)
    system(s3_dl_string,wait=TRUE)
    
    ##Loads fit model.
    print(paste("Loading model..."))
    model <- readRDS(paste(proj_dir,local_path,sep=""))
    spp_stats <- model[[1]]$stats
    spp_sdm <- model[[1]]$final_models
    spp_weights <- spp_stats$ensemble_weight
    spp_obs <- model[[1]]$final_models@models$TR_PRES$glm$'1'@object$data
    spp_pred1 <- model[[1]]$final_models@models$TR_PRES$glm$'1'@evaluation$training@predicted
    spp_pred2 <- model[[1]]$final_models@models$TR_PRES$brt$'2'@evaluation$training@predicted
    spp_pred3 <- model[[1]]$final_models@models$TR_PRES$svm$'3'@evaluation$training@predicted
    spp_pred4 <- model[[1]]$final_models@models$TR_PRES$maxent$'4'@evaluation$training@predicted
    
    ##Preps calibration data.
    print(paste("Prepping calibration data..."))
    pres_focal <- filter(pres_abs, AcceptedNa == spp[i])
    pres_focal <- filter(pres_focal,!duplicated(Loc))
    if(nrow(pres_focal) > 0){
      pres_locs <- unique(pres_focal$Loc)
      pres_focal$PRES <- 1
    }else{
      pres_locs <- NA
      pres_focal$PRES <- integer(0)
    }
    
    focal_gen <- pres_focal$AcceptedGe[1]
    gen_locs <- unique(pres_abs$Loc[pres_abs$AcceptedNa==paste(focal_gen,"NA")])
    abs_focal <- filter(pres_abs, AcceptedNa != spp[i] &
                          !duplicated(Loc) &
                          !(Loc %in% pres_locs) &
                          !(Loc %in% gen_locs))
    abs_focal$PRES <- 0
    pres_abs_focal <- rbind(pres_focal,abs_focal)
    
    ##Holds out 20% of calibration data for testing.
    pres_abs_test <- sample_frac(pres_abs_focal,size=0.5)
    pres_abs_test$train_test <- "test"
    pres_abs_test$added <- FALSE
    pres_abs_focal2 <- filter(pres_abs_focal,!(pres_abs_focal$Loc %in% pres_abs_test$Loc))
    pres_abs_focal3 <- pres_abs_focal2[,c(12,13,2,3,19:52,17)]
    pres_abs_focal3$train_test <- "train"
    pres_abs_focal3$added <- FALSE
    
    ##Adds a few points if there are no calibration presences.
    if(sum(pres_abs_focal3$PRES)<3){
      pres_added <- sample_n(filter(pres,species==spp[i]),size=3)
      pres_added$Loc <- paste("W",round(pres_added$X/200),"N",round(pres_added$Y/200),sep="")
      pres_added <- pres_added[,c(1,2,15:48,51)]
      pres_added$PRES <- 1
      pres_added$PlotAreaHa <- 0.04
      pres_added$train_test <- "train"
      pres_added$added <- TRUE
      names(pres_added) <- names(pres_abs_focal3)
      pres_abs_focal3 <- rbind(pres_added,pres_abs_focal3)
    }
    
    ##Predicts presence-only value for all pres-abs points.
    print(paste("Predicting presence-only values..."))
    pres_abs_focal3 <- as.data.frame(pres_abs_focal3)
    pres_abs_test <- as.data.frame(pres_abs_test)

    ens_pred <- predict(spp_sdm,newdata=pres_abs_focal3,type="response")
    pres_abs_focal3$Ens_pred <- apply(ens_pred,FUN=wfun,MARGIN = 1)
    pres_abs_focal3$Ens_pred[pres_abs_focal3$Ens_pred <= 1e-4] <- 1e-4
    pres_abs_focal3$pred_logit <- logit(pres_abs_focal3$Ens_pred)
    
    ens_pred_test <- predict(spp_sdm,newdata=pres_abs_test,type="response")
    pres_abs_test$Ens_pred <- apply(ens_pred_test,FUN=wfun,MARGIN = 1)
    pres_abs_test$Ens_pred[pres_abs_test$Ens_pred <= 1e-4] <- 1e-4
    pres_abs_test$pred_logit <- logit(pres_abs_test$Ens_pred)
    
    ##Measures baseline performance of model on independent data.
    pres_abs_orig_eval <- evaluate(p=as.numeric(pres_abs_focal3$Ens_pred[pres_abs_focal3$PRES==1]),
                                   a=as.numeric(pres_abs_focal3$Ens_pred[pres_abs_focal3$PRES==0]))
    orig_AUC <- pres_abs_orig_eval@auc
    print(paste("Uncorected AUC",round(orig_AUC,4)))
    
    ##Preps presence-only points for geographic correction
    print(paste("Prepping presence-background data..."))
    pred_all <- cbind(spp_pred1,spp_pred2,spp_pred3,spp_pred4)
    pred_weighted <- apply(pred_all,FUN=wfun,MARGIN = 1)
    pred_weighted[pred_weighted <= 1e-4] <- 1e-4
    
    pred_logit <- logit(pred_weighted)
    corr_data <- model[[1]]$final_models@models$TR_PRES$glm$'1'@object$data
    corr_data$ens_pred <- pred_weighted
    corr_data$pred_logit <- pred_logit
    
    ##Fits the kriging model for geographic correction.
    print(paste("Kriging with presence-background data..."))

    corr_gam <- gam(TR_PRES~s(pred_logit,k=6)+s(PCO_XSC,PCO_YSC,bs="gp"),data=corr_data,
                    family=binomial(link = "logit"))
    
    ##Predicts for heldout data.
    Corr_pred <- predict(corr_gam,newdata=pres_abs_focal3)
    Corr_pred[Corr_pred < -10] <- -10
    Corr_pred_prob <- inv_logit(Corr_pred)
    pres_abs_corr_eval <- evaluate(p=as.numeric(Corr_pred_prob[pres_abs_focal3$PRES==1]),
                                   a=as.numeric(Corr_pred_prob[pres_abs_focal3$PRES==0]))
    base_AUC <- pres_abs_corr_eval@auc
    print(paste("Baseline corrected AUC",round(base_AUC,4)))
    
    
    if(sum(pres_abs_test$PRES,na.rm=TRUE) > 3){
      ks <- c(20,40,60,80)
      corr_gam_test <- as.list(rep(NA,length(ks)))
      new_AUC <- rep(NA,length(ks))
      for(j in 1:length(ks)){
        print(paste("Testing k =",ks[j]))
        corr_gam_test[[j]] <- try(gam(TR_PRES~s(pred_logit,k=7)+s(PCO_XSC,PCO_YSC,bs="gp",k=ks[j]),data=corr_data,
                                  family=binomial(link = "logit")))
        Corr_pred_test <- predict(corr_gam_test[[j]],newdata=pres_abs_focal3)
        Corr_pred_test[Corr_pred_test < -10] <- -10
        Corr_pred_test_prob <- inv_logit(Corr_pred_test)
        Corr_pred_test_eval <- evaluate(p=as.numeric(Corr_pred_test_prob[pres_abs_focal3$PRES==1]),
                                        a=as.numeric(Corr_pred_test_prob[pres_abs_focal3$PRES==0]))
        new_AUC[j] <- Corr_pred_test_eval@auc
      }
      if(any(new_AUC > base_AUC)){
        base_AUC <- max(new_AUC)
        corr_gam <- corr_gam_test[[which.max(new_AUC)]]
      }
    }
    print(paste("Best Corrected AUC",round(base_AUC,digits=4)))
    
    ##Predicts for training data
    corr_data$Corr_pred_prob <- predict(corr_gam,newdata=corr_data,type="response")
    pres_abs_focal3$Corr_pred <- predict(corr_gam,newdata=pres_abs_focal3)
    pres_abs_focal3$Corr_pred[pres_abs_focal3$Corr_pred < -10] <- -10
    pres_abs_focal3$Corr_pred_prob <- inv_logit(pres_abs_focal3$Corr_pred)
    
    ##Predicts for test data
    pres_abs_test$Corr_pred <- predict(corr_gam,newdata=pres_abs_test)
    pres_abs_test$Corr_pred[pres_abs_test$Corr_pred < -10] <- -10
    pres_abs_test$Corr_pred_prob <- inv_logit( pres_abs_test$Corr_pred)

    ##Fits the calibration model
    print(paste("Calibrating with presence-absence data..."))
    
    ##Observed and predicted prevalence
    obs_prev <- mean(pres_abs_focal3$PRES)
    pred_prev <- mean(pres_abs_focal3$Corr_pred_prob) 
    
    ##Checks if a linear model has a positive slope.
    calib_glm0 <- glm(PRES~Corr_pred+PlotAreaHa,data=pres_abs_focal3)
    pred_slope_glm0 <- coef(calib_glm0)[2]
    plot_slope_glm0 <- coef(calib_glm0)[1]
    if(pred_slope_glm0 > 0 & plot_slope_glm0 > 0){
      print(paste("Calibration slope > 0 plotsize slope > 0, Calibrating with plot size."))
      calib_gam <- gam(PRES~s(Corr_pred,k=7)+PlotAreaHa,data=pres_abs_focal3,
                       family=binomial(link = "logit"))
      ##Predicts calibrated values at pres-abs points
      pres_abs_focal3$Calib_pred <- predict(calib_gam,newdata=pres_abs_focal3,
                                            type="response")
      
      ##Predicts calibrated values at heldout points.
      pres_abs_test$Calib_pred <- predict(calib_gam,newdata=pres_abs_test,
                                          type="response")
      calib_type <- "gam_plotsize"
    }else if(pred_slope_glm0 > 0){
      print(paste("Plot Area slope < 0, Calibrating without plot size."))
      calib_gam <- gam(PRES~s(Corr_pred,k=6),data=pres_abs_focal3,
                       family=binomial(link = "logit"))
      ##Predicts calibrated values at pres-abs points
      pres_abs_focal3$Calib_pred <- predict(calib_gam,newdata=pres_abs_focal3,
                                            type="response")
      calib_type <- "gam"
      
      ##Predicts calibrated values at heldout points.
      pres_abs_test$Calib_pred <- predict(calib_gam,newdata=pres_abs_test,
                                          type="response")
    }else{
      print(paste("Calibration slope < 0 with pres-abs data, using manual calibration."))
      calib_gam <- function(x){x*(obs_prev/pred_prev)}
      
      ##Predicts calibrated values at pres-abs points
      pres_abs_focal3$Calib_pred <- calib_gam(pres_abs_focal3$Corr_pred_prob)
      
      ##Predicts calibrated values at heldout points.
      pres_abs_test$Calib_pred <- calib_gam(pres_abs_test$Corr_pred_prob)
      calib_type <- "obs_prevalence"
    }
    
    ##Packages data and models.
    out_list <- list(corr_data=corr_data,
                     corr_model=corr_gam,
                     calib_data=pres_abs_focal3,
                     calib_type=calib_type,
                     obs_prev=obs_prev,
                     pred_prev=pred_prev,
                     calib_model=calib_gam,
                     test_data=pres_abs_test)
    
    # ggplot(out_list$test_data)+
    #   geom_point(aes(x=X,y=Y,color=(Corr_pred_prob - Ens_pred)))+
    #   scale_color_distiller(type="div")+
    #   theme_bw()
    
    ##Writes output to disk.
    out_name <- gsub(".Rdata","_calib.Rdata",model_file,fixed=TRUE)
    out_path <- paste(proj_dir,out_dir,out_name,sep="")
    saveRDS(out_list,file=out_path)
    
    ##Uploads to S3
    upl_string <- paste("~/miniconda2/bin/aws s3 cp",out_path,paste("s3://sdmdata/models_calibrated/",out_name,sep=""))
    system(upl_string,wait=TRUE)
    
    ##Removes local files.
    model_path <- paste(proj_dir,model_dir,model_file,sep="")
    rm_string1 <- paste("rm",model_path)
    rm_string2 <- paste("rm",out_path)
    system(rm_string1,wait=TRUE)
    system(rm_string2,wait=TRUE)
    
    ##Logs completion
    cat(paste("Calibration complete for ",spp[i],"(",i,"of",length(spp),") on",
              Sys.time(),"\n"),file="./scratch/sdm_progress.log",append=TRUE)

    # ##Brings in the original raster prediction.
    # rast_path <- paste(raw_pred_path,gsub(" ","_",spp[i]),"_mosaic.tif",sep="")
    # rast_pred <- raster(rast_path)
    # XY_stack <- brick("WNA_albers_XY.tif")
    # names(XY_stack) <- c('PCO_XSC', 'PCO_YSC')
    # NAvalue(XY_stack) <- 65535
    # 
    # pred_logit <- calc(rast_pred,fun=logit,filename=gsub(".tif","_logit.tif",rast_path,fixed=TRUE),
    #                    overwrite=TRUE)
    # pred_brick <- brick(XY_stack$PCO_XSC,XY_stack$PCO_YSC,pred_logit)
    # names(pred_brick) <- c("PCO_XSC","PCO_YSC","pred_logit")
    # 
    # corr_pred <- raster::predict(pred_brick, model=corr_gam, progress="text",
    #                              overwrite=TRUE,filename="~/GIS/Xerophyllum_tenax_corrected.img")
    # calib_brick <- brick(corr_pred,corr_pred)
    # calib_brick[[2]] <- 0.04
    # names(calib_brick) <- c("Corr_pred","PlotAreaHa")
    # Calib_pred <- raster::predict(calib_brick, model=calib_gam, type="response", progress="text",
    #                               filename="~/GIS/Xerophyllum_tenax_calibrated.img",overwrite=TRUE)
  }
  (out_list)
}
#stopCluster(cl)
