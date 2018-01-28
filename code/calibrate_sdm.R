##Script to calibrate and geographically correct SDMs
library(sdm)
library(readr)
library(dplyr)
library(mgcv)
library(foreach)
library(doParallel)
library(doRNG)
library(raster)
##library(openblasctl)

proj_dir <- "~/code/WesternNA_SDM/"
model_dir <- "./scratch/models/"
out_dir <- "./scratch/models_calib/"
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
pres_abs$Loc <- paste("W",round(pres_abs$X/1000),"N",round(pres_abs$Y/1000),sep="")

##Computes SDM corrections and calibrations for each focal species.
set.seed(37)

cl <- makeCluster(47)
registerDoParallel(cl)
overwrite <- FALSE

foreach(i=1:length(spp),.packages=c()) %dorng% {
  
  ##Defines local utility functions
  logit <- function(x){log(replace(x,x < 1e-5, 1e-5)/(1-replace(x,x < 1e-5, 1e-5)))} 
  inv_logit <- function(x){exp(x)/(1+exp(x))}
  wfun <- function(x){weighted.mean(x,w=spp_weights)}
  
  ##Prevents multithreaded linear algebra library from implicitly parallelizing.
  openblas_set_num_threads(1)
  
  ##Checks to see if file already exists on S3
  fit_file <- paste("sdm_",gsub(" ","_",spp[i]),"_calib.Rdata",sep="")
  file_exists_aws <- system(paste("~/.local/bin/aws s3 ls s3://sdmdata/models_calibrated/",
                                  fit_file))
  
  if(file_exists_aws == 0 & overwrite==FALSE){
    
    cat(paste("Output file",paste(fit_file),"exists, skipping...\n"),
        file="./scratch/sdm_progress.log",append=TRUE)
    
  }else{
    
    cat(paste("Calibrating model for ",spp[i],"(",i,"of",length(spp),") on",
              Sys.time(),"\n"),file="./scratch/sdm_progress.log",append=TRUE)
    
    ##Downloads model from S3
    model_file <- paste("sdm_",gsub(" ","_",spp[i]),".Rdata", sep="")
    s3_path <- paste("s3://sdmdata/models/",model_file, sep="")
    local_path <- paste(model_dir,model_file,sep="")
    s3_dl_string <- paste("~/miniconda2/bin/aws s3 cp",s3_path,local_path)
    system(s3_dl_string)
    
    ##Loads fit model.
    model <- readRDS(local_path)
    spp_stats <- model[[1]]$stats
    spp_sdm <- model[[1]]$final_models
    spp_weights <- spp_stats$ensemble_weight
    spp_obs <- model[[1]]$final_models@models$TR_PRES$glm$'1'@object$data
    spp_pred1 <- model[[1]]$final_models@models$TR_PRES$glm$'1'@evaluation$training@predicted
    spp_pred2 <- model[[1]]$final_models@models$TR_PRES$brt$'2'@evaluation$training@predicted
    spp_pred3 <- model[[1]]$final_models@models$TR_PRES$svm$'3'@evaluation$training@predicted
    spp_pred4 <- model[[1]]$final_models@models$TR_PRES$maxent$'4'@evaluation$training@predicted
    
    ##Preps calibration data.
    pres_focal <- filter(pres_abs, AcceptedNa == spp[i])
    pres_focal <- filter(pres_focal,!duplicated(Loc))
    
    pres_locs <- unique(pres_focal$Loc)
    pres_focal$PRES <- 1
    focal_gen <- pres_focal$AcceptedGe[1]
    gen_locs <- unique(pres_abs$Loc[pres_abs$AcceptedNa==paste(focal_gen,"NA")])
    abs_focal <- filter(pres_abs, AcceptedNa != spp[i] &
                          !duplicated(Loc) &
                          !(Loc %in% pres_locs) &
                          !(Loc %in% gen_locs))
    abs_focal$PRES <- 0
    pres_abs_focal <- rbind(pres_focal,abs_focal)
    
    ##Holds out 20% of calibration data for testing.
    pres_abs_test <- sample_frac(pres_abs_focal,size=0.2)
    pres_abs_test$train_test <- "test"
    pres_abs_test$added <- FALSE
    pres_abs_focal2 <- filter(pres_abs_focal,!(pres_abs_focal$Loc %in% pres_abs_test$Loc))
    pres_abs_focal3 <- pres_abs_focal2[,c(12,13,2,3,19:52,17)]
    pres_abs_test$train_test <- "train"
    pres_abs_focal3$added <- FALSE
    
    ##Adds a few points if there are no calibration presences.
    if(sum(pres_abs_focal2$PRES)<3){
      pres_added <- sample_n(filter(pres,species==spp[i]),size=3)
      pres_added$Loc <- paste("W",round(pres_added$X/1000),"N",round(pres_added$Y/1000),sep="")
      pres_added <- pres_added[,c(1,2,15:48,51)]
      pres_added$PRES <- 1
      pres_added$PlotAreaHa <- 0.04
      pres_added$added <- TRUE
      names(pres_added) <- names(pres_abs_focal3)
      pres_abs_focal3 <- rbind(pres_added,pres_abs_focal3)
    }
    
    ##Predicts presence-only value for all pres-abs points.
    pres_abs_focal3 <- as.data.frame(pres_abs_focal3)
    pres_abs_test <- as.data.frame(pres_abs_test)
    ens_pred <- predict(spp_sdm,newdata=pres_abs_focal3,type="response")
    pres_abs_focal3$Ens_pred <- apply(ens_pred,FUN=wfun,MARGIN = 1)
    pres_abs_focal3$Ens_pred[pres_abs_focal3$Ens_pred <= 0.00001] <- 1e-5
    pres_abs_focal3$pred_logit <- logit(pres_abs_focal3$Ens_pred)
    
    ens_pred_test <- predict(spp_sdm,newdata=pres_abs_test,type="response")
    pres_abs_test$Ens_pred <- apply(ens_pred_test,FUN=wfun,MARGIN = 1)
    pres_abs_test$Ens_pred[pres_abs_test$Ens_pred <= 0.00001] <- 1e-5
    pres_abs_test$pred_logit <- logit(pres_abs_test$Ens_pred)
    
    ##Preps presence-only points for geographic correction
    pred_all <- cbind(spp_pred1,spp_pred2,spp_pred3,spp_pred4)
    pred_weighted <- apply(pred_all,FUN=wfun,MARGIN = 1)
    pred_weighted[pred_weighted <= 0.00001] <- 1e-5
    
    pred_logit <- logit(pred_weighted)
    corr_data <- model[[1]]$final_models@models$TR_PRES$glm$'1'@object$data
    corr_data$pred_logit <- pred_logit
    
    ##Fits the kriging model for geographic correction.
    corr_gam <- gam(TR_PRES~s(pred_logit,k=5)+s(PCO_XSC,PCO_YSC,bs="gp",k=50),data=corr_data,
                    family=binomial(link = "logit"))
    
    ##Predicts for training data
    pres_abs_focal3$Corr_pred <- predict(corr_gam,newdata=pres_abs_focal3)
    pres_abs_focal3$Corr_pred[pres_abs_focal3$Corr_pred < -10] <- -10
    pres_abs_focal3$Corr_pred_prob <- inv_logit(pres_abs_focal3$Corr_pred)
    
    ##Predicts for heldout data.
    pres_abs_test$Corr_pred <- predict(corr_gam,newdata=pres_abs_test)
    pres_abs_test$Corr_pred[pres_abs_test$Corr_pred < -10] <- -10
    pres_abs_test$Corr_pred_prob <- inv_logit(pres_abs_test$Corr_pred)
    
    ##Fits the calibration model
    calib_gam <- gam(PRES~s(Corr_pred,k=4)+PlotAreaHa,data=pres_abs_focal3,
                     family=binomial(link = "logit"))
    
    ##Predicts calibrated values at pres-abs points
    pres_abs_focal3$Calib_pred <- predict(calib_gam,newdata=pres_abs_focal3,
                                          type="response")
    
    ##Predicts calibrated values at heldout points.
    pres_abs_test$Calib_pred <- predict(calib_gam,newdata=pres_abs_test,
                                        type="response")
    
    ##Packages data and models.
    out_list <- list(corr_data=corr_data,
                     corr_model=corr_gam,
                     calib_data=pres_abs_focal3,
                     calib_model=calib_gam,
                     test_data=pres_abs_test)
    
    ##Writes output to disk.
    out_name <- gsub(".Rdata","_calib.Rdata",model_file,fixed=TRUE)
    out_path <- paste(out_dir,out_name,sep="")
    saveRDS(out_list,file=out_path)
    
    ##Uploads to S3
    upl_string <- paste("~/miniconda2/bin/aws s3 cp",out_path,paste("s3://sdmdata/models_calibrated/",out_name,sep=""))
    system(upl_string,wait=TRUE)
    
    ##Removes local files.
    rm_string1 <- paste("rm",model_file)
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
}

