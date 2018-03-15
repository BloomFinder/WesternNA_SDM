##Script to calibrate and geographically correct SDMs
library(readr)
library(dplyr)
library(mgcv)
library(foreach)
library(doParallel)
install.packages("doRNG")
library(doRNG)
library(raster)
library(sdm)
library(dismo)
library(openblasctl)

#proj_dir <- "~/code/WesternNA_SDM/"
proj_dir <- "/home/rstudio/WesternNA_SDM/"
model_dir <- "scratch/models/"
out_dir <- "scratch/models_calib/"
raw_pred_path <- "output/WNA_mosaic/"
calib_pred_path <- "output/WNA_mosaic_calib/"
rast_pred_path <- "~/GIS/"
#gdal_path <- "/Library/Frameworks/GDAL.framework/Programs/"
gdal_path <- "/usr/bin/"
#aws_path <- "~/miniconda2/bin/"
aws_path <- "/home/rstudio/.local/bin/"
setwd(proj_dir)

##Adds custom svm function to sdm package.
#mnames <- names(sdm::getmethodNames())
#if(!("svm4" %in% mnames)){
#  source(paste(proj_dir,"/code/svm4.R",sep=""))
#  sdm::add(methodInfo,w='sdm')
#}
#if(!("gbmstep3" %in% mnames)){
#  source(paste(proj_dir,"/code/gbmstep.R",sep=""))
#  sdm::add(methodInfo,w='sdm')
#}


##Downloads presence point data from Amazon S3 if it doesn't already exist.
if(!file.exists("./data/occurences_final_3_1_2018.csv")){
  aws_dl <- paste(aws_path,"aws s3 cp s3://sdmdata/occurences/occurences_final_3_1_2018.tar.gz ./data/occurences_final_3_1_2018.tar.gz",sep="")
  system(paste("cd",proj_dir,"&&",aws_dl),wait=TRUE)
  tar_dl <- "tar -xf ./data/occurences_final_3_1_2018.tar.gz -C ./data/"
  system(paste("cd",proj_dir,"&&",tar_dl),wait=TRUE)
}

##Downloads presence-absence point data from Amazon S3 if it doesn't already exist.
if(!file.exists("./data/calibration_plotdata_3_1_2018.csv")){
  aws_dl2 <- paste(aws_path,"aws s3 cp s3://sdmdata/occurences/calibration_plotdata_3_1_2018.tar.gz ./data/calibration_plotdata_3_1_2018.tar.gz",sep="")
  system(paste("cd",proj_dir,"&&",aws_dl2),wait=TRUE)
  tar_dl2 <- "tar -xf ./data/calibration_plotdata_3_1_2018.tar.gz -C ./data/"
  system(paste("cd",proj_dir,"&&",tar_dl2),wait=TRUE)
}

##Reads in presence data for species list.
pres <- read_csv(paste(proj_dir,"data/occurences_final_3_1_2018.csv",sep=""))
spp <- unique(pres$species)
gen <- unique(pres$genus)

##Reads in presence-absence data.
pres_abs <- read_csv(paste(proj_dir,"data/calibration_plotdata_3_1_2018.csv",sep=""))
pres_abs <- filter(pres_abs,acceptedge %in% gen)

pres_abs$plotareaha[is.na(pres_abs$plotareaha)] <- 0.04
pres_abs <- pres_abs[which(complete.cases(pres_abs[,20:51])),]
pres_abs$Loc <- paste("W",round(pres_abs$x/200),"N",round(pres_abs$y/200),sep="")

##Computes SDM corrections and calibrations for each focal species.
set.seed(37)

overwrite <- TRUE
# 
# spp <- c("Agastache pallidiflora",
#          "Campanula scabrella",
#          "Asclepias tuberosa",
#          "Agrimonia striata",
#          "Ageratina_herbacea",
#          "Calochortus gunnisonii",
#          "Anemone multifida",
#          "Cassiope mertensiana",
#          "Arnica latifolia",
#          "Aconitum columbianum",
#          "Agastache urticifolia",
#          "Bistorta bistortoides",
#          "Ceanothus velutinus",
#          "Conioselinum scopulorum",
#          "Balsamorhiza sagittata",
#          "Erythronium montanum",
#          "Amelanchier utahensis",
#          "Ligusticum porteri",
#          "Lilium columbianum",
#          "Oplopanax horridus",
#          "Oreostemma alpigenum",
#          "Maianthemum racemosum",
#          "Pedicularis contorta",
#          "Pedicularis attollens",
#          "Silene douglasii",
#          "Rosa woodsii",
#          "Xerophyllum tenax")
#spp <- spp[1:7]

##PROBLEM: sdm predict() function fails when run in parallel with %dopar% or %dorng%
cl <- makeCluster(95,outfile="sdm_messages.log")
registerDoParallel(cl)

cals <- foreach(i=1:length(spp),.packages=c("dplyr")) %dorng% {
  setwd(proj_dir)
  library(gbm)
  library(sdm)

  
  mnames <- names(sdm::getmethodNames())
  if(!("svm4" %in% mnames)){
    source(paste(proj_dir,"/code/svm4.R",sep=""))
    sdm::add(methodInfo,w='sdm')
  }
  if(!("gbmstep3" %in% mnames)){
    source(paste(proj_dir,"/code/gbmstep.R",sep=""))
    sdm::add(methodInfo,w='sdm')
  }
  
  
  ##Defines local utility functions
  logit <- function(x){log(replace(x,x < 1e-4, 1e-4)/(1-replace(x,x < 1e-4, 1e-4)))} 
  inv_logit <- function(x){exp(x)/(1+exp(x))}
  wfun <- function(x){weighted.mean(x,w=spp_weights)}
  
  ##Prevents multithreaded linear algebra library from implicitly parallelizing.
  openblas_set_num_threads(1)
  
  ##Checks to see if file already exists on S3
  fit_file <- paste("sdm_",gsub(" ","_",spp[i]),"_calib.Rdata",sep="")
  file_exists_aws <- system(paste(aws_path, "aws s3 ls ",paste("s3://sdmdata/models_calibrated/",
                                  fit_file,sep=""),sep=""),wait=TRUE)
  
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
    s3_dl_string <- paste(aws_path,"aws s3 cp ",s3_path," ",local_path, sep="")
    system(s3_dl_string,wait=TRUE)
    
    ##Loads fit model.
    print(paste("Loading model..."))
    model <- readRDS(paste(proj_dir,local_path,sep=""))
    spp_stats <- model[[1]]$stats
    spp_sdm <- model[[1]]$final_models
    spp_weights <- spp_stats$ensemble_weight
   
    ##Preps calibration data.
    print(paste("Prepping calibration data..."))
    pres_focal <- filter(pres_abs, acceptedna == spp[i])
    pres_focal <- filter(pres_focal,!duplicated(Loc))
    if(nrow(pres_focal) > 0){
      pres_locs <- unique(pres_focal$Loc)
      pres_focal$PRES <- 1
    }else{
      pres_locs <- NA
      pres_focal$PRES <- integer(0)
    }
    
    focal_gen <- pres_focal$acceptedge[1]
    gen_locs <- unique(pres_abs$Loc[pres_abs$acceptedna==paste(focal_gen,"NA")])
    abs_focal <- filter(pres_abs, acceptedna != spp[i] &
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
    pres_abs_focal3 <- pres_abs_focal2[,c(13,14,3,4,20:53,18)]
    pres_abs_focal3$train_test <- "train"
    pres_abs_focal3$added <- FALSE
    
    ##Adds a few points if there are no calibration presences.
    if(sum(pres_abs_focal3$PRES)<3){
      pres_added <- sample_n(filter(pres,species==spp[i]),size=3)
      pres_added$Loc <- paste("W",round(pres_added$x/200),"N",round(pres_added$y/200),sep="")
      pres_added <- pres_added[,c(2,3,16:49,52)]
      pres_added$PRES <- 1
      pres_added$plotareaha <- 0.04
      pres_added$train_test <- "train"
      pres_added$added <- TRUE
      names(pres_added) <- names(pres_abs_focal3)
      pres_abs_focal3 <- rbind(pres_added,pres_abs_focal3)
    }
    
    ##Predicts presence-only value for all pres-abs points.
    print(paste("Predicting presence-only values..."))
    pres_abs_focal3 <- as.data.frame(pres_abs_focal3)
    pres_abs_test <- as.data.frame(pres_abs_test)

    ens_pred1 <- predict.gbm(object=spp_sdm@models$TR_PRES$gbmstep3$`1`@object,newdata=pres_abs_focal3,type="response",
                                n.trees=spp_sdm@models$TR_PRES$gbmstep3$`1`@object$gbm.call$best.trees)
    ens_pred2 <- kernlab::predict(spp_sdm@models$TR_PRES$svm4$`2`@object,newdata=pres_abs_focal3,type="probabilities")[,2]
    ens_pred3 <- dismo::predict(spp_sdm@models$TR_PRES$maxent$`3`@object,x=pres_abs_focal3,type="response")
    ens_pred <- apply(cbind(ens_pred1,ens_pred2,ens_pred3),MARGIN=1,FUN=wfun)
    pres_abs_focal3$Ens_pred <- ens_pred
    pres_abs_focal3$Ens_pred[pres_abs_focal3$Ens_pred <= 1e-4] <- 1e-4
    pres_abs_focal3$Ens_pred[pres_abs_focal3$Ens_pred > (1 - 1e-4)] <- (1 - 1e-4)
    pres_abs_focal3$pred_logit <- logit(pres_abs_focal3$Ens_pred)
    
    ens_predt1 <- predict.gbm(object=spp_sdm@models$TR_PRES$gbmstep3$`1`@object,newdata=pres_abs_test,type="response",
                             n.trees=spp_sdm@models$TR_PRES$gbmstep3$`1`@object$gbm.call$best.trees)
    ens_predt2 <- kernlab::predict(spp_sdm@models$TR_PRES$svm4$`2`@object,newdata=pres_abs_test,type="probabilities")[,2]
    ens_predt3 <- dismo::predict(spp_sdm@models$TR_PRES$maxent$`3`@object,x=pres_abs_test,type="response")
    ens_pred_test <- apply(cbind(ens_predt1,ens_predt2,ens_predt3),MARGIN=1,FUN=wfun)
    
    pres_abs_test$Ens_pred <- ens_pred_test
    pres_abs_test$Ens_pred[pres_abs_test$Ens_pred <= 1e-4] <- 1e-4
    pres_abs_test$Ens_pred[pres_abs_test$Ens_pred > (1 - 1e-4)] <- (1 - 1e-4)
    pres_abs_test$pred_logit <- logit(pres_abs_test$Ens_pred)
    
    ##Measures baseline performance of model on independent data.
    pres_abs_orig_eval <- evaluate(p=as.numeric(pres_abs_focal3$Ens_pred[pres_abs_focal3$PRES==1]),
                                   a=as.numeric(pres_abs_focal3$Ens_pred[pres_abs_focal3$PRES==0]))
    orig_AUC <- pres_abs_orig_eval@auc
    print(paste("Uncorrected AUC",round(orig_AUC,4)))
    
    ##Preps presence-only points for geographic correction
    print(paste("Prepping presence-background data..."))
    spp_data <- model[[1]]$data
    spp_data <- spp_data[complete.cases(spp_data[,3:35]),]
    
    pred_a1 <- predict.gbm(object=spp_sdm@models$TR_PRES$gbmstep3$`1`@object,newdata=spp_data,type="response",
                              n.trees=spp_sdm@models$TR_PRES$gbmstep3$`1`@object$gbm.call$best.trees)
    pred_a2 <- kernlab::predict(spp_sdm@models$TR_PRES$svm4$`2`@object,newdata=spp_data,type="probabilities")[,2]
    pred_a3 <- dismo::predict(spp_sdm@models$TR_PRES$maxent$`3`@object,x=spp_data,type="response")
    
    
    pred_weighted <- apply(cbind(pred_a1,pred_a2,pred_a3), FUN=wfun,MARGIN = 1)
    pred_weighted[pred_weighted <= 1e-4] <- 1e-4
    pred_weighted[pred_weighted >= (1 - 1e-4)] <- (1 - 1e-4)
    
    
    pred_logit <- logit(pred_weighted)
    corr_data <- spp_data
    corr_data$ens_pred <- pred_weighted
    corr_data$pred_logit <- pred_logit
    
    ##Fits the kriging model for geographic correction.
    print(paste("Kriging with presence-background data..."))

    corr_gam <- gam(TR_PRES~s(pred_logit,k=14)+s(PCO_XSC,PCO_YSC,bs="gp"),data=corr_data,
                    family=binomial(link = "logit"))
    
    ##Predicts for heldout data.
    Corr_pred <- predict(corr_gam,newdata=pres_abs_focal3)
    Corr_pred[Corr_pred < -10] <- -10
    Corr_pred_prob <- inv_logit(Corr_pred)
    pres_abs_corr_eval <- evaluate(p=as.numeric(Corr_pred_prob[pres_abs_focal3$PRES==1]),
                                   a=as.numeric(Corr_pred_prob[pres_abs_focal3$PRES==0]))
    base_AUC <- pres_abs_corr_eval@auc
    corr_k <- NA
    print(paste("Baseline corrected AUC",round(base_AUC,4)))
    
    if(sum(pres_abs_test$PRES,na.rm=TRUE) > 3){
      ks <- c(20,40,60,80,100)
      corr_gam_test <- as.list(rep(NA,length(ks)))
      new_AUC <- rep(NA,length(ks))
      for(j in 1:length(ks)){
        corr_gam_test[[j]] <- try(gam(TR_PRES~s(pred_logit,k=14)+s(PCO_XSC,PCO_YSC,bs="gp",k=ks[j]),data=corr_data,
                                  family=binomial(link = "logit")))
        Corr_pred_test <- predict(corr_gam_test[[j]],newdata=pres_abs_focal3)
        Corr_pred_test[Corr_pred_test < -10] <- -10
        Corr_pred_test_prob <- inv_logit(Corr_pred_test)
        Corr_pred_test_eval <- evaluate(p=as.numeric(Corr_pred_test_prob[pres_abs_focal3$PRES==1]),
                                        a=as.numeric(Corr_pred_test_prob[pres_abs_focal3$PRES==0]))
        new_AUC[j] <- Corr_pred_test_eval@auc
        print(paste("Testing k =",ks[j],", corrected AUC", new_AUC[j]))
      }
      if(any(new_AUC > base_AUC)){
        base_AUC <- max(new_AUC)
        corr_gam <- corr_gam_test[[which.max(new_AUC)]]
        corr_k <- ks[which.max(new_AUC)]
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
    calib_glm0 <- glm(PRES~Corr_pred+plotareaha,data=pres_abs_focal3)
    pred_slope_glm0 <- coef(calib_glm0)[2]
    plot_slope_glm0 <- coef(calib_glm0)[1]
    
    if(pred_slope_glm0 > 0 & plot_slope_glm0 > 0){
      print(paste("Calibration slope > 0 plotsize slope > 0, Calibrating with plot size."))
      calib_gam <- gam(PRES~s(Corr_pred,k=5)+plotareaha,data=pres_abs_focal3,
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
      calib_gam <- gam(PRES~s(Corr_pred,k=5),data=pres_abs_focal3,
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
    out_list <- list(orig_stats=spp_stats,
                     corr_data=corr_data,
                     corr_model=corr_gam,
                     corr_k=corr_k,
                     calib_data=pres_abs_focal3,
                     calib_type=calib_type,
                     obs_prev=obs_prev,
                     pred_prev=pred_prev,
                     calib_prev=mean(pres_abs_focal3$Calib_pred),
                     calib_model=calib_gam,
                     test_data=pres_abs_test)
    
    # library(ggplot2)
    # ggplot(out_list$test_data)+
    #  geom_point(aes(x=x,y=y,color=(Corr_pred_prob - Ens_pred)))+
    #  scale_color_distiller(type="div",limits=c(-1,1))+
    #  theme_bw()
    
    ##Writes output to disk.
    out_name <- gsub(".Rdata","_calib.Rdata",model_file,fixed=TRUE)
    out_path <- paste(proj_dir,out_dir,out_name,sep="")
    saveRDS(out_list,file=out_path)
    
    ##Uploads to S3
    upl_string <- paste(aws_path,"aws s3 cp ",out_path," ",paste("s3://sdmdata/models_calibrated/",out_name,sep=""),sep="")
    system(upl_string,wait=TRUE)
    
    ##Removes local files.
    model_path <- paste(proj_dir,model_dir,model_file,sep="")
    rm_string1 <- paste("rm",model_path)
    rm_string2 <- paste("rm",out_path)
    system(rm_string1,wait=TRUE)
    system(rm_string2,wait=TRUE)

    ##Logs completion
    cat(paste("Calibration fit complete for ",spp[i],"(",i,"of",length(spp),") on",
              Sys.time(),"\n"),file="./scratch/sdm_progress.log",append=TRUE)
    
    # #####Creates calibrated rasters####
    # print(paste("Creating calibrated rasters..."))
    # 
    # ##Brings in the original raster prediction.
    # rast_name <- paste(gsub(" ","_",spp[i]),"_mosaic.tif",sep="")
    # rast_path <- paste(raw_pred_path,rast_name,sep="")
    # rast_pred <- raster(rast_path)
    # test_ext <- extent(-877449,-682475,1732606,1916111)
    # XY_stack <- brick(paste(rast_pred_path,"WNA_XY_albers_crop.tif",sep=""))
    # XY_stack <- crop(XY_stack,test_ext)
    # names(XY_stack) <- c('PCO_XSC', 'PCO_YSC')
    # NAvalue(XY_stack) <- 65535
    # 
    # pred_logit <- calc(rast_pred,fun=logit,
    #                    filename=paste(calib_pred_path,gsub(".tif","_logit.tif",rast_name,fixed=TRUE),sep=""),
    #                    overwrite=TRUE,progress="text")
    # pred_logit <- crop(pred_logit,test_ext)
    # pred_brick <- brick(XY_stack$PCO_XSC,XY_stack$PCO_YSC,pred_logit,overwrite=TRUE,
    #                     filename=paste(calib_pred_path,gsub(" ","_",spp[i]),"_corr_brick.grd",sep=""))
    # names(pred_brick) <- c("PCO_XSC","PCO_YSC","pred_logit")
    # corr_pred <- raster::predict(pred_brick, model=corr_gam, progress="text",ext=test_ext,
    #                              overwrite=TRUE,filename=paste(calib_pred_path,gsub(" ","_",spp[i]),"_corrected.grd",sep=""))
    # 
    # corr_out <- calc(corr_pred,fun=function(x){x*10000},datatype="INT2U",
    #                  overwrite=TRUE,filename=paste(calib_pred_path,gsub(" ","_",spp[i]),"_corrected.tif",sep=""))
    # calib_brick <- brick(corr_pred,corr_pred,overwrite=TRUE,
    #                      filename=paste(calib_pred_path,gsub(" ","_",spp[i]),"_calib_brick.grd",sep=""))
    # calib_brick[[2]] <- 0.04
    # names(calib_brick) <- c("Corr_pred","plotareaha")
    # Calib_pred <- raster::predict(calib_brick, model=calib_gam, type="response", progress="text",
    #                               overwrite=TRUE,filename=paste(calib_pred_path,gsub(" ","_",spp[i]),"_calib.grd",sep=""))
    # outname <- paste(proj_dir,calib_pred_path,gsub(" ","_",spp[i]),"_calibrated.tif",sep="")
    # Calib_out <- calc(Calib_pred,fun=function(x){x*10000},datatype="INT2U", overwrite=TRUE,
    #                   filename=outname)
    # warpname <- gsub("_calibrated.tif","_calibrated_webm.tif",outname,fixed=TRUE)
    # 
    # ##Creates cloud-optimized geotiff in web mercator projection.
    # cogname <- gsub("_calibrated.tif","_calibrated_webmcog.tif",outname,fixed=TRUE)
    # gdalw_call <- paste(gdal_path,"gdalwarp ",outname," ",warpname," -overwrite -t_srs EPSG:3857 -r bilinear -co TILED=YES -co COMPRESS=DEFLATE",
    #                     sep="")
    # system(gdalw_call,wait=TRUE)
    # gdala_call <- paste(gdal_path,"gdaladdo -r average ",warpname," 2 4 8 16 32",
    #                     sep="")
    # system(gdala_call,wait=TRUE)
    # gdalt_call <- paste(gdal_path,"gdal_translate ",warpname," ",cogname, " -ot UInt16 -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=DEFLATE -co BLOCKXSIZE=512 -co BLOCKYSIZE=512 --config GDAL_TIFF_OVR_BLOCKSIZE 512",
    #                      sep="")
    # system(gdalt_call,wait=TRUE)
    # 
    # ##Uploads to Amazon S3
    # rast_upl_string1 <- paste(aws_path, "aws s3 cp ",cogname,
    #                          paste(" s3://bloomfindersdm/cog/",
    #                                gsub("_mosaic.tif","_calibcog.tif",rast_name),sep=""))
    # system(rast_upl_string1,wait=TRUE)
    # rast_upl_string2 <- paste(aws_path,"aws s3 cp ",outname,
    #                           paste(" s3://sdmdata/WNA_mosaic/",
    #                                 gsub("_mosaic.tif","_calib.tif",rast_name),sep=""))
    # system(rast_upl_string2,wait=TRUE)
    # 
    # ##Removes intermediate files.
    # print(paste("Removing intermediate files..."))
    # rm_list1 <- list.files(raw_pred_path,pattern=gsub(" ","_",spp[i]),
    #                        full.names=TRUE)
    # file.remove(rm_list1)
    # rm_list2 <- list.files(calib_pred_path,pattern=gsub(" ","_",spp[i]),
    #                        full.names=TRUE)
    # file.remove(rm_list2)
    # 
    # cat(paste("Calibration prediction complete for ",spp[i],"(",i,"of",length(spp),") on",
    #             Sys.time(),"\n"),file="./scratch/sdm_progress.log",append=TRUE)
    
  }
  (out_list)
}
stopCluster(cl)
