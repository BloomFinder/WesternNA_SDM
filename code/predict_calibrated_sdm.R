####Creates raster predictions of calibrated SDM models.####
library(sdm)
library(readr)
library(dplyr)
library(foreach)
library(doParallel)
library(gdalUtils)
install.packages("doRNG")
library(doRNG)
library(raster)
#devtools::install_github("wrathematics/openblasctl")
library(openblasctl)

proj_dir <- "/home/rstudio/WesternNA_SDM"
setwd(proj_dir)

##Adds custom svm and gbm functions to sdm package.
#mnames <- names(sdm::getmethodNames())
#if(!("svm4" %in% mnames)){
#  source(paste(proj_dir,"/code/svm4.R",sep=""))
#  sdm::add(methodInfo,w='sdm')
#}
#if(!("gbmstep3" %in% mnames)){
#  source(paste(proj_dir,"/code/gbmstep.R",sep=""))
#  sdm::add(methodInfo,w='sdm')
#}


####Creates spatial predictions using the weighted average of the models.####

tile_path <- paste(proj_dir,"/data/SDM_tiles_WNA/",sep="")
raw_model_path <-paste(proj_dir,"/scratch/models/",sep="")
calib_model_path <-paste(proj_dir,"/scratch/models_calib/",sep="")
out_path <- paste(proj_dir,"/output/WNA_tiles/",sep="")
mosaic_path <- paste(proj_dir,"/output/WNA_mosaic_calib/",sep="")
log_path <- paste(proj_dir,"/scratch/sdm_progress.log",sep="")
aws_path <- "/home/rstudio/.local/bin/"
#aws_path <- "~/miniconda2/bin/"
gdal_path <- "/usr/bin/"
#gdal_path <- "/Library/Frameworks/GDAL.framework/Programs/"

#Downloads full raster data from Amazon S3 if it doesn't already exist.
if(!file.exists("./data/SDM_tiles_WNA.tar.gz")){
  aws_dl2 <- "~/.local/bin/aws s3 cp s3://sdmdata/predictors/SDM_tiles_WNA.tar.gz ./data/SDM_tiles_WNA.tar.gz"
  system(paste("cd",proj_dir,"&&",aws_dl2),wait=TRUE)
  tar_dl2 <- "tar -xf ./data/SDM_tiles_WNA.tar.gz -C ./data/"
  system(paste("cd",proj_dir,"&&",tar_dl2),wait=TRUE)
}

##Downloads test raster data from Amazon S3 if it doesn't already exist.
if(!file.exists("./data/SDM_tiles_PNW.tar.gz")){
  aws_dl1 <- paste(aws_path,"aws s3 cp s3://sdmdata/predictors/SDM_tiles_PNW.tar.gz ./data/SDM_tiles_PNW.tar.gz",sep="")
  system(paste("cd",proj_dir,"&&",aws_dl1))
  tar_dl1 <- "tar -xf ./data/SDM_tiles_PNW.tar.gz -C ./data/"
  system(paste("cd",proj_dir,"&&",tar_dl1),wait=TRUE)
}

##Downloads point data from Amazon S3 if it doesn't already exist.
if(!file.exists("./data/occurences_final_3_1_2018.csv")){
  aws_dl3 <- "~/.local/bin/aws s3 cp s3://sdmdata/occurences/occurences_final_3_1_2018.tar.gz ./data/occurences_final_3_1_2018.tar.gz"
  system(paste("cd",proj_dir,"&&",aws_dl3),wait=TRUE)
  tar_dl3 <- "tar -xf ./data/occurences_final_3_1_2018.tar.gz -C ./data/"
  system(paste("cd",proj_dir,"&&",tar_dl3),wait=TRUE)
}

spd <- read_csv("./data/occurences_final_3_1_2018.csv")
all_spp <- unique(spd$species)

exclude_spp <-  c("Agastache urticifolia",    "Angelica arguta",          "Antennaria rosea",        
                  "Castilleja miniata",       "Cistanthe umbellata",      "Erysimum capitatum",      
                  "Ipomopsis aggregata",      "Mimulus guttatus",         "Parnassia fimbriata",     
                  "Senecio integerrimus",     "Stellaria longipes",       "Symphyotrichum foliaceum",
                  "Veronica serpyllifolia",   "Vicia americana",          "Viola macloskeyi",
                  "Valeriana californica",    "Fallugia paradoxa",        "Castilleja austromontana",
                  "Linaria vulgaris")

test_spp <- all_spp[!(all_spp %in% exclude_spp)]
rm(spd)

pred_tiles <- list.files(tile_path, pattern=".tif$", full.names=TRUE)
pred_names <- list.files(tile_path,pattern=".tif$", full.names=FALSE)

model_files <- list.files(raw_model_path,pattern=".Rdata",full.names=TRUE)
overwrite=TRUE

##Sets up cluster.
cl <- makeCluster(70)
registerDoParallel(cl)

##Raster predictions.
foreach (i=1:length(test_spp),.packages=c("raster","sdm","gdalUtils","openblasctl"),
         .errorhandling='remove') %dopar% {
  
  openblas_set_num_threads(1)
  
  ##Adds custom svm and gbm functions to sdm package.
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
  
  ##Checks to see if the mosaic already exists on Amazon S3.
  mos_exists_az <- system(paste(aws_path,"aws s3 ls s3://bloomfindersdm/cog/",gsub(" ","_",test_spp[i]),"_calib_webmcog.tif",sep=""),
                          wait=TRUE)
  if(mos_exists_az==0 & overwrite==FALSE){
    cat(paste("Raster predictions for",test_spp[i],"(",i,"of",length(test_spp),"already exist in S3, skipping...\n"),
        file=log_path,append=TRUE)
    next
  }
    
    ##downloads fit model
    model_dl_cmd <- paste(aws_path,"aws s3 cp s3://sdmdata/models/sdm_", gsub(" ","_",test_spp[i]),".Rdata ",
                          raw_model_path,"sdm_",gsub(" ","_",test_spp[i]),".Rdata",sep="")
    system(model_dl_cmd)
    
    cal_model_dl_cmd <- paste(aws_path,"aws s3 cp s3://sdmdata/models_calibrated/sdm_", gsub(" ","_",test_spp[i]),"_calib.Rdata ",
                          calib_model_path,"sdm_",gsub(" ","_",test_spp[i]),"_calib.Rdata",sep="")
    system(cal_model_dl_cmd)
    
    ##Loads fit models.
    model_file <- paste(raw_model_path,"sdm_",gsub(" ","_",test_spp[i]),".Rdata",sep="")
    calib_model_file <- paste(calib_model_path,"sdm_",gsub(" ","_",test_spp[i]),"_calib.Rdata",sep="")
    out <- readRDS(model_file)
    out_calib <- readRDS(calib_model_file)
    models <- out[[1]]$final_models
    stats <- out[[1]]$stats
    spp <- names(out)
    remove(out)
    
    corr_model <- out_calib$corr_model
    calib_model <- out_calib$calib_model
    
    cat(paste("Raster predictions for",spp,"(",i,"of",length(test_spp),") started on",
              Sys.time(),"\n"),file=log_path,append=TRUE)
    
    ##Creates output directory for each species if it doesn't exist.
    spp_dir <- paste(out_path,gsub(" ","_",spp),"/",sep="")
    if(!dir.exists(spp_dir)){dir.create(spp_dir)}
      
    for(j in 1:length(pred_tiles)) {
    
      outfile <- paste(spp_dir,"sdm_tile_",j,
                       "_",gsub(" ","_",spp),".img",sep="")
      if(file.exists(outfile)){
        print(paste("File exists, skipping..."))
        next
      }
      pred_tile <- readAll(brick(pred_tiles[j]))
      
      if(is.na(minValue(pred_tile[[1]]))){
        next
      }
      
      names(pred_tile) <- c('PCL_MAN', 'PCL_SE1', 'PCL_SE2', 'PCL_SE3', 'PCM_BFP',
                            'PCM_CMD', 'PCM_DD5', 'PCM_MAP', 'PCM_PAS', 'PCM_TD',
                            'PCT_ECO', 'PCT_EFW', 'PLC_HRB', 'PLC_TRE', 'PLC_URB',
                            'PSL_BDR', 'PSL_CAR','PSL_PHO', 'PSL_SND', 'PSL_TUS',
                            'PSL_TWL', 'PTP_ELV', 'PTP_RLV', 'PTP_SLP', 'PTP_WET',
                            'PTP_ASP', 'PTP_SOL', 'PCL_MRA', 'PSW_DIS', 'PSW_OCC',
                            'PCO_XSC', 'PCO_YSC')
      
      raw_pred <- try(ensemble(models,newdata=pred_tile,
                           setting=list(method='weighted',
                                        weights=stats$ensemble_weight[models@run.info$success]),
                           filename=paste(spp_dir,"sdmtemp.img",sep=""),
                           overwrite=TRUE,progress='text'))
      
      raw_pred[raw_pred < (1e-4)] <- 1e-4
      raw_pred[raw_pred > (1 - 1e-4)] <- (1 - 1e-4)
      
      XSC <- pred_tile$PCO_XSC
      YSC <- pred_tile$PCO_YSC
      rm(pred_tile)
      
      pred_logit <- calc(raw_pred,fun=logit,progress='text')
      rm(raw_pred)
      #pred_logit[pred_logit < -4] <- -4
      corr_brick <- brick(XSC,YSC,pred_logit)
      names(corr_brick) <- c("PCO_XSC","PCO_YSC","pred_logit")
      rm(XSC,YSC)
      
      corr_pred <- raster::predict(corr_brick, model=corr_model, progress='text')
      #corr_pred[corr_pred < -10] <- -10
      #corr_pred_prob <- calc(corr_pred,fun=inv_logit,progress='text')
      rm(corr_brick)

      calib_brick <- brick(corr_pred,corr_pred)
      calib_brick[[2]] <- 0.04
      names(calib_brick) <- c("Corr_pred",names(calib_model$coefficients)[2])
      Calib_pred <- raster::predict(calib_brick, model=calib_model, progress="text")
      rm(calib_brick)
      Calib_pred[Calib_pred < -32] <- -32
      Calib_pred[Calib_pred > 32] <- 32
      #Calib_pred_prob <- calc(Calib_pred,fun=inv_logit,progress='text')
      Calib_out <- calc(Calib_pred,fun=function(x){x*1000},datatype="INT2S", overwrite=TRUE,
                        filename=outfile)
    }
    
    ##Removes model files.
    rm_string_mod <- paste("rm",paste(model_file,calib_model_file))
    system(rm_string_mod,wait=TRUE)
    
    ##Removes temporary raster.
    rm_string_temp <- paste("rm ",paste(spp_dir,"sdmtemp.img ",sep=""), paste(spp_dir,"sdmtemp.img.aux.xml",sep=""),sep="")
    system(rm_string_temp,wait=TRUE)
    
    ##Merges output tiles to single raster.
    tiles <- list.files(spp_dir,pattern=".img$")
    setwd(spp_dir)
    outname <- paste(mosaic_path,gsub(" ","_",test_spp[i]),"_calib.tif",sep="")
    mosaic_rasters(gdalfile=tiles,dst_dataset=outname,
                   verbose=TRUE)
    
    ##Removes tiles
    all_tile_files <- list.files(spp_dir)
    rm_string <- paste("rm",paste(paste(spp_dir,all_tile_files,sep=""),collapse=" "))
    system(rm_string,wait=TRUE)
    
    ##Copies raster mosaic to S3 bucket.
    cp_string <- paste(aws_path,"aws s3 cp ",outname,
                       paste(" s3://sdmdata/WNA_mosaic/", gsub(" ","_",test_spp[i]), "_calib.tif",sep=""),sep="")
    system(cp_string,wait=TRUE)
    
    ##Creates a cloud-optimized geotiff in web mercator projection.
    warpname <- gsub("_calib.tif","_calib_webm.tif",outname,fixed=TRUE)
    cogname <- gsub("_calib.tif","_calib_webmcog.tif",outname,fixed=TRUE)
    
    gdalw_call <- paste(gdal_path,"gdalwarp ",outname," ",warpname," -overwrite -t_srs EPSG:3857 -r bilinear -co TILED=YES -co COMPRESS=DEFLATE -co PREDICTOR=2 -co BLOCKXSIZE=512 -co BLOCKYSIZE=512",
                        sep="")
    system(gdalw_call,wait=TRUE)
    gdala_call <- paste(gdal_path,"gdaladdo -r average ",warpname," 2 4 8 16 32 64 128 256 512 1024 2048",
                        sep="")
    system(gdala_call,wait=TRUE)
    gdalt_call <- paste(gdal_path,"gdal_translate ",warpname," ",cogname, " -ot Int16 -co TILED=YES -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=DEFLATE -co PREDICTOR=2 -co BLOCKXSIZE=512 -co BLOCKYSIZE=512 --config GDAL_TIFF_OVR_BLOCKSIZE 512",
                        sep="")
    system(gdalt_call,wait=TRUE)
    
    ##Copies web COG to S3
    cp_string <- paste(aws_path,"aws s3 cp ",cogname,
                       paste(" s3://bloomfindersdm/cog_wna/", gsub(" ","_",test_spp[i]), "_calib_webmcog.tif",sep=""),sep="")
    system(cp_string,wait=TRUE)
    
    ##Removes mosaic files.
    all_mos_files <- list.files(mosaic_path,pattern=gsub(" ","_",test_spp[i]))
    rm_mos_string <- paste("rm",paste(paste(mosaic_path,all_mos_files,sep=""),collapse=" "))
    system(rm_mos_string,wait=TRUE)
    
    cat(paste("Raster predictions for",spp,"(",i,"of",length(test_spp),") completed on",
              Sys.time(),"\n"),file=log_path,append=TRUE)
         }
stopCluster(cl)