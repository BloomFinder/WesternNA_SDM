####Fits SDMs to occupancy data.####
library(sdm)
library(readr)
library(dplyr)
library(foreach)
library(doParallel)
library(gdalUtils)
library(doRNG)
library(raster)
#devtools::install_github("wrathematics/openblasctl")
#library(openblasctl)

proj_dir <- "/Users/Ian/code/WesternNA_SDM"
#aws_path <- "~/.local/bin/"
aws_path <- "/Users/ian/miniconda2/bin/"
setwd(proj_dir)

#mnames <- names(sdm::getmethodNames())
#if(!("svm4" %in% mnames)){
#  source(paste(proj_dir,"/code/svm4.R",sep=""))
#  sdm::add(methodInfo,w='sdm')
#}
#if(!("gbmstep3" %in% mnames)){
#  source(paste(proj_dir,"/code/gbmstep.R",sep=""))
#  sdm::add(methodInfo,w='sdm')
#}

# ##Downloads test raster data from Amazon S3 if it doesn't already exist.
# if(!file.exists("./data/SDM_tiles_PNW.tar.gz")){
#   aws_dl1 <- paste(aws_path,"aws s3 cp s3://sdmdata/predictors/SDM_tiles_PNW.tar.gz ./data/SDM_tiles_PNW.tar.gz",sep="")
#   system(paste("cd",proj_dir,"&&",aws_dl1))
#   tar_dl1 <- "tar -xf ./data/SDM_tiles_PNW.tar.gz -C ./data/"
#   system(paste("cd",proj_dir,"&&",tar_dl1))
# }
# 
# #Downloads full raster data from Amazon S3 if it doesn't already exist.
# if(!file.exists("./data/SDM_tiles_WNA.tar.gz")){
#   aws_dl2 <- paste(aws_path,"aws s3 cp s3://sdmdata/predictors/SDM_tiles_WNA.tar.gz ./data/SDM_tiles_WNA.tar.gz",sep="")
#   system(paste("cd",proj_dir,"&&",aws_dl2),wait=TRUE)
#   tar_dl2 <- "tar -xf ./data/SDM_tiles_WNA.tar.gz -C ./data/"
#   system(paste("cd",proj_dir,"&&",tar_dl2),wait=TRUE)
# }

##Downloads point data from Amazon S3 if it doesn't already exist.
if(!file.exists("./data/occurences_final_3_1_2018.tar.gz")){
  aws_dl3 <- paste(aws_path,"aws s3 cp s3://sdmdata/occurences/occurences_final_3_1_2018.tar.gz ./data/occurences_final_3_1_2018.tar.gz",sep="")
  system(paste("cd",proj_dir,"&&",aws_dl3),wait=TRUE)
  tar_dl3 <- "tar -xf ./data/occurences_final_3_1_2018.tar.gz -C ./data/"
  system(paste("cd",proj_dir,"&&",tar_dl3),wait=TRUE)
}

##Downloads glacier random point data from Amazon S3 if it doesn't already exist.
if(!file.exists("./data/Randolph_glacier_random_points_attrib.tar.gz")){
  aws_dl4 <- paste(aws_path,"aws s3 cp s3://sdmdata/occurences/Randolph_glacier_random_points_attrib.tar.gz ./data/Randolph_glacier_random_points_attrib.tar.gz",sep="")
  system(paste("cd",proj_dir,"&&",aws_dl4),wait=TRUE)
  tar_dl4 <- "tar -xf ./data/Randolph_glacier_random_points_attrib.tar.gz -C ./data/"
  system(paste("cd",proj_dir,"&&",tar_dl4),wait=TRUE)
}

##Reads data in
spd <- read_csv("./data/occurences_final_3_1_2018.csv")
glac <- read_csv("./data/Randolph_glacier_random_points_attrib.csv")
glac <- glac[complete.cases(glac),]

##Species list for analysis
test_spp <- unique(spd$species)
#test_spp <- sample(unique(spd$species),size=50,replace=FALSE)
test_spp <- c("Aquilegia formosa",
               "Ranunculus adoneus",
               "Mimulus guttatus",
               "Vicia americana",
               "Chamerion angustifolium",
               "Maianthemum stellatum",
               "Phacelia heterophylla",
               "Sedum stenopetalum",
               "Ipomopsis aggregata",
               "Claytonia lanceolata",
               "Rudbeckia occidentalis",
               "Veratrum californicum",
               "Agoseris aurantiaca",
               "Sedum lanceolatum",
               "Xerophyllum tenax")

## Model fitting for focal species.
set.seed(38)

cl <- makeCluster(3)
registerDoParallel(cl)
overwrite <- TRUE

all_stats <- foreach(i=1:length(test_spp),.packages=c("dplyr"),
                     .combine="rbind") %dorng% {
                       
                       setwd(proj_dir)
                       library(sdm)
                       #library(openblasctl)
                       
                       mnames <- names(sdm::getmethodNames())
                       if(!("svm4" %in% mnames)){
                         source(paste(proj_dir,"/code/svm4.R",sep=""))
                         sdm::add(methodInfo,w='sdm')
                       }
                       if(!("gbmstep3" %in% mnames)){
                         source(paste(proj_dir,"/code/gbmstep.R",sep=""))
                         sdm::add(methodInfo,w='sdm')
                       }
                       
                       ##Prevents multithreaded linear algebra library from parallelizing.
                       #openblas_set_num_threads(1)
                       
                       ##Checks if model file already exists on AWS.
                       fit_file <- paste(proj_dir,"/scratch/models/sdm_",gsub(" ","_",test_spp[i]),".Rdata",sep="")
                       file_exists_aws <- system(paste(aws_path,"aws s3 ls s3://sdmdata/models/sdm_",gsub(" ","_",test_spp[i]),".Rdata",sep=""))
                       
                       if(file_exists_aws == 0 & overwrite==FALSE){
                         cat(paste("Output file",paste(fit_file),"exists, skipping...\n"),
                             file="./scratch/sdm_progress.log",append=TRUE)
                       }else{
                         cat(paste("Fitting model for ",test_spp[i],"(",i,"of",length(test_spp),") on",
                                   Sys.time(),"\n"),file="./scratch/sdm_progress.log",append=TRUE)
                         print(paste("Fitting model for ",test_spp[i],"(",i,"of",length(test_spp),")"))
                         print(paste("Preparing data..."))
                         tr_pres_abs <- spd
                         colnames(tr_pres_abs)[2:3] <- c("X","Y")
                         tr_pres_abs$TR_PRES <- as.numeric(tr_pres_abs$species==test_spp[i])
                         
                         tr_pres <- as.data.frame(tr_pres_abs[tr_pres_abs$TR_PRES==1,c(2,3,18:ncol(tr_pres_abs))])
                         
                         tr_abs <- tr_pres_abs[tr_pres_abs$TR_PRES==0,c(2,3,18:ncol(tr_pres_abs))]
                         tr_abs <- tr_abs[-which(tr_abs$loc %in% unique(tr_pres$loc)),]
                         tr_abs <- as.data.frame(sample_n(tr_abs,size=max((nrow(tr_pres) * 5),12000)))
                         
                         ##Adds a fixed proportion of random absences on glaciers.
                         n_abs <- nrow(tr_abs)
                         gl_abs <- sample_n(glac,size=n_abs*0.05,weight=glac$PT_DENS)
                         gl_abs <- gl_abs[,-c(3,36)]
                         gl_abs$location <- NA
                         gl_abs$loc <- NA
                         gl_abs$TR_PRES <- 0
                         tr_abs <- rbind(tr_abs,gl_abs)
                         
                         tr_abs <- tr_abs[,-which(colnames(tr_abs)=="loc")]
                         tr_pres <- tr_pres[,-which(colnames(tr_pres)=="loc")]
                         tr_pres_abs2 <- rbind(tr_abs,tr_pres)
                         
                         ##Creates a test dataset which takes data from a subset of 90km blocks
                         ##where a species is recorded as present, and another set where it is absent.
                         blocks <- unique(tr_pres_abs2$location)
                         blocks_pres <- unique(tr_pres_abs2$location[tr_pres_abs2$TR_PRES==1])
                         blocks_abs <- blocks[!(blocks %in% blocks_pres)]
                         n_valid_pblocks <- round(length(blocks_pres)*0.1)
                         sample_pblocks <- sample(blocks_pres,size=n_valid_pblocks+1,replace=FALSE)
                         n_valid_ablocks <- round(length(blocks_abs)*0.1)
                         sample_ablocks <- sample(blocks_abs,size=n_valid_ablocks+1,replace=FALSE)
                         valid_blocks <- c(sample_pblocks,sample_ablocks)
                         
                         test_data <- tr_pres_abs2[tr_pres_abs2$location %in% valid_blocks,-which(colnames(tr_pres_abs2) == "location")]
                         train_data <- tr_pres_abs2[!(tr_pres_abs2$location %in% valid_blocks),-which(colnames(tr_pres_abs2) == "location")]
                         all_data <- rbind(train_data,test_data)
                         
                         all_data$train_test <- c(rep("train",nrow(train_data)),rep("test",nrow(test_data)))
                         
                         #par(mfrow=c(1,2))
                         #plot(train_data$X,train_data$Y,col=train_data$TR_PRES + 1,
                         #     main=paste(test_spp[i],"Training Data"))
                         #plot(test_data$X,test_data$Y,col=test_data$TR_PRES + 1,
                         #     main=paste(test_spp[i],"Test Data"))
                         
                         
                         tr_sdmd <- sdmData(TR_PRES ~ . + f(PCT_EFW) + f(PCT_ECO) + f(PSL_TUS) + f(PSL_TWL),
                                            train=train_data,test=test_data)
                         print(tr_sdmd)
                         
                         ##Fits the models with heldout data to measure performance.
                         print(paste("Initial fitting..."))
                         
                         tr_sdm1 <- sdm(TR_PRES ~ PLC_TRE + PLC_HRB + PCM_CMD + PCM_TD + PCM_PAS + PCM_DD5 + PCM_MAP + PSL_BDR + PSL_SND + PSL_CAR + PSL_PHO + PCL_SE1 + PCL_SE2 + PCL_MRA + PSW_DIS + PTP_RLV + PTP_WET,
                                        data=tr_sdmd,methods=c("gbmstep3","svm","maxent"),
                                        var.selection=FALSE,modelSettings=list(gbmstep3=list(learning.rate=0.02,
                                                                                              n.trees=200,
                                                                                              n.cores=1)),
                                                                               maxent=list(beta=2),
                                                                               svm=list(epsilon=5))
                         
                         print(tr_sdm1)
                         
                         ##Extracts test data performance stats from fit models.
                         gbm_stats1 <- tr_sdm1@models$TR_PRES$gbmstep3$'1'@evaluation$test.indep@threshold_based
                         gbm_stats2 <- tr_sdm1@models$TR_PRES$gbmstep3$'1'@evaluation$test.indep@statistics
                         svm_stats1 <- tr_sdm1@models$TR_PRES$svm$'2'@evaluation$test.indep@threshold_based
                         svm_stats2 <- tr_sdm1@models$TR_PRES$svm$'2'@evaluation$test.indep@statistics
                         mxt_stats1 <- tr_sdm1@models$TR_PRES$maxent$'3'@evaluation$test.indep@threshold_based
                         mxt_stats2 <- tr_sdm1@models$TR_PRES$maxent$'3'@evaluation$test.indep@statistics
                         tr_sdm_stats <- data.frame(rbind(c(gbm_stats1,gbm_stats2),
                                                          c(svm_stats1,svm_stats2),
                                                          c(mxt_stats1,mxt_stats2)))
                         tr_sdm_stats$species <- test_spp[i]
                         tr_sdm_stats$method <- c("gbmstep3","svm","maxent")
                         tr_sdm_stats$ensemble_weight <- (as.numeric(tr_sdm_stats$AUC)-0.5)/0.5
                         
                         ##Re-fits models with the full dataset.
                         print(paste("Final fitting..."))
                         tr_sdmd2 <- sdmData(TR_PRES ~ . + f(PCT_EFW) + f(PCT_ECO) + f(PSL_TUS) + f(PSL_TWL),
                                             train=all_data)
                         sdm_all <- sdm(TR_PRES ~ PLC_TRE + PLC_HRB + PCM_CMD + PCM_TD + PCM_PAS + PCM_DD5 + PCM_MAP + PSL_BDR + PSL_SND + PSL_CAR + PSL_PHO + PCL_SE1 + PCL_SE2 + PCL_MRA + PSW_DIS + PTP_RLV + PTP_WET,
                                        data=tr_sdmd2,methods=c("gbmstep3","svm","maxent"),
                                        var.selection=FALSE,modelSettings=list(gbmstep3=list(learning.rate=0.02,
                                                                                             n.trees=200,
                                                                                             n.cores=1)),
                                                                              maxent=list(beta=4),
                                                                              svm=list(epsilon=8))
                         print(sdm_all)
                         out <- list(list(data=all_data,final_models=sdm_all,stats=tr_sdm_stats))
                         names(out) <- test_spp[i]
                         saveRDS(out,file=fit_file)
                         
                         
                         cat(paste("Model object written to",fit_file,"on",
                                   Sys.time(),"\n"),file="./scratch/sdm_progress.log",
                             append=TRUE)
                         
                         
                         ##Copies model file to S3.
                         #cp_string <- paste(aws_path,"aws s3 cp ",fit_file,
                         #                   paste("s3://sdmdata/models/sdm_", gsub(" ","_",test_spp[i]), ".Rdata",sep=""))
                         #system(cp_string,wait=TRUE)
                         
                         ##Removes file locally if it was successfully uploaded.
                         #file_exists_aws2 <- system(paste(aws_path,"aws s3 ls s3://sdmdata/models/sdm_",gsub(" ","_",test_spp[i]),".Rdata",sep=""))
                         
                         #if(file_exists_aws2==0){
                         #  system(paste("rm",fit_file),wait=TRUE)
                         #}
                         (tr_sdm_stats)
                       }
                     }
stopCluster(cl)
AUC_mean <- mean(as.numeric(all_stats$AUC))
AUC_mean
all_stats_tbl <- as.tbl(all_stats)
save(all_stats_tbl,file="./output/sdm_results.Rdata")

# ####Creates spatial predictions using the weighted average of the models.####
# 
# tile_path <- paste(proj_dir,"/data/SDM_tiles_WNA/",sep="")
# model_path <-paste(proj_dir,"/scratch/sdm_fits/",sep="")
# out_path <- paste(proj_dir,"/output/PNW_tiles/",sep="")
# mosaic_path <- paste(proj_dir,"/output/PNW_mosaic/",sep="")
# log_path <- paste(proj_dir,"/scratch/sdm_progress.log",sep="")
# 
# pred_tiles <- list.files(tile_path, pattern=".tif$", full.names=TRUE)
# pred_names <- list.files(tile_path,pattern=".tif$", full.names=FALSE)
# 
# model_files <- list.files(model_path,pattern=".Rdata",full.names=TRUE)
# overwrite=TRUE
# 
# ##Sets up cluster.
# cl <- makeCluster(9)
# registerDoParallel(cl)
# 
# ##Raster predictions.
# foreach(i=1:length(test_spp),.packages=c("raster","sdm","gdalUtils","openblasctl")) %dorng% {
#   
#   openblas_set_num_threads(1)
#   
#   ##Checks to see if the mosaic exists on Amazon S3.
#   mos_exists_az <- system(paste(aws_path,"aws s3 ls s3://sdmdata/PNW_mosaic/",gsub(" ","_",test_spp[i]),"_mosaic.tif",sep=""),
#                           wait=TRUE)
#   if(mos_exists_az==0 & overwrite==FALSE){
#     cat(paste("Raster predictions for",spp,"(",i,"of",length(test_spp),"already exist in S3, skipping...\n"),
#         file=log_path,append=TRUE)
#   }else{
#     
#     ##downloads fit model
#     model_dl_cmd <- paste(aws_path,"aws s3 cp s3://sdmdata/models/sdm_", gsub(" ","_",test_spp[i]),".Rdata ",
#                           model_path,"sdm_",gsub(" ","_",test_spp[i]),".Rdata",sep="")
#     system(model_dl_cmd)
#     
#     ##Loads fit models.
#     model_file <- paste(model_path,"sdm_",gsub(" ","_",test_spp[i]),".Rdata",sep="")
#     out <- readRDS(model_file)
#     models <- out[[1]]$final_models
#     stats <- out[[1]]$stats
#     spp <- names(out)
#     remove(out)
#     
#     cat(paste("Raster predictions for",spp,"(",i,"of",length(test_spp),") started on",
#               Sys.time(),"\n"),file=log_path,append=TRUE)
#     
#     ##Creates output directory for each species if it doesn't exist.
#     spp_dir <- paste(out_path,gsub(" ","_",spp),"/",sep="")
#     if(!dir.exists(spp_dir)){dir.create(spp_dir)}
#     
#     preds <- list()
#     for(j in 1:length(pred_tiles)){
#       outfile <- paste(spp_dir,"sdm_tile_",j,
#                        "_",gsub(" ","_",spp),".img",sep="")
#       if(file.exists(outfile) & overwrite==FALSE){
#         print(paste("File exists, skipping..."))
#       }else{
#         pred_tile <- readAll(brick(pred_tiles[j]))
#         names(pred_tile) <- c('PCL_MAN', 'PCL_SE1', 'PCL_SE2', 'PCL_SE3', 'PCM_BFP',
#                               'PCM_CMD', 'PCM_DD5', 'PCM_MAP', 'PCM_PAS', 'PCM_TD',
#                               'PCT_ECO', 'PCT_EFW', 'PLC_HRB', 'PLC_TRE', 'PLC_URB',
#                               'PSL_BDR', 'PSL_CAR','PSL_PHO', 'PSL_SND', 'PSL_TUS',
#                               'PSL_TWL', 'PTP_ELV', 'PTP_RLV', 'PTP_SLP', 'PTP_WET',
#                               'PTP_ASP', 'PTP_SOL', 'PCL_MRA', 'PSW_DIS', 'PSW_OCC',
#                               'PCO_XSC', 'PCO_YSC')
#         pred <- try(ensemble(models,newdata=pred_tile,
#                              setting=list(method='weighted',
#                                           weights=stats$ensemble_weight[models@run.info$success]),
#                              filename=outfile,overwrite=TRUE,progress='text'))
#         preds[[i]] <- pred
#       }
#       
#     }
#     
#     ##Merges output tiles to single raster.
#     tiles <- list.files(spp_dir,pattern=".img$")
#     setwd(spp_dir)
#     mosaic_rasters(gdalfile=tiles,dst_dataset=paste(mosaic_path,gsub(" ","_",test_spp[i]),"_mosaic.tif",sep=""),
#                    verbose=TRUE)
#     
#     ##Copies raster mosaic to S3 bucket.
#     cp_string <- paste(aws_path,"aws s3 cp ",paste(mosaic_path,gsub(" ","_",test_spp[i]),"_mosaic.tif",sep="")," ",
#                        paste("s3://sdmdata/WNA_mosaic/", gsub(" ","_",test_spp[i]), "_mosaic.tif",sep=""),sep="")
#     system(cp_string,wait=TRUE)
#     
#     ##Removes tiles, mosaic, and model to save disk space.
#     all_tile_files <- list.files(spp_dir)
#     rm_string <- paste("rm",paste(paste(spp_dir,all_tile_files,sep=""),collapse=" "))
#     system(rm_string)
#     
#     rm_string_m <- paste("rm ",mosaic_path,gsub(" ","_",test_spp[i]),"_mosaic.tif",sep="")
#     system(rm_string_m)
#     
#     rm_string_mod <- paste("rm",model_file)
#     system(rm_string_mod)
#     
#     cat(paste("Raster predictions for",spp,"(",i,"of",length(test_spp),") completed on",
#               Sys.time(),"\n"),file=log_path,append=TRUE)
#   }
# }
# stopCluster(cl)
