# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# The Aim of this workflow is to create a spatial stratification based on subsampling elevation
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Diego and I just talked. I suggested also trying elevationally (and potentially latitudinally) biased subsampling (along the lines we did in GEB spatial prior paper). So drop 10, 25, 50, 75% of points starting from lowest or highest elevation and see how that delta AUC hopefully increases. Make sense, Cory?

# [1] Make a species.csv file where we filter 25% of the lowest level elevation.
# 
# Wr source the function get_point_by_elev_quantiles_for_this
source('/Users/diegoellis/Dropbox/20171107_Hummingbird/Check_AUC_subsamples/subset_points_by_elev.R')
require(raster)
points_dir = '/Users/diegoellis/Dropbox/Diego_Cory_Walter/Hummingbirds_AccuratePts_4_9_18/Inputs/Only_accurate_points/Points/'
speciesName = 'Heliangelus_exortis'
elevationstack_path = '/Users/diegoellis/Dropbox/Diego_Cory_Walter/Hummingbirds_AccuratePts_4_9_18/Inputs/Misc/'
desired_quantile_removal = .25
outdir = '/Users/diegoellis/Desktop/Hummingbirds/v39_elevation_subsampling/nsamples_25/Points/'

get_points_by_elev_quantiles(points_dir = points_dir,
                             speciesName = 'Heliangelus_exortis',
                             elevationstack_path =  elevationstack_path,
                             desired_quantile_removal = .25,
                             outdir = outdir)

if(desired_quantile_removal == .25){ sampleDir =  "/Users/diegoellis/Desktop/Hummingbirds/v39_elevation_subsampling/nsamples_25/" }
if(desired_quantile_removal == .50){ sampleDir =  "/Users/diegoellis/Desktop/Hummingbirds/v39_elevation_subsampling/nsamples_50/" }
if(desired_quantile_removal == .75){ sampleDir =  "/Users/diegoellis/Desktop/Hummingbirds/v39_elevation_subsampling/nsamples_75/" }
if(!file.exists(paste0(sampleDir, '/Points/'  ))){dir.create(paste0(sampleDir, '/Points/'  ))}
if(!file.exists(paste0(sampleDir, '/Shapefile/'  ))){dir.create(paste0(sampleDir, '/Shapefile/'  ))}
if(!file.exists(paste0(sampleDir, '/Raster/'  ))){dir.create(paste0(sampleDir, '/Raster/'  ))}

# [2] 
# Now make 10 copies of this and prepare it for workflowMOL
# pres is a data frame containing lat long of species (in our case a speciesCSV filtered by elevation)

shp_indir = '/Users/diegoellis/Dropbox/Diego_Cory_Walter/Hummingbirds_AccuratePts_4_9_18/Inputs/Only_accurate_points/Shapefile/' 
raster_indir = '/Users/diegoellis/Dropbox/Diego_Cory_Walter/Hummingbirds_AccuratePts_4_9_18/Inputs/Only_accurate_points/CorrectedRasters/'
indir =  '/Users/diegoellis/Desktop/Hummingbirds/v39_elevation_subsampling/'
speciesName = 'Heliangelus_exortis'
desired_quantile_removal = .25

copy_10_replicates_points_shp_rasters(indir = indir, speciesName = speciesName, desired_quantile_removal = desired_quantile_removal,
                                      shp_indir = shp_indir, raster_indir = raster_indir)

# [3] Run workflowMOL
library(cmsdm)
library(MOLSDM)
library(trinaryMaps)
runName=paste0('Elev_Hummingbirds_elev_25_percent')
#objects_for_sdm <- 
startworkflowMOL(desired_quantile_removal = desired_quantile_removal ,runName = runName)

# [5] Make delta AIC plots:

# Geom_boxplot:
folders <- c(10, 25, 40, 60)
x <- data.frame(matrix(ncol = 4, nrow = 40)) # 40 
path_to_folders <- '/Users/diegoellis/Desktop/Hummingbirds/TEST/'
folders <- c(10, 25,40, 75)
speciesName = 'Heliangelus_exortis'
x <- data.frame(matrix(ncol = 5, nrow = 40)) # 40 
folders <- c(10, 25,40, 50, 75)
speciesName = 'Oreotrochilus_leucopleurus'
x <- data.frame(matrix(ncol = 5, nrow = 50)) # 40 
names(x) <- c('Species','Number_samples','Model_name' ,'Delta_AIC_EM', 'Delta_AIC_NOEM')
species_I_want_list <- list()
for(nsample in folders){
  i <- which(folders %in% nsample)
  print(nsample)
  # nsample <- 10
  stats_folder <- paste0(path_to_folders, 'nsamples_', nsample, '/Elevation_influence_sample_size_', nsample, '_outputs/Statistics/') 
  # Make delta AUC for each sample size
  # file_names <- list.files('/Users/diegoellis/Desktop/Hummingbirds/v39_subsampling_points/nsamples_25/Elevation_influence_sample_size_25_outputs/Statistics/', full.names = T)
  # file.names = list.files(stats_folder, full.names = T)
  # for(i in 1:length(file_names)){
  species_I_want <- list.files(stats_folder, pattern = speciesName, full.names = T)
  species_I_want_list[[i]] <- species_I_want
# return(species_I_want_list)
  }
# big_data = do.call(rbind, species_I_want_list)
species_I_want <- unlist(species_I_want_list)

  for(i in 1:length(species_I_want)){
    # i <- 1
    a <- read.csv(species_I_want[i])
    if(nrow(a[a$elevationSet == 'elevation' & a$probSet != 1e-06
              ,]) !=0){ # If there are SDM with expert range maps. 
    max_auc_elev_EM <- max(a[a$elevationSet == 'elevation' & a$probSet != 1e-06
                             ,]$cv.mean.test.auc)
    max_auc_NO_elev_EM <- max(a[a$elevationSet == 'none'  & a$probSet != 1e-06,]$cv.mean.test.auc)
    delta_auc_EM <- max_auc_elev_EM - max_auc_NO_elev_EM
    x[i,]$Delta_AIC_EM <- delta_auc_EM
    }
    max_auc_elev_NOEM <- max(a[a$elevationSet == 'elevation' & a$probSet == 1e-06
                               ,]$cv.mean.test.auc)
    max_auc_NO_elev_NOEM <- max(a[a$elevationSet == 'none'  & a$probSet == 1e-06,]$cv.mean.test.auc)
    
    delta_auc_NOEM <- max_auc_elev_NOEM - max_auc_NO_elev_NOEM
    # x[i,]$Number_samples <- nsample  
    x[i,]$Model_name <- gsub('.csv', '',basename(species_I_want[i]))
    x[i,]$Delta_AIC_NOEM <- delta_auc_NOEM
    x[i,]$Species <- speciesName
  }
# Add sample size column:
x[1:10,]$Number_samples <- '10' 
x[11:20,]$Number_samples <- '25' 
x[21:30,]$Number_samples <- '40' 
x[31:40,]$Number_samples <- '50' 
x[41:50,]$Number_samples <- '75' 
require(gridExtra)
require(ggplot2)
EM <- ggplot(aes(y = log(Delta_AIC_EM), x = Number_samples), data = x) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +theme_bw()
NOEM <- ggplot(aes(y = log(Delta_AIC_NOEM), x = Number_samples), data = x) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +theme_bw()

p <- ggplot(x[c(1:20, 31:40),], aes(factor(Number_samples), log(Delta_AIC_NOEM))) +theme_bw() +
  xlab('Biased subsampling: Percentage of points removed\nelevationally from lowest to highest elevation')+ ylab(('AUC SDM elevation - AUC SDM no elevation'))+ggtitle(gsub('_', ' ', speciesName))+
# +xlab(('Number of spatially stratified points')) 
  theme(plot.title = element_text(size = 18, face = "italic"),
        element_line(colour = "black", size = 1, linetype = "solid"),
       axis.title.y = element_text(face = 'bold'),
       axis.title.x = element_text(face = 'bold'),
       axis.text = element_text(face = 'bold')) +
  scale_fill_discrete(name="Number of\nPoints")# face = 'bold'

p <- ggplot(x, aes(factor(Number_samples), (Delta_AIC_NOEM))) +theme_bw() +
  xlab('Biased subsampling: Percentage of points removed\nelevationally from lowest to highest elevation')+ ylab(('AUC SDM elevation - AUC SDM no elevation'))+ggtitle(gsub('_', ' ', speciesName))+
  # +xlab(('Number of spatially stratified points')) 
  theme(plot.title = element_text(size = 18, face = "italic"),
        element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(face = 'bold'),
        axis.title.x = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold')) +
  scale_fill_discrete(name="Number of\nPoints")# face = 'bold'

p + geom_violin() + geom_jitter(height = 0, width = 0.1)
p + geom_violin(aes(fill = factor(Number_samples)))+geom_jitter(height = 0, width = 0.1)
p + geom_boxplot(aes(fill = factor(Number_samples)))+geom_jitter(height = 0, width = 0.1)

# Add mean and SD:
# p <- ggplot(x, aes(x=Number_samples, y=Delta_AIC_NOEM)) + 
#   geom_violin(trim=FALSE)
# p + stat_summary(fun.data="mean_sdl", mult=1, 
#                  geom="crossbar", width=0.2 )
# p + stat_summary(fun.data=mean_sdl, mult=1, 
#                  geom="pointrange", color="red")
# p + stat_summary(fun.data=mean_sdl, mult=1, 
#                  geom="pointrange", color="red", notch = TRUE)

ggplot(x, aes(x=Number_samples, y=Delta_AIC_NOEM)) +
  geom_boxplot(
    notch = TRUE, aes(fil=factor(Number_samples)))+
  geom_jitter(height = 0, width = 0.1)
# p +geom_boxplot(size = 1, notch = TRUE) 
grid.arrange(EM, NOEM)


#mm = melt(df, id=c('id','factor.col'))
#ggplot(mm)+geom_boxplot(aes(x=paste(variable,factor.col,sep="_"), y=value))




#mean(x$Delta_AIC_EM)
#mean(x$Delta_AIC_NOEM)



# 
# 
# # Filter top 10 percentile of a column using dplyr
# library(data.table)
# range(pdf_df$max_elev)
# nrow(pdf_df)
# 
# 
# sub <- pdf_df[quantile(pdf_df$max_elev,.50) <= pdf_df$max_elev,] 
# 
# range(sub$max_elev)
# nrow(sub)
# 
# # a <- setDT(pdf_df)[,.SD[quantile(max_elev, 0.75) > max_elev]]
# 
# 
# 
# 
# 
# 
# # === === === === === === === === === === === === === === ===
# # Check sample size and AIC influence
# #
# #
# #
# # === === === === === === === === === === === === === === ===
# 
# message('Get the standard deviation of AIC for the 100 replicates for each sample size')
# 
# 
# message('Get the mean AIC for the  the 100 replicates for each sample size')
# 
# message('Make a barplot for this')
# # === === === === === === === === === === === === === === ===
# require(reshape)
# library(reshape2)
# require(ggplot2)
# 
# # Make a fake data frame to prepare for this
# x <- data.frame(matrix(ncol = 2, nrow = 300))
# names(x) <- c('Number_samples', 'AIC')
# x[,c('Number_samples')] <- rep(c(20,50,100),100 ) #c
# x[,c('AIC')] <- rnorm(300)
# x$Number_samples <- as.factor(x$Number_samples)
# head(x)
# 
# ggplot(aes(y = AIC, x = Number_samples, color = Number_samples), data = x) +
#   geom_boxplot(outlier.colour="black", outlier.shape=16,
#                outlier.size=2, notch=FALSE) +
#   theme_bw()#  + geom_jitter(height = 0, width = 0.1)
# 
# ggplot(aes(y = AIC, x = Number_samples, color = Number_samples), data = x) +
#   geom_boxplot(outlier.colour="black", outlier.shape=16,
#                outlier.size=2, notch=FALSE) +
#   theme_bw()  + geom_jitter(height = 0, width = 0.1)
# 
# ggplot(aes(y = AIC, x = Number_samples, color = Number_samples), data = x) + geom_violin() + geom_jitter(height = 0, width = 0.1) + theme_bw()
# 
# # ggplot(aes(y = AIC, x = Number_samples, color = Number_samples, fill = factor(Number_samples)), data = x) + geom_boxplot() + theme_bw() + geom_jitter(height = 0, width = 0.1)
# 
# geom_boxplot(outlier.colour="black", outlier.shape=16,
#              outlier.size=2, notch=FALSE)
# 
# # fill = AIC
# require(plyr)
# ddply(x, 'Number_samples', function(x){ mean(x$AIC, na.rm = TRUE)})
# 
# 
# df_long <- melt(df, id=c("id","factor.col"))
# 
# mm = melt(x, id=c('Number_samples','AIC'))
# 
# ggplot(mm)+geom_boxplot(aes(x=paste(variable,factor.col,sep="_"), y=value))
# 


# TO DO: DO WITH ALL THE POINTS! not only accurate points ! 