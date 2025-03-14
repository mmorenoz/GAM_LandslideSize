# The code illustrates the procedure described in the manuscript--
# Modeling the size of co-seismic landslides via data-driven models: the Kaikōura example
# by authors Mateo Moreno, Stefan Steger, Hakan Tanyas and Luigi Lombardo


# INITIAL SETTINGS --------------------------------------------------------

# clean environment
rm(list = ls())

# install package renv for reproductibility
# install.packages("renv")

# install packages with respective version
renv::restore()

# installing necessary libraries
list.packages = c("tidyverse", "sf", "mapview", "mgcv", "sperrorest",
                  "ggplot2", "reshape", "RColorBrewer")
vapply(list.packages, library, logical(1), character.only = TRUE, logical.return = TRUE, quietly = TRUE)
rm(list.packages)

# data loading
SU = sf::st_read("./dat/SU.gpkg")
SU = as.data.frame(SU)

# FUNCTIONS ---------------------------------------------------------------
# log transform
my.log = function(x){
  a = log(x)
  a[a == -Inf] = NA
  return(a)}

# To limit the road distance
my.trafo = function(x) {
  x$Roads_mean[x$Roads_mean > 500 ] = 500
  x$Roads_stdev[x$Roads_stdev > 500 ] = 500
  return(x)}

# Performance metric - RMSE
my.RMSE = function(m, o, na.rm){
  sqrt(mean((m - o)^2, na.rm = na.rm))
}

# Performance metric - MAE
my.MAE = function(m, o, na.rm){
  mean(abs(m - o), na.rm = na.rm)
}

# error function for sperrorest package
my.error = function (obs, pred){
  mmin <- function(x, ...) {
    if (length(x) == 0) {
      x <- Inf
    }
    return(min(x, ...))
  }
  mmax <- function(x, ...) {
    if (length(x) == 0) {
      x <- -Inf
    }
    return(max(x, ...))
  }
  if (is.logical(obs)) {
    obs <- factor(obs, levels = c("FALSE", "TRUE"))
  }
  if (is.logical(pred)) {
    pred <- factor(pred, levels = c("FALSE", "TRUE"))
  }
  if (any(is.na(pred))) {
    index_na <- which(pred %in% NA)
    obs <- obs[-index_na]
    pred <- pred[-index_na]
  }
  if (any(is.na(pred))) {
    index_na <- which(pred %in% NA)
    obs <- obs[-index_na]
    pred <- pred[-index_na]
  }
  if (is.factor(obs)) {
    if (is.factor(pred)) {
      pred <- as.character(pred)
      pred <- factor(pred, levels = levels(obs))
      err <- list(error = mean(obs != pred), accuracy = mean(obs == 
                                                               pred))
      if (nlevels(obs) == 2) {
        npos <- sum(obs == levels(obs)[2])
        nneg <- sum(obs == levels(obs)[1])
        ntruepos <- sum((obs == levels(obs)[2]) & (pred == 
                                                     levels(obs)[2]))
        ntrueneg <- sum((obs == levels(obs)[1]) & (pred == 
                                                     levels(obs)[1]))
        err$sensitivity <- ntruepos/npos
        err$specificity <- ntrueneg/nneg
        npospred <- sum(pred == levels(obs)[2])
        nnegpred <- sum(pred == levels(obs)[1])
        err$ppv <- ntruepos/npospred
        err$npv <- ntrueneg/nnegpred
        n <- length(obs)
        pexp <- (npos/n) * (npospred/n) + (nneg/n) * 
          (nnegpred/n)
        if (pexp == 1) {
          err$kappa <- NA
        }
        else {
          err$kappa <- (err$accuracy - pexp)/(1 - pexp)
        }
      }
    }
    else {
      if (!is.vector(pred)) {
        pred <- as.numeric(pred)
      }
      if (!is.vector(pred) && !is.matrix(pred)) {
        pred <- as.numeric(pred)
      }
      predobj <- ROCR::prediction(pred, obs)
      auroc <- ROCR::performance(predobj, measure = "auc")@y.values[[1]]
      err <- list(auroc = auroc)
      pos <- levels(obs)[2]
      err$error <- mean((obs == pos) != (pred >= 0.5), 
                        na.rm = TRUE)
      err$accuracy <- 1 - err$error
      tpr <- function(o, p, np, t) sum(o & (p >= t))/np
      fpr <- function(o, p, nn, t) sum(!o & (p >= t))/nn
      npos <- sum(obs == pos)
      nneg <- sum(obs != pos)
      err$sensitivity <- tpr(obs == pos, pred, npos, 0.5)
      err$specificity <- 1 - fpr(obs == pos, pred, nneg, 
                                 0.5)
      thrs <- unique(pred)
      if (length(thrs) > 500) {
        thrs <- seq(min(pred) + 1e-04, max(pred) + 1e-04, 
                    length = 500)
      }
      thr <- mmax(thrs[sapply(thrs, function(x) {
        tpr(obs == pos, pred, npos, x) >= 0.7
      })])
      err$fpr70 <- fpr(obs == pos, pred, nneg, thr)
      thr <- mmax(thrs[sapply(thrs, function(x) {
        tpr(obs == pos, pred, npos, x) >= 0.8
      })])
      err$fpr80 <- fpr(obs == pos, pred, nneg, thr)
      thr <- mmax(thrs[sapply(thrs, function(x) {
        tpr(obs == pos, pred, npos, x) >= 0.9
      })])
      err$fpr90 <- fpr(obs == pos, pred, nneg, thr)
      thr <- mmin(thrs[sapply(thrs, function(x) {
        fpr(obs == pos, pred, nneg, x) <= 0.2
      })])
      err$tpr80 <- tpr(obs == pos, pred, npos, thr)
      thr <- mmin(thrs[sapply(thrs, function(x) {
        fpr(obs == pos, pred, nneg, x) <= 0.1
      })])
      err$tpr90 <- tpr(obs == pos, pred, npos, thr)
      thr <- mmin(thrs[sapply(thrs, function(x) {
        fpr(obs == pos, pred, nneg, x) <= 0.05
      })])
      err$tpr95 <- tpr(obs == pos, pred, npos, thr)
    }
    err$events <- sum(obs == levels(obs)[2])
  }
  else {
    err <- list(bias = mean(obs - pred, na.rm=T),
                stddev = sd(obs - pred, na.rm=T), 
                mae = mean(abs(obs - pred), na.rm=T),
                rpearson = cor(obs, pred, use="complete.obs"),
                rsquared = (cor(obs, pred, use="complete.obs")^2),
                rmse = sqrt(mean((obs - pred)^2, na.rm=T)), 
                mad = mad(obs - pred, na.rm=T),
                median = median(obs - pred, na.rm=T),
                iqr = IQR(obs - pred, na.rm=TRUE))
  }
  err$count <- length(obs)
  return(err)
}

# DATA PRE-PROCESSING -----------------------------------------------------
# exclude covariates based on the variable selecction procedure
names(SU)
SU = SU %>%
  dplyr::select(-c(TC_stdev, TPI_stdev, TRI_mean, TRI_stdev, Aspect_mean,
                   Aspect_stdev, TC_stdev, Geology, Landcover,
                   MMI_mean, MMI_stdev, PGV_mean, SUM_AREA_M2,
                   MAX_AREA_M2, MAX_TRUEAREA_M2, AREA_SU_M2)) %>%
  dplyr::rename(Id = SU3_ID,
                Bin = Landslide,
                Area_sum_true = SUM_TRUEAREA_M2,
                Area_SU = TRUEAREA_SU_M2) %>%
  mutate(Relief_range = as.numeric(Relief_range)) %>%
  mutate(Soildepth_majority = as.factor(Soildepth_majority))

# generate log(area)
SU = SU %>%
  mutate(Log_area_sum_true = my.log(Area_sum_true)) %>%
  relocate(Log_area_sum_true, .after=Area_sum_true)
names(SU)

# coping Slope Units and removing na values
SU_whole = SU
SU = SU %>% tidyr::drop_na()

# MODELING ----------------------------------------------------------------

# limit road distance 
SU = my.trafo(SU)

# define formula for GAM
smooth_term = 4

Formula_truesum = Log_area_sum_true ~ s(Slope_mean, k = smooth_term) + #topography
  s(Slope_stdev, k = smooth_term) + s(Northerness_mean, k = smooth_term) + #topography
  s(Easterness_mean, k = smooth_term) + s(Easterness_stdev, k = smooth_term) + #topography
  s(Relief_range, k = smooth_term) + s(PRC_stdev, k = smooth_term) + s(TC_mean, k = smooth_term) + #topography
  Soildepth_majority + #soil
  s(Roads_mean, k = smooth_term) + s(Roads_stdev, k = smooth_term) + #roads
  s(PGA_mean, k = 3) + s(PGA_stdev, k=smooth_term) #earthquake

# fitting GAM
Mod_truesum = mgcv::gam(Formula_truesum, family = gaussian, data = SU)
summary(Mod_truesum)

# component smooth function plots
plot(Mod_truesum, all.terms=T, residuals=F, shade=T, shade.col="#B0C4DE",
     sewithMean=T, se=T)

# fitting performance
SU = SU %>%
  mutate(Pred = mgcv::predict.gam(Mod_truesum, SU)) %>%
  relocate(Pred, .after = Log_area_sum_true)

# fitting performance metrics 
cor(SU$Log_area_sum_true, SU$Pred, use="complete.obs")
my.RMSE(SU$Log_area_sum_true, SU$Pred, na.rm = T)
my.MAE(SU$Log_area_sum_true, SU$Pred, na.rm = T)

# plot agreement between observed and fitted data
ggplot(SU, aes(Log_area_sum_true, Pred)) + geom_hex(bins = 58) + 
  geom_abline(slope = 1, intercept =  0, col = "red", size = 0.7, linetype = "dashed") +
  expand_limits(x = 0, y= 0) +
  scale_x_continuous(expand = c(0, 0), limits = c(5, 15)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(5, 15)) +
  coord_fixed() +
  scale_fill_continuous(type = "viridis", direction = 1) +
  labs(title = "FIT", x = bquote("\nObserved landslide area"~log(m^2)),
       y = bquote("Estimated landslide area "~log(m^2))) +
  theme_linedraw() + theme_bw() + 
  theme(text = element_text(size = 14, family = "LM Roman 10"), plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.grid = element_blank(), axis.text = element_text(size = 12, colour = "black"),
        legend.title = element_blank())


# NON-SPATIAL CROSS-VALIDATION --------------------------------------------
# exploring random partition
resamp = partition_cv(SU, nfold = 3, repetition = 1, seed1= 1, coords = c("X", "Y")) # example with 3 folds for SUs with landslides
plot(resamp, SU, coords = c("X", "Y"), cex = 0.01, pch = 19)

# performing non-spatial cross-validation, 10 rep and 10 folds. This step takes a while
pred_cv = sperrorest(data = SU_whole, formula = Formula_truesum,
                     model_fun = mgcv::gam,
                     coords = c("X", "Y"),
                     model_args = list(family = gaussian),
                     pred_fun = predict,
                     progress = 2,
                     err_fun = my.error,
                     smp_fun = partition_cv,
                     smp_args = list(repetition = 1:10, nfold = 10, seed1= 1))

# check the results
summary(pred_cv$error_rep)
summary(pred_cv$error_fold)
summary(pred_cv$represampling)
pred_cv$error_rep

# generation of predictions
SU_whole$Pred_cv = NA
resamp = partition_cv(SU_whole, nfold = 10, repetition = 1, seed1= 1, coords = c("X", "Y")) 

# looping training and prediction 
# creating training and test dataframes
for (i in 1:10){
  cat("Running classification iteration", i, "\n")
  id.holdout = resamp[["1"]][[i]]$test
  SU_whole.train = SU_whole[-id.holdout, ]
  SU_whole.test = SU_whole[id.holdout, ]
  
  # training model
  Mod1_truesum = mgcv::gam(Formula_truesum, family = gaussian, data = SU_whole.train)
  
  # predicting model
  SU_whole$Pred_cv[id.holdout] = mgcv::predict.gam(Mod1_truesum, SU_whole.test)
}

# plot agreement between observed and estimated data (RCV)
ggplot(SU_whole, aes(Log_area_sum_true, Pred_cv)) + geom_hex(bins = 60) + 
  geom_abline(slope = 1, intercept =  0, col = "red", size = 0.7, linetype = "dashed") +
  expand_limits(x = 0, y= 0) +
  scale_x_continuous(expand = c(0, 0), limits = c(5, 15)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(5, 15)) +
  coord_fixed() + scale_fill_continuous(type = "viridis") +
  labs(title = "CV", x = bquote("\nResponse  - landslide area"~log(m^2)),
       y = bquote("Fitted values - landslide area "~log(m^2))) +
  theme_linedraw() + theme_bw() + 
  theme(text = element_text(size = 14, family = "LM Roman 10"), plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.grid = element_blank(), axis.text = element_text(size = 12, colour = "black"),
        legend.title = element_blank())


# SPATIAL CROSS-VALIDATION ------------------------------------------------

# exploring the data partition
resamp_s = partition_kmeans(SU, nfold = 5, repetition = 1, seed1= 1, coords = c("X", "Y")) #example with 3 folds for Su with landslide data
plot(resamp_s, SU, coords = c("X", "Y"), pch = 19, cex = 0.001)


# performing non-spatial cross-validation, 10 rep and 10 folds. This step takes a while
pred_scv = sperrorest(data = SU_whole,
                      coords = c("X", "Y"),
                      formula = Formula_truesum,
                      model_fun = mgcv::gam,
                      model_args = list(family = gaussian),
                      pred_fun = predict,
                      progress = 2,
                      smp_fun = partition_kmeans,
                      err_fun = my.error,
                      smp_args = list(repetition = 1:10, nfold = 10, seed1= 1))

# check the results
summary(pred_scv$error_rep)
summary(pred_scv$error_fold)
summary(pred_scv$represampling)
pred_scv$error_rep

# Generation of predictions
SU_whole$Pred_scv = NA
resamp_s = partition_kmeans(SU_whole, nfold = 10, repetition = 1, seed1= 1, coords = c("X", "Y")) 

# looping training and prediction
# creating training and test dataframes- This should be the correct way.
for (i in 1:10){
  cat("Running classification iteration", i, "\n")
  id.holdout.s = resamp_s[["1"]][[i]]$test
  SU_whole.train.s = SU_whole[-id.holdout.s, ]
  SU_whole.test.s = SU_whole[id.holdout.s, ]
  
  ## training model
  Mod1_truesum.s = mgcv::gam(Formula_truesum, family = gaussian, data = SU_whole.train.s)
  
  ## predicting model
  SU_whole$Pred_scv[id.holdout.s] = mgcv::predict.gam(Mod1_truesum.s, SU_whole.test.s)
}

# plot agreement between observed and estimated data (RCV)
ggplot(SU_whole, aes(Log_area_sum_true, Pred_scv)) + geom_hex(bins = 52) + 
  geom_abline(slope = 1, intercept =  0, col = "red", size = 0.7, linetype = "dashed") +
  expand_limits(x = 0, y= 0) + coord_fixed() +
  scale_x_continuous(expand = c(0, 0), limits = c(5, 15)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(5, 15)) +
  scale_fill_continuous(type = "viridis") +
  labs(title = "SCV", x = bquote("\nObserved  - landslide area"~log(m^2)),
       y = bquote("Fitted values - landslide area "~log(m^2))) +
  theme_linedraw() + theme_bw() + 
  theme(text = element_text(size = 14, family = "LM Roman 10"), plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.grid = element_blank(), axis.text = element_text(size = 12, colour = "black"),
        legend.title = element_blank())

# PLOT RESULTS SPERROREST -------------------------------------------------

# creating dataframe with performance metrics (Pearson, MAE, RMSE)
validation = data.frame(CV_r = unlist(summary(pred_cv$error_rep, level = 1)$test_rpearson),
                        CV_rsq = unlist(summary(pred_cv$error_rep, level = 1)$test_rsquared),
                        CV_mae = unlist(summary(pred_cv$error_rep, level = 1)$test_mae),
                        CV_rmse = unlist(summary(pred_cv$error_rep, level = 1)$test_rmse),
                        SCV_r = unlist(summary(pred_scv$error_rep, level = 1)$test_rpearson),
                        SCV_rsq = unlist(summary(pred_scv$error_rep, level = 1)$test_rsquared),
                        SCV_rmse = unlist(summary(pred_scv$error_rep, level = 1)$test_rmse),
                        SCV_mae = unlist(summary(pred_scv$error_rep, level = 1)$test_mae))

# reshaping data
validation = data.frame(cbind(ID = rownames(validation), validation))
validation_melt = melt(validation, id.vars="ID", measure.vars=c("CV_r", "CV_rsq", 
                                                                "CV_mae", "CV_rmse",
                                                                "SCV_r", "SCV_rsq",
                                                                "SCV_rmse", "SCV_mae"))

# dividing boxplor per metric
validationmelt_rmse = validation_melt %>% 
  dplyr::filter(variable == "CV_rmse"|variable == "SCV_rmse") %>%
  dplyr::mutate(variable = recode(variable, "CV_rmse" = "CV", "SCV_rmse" = "SCV"))
validationmelt_mae = validation_melt %>% 
  dplyr::filter(variable == "CV_mae"|variable == "SCV_mae") %>%
  dplyr::mutate(variable = recode(variable, "CV_mae" = "CV", "SCV_mae" = "SCV"))
validationmelt_r = validation_melt %>% 
  dplyr::filter(variable == "CV_r"|variable == "SCV_r") %>%
  dplyr::mutate(variable = recode(variable, "CV_r" = "CV", "SCV_r" = "SCV"))
validationmelt_rsq = validation_melt %>% 
  dplyr::filter(variable == "CV_rsq"|variable == "SCV_rsq") %>%
  dplyr::mutate(variable = recode(variable, "CV_rsq" = "CV", "SCV_rsq" = "SCV"))

# boxplots performance
# RMSE
ggplot(data = validationmelt_rmse, aes(x = names, fill = variable, color = variable)) + 
  geom_boxplot(width = 0.5, aes(x=variable, y = value)) +
  labs(title = NULL, x = NULL, y = bquote("RMSE " ~ log(m^2))) +
  stat_summary(mapping = aes(x= variable, y=value), fun = "mean", geom = "point",
               shape = 21, size = 2, color = "black", fill = "yellow1") +
  scale_fill_manual(name = "variable", values = c("#E53232", "steelblue3")) +
  scale_color_manual(name = "variable", values = c("black", "black")) +
  theme_linedraw() + theme_bw() +
  theme(text = element_text(size = 14, family = "LM Roman 10"), plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.grid = element_blank(), axis.text = element_text(size = 12, colour = "black"),
        legend.title = element_blank(), legend.position = "none")

# MAE
ggplot(data = validationmelt_mae, aes(x = names, fill = variable, color = variable)) + 
  geom_boxplot(width = 0.5, aes(x=variable, y = value)) +
  labs(title = NULL, x = NULL, y = bquote("MAE " ~ log(m^2))) +
  stat_summary(mapping = aes(x= variable, y=value), fun = "mean", geom = "point",
               shape = 21, size = 2, color = "black", fill = "yellow1") +
  scale_fill_manual(name = "variable", values = c("#E53232", "steelblue3")) +
  scale_color_manual(name = "variable", values = c("black", "black")) +
  theme_linedraw() + theme_bw() +
  theme(text = element_text(size = 14, family = "LM Roman 10"), plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.grid = element_blank(), axis.text = element_text(size = 12, colour = "black"),
        legend.title = element_blank(), legend.position = "none")
  
# Pearson
ggplot(data = validationmelt_r, aes(x = names, fill = variable, color = variable)) + 
    geom_boxplot(width = 0.5, aes(x=variable, y = value)) +
    labs(title = NULL, x = NULL, y = bquote("R-Pearson ")) +
    stat_summary(mapping = aes(x= variable, y=value), fun = "mean", geom = "point",
                 shape = 21, size = 2, color = "black", fill = "yellow1") +
    scale_fill_manual(name = "variable", values = c("#E53232", "steelblue3")) +
    scale_color_manual(name = "variable", values = c("black", "black")) +
    theme_linedraw() + theme_bw() +
    theme(text = element_text(size = 14, family = "LM Roman 10"), plot.margin = unit(c(1, 1, 1, 1), "cm"),
          panel.grid = element_blank(), axis.text = element_text(size = 12, colour = "black"),
          legend.title = element_blank(), legend.position = "none")

# compiling results in single file
fit = data.frame(cbind(Id=SU$Id, Pred=SU$Pred))
SU_whole = left_join(SU_whole, fit, by="Id")
SU_whole = sf::st_as_sf(SU_whole)

# building color palette
col.na = rgb(120, 120, 120, max = 255, alpha = 0)
my.color = c("#E0E0E0", "#87BBFA", "#1FC44E", "#8BD100", "#FFFF00", "#FA8900", "#E60000")

# checking results visually
mapview(SU_whole, zcol="Log_area_sum_true", col.regions = my.color,
        na.color = col.na, map.types = "Esri.WorldImagery", layer.name = "Observed LA",
        at = c(5.4, 6.8, 7.7, 8.6, 9.5, 10.4, 11.4, 14)) # observed  
mapview(SU_whole, zcol="Pred_cv", col.regions = my.color,
        na.color = col.na, map.types = "Esri.WorldImagery", layer.name = "CV preditions",
        at = c(5.4, 6.8, 7.7, 8.6, 9.5, 10.4, 11.4, 14)) # cv predictions
mapview(SU_whole, zcol="Pred_scv", col.regions = my.color,
        na.color = col.na, map.types = "Esri.WorldImagery", layer.name = "SCV predictions",
        at = c(5.4, 6.8, 7.7, 8.6, 9.5, 10.4, 11.4, 14)) # scv predictions











