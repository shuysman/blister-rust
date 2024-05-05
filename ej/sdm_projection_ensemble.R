library(terra)
library(sdmTMB)
library(tidyverse)
library(glue)
library(tidyterra)
library(pbmcapply)

cores <- parallel::detectCores() - 4

sdm_model <- readRDS("sdm_t_2_x_p.RDS")

##maca_dir <- "/media/smithers/shuysman/data/MACA/wbp_range/monthly/"
maca_dir <- "~/data/MACA/wbp_range/monthly/"

##out_dir <- "/media/smithers/shuysman/data/out/blister/sdm_v2/"
out_dir <- "~/out/blister/sdm_v2/"

models <- c("NorESM1-M", "MRI-CGCM3", "MIROC-ESM-CHEM", "MIROC5", "IPSL-CM5A-LR", "inmcm4", "HadGEM2-CC365", "CSIRO-Mk3-6-0", "CNRM-CM5", "CanESM2", "BNU-ESM", "GFDL-ESM2G")
### CCSM4 tmmn data missing
##models <- c("NorESM1-M")
scenarios <- c('rcp85', 'rcp45')
##scenarios <- c("rcp45")

current <- seq(2010, 2025)
near <- seq(2026, 2050)
mid <- seq(2051, 2075)
end <- seq(2076, 2099)

##log_dbh <- c(-2, -1, 0, 1, 2, 3.5, 4)

log_dbh <- log(c(1,5,10,20,40,80)) %>% round(1)

crs_lcc <- "+proj=lcc +lon_0=-117 +lat_1=36 +lat_2=49 +units=km"

##wbp_existing <- project(rast("/home/steve/OneDrive/whitebark/data/existing/wbp_existing.tif"), crs(crs_lcc))

##wbp_existing <- wbp_existing == 1

process_gcm <- function(options) {
    scenario <- options$scenario
    pred_range <- options$pred_range
    log_dbh <- options$log_dbh

    tmin_filename <- paste0(maca_dir, "tmmn_", models, "_", scenario, "_2006-2099_monthly_gye.nc")
    tmin <- rast(tmin_filename) %>%
        subset(month(time(.)) == 8 | month(time(.)) == 9) %>%
        subset(year(time(.)) %in% pred_range) %>%
        terra::tapp(index = "days", fun = mean)
    tmax_filename <- paste0(maca_dir, "tmmx_", models, "_", scenario, "_2006-2099_monthly_gye.nc")
    tmax <- rast(tmax_filename) %>%
        subset(month(time(.)) == 8 | month(time(.)) == 9) %>%
        subset(year(time(.)) %in% pred_range) %>%
        terra::tapp(index = "days", fun = mean)
    pr_filename <- paste0(maca_dir, "pr_", models, "_", scenario, "_2006-2099_monthly_gye.nc")
    pr <- rast(pr_filename) %>%
        subset(month(time(.)) == 8 | month(time(.)) == 9) %>%
        subset(year(time(.)) %in% pred_range) %>%
        terra::tapp(index = "days", fun = mean)

    tminmax <- c(tmin, tmax)
    tavg <- tapp(tminmax, fun = mean, index = "years")
    names(tavg) <- rep("asT", nlyr(tavg))
    tavg <- tavg - 273.15 ## convert K to C

    tavg_global <- mean(tavg)
    names(tavg_global) <- "asT"

    pr_sum <- tapp(pr, fun = sum, index = "years")
    names(pr_sum) <- rep("asP", nlyr(pr_sum))

    pr_global <- mean(pr_sum)
    names(pr_global) <- "asP"
    
    out_file <- file.path(out_dir, glue("ensemble_{scenario}_{head(pred_range, 1)}-{tail(pred_range, 1)}_{log_dbh}-logdbh_brinfpr.nc"))
    if (file.exists(out_file)) {
        next
    }
    


    as_t_p <- c(tavg_global, pr_global)

    as_t_p_lcc <- project(as_t_p, crs(crs_lcc))
    
    as_t_p_lcc_df <- as.data.frame(as_t_p_lcc, xy = TRUE, na.rm = TRUE) %>%
        rename(X = x, Y = y)

    as_t_p_lcc_df$log_dbh_cm = log_dbh

    predicted <- predict(sdm_model,
                         as_t_p_lcc_df,
                         type = "response"
                         )

    p_rast <- rast(predicted)

    crs(p_rast) <- crs(crs_lcc)

    writeCDF(p_rast,
             filename = out_file,
             split = TRUE
             ## varname = glue("{pred_range[1]}-{pred_range[2]}_brinfpr"),
             ## longname = glue("{pred_range[1]}-{pred_range[2]} Probability of Blister Rust infection")
             ) 
}


options <- expand.grid(scenario = scenarios, pred_range = list(current, near, mid, end), log_dbh = log_dbh) %>%
    t() %>%
    data.frame()  ### How to use (mc)lapply with expand.grid https://gist.github.com/ihrke/c54e8d7c82cd91ae51d1d62f8d29b936

pbmclapply(options, 
         FUN = process_gcm,
         mc.cores = cores)
