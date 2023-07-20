library(terra)
library(tidyverse)
library(ggeffects)
library(sf)

maca_dir <- "/media/smithers/shuysman/data/MACA/"
out_dir <- "/media/smithers/shuysman/data/out/blister/"

m17.5 <- readRDS("m17_5.RDS")

recentdbh <- 35
## I was not able to download some of the monthly aggregrated
## forecasts for t/rh from MACA thredds as of <2023-07-19 Wed>
## NorESM1-M
## MIROC-ESM-CHEM
## CSIRO-Mk3-6-0
## CanESM2
## CCSM4
models <- c("NorESM1-M", "MRI-CGCM3", "MIROC-ESM-CHEM", "MIROC5", "IPSL-CM5A-LR", "inmcm4", "HadGEM2-CC365", "CSIRO-Mk3-6-0", "CNRM-CM5", "CanESM2", "BNU-ESM", "CCSM4", "GFDL-ESM2G")
##models <- c("NorESM1-M")
scenarios <- c('rcp85', 'rcp45')

gye_boundary <- st_read("/home/steve/OneDrive/whitebark/gyeboundary/GYE_boundary_dd.shp")

near_century <- seq(2024, 2049)
mid_century <- seq(2050, 2075)
end_century <- seq(2075, 2099)

##rmax_MRI-CGCM3_rcp85_2006-2099_monthly_gye.nc
for (model in models) {
    for (scenario in scenarios) {
        print(paste(model, scenario))

        ## model <- "BNU-ESM"
        ## scenario <- "rcp45"

        skip_to_next <- FALSE
        tryCatch(
            {
                rmin_filename <- paste0(maca_dir, "rmin_", model, "_", scenario, "_2006-2099_monthly_gye.nc")
                rmin <- rast(rmin_filename) %>% subset(month(time(.)) == 8 | month(time(.)) == 9)
                names(rmin) <- rep("rmin", nlyr(rmin))
                rmax_filename <- paste0(maca_dir, "rmax_", model, "_", scenario, "_2006-2099_monthly_gye.nc")
                rmax <- rast(rmax_filename) %>% subset(month(time(.)) == 8 | month(time(.)) == 9)
                names(rmax) <- rep("rmax", nlyr(rmin))

                tmin_filename <- paste0(maca_dir, "tmmn_", model, "_", scenario, "_2006-2099_monthly_gye.nc")
                tmin <- rast(tmin_filename) %>% subset(month(time(.)) == 8 | month(time(.)) == 9)
                tmax_filename <- paste0(maca_dir, "tmmx_", model, "_", scenario, "_2006-2099_monthly_gye.nc")
                tmax <- rast(tmax_filename) %>% subset(month(time(.)) == 8 | month(time(.)) == 9)
            },
            error = function(e) {
                skip_to_next <<- TRUE
            }
        )

        if (skip_to_next == TRUE) {
            print(paste(model, scenario, "file does not exist."))
            break
        }

        rhminmax <- c(rmin, rmax)

        rhavg <- tapp(rhminmax, fun = mean, index = "years")
        
        tminmax <- c(tmin, tmax)

        tavg <- tapp(tminmax, fun = mean, index = "years")

        names(rhavg) <- rep("asRH", nlyr(rhavg))
        names(tavg) <- rep("asTEMP", nlyr(tavg))

        tavg <- tavg - 273.15 ## convert K to C

        as_rh_temp <- c(rhavg, tavg)

        for (year_set in list(near_century, mid_century, end_century)) {
            br_prob <- rast()
            for (year in year_set) {
                r <- subset(as_rh_temp, time(as_rh_temp) == year)
                predicted <- predict(r,
                    model = m17.5,
                    const = data.frame(recentdbh = recentdbh),
                    re.form = ~0,
                    type = "response"
                )
                time(predicted) <- year
                br_prob <- c(br_prob, predicted)
            }
            br_prob_mean <- mean(br_prob)
            writeCDF(br_prob_mean,
                     paste0(out_dir, "brinfpr_", model, "_", scenario, "_", year_set[1], "-", tail(year_set, 1), "_", recentdbh, "cm.nc"),
                     overwrite = TRUE)
        }
    }
}


## Ensembles
nc45_files <- list.files(path = out_dir, pattern = paste0("brinfpr_.*rcp45_2024-2049_", recentdbh, "cm.nc"))
ensemble_nc45 <- mean(rast(paste0(out_dir, nc45_files)))
ggplot(as.data.frame(ensemble_nc45, xy = TRUE)) +
    geom_raster(mapping = aes(x = x, y = y, fill = mean)) +
    geom_sf(data = gye_boundary, alpha = 0.15, fill = NA, lwd = 2) +
    theme_bw() +
    scico::scale_fill_scico(palette = "bilbao", begin = 0, limits = c(0,1)) +
    ggtitle(paste("RCP45", near_century[1], "-", tail(near_century, 1), "DBH", recentdbh, "cm Ensemble Mean P(infection)")) +
    scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
    scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
ggsave(filename = paste0("img/", "RCP45_", near_century[1], "-", tail(near_century, 1), "_", recentdbh, "cm.png"), width = 8, height = 8) 

nc85_files <- list.files(path = out_dir, pattern = paste0("brinfpr_.*rcp85_2024-2049_", recentdbh, "cm.nc"))
ensemble_nc85 <- mean(rast(paste0(out_dir, nc85_files)))
ggplot(as.data.frame(ensemble_nc85, xy = TRUE)) +
    geom_raster(mapping = aes(x = x, y = y, fill = mean)) +
    geom_sf(data = gye_boundary, alpha = 0.15, fill = NA, lwd = 2) +
    theme_bw() +
    scico::scale_fill_scico(palette = "bilbao", limits = c(0,1)) +
    ggtitle(paste("RCP85", near_century[1], "-", tail(near_century, 1), "DBH", recentdbh, "cm Ensemble Mean P(infection)")) +
    scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
    scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
ggsave(filename = paste0("img/", "RCP85_", near_century[1], "-", tail(near_century, 1), "_", recentdbh, "cm.png"), width = 8, height = 8)

mc45_files <- list.files(path = out_dir, paste0("brinfpr_.*rcp45_2050-2075_", recentdbh, "cm.nc"))
ensemble_mc45 <- mean(rast(paste0(out_dir, mc45_files)))
ggplot(as.data.frame(ensemble_mc45, xy = TRUE)) +
    geom_raster(mapping = aes(x = x, y = y, fill = mean)) +
    geom_sf(data = gye_boundary, alpha = 0.15, fill = NA, lwd = 2) +
    theme_bw() +
    scico::scale_fill_scico(palette = "bilbao", limits = c(0,1)) +
    ggtitle(paste("RCP45", mid_century[1], "-", tail(mid_century, 1), "DBH", recentdbh, "cm Ensemble Mean P(infection)")) +
    scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
    scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
ggsave(filename = paste0("img/", "RCP45_", mid_century[1], "-", tail(mid_century, 1), "_", recentdbh, "cm.png"), width = 8, height = 8)


mc85_files <- list.files(path = out_dir, paste0("brinfpr_.*rcp85_2050-2075_", recentdbh, "cm.nc"))
ensemble_mc85 <- mean(rast(paste0(out_dir, mc85_files)))
ggplot(as.data.frame(ensemble_mc85, xy = TRUE)) +
    geom_raster(mapping = aes(x = x, y = y, fill = mean)) +
    geom_sf(data = gye_boundary, alpha = 0.15, fill = NA, lwd = 2) +
    theme_bw() +
    scico::scale_fill_scico(palette = "bilbao", limits = c(0,1)) +
    ggtitle(paste("RCP85", mid_century[1], "-", tail(mid_century, 1), "DBH", recentdbh, "cm Ensemble Mean P(infection)")) +
    scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
    scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
ggsave(filename = paste0("img/", "RCP85_", mid_century[1], "-", tail(mid_century, 1), "_", recentdbh, "cm.png"), width = 8, height = 8)

ec45_files <- list.files(path = out_dir, paste0("brinfpr_.*rcp45_2075-2099_", recentdbh, "cm.nc"))
ensemble_ec45 <- mean(rast(paste0(out_dir, ec45_files)))
ggplot(as.data.frame(ensemble_ec45, xy = TRUE)) +
    geom_raster(mapping = aes(x = x, y = y, fill = mean)) +
    geom_sf(data = gye_boundary, alpha = 0.15, fill = NA, lwd = 2) +
    theme_bw() +
    scico::scale_fill_scico(palette = "bilbao", limits = c(0,1)) +
    ggtitle(paste("RCP45", end_century[1], "-", tail(end_century, 1), "DBH", recentdbh, "cm Ensemble Mean P(infection)")) +
    scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
    scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
ggsave(filename = paste0("img/", "RCP45_", end_century[1], "-", tail(end_century, 1), "_", recentdbh, "cm.png"), width = 8, height = 8)


ec85_files <- list.files(path = out_dir, paste0("brinfpr_.*rcp85_2075-2099_", recentdbh, "cm.nc"))
ensemble_ec85 <- mean(rast(paste0(out_dir, ec85_files)))
ggplot(as.data.frame(ensemble_ec85, xy = TRUE)) +
    geom_raster(mapping = aes(x = x, y = y, fill = mean)) +
    geom_sf(data = gye_boundary, alpha = 0.15, fill = NA, lwd = 2) +
    theme_bw() +
    scico::scale_fill_scico(palette = "bilbao", limits = c(0,1)) +
    ggtitle(paste("RCP85", end_century[1], "-", tail(end_century, 1), "DBH", recentdbh, "cm Ensemble Mean P(infection)")) +
    scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
    scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
ggsave(filename = paste0("img/", "RCP85_", end_century[1], "-", tail(end_century, 1), "_", recentdbh, "cm.png"), width = 8, height = 8)
