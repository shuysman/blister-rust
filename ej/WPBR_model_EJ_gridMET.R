# ------------------------------------------
# SIP Project: development and validation of a WPBR sensitivity model
# 7/24/34
# EKJ
# ------------------------------------------
#

# Script description ------------------------------------------------------
# LINE 42: Section 1- Load and format data
# LINE 175: Section 2- Model development
# LINE 311: Section 3- Model selection
# LINE 455: Section 4- Model validation
# LINE 517: Section 5- Modeling with VPD
# LINE 938: Setion 6- Model figures

# Load packages -----------------------------------------------------------
library(dplyr)
library(lme4)
library(scales)
library(reshape2)
library(ggplot2)
library(plotrix)
library(plotly)
library(naniar)
library(lmerTest)
library(ggrepel)
library(ggthemes)
library(ggpattern)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(purrr)
library(lubridate)

# model evaluation
library(MuMIn)
library(caret)
library(car)
library(Hmisc)
library(sjPlot) # tabl_model
library(report) # report
library(shiny)
library(DescTools) # pseudo R2 value


###############################################################
# Sections 1-7 below create the following csv file:
mod_data <- read.csv(file = "./gridMET_climate_data_SUMMARY.csv") %>% dplyr::filter(park != "GRBA" & park != "GRSA" & park != "CRMO_CRMP" & park != "ROMO")
head(mod_data)
tail(mod_data)
###############################################################


###############################################################
# SECTION 1: LOAD AND FORMAT DATA ----------------------------------------------------
###############################################################
# Copied from WPBR_sensitivity_model_development.R

# 1. Tree data from IMD high-elevation five needle pine monitoring sites
#   - csv created in "clean_site_monitoring_dataset.R"

# Load tree data
tree_data <- read.csv("./Complete_Site_Monitoring_Dataset.csv")
head(tree_data)
length(tree_data$transect_id)

# add "park"
park_lookup <- read.csv("./SITE_LOCATIONS.csv")
park_lookup <- dplyr::select(park_lookup, c("park", "transect_id"))
head(park_lookup)
tree_data <- merge(park_lookup, tree_data, by = c("transect_id"))
head(tree_data)
which(is.na(tree_data$park))

length(tree_data$transect_id)
table(tree_data$br_status, by = tree_data$park)

# 2. Extract site infection status

# count number of infected trees at the transect then
# assign blister rust infection status:
#   0 = 0 trees at site are infected
#   1 = 1+ trees at site are infected
site_status <- tree_data %>%
  group_by(park, transect_id) %>%
  dplyr::summarize(inf_trees = sum(br_status),
                   recent_br_status = cut(inf_trees, breaks = c(-1,1,200),
                                          labels = c(0, 1))) %>%
  dplyr::select(c("park", "transect_id", "recent_br_status")) %>% 
  as.data.frame()

head(site_status); tail(site_status)
length(site_status$transect_id)
unique(site_status$park)

# explore blister rust status by park
table(site_status$recent_br_status, by = site_status$park)

# 3. Extract and summarize gridMET data for each transcect

# Combining csv files

#read gridMET cvs files into a list
files <- list.files(path = "out/elev_corrected",  # Identify all CSV files
                    pattern = "*.csv", full.names = TRUE)

#extract data from csv files and combine into a single dataframe
all_data <- do.call(rbind, lapply(files, function(x) 
  transform(read.csv(x), transect_id  = basename(x))))
head(all_data); tail(all_data)

#remove ".csv" from transect_id name
all_data$transect_id = substr(all_data$transect_id,1,nchar(all_data$transect_id)-4)
head(all_data); tail(all_data)

#remove column X (row number from csv)
all_data <- all_data[,!names(all_data) %in% c("X")]
head(all_data); tail(all_data)

# move transect_id to first column
all_data <- all_data %>% relocate(transect_id)
head(all_data); tail(all_data)

#change Date column from character to date
str(all_data)
all_data$Date <- ymd(all_data$Date)
str(all_data)

#subset to August and September
gridmet <- subset(all_data, format.Date(Date, "%m")=="08" | format.Date(Date, "%m")=="09")
head(gridmet); tail(gridmet)
gridmet

#subset to 2000-2020 (years monitoring data was collected)
gridmet <- gridmet[gridmet$Date > "2000-01-01" &    # Extract data frame subset
                     gridmet$Date < "2021-01-01", ]
head(gridmet); tail(gridmet)

#export data as csv
write.csv(gridmet, file = "Output/Data/gridMET_climate_data/gridMET_climate_data.csv", row.names = FALSE)

# 4. Calculate 2000-2020 mean at each transect
trans_mean <- gridmet %>% 
  group_by(transect_id) %>%
  dplyr::summarize(asTEMP = mean(T_cor, na.rm = TRUE),
                   asRH = mean(RH, na.rm = TRUE),
                   min_t = min(Tmin_cor),
                   max_t = max(Tmax_cor),
                   se_t = std.error(T_cor),
                   min_rh = min(RHmin),
                   max_rh = max(RHmax),
                   se_rh = std.error(RH)) %>% 
  as.data.frame()

head(trans_mean);tail(trans_mean)

trans_mean$period = "2000-2020"

head(trans_mean);tail(trans_mean)
length(trans_mean$transect_id)

# 5. Merge transect-level infection data with climatic correlate data
mod_data_full <- merge(site_status, trans_mean, by = c("transect_id"))
length(mod_data_full$transect_id)
head(mod_data_full)

# reorder mod_data
mod_data_full <- mod_data_full[with(mod_data_full, order(park, transect_id)),]
head(mod_data_full); tail(mod_data_full)

table(mod_data_full$recent_br_status, by = mod_data_full$park)

# format mod_data_full
mod_data_full
head(mod_data_full)
str(mod_data_full)

mod_data_full$recent_br_status <- dplyr::recode(mod_data_full$recent_br_status,
                                                "0" = 0,
                                                "1" = 1)
str(mod_data_full)
table(mod_data_full$recent_br_status, by = mod_data_full$park)


# 6. Format and summarize mod_data_full

# check for NA
which(is.na(mod_data_full))

# summarize mod_data
head(mod_data_full)
summary(mod_data_full$asTEMP)
summary(mod_data_full$asRH)
table(mod_data_full$recent_br_status, by = mod_data$park)
table(mod_data_full$park)
table(mod_data_full$recent_br_status)
length(mod_data_full$transect_id)
unique(mod_data_full$transect_id)


# 7. CACLULATE VPD -----------------------------------------------------------
mod_data_full$svp <- 610.8*exp((17.27*mod_data_full$asTEMP/(237.3+mod_data_full$asTEMP)))
head(mod_data_full)

mod_data_full$asRHcalc <- mod_data_full$asRH/100
head(mod_data_full)

mod_data_full$vpd <- (((100-mod_data_full$asRHcalc*100)/100)*mod_data_full$svp)
head(mod_data_full); tail(mod_data_full)

length(mod_data_full$recent_br_status)

mod_data <- mod_data_full

### EXPORT CVS ####
write.csv(mod_data, file = "Output/Data/gridMET_climate_data/gridMET_climate_data_SUMMARY.csv", row.names = FALSE)


#### Data Vis
lm_min_temp_month2 <- readRDS("lm_min_temp_month2.RDS")
lm_max_temp_month2 <- readRDS("lm_max_temp_month2.RDS")


###############################################################
# SECTION 2: MODEL DEVELOPMENT -----------------------------------------------
###############################################################

# Transect-level model
# logistic regression that predicts probability of transect level infection

# DETERMINE IF MIXED EFFECTS IS JUSTIFIED -------------------------

# https://slcladal.github.io/regression.html#Mixed-Effects_Binomial_Logistic_Regression

# BASELINE FIXED EFFECT
b_glm <- glm(data = mod_data, recent_br_status ~ 1, family = binomial(link = "logit"))

summary(b_glm)
logLik(b_glm)
aic_glm <- AIC(logLik(b_glm)); aic_glm

# BASELINE MIXED EFFECT
b_glmer <- glmer(data = mod_data, recent_br_status ~ (1|park), family = binomial(link = "logit"))

summary(b_glmer)
logLik(b_glmer)
aic_glmer <- AIC(logLik(b_glmer)); aic_glmer

# Model Lilkelihood Ratio Test
# tests whether b_glmer is significantly better than b_glm
null.id = -2 * logLik(b_glm) + 2 * logLik(b_glmer)
pchisq(as.numeric(null.id), df=1, lower.tail=F) # p-value < 0.05
# result: b_glmer is better than b_glm
#         mixed-effects model could be justified


# GLMER MODELS ------------------------------------------------------------

# Rescale variables for mixed effects
mod_data$re_asRH <- scales::rescale(mod_data$asRH, newrange = c(0,1))
mod_data$re_asTEMP <- scales::rescale(mod_data$asTEMP, newrange = c(0,1))
head(mod_data)

# GLMER_REGION
# random effect = region
# fixed effect = interaction
glmer_region <- glmer(data = mod_data, recent_br_status ~ (1|region) +
                        re_asRH + re_asTEMP + re_asRH * re_asTEMP,
                      family=binomial(link = "logit"))

summary(glmer_region)

# dredge function to identify other model candidates
options(na.action = na.fail)
dredge(glmer_region)
# result: model with interaction has lowest AIC

# GLMER_PARK
# random effect = park
# fixed effect = interaction
glmer_park <- glmer(data = mod_data, recent_br_status ~ (1|park) +
                      re_asRH + re_asTEMP + re_asRH * re_asTEMP,
                    family=binomial(link = "logit"))

summary(glmer_park)

# dredge function to identify other model candidates
options(na.action = na.fail)
dredge(glmer_park)
# RESULT: following models are equivalent based on AIC
#         - br ~ re_asRH
#         - br ~ re_asRH + re_asTEMP
#         - br ~ re_asRH * re_asTEMP

# GLMER_REGION VS. GLMER_PARK
aic_region <- AIC(logLik(glmer_region)); aic_region
aic_park <- AIC(logLik(glmer_park)); aic_park

# Model Lilkelihood Ratio Test
# tests whether b_glmer is significantly better than b_glm
null.id = -2 * logLik(glmer_region) + 2 * logLik(glmer_park)
pchisq(as.numeric(null.id), df=1, lower.tail=F) # p-value < 0.05
# RESULT: glmer_park is significantly better than glmer_region

# BEST MODEL WITH RANDOM EFFECT
glmer_park <- glmer(data = mod_data, recent_br_status ~ (1|park) +
                      asRH,
                    family=binomial(link = "logit"))

summary(glmer_park)

# GLM MODELS --------------------------------------------------------------

# GLM_FULL
# t, rh, and interaction
glm_full <- glm(data = mod_data, recent_br_status ~ asRH + asTEMP + asRH * asTEMP,
                family=binomial(link = "logit"))
summary(glm_full)

# dredge function to identify other model candidates
options(na.action = na.fail)
dredge(glm_full)
# RESULT: based on AIC, the optimal model includes asRH, asTEMP, and their interaction

# GLM_FULL_RE
# rescaled t, rh, and interaction
glm_full_re <- glm(data = mod_data, recent_br_status ~ re_asRH + re_asTEMP + re_asRH * re_asTEMP,
                   family=binomial(link = "logit"))
summary(glm_full_re)

# dredge function to identify other model candidates
options(na.action = na.fail)
dredge(glm_full_re)
# RESULT: based on AIC, the optimal model includes re_asRH, re_asTEMP, and their interaction

# GLM_FULL VS. GLM_FULL_RE
# odds ratios for rescaled variables aren't interpretable
# select GLM_FULL

# EXPLORE GLMS WITHOUT INTERACTIONS

# GLM without an interaction
# r and t, no interaction
glm_no_int <- glm(data = mod_data, recent_br_status ~ asRH + asTEMP,
                  family=binomial(link = "logit"))
summary(glm_no_int)

# dredge function to identify other model candidates
options(na.action = na.fail)
dredge(glm_no_int)
# RESULT: asRH is very important when there are no interactions

# GLM_RH
# br predicted only with asRH
glm_rh <- glm(data = mod_data, recent_br_status ~ asRH,
              family=binomial(link = "logit"))
summary(glm_rh)

###############################################################
# SECTION 3: MODEL SELECTION --------------------------------------------------------
###############################################################

# CANDIDATE MODEL SPECS ------------------------------------------

# GLM PARK
tab_model(glmer_park)
# RESULT:
#   - R2=0.12, conditional R2=0.60
#   - asRH (only fixed effect) is significant
report(glmer_park)
# RESULT: substantial explanatory power

# GLM_FULL
tab_model(glm_full)
# RESULT:
#   - R2 = 0.43
#   - all terms are significant
report(glm_full)
# RESULT: substantial explanatory power

# GLM_FULL_RE
tab_model(glm_full_re)
# RESULT:
#   - R2 = 0.43
#   - rh is not significant
report(glm_full_re)
# RESULT: substantial explanatory power

# GLM_NO_INT
tab_model(glm_no_int)
# RESULT:
#   - R2 = 0.40
#   - temp is not significant
report(glm_no_int)

# GLM_RH
tab_model(glm_rh)
# RESULT:
#   - R2 = 0.40
#   - asRH (only fixed effect) is significant
report(glm_rh)
# RESULT: substantial explanatory power
###############################################################

# K-FOLD CROSS VALIDATION -------------------------------------------------

# https://www.r-bloggers.com/2015/08/evaluating-logistic-regression-models/
# https://topepo.github.io/caret/measuring-performance.html
# definitions of confusionMatrix outputs: https://m-clark.github.io/confusionMatrix/reference/calc_stats.html

# K-fold cross validation: 10-fold cross-validation is most common

# set seed for reproducible results
set.seed(10)

# set parameters for train function
ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = TRUE, repeats = 5)

head(mod_data)

# change br_status into a factor for classifying
mod_data$br_fac <- as.character(mod_data$recent_br_status)
str(mod_data)

# train function for different versions of the model
# Can't evaluate GLMER (because of random effect)
set.seed(10)
glm_full_train <- train(br_fac ~ asRH + asTEMP + asRH * asTEMP, data = mod_data, 
                        method = "glm", family = "binomial", 
                        trControl = ctrl, tuneLength = 5)

set.seed(10)
glm_full_re_train <- train(br_fac ~ re_asRH + re_asTEMP + re_asRH * re_asTEMP, data = mod_data, 
                           method = "glm", family = "binomial", 
                           trControl = ctrl, tuneLength = 5)

set.seed(10)
glm_no_int_train <- train(br_fac ~ asRH + asTEMP, data = mod_data, 
                          method = "glm", family = "binomial", 
                          trControl = ctrl, tuneLength = 5)

set.seed(10)
glm_rh_train <- train(br_fac ~ asRH, data = mod_data, 
                      method = "glm", family = "binomial", 
                      trControl = ctrl, tuneLength = 5)

set.seed(10)
glm_rh_re_train <- train(br_fac ~ re_asRH, data = mod_data, 
                         method = "glm", family = "binomial", 
                         trControl = ctrl, tuneLength = 5)

# extract Accuracy and Kappa
print(glm_full_train)
print(glm_full_re_train)
print(glm_no_int_train)
print(glm_rh_train)
print(glm_rh_re_train)

# get confusion matrix
confusionMatrix(glm_full_train, norm = "none")
confusionMatrix(glm_full_re_train, norm = "none")
confusionMatrix(glm_no_int_train, norm = "none")
confusionMatrix(glm_rh_train, norm = "none")
confusionMatrix(glm_rh_re_train, norm = "none")

###############################################################

# VISUALIZE MODEL PREDICTIONS ---------------------------------------------

# visualize predicted probability at each park
mod_data_pred <- mod_data
mod_data_pred$pred <- predict(glm_full, mod_data, type = "response")
head(mod_data_pred)

ggplot(mod_data_pred, aes(x = park, y = pred)) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) +
  theme_set(theme_bw(base_size = 10)) +
  theme(legend.position = "top") +
  ylim(0, 1) +
  labs(x = "", y = "Predicted probabilty of blsiter rust") +
  scale_color_manual(values = c("gray20", "gray70"))

plot(glm_full, pch = 20, col = "black", lty = "dotted")

###############################################################

# MODEL SELECTION ---------------------------------------------------------

# WPBR sensitivity model: site_mod_int
mod <- glm_full
summary(mod)

# Save model as RDS file
#saveRDS(mod, file = "Output/Models/WPBR_sensitivity_model.rds")

# Test model loading
mod_load <- readRDS("Output/Models/WPBR_sensitivity_model.rds")
summary(mod_load)


###############################################################
# SECTION 4: MODEL VALIDATION ---------------------------------------------
###############################################################

# Region specific cross validations ---------------------------------------

# set seed for reproducible results
set.seed(10)

# set parameters for train function
ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = TRUE, repeats = 5)

head(mod_data)
unique(mod_data$park)
# subset data to park you want to look at
park_data <- subset(mod_data, park == "GYE")
head(park_data)

# train function for each park
set.seed(10)
park_train <- train(mod, data = park_data, 
                    method = "glm", family = "binomial", 
                    trControl = ctrl, tuneLength = 5)

confusionMatrix(park_train, norm = "none")


# RESULTS -----------------------------------------------------------------

# CAN: 93%
# CRLA: 90%
# GLAC: 93%
# GYE: 85%
#


# OTHER METHOD
head(mod_data)
park_data <- mod_data
head(park_data); tail(park_data)

park_pred <- cbind(park_data, pred = (predict(newdata = park_data, mod, 
                                              allow.new.levels = TRUE,
                                              type = "response")))
head(park_pred)
park_pred$pred_status <- round(park_pred$pred, digits = 0)
head(park_pred)

park_pred$eval <- ifelse(park_pred$recent_br_status > park_pred$pred_status, "FN",
                         ifelse(park_pred$recent_br_status < park_pred$pred_status, "FP",
                                ifelse(park_pred$pred_status == 1, "TP", "TN")))
head(park_pred)

lavo <- subset(park_pred, park == "LAVO")
lavo

park_pred_summary <- table(park_pred$eval, by = park_pred$park)
park_pred_summary
length(park_pred$transect_id)

table(park_pred$recent_br_status, by = park_pred$park)

###############################################################
# SECTION 5: MODELING WITH VPD --------------------------------------------
##############################################################

# Since average temperature (through SVP) and average relative humidity are used to calculate VPD,
# there are some strange results when you include interactions between the three terms
#     check out the dredge function and EXPLORE VPD MODELS section to see this

# SUMMARY:
#   BR ~ asRH * VPD: best model based on AIC. Better at predicting uninfected sites
#                   unsure how to interpret ecologically
#                   VPD seems to drive model
#   BR ~ asTEMP * asRH (original model): next best model based on AIC. A bit worse at predicting uninfected sites
#                   makes more sense ecologically?

head(mod_data)
tail(mod_data)
vpd_mod <- glm(data = mod_data, recent_br_status ~ vpd + asTEMP + asRH + asTEMP*asRH +asTEMP*vpd +asRH*vpd,
               family=binomial(link = "logit"))
summary(vpd_mod)

options(na.action = na.fail)
dredge(vpd_mod)
dredge_mod <- dredge(vpd_mod)

dredge_mod <- as.data.frame(dredge_mod)
dredge_mod

# original model: BR ~ temp and RH
org_mod <- glm(data = mod_data, recent_br_status ~ asTEMP * asRH,
               family=binomial(link = "logit"))
summary(org_mod)

# original model: BR ~ temp and RH
org_mod2 <- glm(data = mod_data, recent_br_status ~ I(asTEMP ^ 2) * asRH  + asTEMP * asRH,
               family=binomial(link = "logit"))
summary(org_mod2)

## temperature and humidity with interactions to compete against VPD
## try to improve confusion matrix
## transmission occurs on warm humid nights
## do seasonal phenology of august/september change in the future?
## Extract t and rh that is ideal for BR infection from historical model, check surrounding months

# VPD model: BR ~ VPD
# significantly higher AIC than mod
vpd_mod <- glm(data = mod_data, recent_br_status ~ vpd,
               family=binomial(link = "logit"))
summary(vpd_mod)

# full model: BR ~ temp + RH + VPD
# equivilant to mod
full_mod <- glm(data = mod_data, recent_br_status ~ vpd + asTEMP + asRH,
                family=binomial(link = "logit"))
summary(full_mod)

 # BR ~ temp + VPD
# slightly higher AIC than mod (but significant)
vpd_temp_mod <- glm(data = mod_data, recent_br_status ~ vpd + asTEMP,
                family=binomial(link = "logit"))
summary(vpd_temp_mod)

# VARIATIONS OF VPD MODELS ------------------------------------------------

# vpd
vpd_mod1 <- glm(data = mod_data, recent_br_status ~ vpd,
                family=binomial(link = "logit"))
summary(vpd_mod1)

tab_model(vpd_mod1)

# vpd*asRH
vpd_mod2 <- glm(data = mod_data, recent_br_status ~ asRH*vpd,
                family=binomial(link = "logit"))
summary(vpd_mod2)

tab_model(vpd_mod2)

# vpd + rh
vpd_mod3 <- glm(data = mod_data, recent_br_status ~ asRH+vpd,
                family=binomial(link = "logit"))
summary(vpd_mod3)

tab_model(vpd_mod3)

#######################################
# 10-FOLD CROSS VALIDATION ------------------------------------------------
#######################################
# K-fold cross validation: 10-fold cross-validation is most common

# set seed for reproducible results
set.seed(10)

# set parameters for train function
ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = TRUE, repeats = 5)

head(mod_data)

# change br_status into a factor for classifying
mod_data$br_fac <- as.character(mod_data$recent_br_status)
str(mod_data)

# train function for different versions of the model
# Can't evaluate GLMER (because of random effect)
set.seed(10)
glm_vpd1_train <- train(br_fac ~ vpd, data = mod_data, 
                        method = "glm", family = "binomial", 
                        trControl = ctrl, tuneLength = 5)

glm_vpd2_train <- train(br_fac ~ vpd + asRH*vpd, data = mod_data, 
                        method = "glm", family = "binomial", 
                        trControl = ctrl, tuneLength = 5)

glm_org_train <- train(br_fac ~ asRH * asTEMP, data = mod_data, 
                        method = "glm", family = "binomial", 
                        trControl = ctrl, tuneLength = 5)

glm_vpd_temp_train <- train(br_fac ~ asTEMP + vpd, data = mod_data, 
                            method = "glm", family = "binomial", 
                            trControl = ctrl, tuneLength = 5)

glm_full_mod_train <- train(br_fac ~ asRH + asTEMP + vpd, data = mod_data, 
                       method = "glm", family = "binomial", 
                       trControl = ctrl, tuneLength = 5)


# extract Accuracy and Kappa
print(glm_vpd1_train)
print(glm_vpd2_train)
print(glm_org_train)
print(glm_vpd_temp_train)
print(glm_full_mod_train)

# get confusion matrix
confusionMatrix(glm_vpd1_train, norm = "none")
confusionMatrix(glm_vpd2_train, norm = "none")
confusionMatrix(glm_org_train, norm = "none")
confusionMatrix(glm_vpd_temp_train, norm = "none")
confusionMatrix(glm_full_mod_train, norm = "none")

###############################################################
# EXPLORING VPD MODELS ----------------------------------------------------

# MODEL BY PARK -----------------------------------------------------------

unique(mod_data$park)
length(mod_data$park)

# CAN
CAN <- subset(mod_data, park == "CAN")
length(CAN$park)
head(CAN)

CAN_mod <- glm(data = CAN, recent_br_status ~ vpd + asTEMP + asRH,
               family=binomial(link = "logit"))
summary(CAN_mod)

# CRLA
CRLA <- subset(mod_data, park == "CRLA")
length(CRLA$park)
head(CRLA)

CRLA_mod <- glm(data = CRLA, recent_br_status ~ vpd + asTEMP + asRH,
                family=binomial(link = "logit"))
summary(CRLA_mod)

# GLAC
GLAC <- subset(mod_data, park == "GLAC")
length(GLAC$park)
head(GLAC)

GLAC_mod <- glm(data = GLAC, recent_br_status ~ vpd + asTEMP + asRH,
                family=binomial(link = "logit"))
summary(GLAC_mod)

# GYE
GYE <- subset(mod_data, park == "GYE")
length(GYE$park)
head(GYE)

GYE_mod <- glm(data = GYE, recent_br_status ~ vpd + asTEMP + asRH,
               family=binomial(link = "logit"))
summary(GYE_mod)

# KICA
KICA <- subset(mod_data, park == "KICA")
length(KICA$park)
head(KICA)

KICA_mod <- glm(data = KICA, recent_br_status ~ vpd + asTEMP + asRH,
                family=binomial(link = "logit"))
summary(KICA_mod)

# LAVO
LAVO <- subset(mod_data, park == "LAVO")
length(LAVO$park)
head(LAVO)

LAVO_mod <- glm(data = LAVO, recent_br_status ~ vpd + asTEMP + asRH,
                family=binomial(link = "logit"))
summary(LAVO_mod)

# MORA
MORA <- subset(mod_data, park == "MORA")
length(MORA$park)
head(MORA)

MORA_mod <- glm(data = MORA, recent_br_status ~ vpd + asTEMP + asRH,
                family=binomial(link = "logit"))
summary(MORA_mod)

# NOCA_LACH
NOCA_LACH <- subset(mod_data, park == "NOCA_LACH")
length(NOCA_LACH$park)
head(NOCA_LACH)

NOCA_LACH_mod <- glm(data = NOCA_LACH, recent_br_status ~ vpd + asTEMP + asRH,
                     family=binomial(link = "logit"))
summary(NOCA_LACH_mod)

# SEQU
SEQU <- subset(mod_data, park == "SEQU")
length(SEQU$park)
head(SEQU)

SEQU_mod <- glm(data = SEQU, recent_br_status ~ vpd + asTEMP + asRH,
                family=binomial(link = "logit"))
summary(SEQU_mod)

# YOSE
YOSE <- subset(mod_data, park == "YOSE")
length(YOSE$park)
head(YOSE)

YOSE_mod <- glm(data = YOSE, recent_br_status ~ vpd + asTEMP + asRH,
                family=binomial(link = "logit"))
summary(YOSE_mod)

###############################################################

#################################################
# VISUALIZE PREDICTIONS ----------------------------------------------------
################################################

summary(vpd_mod1)

# 2000-2020 mean rh and t and pred at each transect
# for points for transect predictions
head(mod_data)
mod_data_pred <- cbind(mod_data, pred = (predict(newdata = mod_data,
                                                 vpd_mod1, 
                                                 allow.new.levels = TRUE, 
                                                 type = "response")))
head(mod_data_pred)

# 2000-2020 mean rh, t, and pred at each park
# for park labels
park_clim_pred <- mod_data_pred %>% 
  group_by(park) %>% 
  dplyr::summarize(park_vpd = mean(vpd, na.rm = TRUE),
                   park_sevpd = std.error(vpd, na.rm = TRUE),
                   park_asRH = mean(asRH, na.rm = TRUE),
                   park_seRH = std.error(asRH, na.rm = TRUE),
                   park_mean_pred = mean(pred, na.rm = TRUE),
                   park_sepred = std.error(pred, na.rm = TRUE)) %>% 
  as.data.frame()

park_clim_pred

# generate data

# Create grid -------------------------------------------------------------

# model predictions across all combination of model inputs
summary(mod_data$vpd)

grid <- as.data.frame(expand.grid(vpd = seq(from = 296, to = 1395, by = 1),
                                  asRH = seq(from = 0.22, to = 0.74, by = 0.1)))
head(grid); tail(grid)

park_clim <- mod_data %>% 
  group_by(park) %>% 
  dplyr::summarize(vpd = mean(vpd, na.rm = TRUE),
                   asRH = mean(asRH, na.rm = TRUE)) %>% 
  as.data.frame()

park_clim


# Create predictions ------------------------------------------------------

# PREDICTOR = VPD #
pred_vpd_mod1 <- cbind(park_clim, pred = (predict(newdata = park_clim, vpd_mod1, 
                                                  allow.new.levels = TRUE, 
                                                  type = "response")))
pred_vpd_mod1

summary(park_clim$vpd)

grid_vpd_mod1 <- cbind(grid, pred = (predict(newdata = grid, vpd_mod1, 
                                             allow.new.levels = TRUE, 
                                             type = "response")))
head(grid_vpd_mod1)
head(mod_data_pred)
head(park_clim_pred)

park_cols <- c("GYE" = "#004586",
               "GLAC" = "#aecf00",
               "CAN" = "#83caff",
               "SEQU" = "#579d1c",
               "KICA" = "#ff420e",
               "YOSE" = "#0084d1",
               "LAVO" = "#314004",
               "CRLA" = "#7e0021",
               "MORA" = "#4b1f6f",
               "NOCA_LACH" = "#ff950e")
# curve
glm_vpd_curve <- ggplot() +
  geom_line(data = grid_vpd_mod1, aes(x = vpd, y = pred), size = 1, color = "black", alpha = 0.7) +
  geom_jitter(data = mod_data_pred, aes(x = vpd, y = pred, color = park),
              width = 3,
              size = 2, alpha = 0.8) +
  geom_label_repel(data = park_clim_pred, aes(x = park_vpd, y = park_mean_pred, label = park, color = park), 
                   size = 3.5, box.padding = 0.7, min.segment.length = 0, 
                   label.r = 0, label.size = 0)+
  scale_color_manual(values = park_cols, name = "Park or region", guide = "none") +
  xlab("Vapor Pressure Deficit") +
  ylab("Probability of WPBR infection") +
  theme_bw() +
  theme(text = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15, face = "plain"),
        axis.title.x = element_text(size = 15, vjust = -1.8, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(0.4, 0.3, 0.5, 0.5), "cm"),
        panel.background = element_rect(fill = "gray99"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = c(0,0), lim = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(expand = c(0,0), lim = c(290, 1400))

print(glm_vpd_curve)
#####

# PREDICTOR = VPD*RH #
pred_vpd_mod2 <- cbind(park_clim, pred = (predict(newdata = park_clim, vpd_mod2, 
                                                  allow.new.levels = TRUE, 
                                                  type = "response")))
pred_vpd_mod2
head(mod_data)

summary(park_clim$vpd)

grid_vpd_mod2 <- cbind(grid, pred = (predict(newdata = grid, vpd_mod2, 
                                             allow.new.levels = TRUE, 
                                             type = "response")))
head(grid_vpd_mod2)
park_clim$asRHcalc <- park_clim$asRH/100
park_clim

# heat map
htmap <- ggplot() +
  geom_tile(data = grid_vpd_mod2, aes(x = vpd, y = asRH, fill = pred)) +
  #annotate("rect", xmin = 36, xmax = 71, ymin = 7, ymax = 12, color = "gray30", 
  #alpha = 0, size = 1) +
  geom_point(data = mod_data, aes(x = vpd, y = asRHcalc, color = park), size = 3) +
  geom_label_repel(data = park_clim, aes(x = vpd, y = asRHcalc, label = park, color = park), 
                   size = 5, box.padding = 0.7, min.segment.length = 10, 
                   label.r = 0, label.size = 0)+
  scale_color_manual(values = park_cols, name = "", guide = "none") +
  scale_fill_gradient(name = "WPBR potential impact", low = "lightyellow1", high = "brown",
                      breaks = c(),
                      limits = c(),
                      labels = c()) +
  #theme_light() +
  xlab("VPD") +
  ylab("RH") +
  theme(text = element_text(size = 15, face = "plain"),
        legend.text = element_text(size = 15, face = "plain"),
        axis.title.x = element_text(size = 15, vjust = -1.8, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        plot.margin = unit(c(0.4, 0.3, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#scale_y_continuous(expand = c(0,0), breaks = c()) +
#scale_x_continuous(expand = c(0,0), breaks = c())

print(htmap)
############

# PREDICTOR = VPD + RH #
pred_vpd_mod3 <- cbind(park_clim, pred = (predict(newdata = park_clim, vpd_mod3, 
                                                  allow.new.levels = TRUE, 
                                                  type = "response")))
pred_vpd_mod3
head(mod_data)

summary(park_clim$vpd)

grid_vpd_mod3 <- cbind(grid, pred = (predict(newdata = grid, vpd_mod3, 
                                             allow.new.levels = TRUE, 
                                             type = "response")))
head(grid_vpd_mod3)
park_clim$asRHcalc <- park_clim$asRH/100
park_clim

# heat map
htmap <- ggplot() +
  geom_tile(data = grid_vpd_mod3, aes(x = vpd, y = asRH, fill = pred)) +
  #annotate("rect", xmin = 36, xmax = 71, ymin = 7, ymax = 12, color = "gray30", 
  #alpha = 0, size = 1) +
  geom_point(data = mod_data, aes(x = vpd, y = asRHcalc, color = park), size = 3) +
  geom_label_repel(data = park_clim, aes(x = vpd, y = asRHcalc, label = park, color = park), 
                   size = 5, box.padding = 0.7, min.segment.length = 10, 
                   label.r = 0, label.size = 0)+
  scale_color_manual(values = park_cols, name = "", guide = "none") +
  scale_fill_gradient(name = "WPBR potential impact", low = "lightyellow1", high = "brown",
                      breaks = c(),
                      limits = c(),
                      labels = c()) +
  #theme_light() +
  xlab("VPD") +
  ylab("RH") +
  theme(text = element_text(size = 15, face = "plain"),
        legend.text = element_text(size = 15, face = "plain"),
        axis.title.x = element_text(size = 15, vjust = -1.8, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        plot.margin = unit(c(0.4, 0.3, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#scale_y_continuous(expand = c(0,0), breaks = c()) +
#scale_x_continuous(expand = c(0,0), breaks = c())

print(htmap)

###############################################################

###############################################################
# SECTION 6: MODEL FIGURES -----------------------------------------------------------------
###############################################################

# GLM_FULL and GLM_NO_INT -------------------------------------------------

### HEAT MAP ### 

# model predictions across all combination of model inputs
summary(mod_data$asRH)
summary(mod_data$asTEMP)

grid <- as.data.frame(expand.grid(asRH = seq(from = 22.72, to = 74.79, by = 1), 
                                  asTEMP = seq(from = 5.2, to = 15.9, by = 0.5)))
head(grid); tail(grid)

# rescaled values
#grid$re_asRH <- scales::rescale(grid$asRH, newrange = c(0,1))
#grid$re_asTEMP <- scales::rescale(grid$asTEMP, newrange = c(0,1))

# GLM_FULL
grid <- cbind(grid, pred = (predict(newdata = grid, glm_full, 
                                    allow.new.levels = TRUE, 
                                    type = "response")))
head(grid)

# OR: GLM_NO_INT
#summary(glm_no_int)
#grid <- cbind(grid, pred = (predict(newdata = grid, glm_no_int, 
#                                    allow.new.levels = TRUE, 
#                                    type = "response")))
#head(grid)

# 2000-2020 mean across all transects in each park
# for park label in heat map
head(mod_data)

park_clim <- mod_data %>% 
  group_by(park) %>% 
  dplyr::summarize(asTEMP = mean(asTEMP, na.rm = TRUE),
                   asRH = mean(asRH, na.rm = TRUE)) %>% 
  as.data.frame()

park_clim
park_clim_pred <- cbind(park_clim, pred = (predict(newdata = park_clim, glm_full, 
                                                   allow.new.levels = TRUE, 
                                                   type = "response")))
park_clim_pred <- arrange(park_clim_pred, pred)
park_clim_pred
summary(park_clim$asTEMP)
summary(park_clim$asRH)

# PLOTTING
#ggthemes_data$calc
#show_col(calc_pal()(12))

park_cols <- c("GYE" = "#004586",
               "GLAC" = "#aecf00",
               "CAN" = "#83caff",
               "SEQU" = "#579d1c",
               "KICA" = "#ff420e",
               "YOSE" = "#0084d1",
               "LAVO" = "#314004",
               "CRLA" = "#7e0021",
               "MORA" = "#4b1f6f",
               "NOCA_LACH" = "#ff950e")

# heat map
htmap <- ggplot() +
  geom_tile(data = grid, aes(x = asRH, y = asTEMP, fill = pred)) +
  #annotate("rect", xmin = 36, xmax = 71, ymin = 7, ymax = 12, color = "gray30", 
  #alpha = 0, size = 1) +
  geom_point(data = mod_data, aes(x = asRH, y = asTEMP, color = park), size = 3) +
  geom_label_repel(data = park_clim, aes(x = asRH, y = asTEMP, label = park, color = park), 
                   size = 5, box.padding = 0.7, min.segment.length = 10, 
                   label.r = 0, label.size = 0)+
  scale_color_manual(values = park_cols, name = "", guide = "none") +
  scale_fill_gradient(name = "WPBR potential impact", low = "lightyellow1", high = "brown",
                      breaks = c(0, 0.5, 1),
                      limits = c(0,1),
                      labels = c(0,0.5,1)) +
  #theme_light() +
  xlab("Aug - Sept mean relative humidity (%)") +
  ylab("Aug - Sept mean temperature (°C)") +
  theme(text = element_text(size = 15, face = "plain"),
        legend.text = element_text(size = 15, face = "plain"),
        axis.title.x = element_text(size = 15, vjust = -1.8, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        plot.margin = unit(c(0.4, 0.3, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = c(0,0), breaks = c(5, 10, 15)) +
  scale_x_continuous(expand = c(0,0), breaks = c(25, 50, 75))

print(htmap)

## 3-D PLOT
head(grid)
RH <- grid[,1]
TEMP <- grid[,2]
Impact <- grid[,3]

p <- plot_ly(x = ~RH, y = ~TEMP, z = ~Impact, type = 'mesh3d',
             intensity = Impact, 
             #colorscale = 'Electric',
             colorscale = list(c(1, 'rgb(0, 0, 255)'), c(0.5, 'rgb(0, 255, 0)'), c(0, 'rgb(0, 0, 300)')), 
             showscale =TRUE)%>%
  layout(scene = list(
    xaxis = list(title = "Relative humitidy (%)"),
    yaxis = list(title = "Temperature (°C)"),
    zaxis = list(title = "WPBR potential impact"))#,
    #xaxis = list(tickfont = list(size = 15)), 
    #yaxis = list(titlefont = list(size = 25), title = "test")
  )

p
###############################################################

# GLM_RH ------------------------------------------------------------------
### Prediction Curve ###
### BR ~ RH ###
summary(glm_rh)

# 2000-2020 mean rh and t and pred at each transect
# for points for transect predictions
head(mod_data)
mod_data_pred <- cbind(mod_data, pred = (predict(newdata = mod_data,
                                                 glm_rh, 
                                                 allow.new.levels = TRUE, 
                                                 type = "response")))
head(mod_data_pred)

# 2000-2020 mean rh, t, and pred at each park
# for park labels
park_clim_pred <- mod_data_pred %>% 
  group_by(park) %>% 
  dplyr::summarize(park_asTEMP = mean(asTEMP, na.rm = TRUE),
                   park_seTEMP = std.error(asTEMP, na.rm = TRUE),
                   park_asRH = mean(asRH, na.rm = TRUE),
                   park_seRH = std.error(asRH, na.rm = TRUE),
                   park_mean_pred = mean(pred, na.rm = TRUE),
                   park_sepred = std.error(pred, na.rm = TRUE)) %>% 
  as.data.frame()

park_clim_pred

# generate data

# range used to build model
summary(mod_data$asRH)
summary(mod_data$asTEMP)

# grid
grid <- as.data.frame(expand.grid(asRH = seq(from = 22, to = 75, by = 0.5), 
                                  asTEMP = seq(from = 5, to = 16, by = 0.5)))
head(grid); tail(grid)


# apply trans_mod_rh
grid <- cbind(grid, pred = (predict(newdata = grid, glm_rh, 
                                    allow.new.levels = TRUE, 
                                    type = "response")))
head(grid)


# line graph
glm_rh_curve <- ggplot() +
  geom_line(data = grid, aes(x = asRH, y = pred), size = 1, color = "black", alpha = 0.7) +
  geom_jitter(data = mod_data_pred, aes(x = asRH, y = pred, color = park),
              width = 3,
              size = 2, alpha = 0.8) +
  geom_label_repel(data = park_clim_pred, aes(x = park_asRH, y = park_mean_pred, label = park, color = park), 
                   size = 3.5, box.padding = 0.7, min.segment.length = 0, 
                   label.r = 0, label.size = 0)+
  scale_color_manual(values = park_cols, name = "Park or region", guide = "none") +
  xlab("Aug - Sept mean relative humidity (%)") +
  ylab("WPBR potential impact") +
  theme_bw() +
  theme(text = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15, face = "plain"),
        axis.title.x = element_text(size = 15, vjust = -1.8, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(0.4, 0.3, 0.5, 0.5), "cm"),
        panel.background = element_rect(fill = "gray99"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = c(0,0), lim = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(expand = c(0,0), lim = c(20,80))

print(glm_rh_curve)

###############################################################
# ADDITIONAL FIGURES ------------------------------------------------------
###############################################################
# GLMER_PARK --------------------------------------------------------------
### PREDICTION CURVES ###

summary(mod_data$asTEMP)
summary(mod_data$asRH)

park_names <- as.vector(unique(mod_data$park))
park_names

grid <- as.data.frame(expand.grid(asRH = seq(from = 22, to = 75, by = 1), 
                                  asTEMP = seq(from = 5, to = 16, by = 1),
                                  park = rep(park_names, each = 648)))

grid$re_asRH <- scales::rescale(grid$asRH, newrange = c(0,1))
grid$re_asTEMP <- scales::rescale(grid$asTEMP, newrange = c(0,1))
head(grid);tail(grid)

# GLMER_PARK
grid <- cbind(grid, pred = (predict(newdata = grid, glmer_park, 
                                    allow.new.levels = TRUE, 
                                    type = "response")))
head(grid);tail(grid)

park_cols <- c("LAVO" = "#314004",
               "CRLA" = "#7e0021",
               "NOCA_LACH" = "#ff950e",
               "MORA" = "#4b1f6f",
               "CAN" = "#83caff",
               "GLAC" = "#aecf00",
               "GYE" = "#004586",
               "KICA" = "#ff420e",
               "YOSE" = "#0084d1",
               "SEQU" = "#579d1c")

# line graph
glmer_park_curve <- ggplot() +
  geom_line(data = grid, aes(x = asRH, y = pred, color = park), size = 1.1) +
  scale_color_manual(values = park_cols, name = "Park or region") +
  xlab("Aug - Sept mean relative humidity (%)") +
  ylab("WPBR infection probability") +
  theme_bw() +
  theme(text = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15, face = "plain"),
        axis.title.x = element_text(size = 15, vjust = -1.8, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(0.4, 0.3, 0.5, 0.5), "cm"),
        panel.background = element_rect(fill = "gray99"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = c(0,0), lim = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(expand = c(0,0), lim = c(20,80))

print(glmer_park_curve)



# COMPARING PREDICTIONS WITH REALITY --------------------------------------

### GROUPED BAR CHART ###
# PROPORTION OF TRANSECTS AND TREES INFECTED AND PREDICTED PROBABILITY
# paired bar chart for proportion of transects infected
# can also change to proportion of trees infected

# apply model to each transect
head(mod_data)
mod_data_pred <- cbind(mod_data, pred = (predict(newdata = mod_data,
                                                 glm_rh, 
                                                 allow.new.levels = TRUE, 
                                                 type = "response")))

head(mod_data_pred)

# 2000-2020 mean rh, t, and pred at each park
park_clim_pred <- mod_data_pred %>% 
  group_by(park) %>% 
  dplyr::summarize(park_asTEMP = mean(asTEMP, na.rm = TRUE),
                   park_seTEMP = std.error(asTEMP, na.rm = TRUE),
                   park_asRH = mean(asRH, na.rm = TRUE),
                   park_seRH = std.error(asRH, na.rm = TRUE),
                   park_mean_pred = mean(pred, na.rm = TRUE),
                   park_sepred = std.error(pred, na.rm = TRUE)) %>% 
  as.data.frame()

park_clim_pred

# calculate proportion of sites infected
head(mod_data)

prop_site_inf <- mod_data %>%
  group_by(park) %>%
  dplyr::summarize(site_count = n(),
                   inf_count = sum(recent_br_status),
                   prop_site_inf = inf_count/site_count) %>%
  as.data.frame()

prop_site_inf

prop_site_inf <- select(prop_site_inf, c("park", "prop_site_inf"))
prop_site_inf$se_site <- 0
prop_site_inf

# extract proportion of trees infected
# summarize blister rust severity in each park
head(tree_data)

prop_tree_inf <- tree_data %>%
  group_by(park) %>%
  dplyr::summarize(transect_count= length(unique(transect_id)),
                   tree_count = n(),
                   species_count = length(unique(tree_species)),
                   tree_inf_count = sum(br_status),
                   prop_tree_inf = tree_inf_count/tree_count) %>%
  as.data.frame()
print(prop_tree_inf)

prop_tree_inf <- select(prop_tree_inf, c("park", "prop_tree_inf"))
prop_tree_inf

# merge prop of transects infected with br predictions
prop_site_inf
prop_tree_inf
park_pred <- select(park_clim_pred, c("park", "park_mean_pred", "park_sepred"))
park_pred

br <- merge(park_pred, prop_site_inf, by = "park")
br

br <- select(br, c("park", "park_mean_pred", "prop_site_inf"))
br

# melt to format for grouped bar chart
br_melt <- melt(br, id = "park")
br_melt

# add se data
prop_site_inf
se <- as.vector(c(park_pred$park_sepred, prop_site_inf$se_site))
se

br_melt$se <- se
br_melt

# add prop trees infected
prop_tree_inf
prop_tree_inf <- subset(prop_tree_inf, park != "CRMO_CRMP")
prop_tree_inf <- subset(prop_tree_inf, park != "ROMO")
prop_tree_inf <- subset(prop_tree_inf, park != "GRBA")
prop_tree_inf <- subset(prop_tree_inf, park != "GRSA")
prop_tree_inf

br_melt$prop_tree_inf <- c(rep(0, 10), prop_tree_inf$prop_tree_inf)
br_melt
br_melt$pattern <- "pattern"
br_melt

br_comp <- ggplot(data = br_melt, aes(x = reorder(park, value), y = value, 
                                      fill = variable, pattern = pattern)) +
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  geom_errorbar(aes(ymin = value - se, ymax = value + se), 
                position = "dodge", size = 0.5) +
  scale_y_continuous(expand = c(0,0), lim = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_bar_pattern(stat = "identity", position = "dodge",
                   aes(x = reorder(park, value), y = prop_tree_inf),
                   color = "black", pattern_fill = "black", pattern_angle = 45,
                   pattern_density = 0.3)+scale_fill_manual(name = "Value", values = c("park_mean_pred" = "darkslategray4",
                                                                                       "prop_site_inf" = "firebrick4"),
                                                            labels = c("potential impact", "proportion of transects infected")) +
  scale_pattern_manual(name = "", values = c(pattern = "stripe"), labels = c("proportion of trees infected")) +
  xlab("Park or region") +
  ylab("Value")+
  #labs(fill = "Value")+
  theme_bw()+
  theme(text = element_text(size = 15, face = "bold"),
        #plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "right",
        #legend.title = element_text(size =),
        #legend.key = element_rect(fill = "black"),
        legend.key.size = unit(0.6, "cm"),
        legend.margin = margin(r = 8, l = 8, t = 0, b = 8),
        legend.text = element_text(face = "plain"),
        axis.title.x = element_text(size = 15, vjust = -1, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(0.75, 0, 0.5, 0.75), "cm"),
        panel.background = element_rect(fill = "gray99"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


br_comp

# -------------------------------------------------------------------------

### CLIMATE AT TRANSECTS ###
head(br_summary_trans);tail(br_summary_trans)
min(br_summary_trans$mean_t, na.rm = TRUE);max(br_summary_trans$mean_t, na.rm = TRUE)
min(br_summary_trans$mean_rh, na.rm = TRUE); max(br_summary_trans$mean_rh, na.rm = TRUE)
head(br_summary_trans)

br_summary_trans$re_asRH <- rescale(br_summary_trans$mean_t, newrange = c(0,1))
head(br_summary_trans)

br_summary_trans$pred <- predict(trans_mod, br_summary_trans, type = "response", 
                                 allow.new.levels = TRUE)
head(br_summary_trans)

summary(trans_mod)

cc_trans <- ggplot(br_summary_trans, aes(x = park, y = mean_rh, color = park)) +
  geom_point(size = 3, alpha = 0.2)+
  stat_ellipse(aes(color = park), size = 1)+
  labs(title = "BR climate correlates at WBP transects calculated from 1980-2020 mean")+
  theme_bw()+
  theme(text = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.3, vjust = 5, size = 12),
        legend.position = "right",
        legend.text = (element_text(size = 15, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 15, vjust = -1.8, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0.75, 0.75, 0.75, 0.75), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("rect", xmin = 7, xmax = 11, ymin = 50, ymax = 100, fill = "gray50", alpha = 0.3) +
  #scale_color_manual(name = "Park", values = park_cols) +
  scale_color_calc(name = "Park") +
  scale_y_continuous(name = "Aug - Sept mean relative humidity (%)", 
                     expand = c(0,0), lim = c(0,100), breaks = c(0, 25, 50, 75, 100))

cc_trans

### CLIMATE BY YEAR ###
head(br_summary_year);tail(br_summary_year)
min(br_summary_year$mean_t, na.rm = TRUE);max(br_summary_year$mean_t, na.rm = TRUE)
min(br_summary_year$mean_rh, na.rm = TRUE); max(br_summary_year$mean_rh, na.rm = TRUE)
head(br_summary_year)

br_summary_year$re_asRH <- rescale(br_summary_year$mean_t, newrange = c(0,1))
head(br_summary_year)

br_summary_year$pred <- predict(trans_mod, br_summary_year, type = "response", 
                                allow.new.levels = TRUE)
head(br_summary_year)

cc_year <- ggplot(br_summary_year, aes(x = year, y = pred, color = park)) +
  geom_line(size = 1.2, alpha = 0.8)+
  #stat_ellipse(aes(color = park), size = 1)+
  labs(title = "BR vulnerability by year")+
  theme_bw()+
  theme(text = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.3, vjust = 5, size = 12),
        legend.position = "right",
        legend.text = (element_text(size = 15, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 15, vjust = -1.8, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0.75, 0.75, 0.75, 0.75), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #annotate("rect", xmin = 7, xmax = 11, ymin = 50, ymax = 100, fill = "gray50", alpha = 0.3) +
  #scale_color_manual(name = "Park", values = park_cols) +
  scale_color_calc(name = "Park") +
  scale_y_continuous(name = "Blister rust vulnerability", 
                     expand = c(0,0), lim = c(0,1), breaks = c(0, 0.25, 0.50, 0.75, 0.100)) +
  scale_x_continuous(expand = c(0,0), lim = c(1979, 2021))

print(cc_year)


# -------------------------------------------------------------------------

# BR BY YEAR --------------------------------------------------------------

correl_trans$period <- cut(correl_trans$year, breaks = c(1979, 1989, 1999, 2009, 2019, 2020),
                           labels = c("1980-1989", "1990-1999", "2000-2009","2010-2019", "2020"))
head(correl_trans)

# calculate BR across year by each park
correl_trans$re_asRH <- rescale(correl_trans$rh, newrange = c(0,1))
correl_trans$re_asTEMP <- rescale(correl_trans$t, newrange = c(0,1))

correl_trans$pred <- predict(trans_mod_full, correl_trans, type = "response")
head(correl_trans)

br_change <- subset(correl_trans, period != "2020")
head(br_change)
br_change

br_change <- br_change %>% 
  group_by(park, period) %>% 
  dplyr::summarize(t = mean(t, na.rm = TRUE),
                   rh = mean(rh, na.rm = TRUE),
                   vul = mean(pred, na.rm = TRUE),
                   se_vul = std.error(pred)) %>% 
  as.data.frame()

head(br_change)

# plot: change from 1980-1989 to 2010-2019 means- arrows
br_year <- ggplot() +
  #stat_ellipse(data = clim_year, aes(x = t, y = rh, color = park), size = 1, alpha = 0.5)+
  #geom_point(data = clim_year, aes(x = t, y = rh, fill = park), 
  #shape = 21, color = "black", size = 2, alpha = 0.6)+
  geom_line(data = br_change, aes(x = period, y = vul, group = park, color = park),
            size = 1.5)+
  #labs(title = "Change in BR climate correlates calculated from 1980-1989 to 2010-2019 means")+
  theme_bw()+
  theme(text = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.1, vjust = 5, size = 12),
        legend.position = "right",
        legend.text = (element_text(size = 15, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 15, vjust = -1.8, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0.75, 0.75, 0.75, 0.75), "cm"),
        panel.background = element_rect(fill = "gray99"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #annotate("rect", xmin = 7, xmax = 11, ymin = 50, ymax = 100, fill = "red", alpha = 0.3) +
  scale_color_manual(values = park_cols, name = "Park or region") +
  #scale_fill_manual(values = park_cols, guide = "none") +
  scale_y_continuous(name = "Blister rust vulnerability", 
                     expand = c(0,0), lim = c(0,1), breaks = c(0, 0.25, 0.50, 0.75, 1))

br_year

length(unique(correl_trans$transect_id))

###############################################################

# -------------------------------------------------------------------------
### CLIMATE AND BR BY PARK ###

# FORMAT DATA
head(correl_trans)

correl_trans$re_asRH <- rescale(correl_trans$rh, newrange = c(0,1))
head(correl_trans)

correl_trans$pred <- predict(trans_mod, correl_trans, type = "response", 
                             allow.new.levels = TRUE)

head(correl_trans)

# summarize pred by year
pred_year <- correl_trans %>% 
  group_by(park, year) %>% 
  summarize(num_trans = n(),
            mean_pred = mean(pred, na.rm = TRUE),
            min_pred = min(pred, na.rm = TRUE),
            max_pred = max(pred, na.rm = TRUE),
            se_pred = std.error(pred, na.rm = TRUE),
            asTEMP = mean(t, na.rm = TRUE),
            min_asTEMP = min(t, na.rm = TRUE),
            max_asTEMP = max(t, na.rm = TRUE),
            se_asTEMP = std.error(t, na.rm = TRUE),
            asRH = mean(rh, na.rm = TRUE),
            min_asRH = min(rh, na.rm = TRUE),
            max_asRH = max(rh, na.rm = TRUE),
            se_asRH = std.error(rh, na.rm = TRUE)) %>% 
  as.data.frame()

head(pred_year)

coeff_rh <- 100

br_prop_trans


# FIGURES

# GRYN #
GRYN <- ggplot(data = subset(pred_year, park == "GRYN"), 
               aes(x = year)) +
  geom_boxplot(data = subset(correl_trans, park == "GRYN"), 
               aes(x = year, y = rh/coeff_rh, group = year), color = "#80B1D3", 
               alpha = 0) +
  geom_smooth(aes(y = mean_pred), color = "#B3B3B3", fill = "#B3B3B3", 
              method = lm, se = TRUE)+
  geom_errorbar(aes(ymin = mean_pred + se_pred, 
                    ymax = mean_pred - se_pred), color = "black", 
                size = 0.7, width = 0.7) +
  geom_point(aes(y = mean_pred), fill = "#B3B3B3", color = "black", 
             size = 3.3, pch = 21)+
  labs(title = "A. GRYN (n = 184)")+
  theme_classic()+
  theme(text = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.03, vjust = -2, size = 14),
        legend.position = "none",
        legend.text = (element_text(size = 12, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 12, vjust = -2.5, face = "bold"),
        axis.title.y = element_text(size = 12, vjust = 4, face = "bold", color = "gray50"),
        axis.title.y.right = element_text(size = 12, vjust = 1.5, face = "bold", color = "#80B1D3"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "",
                     sec.axis = sec_axis(~.*coeff_rh, name = "", breaks = c(0,25,50,75, 100)),
                     expand = c(0,0), lim = c(-0.02,1.02), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(name = "", expand = c(0,0), lim = c(1979, 2021))

print(GRYN)


# GLAC #
GLAC <- ggplot(data = subset(pred_year, park == "GLAC"), 
               aes(x = year)) +
  geom_boxplot(data = subset(correl_trans, park == "GLAC"), 
               aes(x = year, y = rh/coeff_rh, group = year), color = "#80B1D3", 
               alpha = 0) +
  geom_smooth(aes(y = mean_pred), color = "#B3B3B3", fill = "#B3B3B3", 
              method = lm, se = TRUE)+
  geom_errorbar(aes(ymin = mean_pred + se_pred, 
                    ymax = mean_pred - se_pred), color = "black", 
                size = 0.7, width = 0.7) +
  geom_point(aes(y = mean_pred), fill = "#B3B3B3", color = "black", 
             size = 3.3, pch = 21)+
  labs(title = "A. GLAC (n = 66)")+
  theme_classic()+
  theme(text = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.03, vjust = -2, size = 14),
        legend.position = "none",
        legend.text = (element_text(size = 12, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 12, vjust = -2.5, face = "bold"),
        axis.title.y = element_text(size = 12, vjust = 4, face = "bold", color = "gray50"),
        axis.title.y.right = element_text(size = 12, vjust = 1.5, face = "bold", color = "#80B1D3"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "",
                     sec.axis = sec_axis(~.*coeff_rh, name = "", breaks = c(0,25,50,75, 100)),
                     expand = c(0,0), lim = c(-0.02,1.02), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(name = "", expand = c(0,0), lim = c(1979, 2021))

print(GLAC)


# CAN #
CAN <- ggplot(data = subset(pred_year, park == "CAN"), 
              aes(x = year)) +
  geom_boxplot(data = subset(correl_trans, park == "CAN"), 
               aes(x = year, y = rh/coeff_rh, group = year), color = "#80B1D3", 
               alpha = 0) +
  geom_smooth(aes(y = mean_pred), color = "#B3B3B3", fill = "#B3B3B3", 
              method = lm, se = TRUE)+
  geom_errorbar(aes(ymin = mean_pred + se_pred, 
                    ymax = mean_pred - se_pred), color = "black", 
                size = 0.7, width = 0.7) +
  geom_point(aes(y = mean_pred), fill = "#B3B3B3", color = "black", 
             size = 3.3, pch = 21)+
  labs(title = "A. CAN (n = 112)")+
  theme_classic()+
  theme(text = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.03, vjust = -2, size = 14),
        legend.position = "none",
        legend.text = (element_text(size = 12, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 12, vjust = -2.5, face = "bold"),
        axis.title.y = element_text(size = 12, vjust = 4, face = "bold", color = "gray50"),
        axis.title.y.right = element_text(size = 12, vjust = 1.5, face = "bold", color = "#80B1D3"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "",
                     sec.axis = sec_axis(~.*coeff_rh, name = "", breaks = c(0,25,50,75, 100)),
                     expand = c(0,0), lim = c(-0.02,1.02), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(name = "", expand = c(0,0), lim = c(1979, 2021))

print(CAN)


# NOCA_LACH #
NOCA_LACH <- ggplot(data = subset(pred_year, park == "NOCA_LACH"), 
                    aes(x = year)) +
  geom_boxplot(data = subset(correl_trans, park == "NOCA_LACH"), 
               aes(x = year, y = rh/coeff_rh, group = year), color = "#80B1D3", 
               alpha = 0) +
  geom_smooth(aes(y = mean_pred), color = "#B3B3B3", fill = "#B3B3B3", 
              method = lm, se = TRUE)+
  geom_errorbar(aes(ymin = mean_pred + se_pred, 
                    ymax = mean_pred - se_pred), color = "black", 
                size = 0.7, width = 0.7) +
  geom_point(aes(y = mean_pred), fill = "#B3B3B3", color = "black", 
             size = 3.3, pch = 21)+
  labs(title = "A. NOCA_LACH (n = 51)")+
  theme_classic()+
  theme(text = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.03, vjust = -2, size = 14),
        legend.position = "none",
        legend.text = (element_text(size = 12, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 12, vjust = -2.5, face = "bold"),
        axis.title.y = element_text(size = 12, vjust = 4, face = "bold", color = "gray50"),
        axis.title.y.right = element_text(size = 12, vjust = 1.5, face = "bold", color = "#80B1D3"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "",
                     sec.axis = sec_axis(~.*coeff_rh, name = "", breaks = c(0,25,50,75, 100)),
                     expand = c(0,0), lim = c(-0.02,1.02), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(name = "", expand = c(0,0), lim = c(1979, 2021))

print(NOCA_LACH)


# MORA #
MORA <- ggplot(data = subset(pred_year, park == "MORA"), 
               aes(x = year)) +
  geom_boxplot(data = subset(correl_trans, park == "MORA"), 
               aes(x = year, y = rh/coeff_rh, group = year), color = "#80B1D3", 
               alpha = 0) +
  geom_smooth(aes(y = mean_pred), color = "#B3B3B3", fill = "#B3B3B3", 
              method = lm, se = TRUE)+
  geom_errorbar(aes(ymin = mean_pred + se_pred, 
                    ymax = mean_pred - se_pred), color = "black", 
                size = 0.7, width = 0.7) +
  geom_point(aes(y = mean_pred), fill = "#B3B3B3", color = "black", 
             size = 3.3, pch = 21)+
  labs(title = "A. MORA (n = 31)")+
  theme_classic()+
  theme(text = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.03, vjust = -2, size = 14),
        legend.position = "none",
        legend.text = (element_text(size = 12, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 12, vjust = -2.5, face = "bold"),
        axis.title.y = element_text(size = 12, vjust = 4, face = "bold", color = "gray50"),
        axis.title.y.right = element_text(size = 12, vjust = 1.5, face = "bold", color = "#80B1D3"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "",
                     sec.axis = sec_axis(~.*coeff_rh, name = "", breaks = c(0,25,50,75, 100)),
                     expand = c(0,0), lim = c(-0.02,1.02), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(name = "", expand = c(0,0), lim = c(1979, 2021))

print(MORA)


# CRLA #
CRLA <- ggplot(data = subset(pred_year, park == "CRLA"), 
               aes(x = year)) +
  geom_boxplot(data = subset(correl_trans, park == "CRLA"), 
               aes(x = year, y = rh/coeff_rh, group = year), color = "#80B1D3", 
               alpha = 0) +
  geom_smooth(aes(y = mean_pred), color = "#B3B3B3", fill = "#B3B3B3", 
              method = lm, se = TRUE)+
  geom_errorbar(aes(ymin = mean_pred + se_pred, 
                    ymax = mean_pred - se_pred), color = "black", 
                size = 0.7, width = 0.7) +
  geom_point(aes(y = mean_pred), fill = "#B3B3B3", color = "black", 
             size = 3.3, pch = 21)+
  labs(title = "A. CRLA (n = 30)")+
  theme_classic()+
  theme(text = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.03, vjust = -2, size = 14),
        legend.position = "none",
        legend.text = (element_text(size = 12, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 12, vjust = -2.5, face = "bold"),
        axis.title.y = element_text(size = 12, vjust = 4, face = "bold", color = "gray50"),
        axis.title.y.right = element_text(size = 12, vjust = 1.5, face = "bold", color = "#80B1D3"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "",
                     sec.axis = sec_axis(~.*coeff_rh, name = "", breaks = c(0,25,50,75, 100)),
                     expand = c(0,0), lim = c(-0.02,1.02), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(name = "", expand = c(0,0), lim = c(1979, 2021))

print(CRLA)


# SEQU #
SEQU <- ggplot(data = subset(pred_year, park == "SEQU"), 
               aes(x = year)) +
  geom_boxplot(data = subset(correl_trans, park == "SEQU"), 
               aes(x = year, y = rh/coeff_rh, group = year), color = "#80B1D3", 
               alpha = 0) +
  geom_smooth(aes(y = mean_pred), color = "#B3B3B3", fill = "#B3B3B3", 
              method = lm, se = TRUE)+
  geom_errorbar(aes(ymin = mean_pred + se_pred, 
                    ymax = mean_pred - se_pred), color = "black", 
                size = 0.7, width = 0.7) +
  geom_point(aes(y = mean_pred), fill = "#B3B3B3", color = "black", 
             size = 3.3, pch = 21)+
  labs(title = "A. SEQU (n = 30)")+
  theme_classic()+
  theme(text = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.03, vjust = -2, size = 14),
        legend.position = "none",
        legend.text = (element_text(size = 12, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 12, vjust = -2.5, face = "bold"),
        axis.title.y = element_text(size = 12, vjust = 4, face = "bold", color = "gray50"),
        axis.title.y.right = element_text(size = 12, vjust = 1.5, face = "bold", color = "#80B1D3"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "",
                     sec.axis = sec_axis(~.*coeff_rh, name = "", breaks = c(0,25,50,75, 100)),
                     expand = c(0,0), lim = c(-0.02,1.02), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(name = "", expand = c(0,0), lim = c(1979, 2021))

print(SEQU)


# KICA #
KICA <- ggplot(data = subset(pred_year, park == "KICA"), 
               aes(x = year)) +
  geom_boxplot(data = subset(correl_trans, park == "KICA"), 
               aes(x = year, y = rh/coeff_rh, group = year), color = "#80B1D3", 
               alpha = 0) +
  geom_smooth(aes(y = mean_pred), color = "#B3B3B3", fill = "#B3B3B3", 
              method = lm, se = TRUE)+
  geom_errorbar(aes(ymin = mean_pred + se_pred, 
                    ymax = mean_pred - se_pred), color = "black", 
                size = 0.7, width = 0.7) +
  geom_point(aes(y = mean_pred), fill = "#B3B3B3", color = "black", 
             size = 3.3, pch = 21)+
  labs(title = "A. KICA (n = 35)")+
  theme_classic()+
  theme(text = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.03, vjust = -2, size = 14),
        legend.position = "none",
        legend.text = (element_text(size = 12, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 12, vjust = -2.5, face = "bold"),
        axis.title.y = element_text(size = 12, vjust = 4, face = "bold", color = "gray50"),
        axis.title.y.right = element_text(size = 12, vjust = 1.5, face = "bold", color = "#80B1D3"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "",
                     sec.axis = sec_axis(~.*coeff_rh, name = "", breaks = c(0,25,50,75, 100)),
                     expand = c(0,0), lim = c(-0.02,1.02), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(name = "", expand = c(0,0), lim = c(1979, 2021))

print(KICA)


# YOSE #
YOSE <- ggplot(data = subset(pred_year, park == "YOSE"), 
               aes(x = year)) +
  geom_boxplot(data = subset(correl_trans, park == "YOSE"), 
               aes(x = year, y = rh/coeff_rh, group = year), color = "#80B1D3", 
               alpha = 0) +
  geom_smooth(aes(y = mean_pred), color = "#B3B3B3", fill = "#B3B3B3", 
              method = lm, se = TRUE)+
  geom_errorbar(aes(ymin = mean_pred + se_pred, 
                    ymax = mean_pred - se_pred), color = "black", 
                size = 0.7, width = 0.7) +
  geom_point(aes(y = mean_pred), fill = "#B3B3B3", color = "black", 
             size = 3.3, pch = 21)+
  labs(title = "A. YOSE (n = 35)")+
  theme_classic()+
  theme(text = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.03, vjust = -2, size = 14),
        legend.position = "none",
        legend.text = (element_text(size = 12, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 12, vjust = -2.5, face = "bold"),
        axis.title.y = element_text(size = 12, vjust = 4, face = "bold", color = "gray50"),
        axis.title.y.right = element_text(size = 12, vjust = 1.5, face = "bold", color = "#80B1D3"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "",
                     sec.axis = sec_axis(~.*coeff_rh, name = "", breaks = c(0,25,50,75, 100)),
                     expand = c(0,0), lim = c(-0.02,1.02), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(name = "", expand = c(0,0), lim = c(1979, 2021))

print(YOSE)


# LAVO #
LAVO <- ggplot(data = subset(pred_year, park == "LAVO"), 
               aes(x = year)) +
  geom_boxplot(data = subset(correl_trans, park == "LAVO"), 
               aes(x = year, y = rh/coeff_rh, group = year), color = "#80B1D3", 
               alpha = 0) +
  geom_smooth(aes(y = mean_pred), color = "#B3B3B3", fill = "#B3B3B3", 
              method = lm, se = TRUE)+
  geom_errorbar(aes(ymin = mean_pred + se_pred, 
                    ymax = mean_pred - se_pred), color = "black", 
                size = 0.7, width = 0.7) +
  geom_point(aes(y = mean_pred), fill = "#B3B3B3", color = "black", 
             size = 3.3, pch = 21)+
  labs(title = "A. LAVO (n = 30)")+
  theme_classic()+
  theme(text = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.03, vjust = -2, size = 14),
        legend.position = "none",
        legend.text = (element_text(size = 12, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 12, vjust = -2.5, face = "bold"),
        axis.title.y = element_text(size = 12, vjust = 4, face = "bold", color = "gray50"),
        axis.title.y.right = element_text(size = 12, vjust = 1.5, face = "bold", color = "#80B1D3"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "",
                     sec.axis = sec_axis(~.*coeff_rh, name = "", breaks = c(0,25,50,75, 100)),
                     expand = c(0,0), lim = c(-0.02,1.02), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(name = "", expand = c(0,0), lim = c(1979, 2021))

print(LAVO)


# COMBINE FIGURES
park_fig <- ggarrange(CAN, GLAC, GRYN, CRLA, SEQU, KICA, YOSE, LAVO, MORA, NOCA_LACH, 
                      nrow = 5, ncol = 2)
print(park_fig)


park_fig <- annotate_figure(park_fig, left = text_grob("Mean blister rust vulnerability", 
                                                       color = "gray50", face = "bold", size = 15, rot = 90),
                            right = text_grob("Relative humidity (%)", 
                                              color = "#80B1D3", face = "bold", size = 15, rot = -90))
print(park_fig)


###############################################################


# CLIMATE FIGURES ---------------------------------------------------------

### CHANGE IN CLIMATE CORRELATES: 1980-89 TO 2010-2019 MEANS
head(correl_trans)

# calculate 1980-1989 and 2010-2019 means for each park
correl_trans$period <- cut(correl_trans$year, breaks = c(1979, 1989, 1999, 2009, 2019, 2020),
                           labels = c("1980-1989", "1990-1999", "2000-2009","2010-2019", "2020"))
head(correl_trans)

# for 1980-98 -> 2010-19 arrows
clim_change <- subset(correl_trans, period == "1980-1989" | period == "2010-2019")

head(clim_change)

clim_change <- clim_change %>% 
  group_by(park, period) %>% 
  dplyr::summarize(t = mean(t, na.rm = TRUE),
                   rh = mean(rh, na.rm = TRUE)) %>% 
  as.data.frame()

clim_change

# for park each year point
clim_year <- correl_trans %>% 
  group_by(park, year) %>% 
  dplyr::summarize(t = mean(t, na.rm = TRUE),
                   rh = mean(rh, na.rm = TRUE)) %>% 
  as.data.frame()

clim_year
length(clim_year$park)


# for transect point
clim_trans <- correl_trans %>% 
  group_by(park, transect_id) %>% 
  dplyr::summarize(t = mean(t, na.rm = TRUE),
                   rh = mean(rh, na.rm = TRUE)) %>% 
  as.data.frame()

clim_trans

# ARROW PLOT
# colors from calc ggtheme
park_cols <- c("MORA" = "#4b1f6f",
               "GLAC" = "#aecf00",
               "CAN" = "#83caff",
               "NOCA_LACH" = "#ff950e",
               "GYE" = "#004586",
               #"BLM_MT" = "gray40",
               "CRLA" = "#7e0021",
               "LAVO" = "#314004",
               "YOSE" = "#0084d1",
               "KICA" = "#ff420e",
               "SEQU" = "#579d1c")

# remove BLM MT
clim_change <- subset(clim_change, clim_change$park != "BLM_MT")

# rename GRYN to GYE
clim_change$park <- recode(clim_change$park, "GRYN" = "GYE")

# plot: change from 1980-1989 to 2010-2019 means- arrows
cc_arrow <- ggplot() +
  #stat_ellipse(data = clim_year, aes(x = t, y = rh, color = park), size = 1, alpha = 0.5)+
  #geom_point(data = clim_year, aes(x = t, y = rh, fill = park), 
  #shape = 21, color = "black", size = 2, alpha = 0.6)+
  geom_line(data = clim_change, aes(x = t, y = rh, color = park),
            size = 1.5,
            arrow = arrow(length=unit(0.35,"cm"), ends = "last", type = "open"))+
  #labs(title = "Change in BR climate correlates calculated from 1980-1989 to 2010-2019 means")+
  theme_bw()+
  theme(text = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.1, vjust = 5, size = 12),
        legend.position = "right",
        legend.text = (element_text(size = 15, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 15, vjust = -1.8, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0.75, 0.75, 0.75, 0.75), "cm"),
        panel.background = element_rect(fill = "gray99"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #annotate("rect", xmin = 7, xmax = 11, ymin = 50, ymax = 100, fill = "red", alpha = 0.3) +
  scale_color_manual(values = park_cols, name = "Park or region") +
  scale_fill_manual(values = park_cols, guide = "none") +
  scale_y_continuous(name = "Aug - Sept mean relative humidity (%)", 
                     expand = c(0,0), lim = c(0,100), breaks = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(name = "Aug - Sept mean temperature (°C)", 
                     expand = c(0,0), lim = c(5,15), breaks = c(5, 10, 15))

cc_arrow

# all info
trans_year <- ggplot() +
  #stat_ellipse(data = clim_year, aes(x = t, y = rh, color = park), size = 1, alpha = 0.5)+
  geom_point(data = correl_trans, aes(x = t, y = rh, color = park), size = 3, alpha = 0.3)+
  geom_line(data = clim_change, aes(x = t, y = rh, color = park),
            size = 1.5,
            arrow = arrow(length=unit(0.35,"cm"), ends="last", type = "closed"))+
  labs(title = "Change in BR climate correlates calculated from 1980-1989 to 2010-2019 means")+
  theme_bw()+
  theme(text = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.1, vjust = 5, size = 12),
        legend.position = "right",
        legend.text = (element_text(size = 15, face = "plain")),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 15, vjust = -1.8, face = "bold"),
        axis.title.y = element_text(size = 15, vjust = 4, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        panel.background = element_rect(fill = "gray99"),
        plot.margin = unit(c(0.75, 0.75, 0.75, 0.75), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #annotate("rect", xmin = 7, xmax = 11, ymin = 50, ymax = 100, fill = "red", alpha = 0.3) +
  scale_color_manual(values = park_cols, name = "Park or region") +
  scale_y_continuous(name = "Aug - Sept mean relative humidity (%)", 
                     expand = c(0,0), lim = c(0,100), breaks = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(name = "Aug - Sept mean temperature (°C)", 
                     expand = c(0,0), lim = c(0,30))

trans_year
