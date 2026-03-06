 
#================================================================#
#  Script created by JW Gaeta, IEP/CDFW on Feb 4, 2020           #
#  Last modified by Shruti Khanna in February 2025               #
#================================================================#

# rm(list=ls())
# graphics.off()
# cat("\f")

# ================================================================
# Attributes of the input file
# ================================================================
# FID         - site ID
# Year        - Year of observation
# CY          - continuous years of herbicide application upto Year in this polygon
# pixFAV      - number of pixels of FAV   at the site for that year
# pixSAV      - number of pixels of SAV   at the site for that year
# pixwater    - number of pixels of water at the site for that year
# percFAV     - fractional cover of FAV   at the site for that year
# percSAV     - fractional cover of SAV   at the site for that year
# percwater   - fractional cover of water at the site for that year; FAV + SAV + water = 1
# spr1cy      - polygon was sprayed 1 year  ago and not this year  (CY from last year)
# spr2cy      - polygon was sprayed 2 years ago and not since then (CY from 2 yrs ago)
# spr3cy      - polygon was sprayed 3 years ago and not since then (CY from 3 yrs ago)
# uniqueID    - FID_yyyy - unique ID for each row in the table
# SprTtlYrs   - total number of years since 2003 that a site was sprayed
# Orig_FID    - unique per polygon (multiple year observations per polygon)
# Habitat     - Channel; SlowShallows (between marsh remnants); DeadendSlough; FloodedIsland
# Stratum     - Region of the Delta (regional differences)
# index       - just an index - ignore
# SprayPolyAreaM2 - Area of polygon (site) calculated in square meters
# Latitude    - Latitude  for centroid of each site
# Longitude   - Longitude for centroid of each site
# RID         - DBW reach ID - matches DBW's numbers for their site IDs
# DBWsiteAreaM2 - Area of DBW reach in square meters
# yr_corr     - the correct year for the spray data (2004-2008 imagery is pre-spray so "yr-1")
# yr_rid      - yyyy_DBW reach ID
# ttlX24D     - amount of 2,4-D sprayed in that DBW reach
# ttlGlyph    - amount of glyphosate sprayed in that DBW reach
# cnt.rws     - number of times each DBW reach was sprayed in one year
# Year.y      - same at yr_corr but if polygon not sprayed then "NA"
# AggHabitat  - DeadEndSlough changed to SlowShallows to aggregate to 3 habitats
# AggStratum  - six strata aggregated to 4 regions
# non_FAV_pixel_count - total SAV + Water pixels in the polygon (not FAV)
# FAV_pixel_count - total FAV pixels in the polygon
# FAVareaM2   - calculated here for getting to rate of application
# ================================================================
# calculated rate of application of 2,4-D and glyphosate
# ttlFAVareaM2 - total area of all sprayed polygons within a DBW reach
# X24D_per_area_FAV   - 2,4-D      applied per 1 sq meter of FAV per DBW reach
# glypho_per_area_FAV - glyphosate applied per 1 sq meter of FAV per DBW reach
# =================================================================

#================
#~ load packages
#================

library('lme4')
library('arm')
library('gamm4')
library('ggeffects')
library('mgcv')
library('tidyverse')
library('effects')
library('readr')
library('performance')
library("performance")
library("blmeco")

################ These paths will change ##################

# workdir  = "X:/delta_sav/raster/analysis/ControlEfficacy/FAV_efficacy/mixed_effects/LastCopies/"
# setwd(workdir)
filename = paste0("figures/", "LongFormFAVdata_Agg_corrCY_Cnt_corrTrt_2511.csv")

###########################################################

# read in the data, create columns with the 2 response variables
dat = read_csv(filename) %>%
      mutate(non_FAV_pixel_count = pixSAV + pixWat) %>%
      mutate(FAV_pixel_count = pixFAV) %>% rename(Year = Year.x) %>%
      # calculate area in sq. meters for FAV only
      # 2004-2008 = 3x3m pixels; 2014-2017 = 2.2x2.5m pixels; 2018 = 1.7x1.7m pixels
      mutate(FAVareaM2 = ifelse(Year == 2018, pixFAV*1.7*1.7, ifelse(Year <= 2008, pixFAV*3*3, pixFAV*2.5*2.5))) %>%
      mutate(TimePeriod = ifelse(Year < 2009, "first", "second"))

# convert year to numeric
dat$Year = as.numeric(dat$Year)

numrows = nrow(dat)
names(dat)

# get the names for all columns
col.names = colnames(dat)
col.yr = match("Year",   col.names)
col.id = match("FID",    col.names)
col.CY = match("CY",     col.names)
col.sp = match("spr1cy", col.names)

# change polygon unique ID to a factor and assign to SID column
dat$SID = as.factor(dat$FID)

ttlx = table(dat$CY)

# for all 10 years of data: 2004-2008; 2014-2018
# 0     1    2    3    4    5    6    7     8     9    10   11   12   13  14   15   16   17   18
# 1228  469  310  367  279  211  175  124   10    9    7    4    2    3   21   19   13   14   12

# order by SID and year
dat = dat[order(dat$SID, dat$Year),]
dat.tab = dat %>% mutate(CY7 = ifelse(CY <= 6, CY, 7))

#=======================================
#~ Check response variable for normality
#=======================================
par(mfrow = c(1,1))
hist(dat$FAV_pixel_count, breaks = 30)
hist(dat$non_FAV_pixel_count, breaks = 30)

#=====================================================
#~ Hypothesis driven model
#-----------------------------------------------------
# non_FAV_pixel_count - non-FAV pixels in site
# FAV_pixel_count - FAV pixels in site (response)
# tensor product predictors:
# ctr-year - center year (-7:7)
# CY - consecutive years of treatment
# Spr1cy - CY if sprayed 1 year  ago | not sprayed since
# Spr2cy - CY if sprayed 2 years ago | not sprayed since
# Spr3cy - CY if sprayed 3 years ago | not sprayed since
# random effect:
# SiteID - multiple observations per site
#=====================================================

# constant conversion from gallons to liters
# doing this to make the range have fewer decimals
# so the unit of herbicide concentration is "milliliters sprayed per square meters"
GtoML = 3785.41

# only for sites that were sprayed, group by Reach ID and calculate total FAV area per Reach
calc.spr.rate = dat %>% filter(CY != 0) %>%
                group_by(yr_rid) %>%
                summarise(ttlFAVareaM2 = sum(FAVareaM2), ttlX24Dr = mean(ttlX24D), ttlGlyphr = mean(ttlGlyph)) %>%
                mutate(X24D_perFAVareaM2 = ifelse(ttlFAVareaM2 != 0,  ttlX24Dr*GtoML/ttlFAVareaM2, 0)) %>%
                mutate(glph_perFAVareaM2 = ifelse(ttlFAVareaM2 != 0, ttlGlyphr*GtoML/ttlFAVareaM2, 0))

# join this newly calculated spray rate to original data by DBW RID and Year
# keep in mind this will also join NON-SPRAYED sites to a non-zero value within that reach
dat.rate = left_join(dat, calc.spr.rate, by = "yr_rid") %>%
           mutate(X24D_perFAVareaM2 = ifelse(is.na(X24D_perFAVareaM2), 0, X24D_perFAVareaM2)) %>%
           mutate(glph_perFAVareaM2 = ifelse(is.na(glph_perFAVareaM2), 0, glph_perFAVareaM2))

# amount of 2,4-D sprayed per unit area of FAV
hist(dat.rate$X24D_perFAVareaM2, breaks = 80)
# amount of glyphosate sprayed per unit area of FAV
hist(dat.rate$glph_perFAVareaM2, breaks = 80)

#===================================================#
# Originally Jereme's section setting up the models #
#         Now greatly modified by Shruti            #
#     Bifurcates here for the two time periods      #
#===================================================#

# dataset dealing with the first time period
dat.first  = dat.rate[dat.rate$TimePeriod == "first",]
# dataset dealing with the second time period
dat.second = dat.rate[dat.rate$TimePeriod == "second",]

# currently doing first time period
dat.cur = dat.first

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# currently doing second time period
#dat.cur = dat.second

#    start here for analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# which rows had value of CY = 0?
zero_index = which(dat.cur$CY == 0)

# replace all spray rate values associated with control sites by 0's
dat.cur$X24D_perFAVareaM2[zero_index] = 0
dat.cur$glph_perFAVareaM2[zero_index] = 0
dat.cur$cnt.rws[zero_index] = 0

# Scale the CY distribution (mean = 0)
dat.scl = dat.cur %>%
          mutate(cbrt_sprfreq = cnt.rws^(1/3)) %>%
          # cube-root of z-scaled amt/unit area calculated to normalize distribution
          mutate(cbrt_X24D = sign(X24D_perFAVareaM2)*abs(X24D_perFAVareaM2)^(1/3)) %>%
          mutate(cbrt_glph = sign(glph_perFAVareaM2)*abs(glph_perFAVareaM2)^(1/3)) %>%
          # natural log of amt/unit area calculated to normalize distribution
          # add 1 because log(0) is undefined
          mutate(logn_X24D = log(1 + X24D_perFAVareaM2)) %>%
          mutate(logn_glph = log(1 + glph_perFAVareaM2)) %>% 
          # log to base 10 of amt/unit area calculated to normalize distribution
          # add 1 because log(0) is undefined
          mutate(log10_X24D = log10(1 + X24D_perFAVareaM2)) %>%
          mutate(log10_glph = log10(1 + glph_perFAVareaM2)) %>% 
          mutate(CY7 = ifelse(CY <= 6, CY, 7), spr1CY7 = ifelse(spr1cy <= 4, spr1cy, 5), spr2CY7 = ifelse(spr2cy <= 4, spr2cy, 5)) %>% 
          mutate(CY2yrs = ifelse(CY <= 6, (ceiling(CY/2)), 4)) %>% 
          mutate(CYscale = (CY7/4), spr1CYscl = (spr1CY7/5), spr2CYscl = (spr2CY7/5)) %>% 
          mutate(flgspr = ifelse(CY != 0, NA, ifelse(spr3cy != 0, 3, ifelse(spr2cy != 0, 2, ifelse(spr1cy != 0, 1, 0))))) %>% 
          mutate(spr1CY2yrs = ifelse(spr1CY7 <= 6, (ceiling(spr1CY7/2)), 4)) %>% 
          mutate(log10_pszsav = log10(1 + pszsav)) %>% 
          mutate(log10_pszfav = log10(1 + pszfav))

table(dat.scl$CY7)
table(dat.scl$spr1CY7)
table(dat.scl$spr2CY7)

# check histograms for normality
#-------------------------------------------
# temp = filter(dat.scl, cbrt_X24D != 0)
# hist(temp$cbrt_X24D,  breaks = 30)
# hist(temp$logn_X24D,  breaks = 30)
# hist(temp$log10_X24D, breaks = 30)
# temp = filter(dat.scl, cbrt_glph != 0)
# hist(temp$cbrt_glph,  breaks = 30)
# hist(temp$logn_glph,  breaks = 50)
# hist(temp$log10_glph, breaks = 50)
# hist(temp$cnt.rws, breaks = 30)
# unique(temp$cnt.rws)
# hist(temp$cbrt_sprfreq, breaks = 30)
# temp = filter(dat.scl, CY != 0)
# hist(temp$CY, breaks = 30)
# hist(temp$CYscale, breaks = 30)
# 
# hist(dat.scl$pszfav, breaks = 30)
# hist(log(dat.scl$pszfav), breaks = 30)
# hist(log10(dat.scl$pszfav), breaks = 30)
# 
# hist(dat.scl$pszsav, breaks = 30)
# hist(log(dat.scl$pszsav), breaks = 30)
# hist(log10(dat.scl$pszsav), breaks = 30)
#-------------------------------------------

# declare CY (continuous years of treatment) and Year as factor
dat.scl$fCY    = as.factor(dat.scl$CY)
dat.scl$fCY2yr = as.factor(dat.scl$CY2yrs)
dat.scl$fYear  = as.factor(dat.scl$Year)
dat.scl$fflgspr= as.factor(dat.scl$flgspr)
dat.scl$fspr1CY2yrs = as.factor(dat.scl$spr1CY2yrs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://stackoverflow.com/questions/33670628/solution-to-the-warning-message-using-glmer
# https://stackoverflow.com/questions/60564119/convergence-issues-glmer-how-to-interpret-allfit-outcomes-and-comparing-models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#=====================================================================================#
# test if year, consecutive years of treatment are related to fractional cover of FAV #
#     LOGISTIC MODEL WITH FAMILY = BINOMIAL; separate models for each habitat         #
#=====================================================================================#

# separate out the 3 habitats as 3 data frames
chan_dat  = subset(dat.scl, dat.scl$AggHabitat == "Channel")       %>% filter(glph_perFAVareaM2 == 0 | CY == 0)
shal_dat  = subset(dat.scl, dat.scl$AggHabitat == "SlowShallows")  %>% filter(glph_perFAVareaM2 == 0 | CY == 0)
flood_dat = subset(dat.scl, dat.scl$AggHabitat == "FloodedIsland") %>% filter(glph_perFAVareaM2 == 0 | CY == 0)

table( chan_dat$CY2yrs)
table( shal_dat$CY2yrs)
table(flood_dat$CY2yrs)
table( chan_dat$flgspr)
table( shal_dat$flgspr)
table(flood_dat$flgspr)

#==================================================================================#
# FINAL MODELS FOR ALL 3 HABITATS: additive or interactive between both herbicides #
#==================================================================================#

chan_mod  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~ log10_X24D + cbrt_sprfreq + fCY2yr +
                   (1|fYear) + (1|SID), data = chan_dat, family = binomial)
shal_mod  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~ log10_X24D + cbrt_sprfreq + fCY2yr +
                   (1|fYear) + (1|SID), data = shal_dat, family = binomial)
flood_mod = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~ log10_X24D + cbrt_sprfreq + fCY2yr +
                   (1|fYear) + (1|SID), data = flood_dat, family = binomial)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### formatting and writing the statistics figures to be reported for model
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cm.sum = summary(chan_mod)
sm.sum = summary(shal_mod)
fm.sum = summary(flood_mod)
cm.coeff = cm.sum$coefficients
sm.coeff = sm.sum$coefficients
fm.coeff = fm.sum$coefficients
am.coeff = data.frame(rbind(cm.coeff, sm.coeff, fm.coeff))
am.coeff$varname = row.names(am.coeff)
colnames(am.coeff) = c("Estimate", "StdError", "Zvalue", "Pvalue", "VarName")
am.coeff.fm <- am.coeff %>% mutate(
  # Columns 1, 2, 3 as decimal with 4 digits
  Estimate = sprintf("%.3f", Estimate),
  StdError = sprintf("%.3f", StdError),
  Zvalue   = sprintf("%.2f", Zvalue),
  Pvalue   = sprintf("%.4e", Pvalue)
  # ".2e" means 2 digits after the decimal in scientific notation (e.g., 9.88e+06)
) %>% mutate(EstSE = paste0(Estimate, " (", StdError, ")"))
rand.eff = rbind(unlist(cbind(cm.sum$devcomp$dims[1], cm.sum$ngrps, cm.sum$varcor)),
                 unlist(cbind(sm.sum$devcomp$dims[1], sm.sum$ngrps, sm.sum$varcor)),
                 unlist(cbind(fm.sum$devcomp$dims[1], fm.sum$ngrps, fm.sum$varcor)))
colnames(rand.eff) = c("nObs", "nDud", "nSID", "nYr", "seSID", "seYr")
df.reff = data.frame(rand.eff) %>% mutate(
  # Columns 1, 2, 3 as decimal with 4 digits
  nObs  = sprintf("%d", nObs),
  nDud  = sprintf("%d", nDud),
  nSID  = sprintf("%d", nSID),
  nYr   = sprintf("%d", nYr),
  seSID = sprintf("%.3f", seSID),
  seYr  = sprintf("%.3f", seYr))

write_csv(am.coeff.fm, paste0(workdir, "fixeff_allhab_glmm_0408_noGlyph.csv"))
write_csv(df.reff,     paste0(workdir, "rndeff_allhab_glmm_0408_noGlyph.csv"))

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get model residuals
 chan.rsdl = resid(chan_mod,  type = "pearson")
 shal.rsdl = resid(shal_mod,  type = "pearson")
flood.rsdl = resid(flood_mod, type = "pearson")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MODEL VALIDATION PLOTS FOR THE THREE FINAL HABITAT MODELS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CHANNEL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# validate the final model by plotting residuals and save to a tiff
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filename, resolution, size, type, etc.
tiff(filename = paste0("figures/", "channel_glmm_validation_0408_noGlph.tiff"), width = 10, height = 6, units = "in", res = 300)
# allows adjusting the first 2 plots to take the entire width
# makes a 2 row - 6 column layout
# 3 columns each to the first row (2 plots)
# 2 columns each to second row (3 plots)
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 2, 6, byrow = TRUE), respect = FALSE)
# global margins
par(mar = c(4, 4, 1, 1) + 0.1, oma = rep(0.5, 4), cex = 1.02)
#op = par(mfrow = c(2, 3), mar = c(4,4,1,1), oma = c(0.5,0.5,1,0), cex = 1.2)
# Plot 1: Residuals vs. Fitted values; should be centered around 0
plot(chan.rsdl ~ fitted(chan_mod), pch = 20, cex.lab = 1.15, las = 1,
     col = gray(0.3,0.3), xlab = "Model Residuals", ylab = "Predicted")
abline(h = 0, lwd = 2, col = gray(0.1, 0.4))
mtext("(a)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plot 2: histogram of the residuals; should be centered around 0
hist(chan.rsdl, las = 1, yaxs = 'i', cex.lab = 1.15, breaks = 30,
     xlab = 'Model Residuals', ylab = "Frequency", main = "")
mtext("(b)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plots 4,5,6: the Residuals vs. all the predictors; should be centered around 0
plot(chan_dat$log10_X24D, chan.rsdl, xlab = bquote(log[10]*"herbicide"~"(ml"%.%m^-2*")"), ylab = "Model residuals",
     pch = 20, cex.lab = 1.15, col = gray(0.3,0.3), las = 1)
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(c)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
plot(chan_dat$cbrt_sprfreq, chan.rsdl, xlab = bquote("sprFreq"^frac(1,3)), ylab = "",
     pch = 20, cex.lab = 1.15, col = gray(0.3,0.3), las = 1)
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(d)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
boxplot(chan.rsdl ~ chan_dat$CY2yrs, xlab = "CY2yrs", ylab = "",
     pch = 20, cex.lab = 1.15, col = gray(0.3,0.3), las = 1)
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(e)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# par(op)
dev.off()

# SLOW SHALLOWS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filename, resolution, size, type, etc.
tiff(filename = paste0("figures/", "slowsh_glmm_validation_0408_noGlph.tiff"), width = 10, height = 6, units = "in", res = 300)
# allows adjusting the first 2 plots to take the entire width
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 2, 6, byrow = TRUE), respect = FALSE)
# global margins
par(mar = c(4, 4, 1, 1) + 0.1, oma = rep(0.5, 4), cex = 1.02)
#op = par(mfrow = c(2, 3), mar = c(4,4,1,1), oma = c(0.5,0.5,1,0), cex = 1.2)
# Plot 1: Residuals vs. Fitted values; should be centered around 0
plot(shal.rsdl ~ fitted(shal_mod), pch = 20, cex.lab = 1.15, las = 1,
     col = gray(0.3,0.3), xlab = "Model Residuals", ylab = "Predicted")
abline(h = 0, lwd = 2, col = gray(0.1, 0.4))
mtext("(a)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plot 2: histogram of the residuals; should be centered around 0
hist(shal.rsdl, las = 1, yaxs = 'i', cex.lab = 1.15, breaks = 30,
     xlab = 'Model Residuals', ylab = "Frequency", main = "")
mtext("(b)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plots 4,5,6: the Residuals vs. all the predictors; should be centered around 0
plot(shal_dat$log10_X24D, shal.rsdl, xlab = bquote(log[10]*"herbicide"~"(ml"%.%m^-2*")"), ylab = "Model residuals",
     pch = 20, cex.lab = 1.15, col = gray(0.3,0.3), las = 1)
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(c)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
plot(shal_dat$cbrt_sprfreq, shal.rsdl, xlab = bquote("sprFreq"^frac(1,3)), ylab = "",
     pch = 20, cex.lab = 1.15, col = gray(0.3,0.3), las = 1)
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(d)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
boxplot(shal.rsdl ~ shal_dat$CY2yrs, xlab = "CY2yrs", ylab = "",
        pch = 20, cex.lab = 1.15, col = gray(0.3,0.3), las = 1)
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(e)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# par(op)
dev.off()

# FLOODED ISLANDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# filename, resolution, size, type, etc.
tiff(filename = paste0("figures/", "floodis_glmm_validation_0408_noGlph.tiff"), width = 10, height = 6, units = "in", res = 300)
# allows adjusting the first 2 plots to take the entire width
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 2, 6, byrow = TRUE), respect = FALSE)
# global margins
par(mar = c(4, 4, 1, 1) + 0.1, oma = rep(0.5, 4), cex = 1.02)
#op = par(mfrow = c(2, 3), mar = c(4,4,1,1), oma = c(0.5,0.5,1,0), cex = 1.2)
# Plot 1: Residuals vs. Fitted values; should be centered around 0
plot(flood.rsdl ~ fitted(flood_mod), pch = 20, cex.lab = 1.15, las = 1,
     col = gray(0.3,0.3), xlab = "Model Residuals", ylab = "Predicted")
abline(h = 0, lwd = 2, col = gray(0.1, 0.4))
mtext("(a)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plot 2: histogram of the residuals; should be centered around 0
hist(flood.rsdl, las = 1, yaxs = 'i', cex.lab = 1.15, breaks = 30,
     xlab = 'Model Residuals', ylab = "Frequency", main = "")
mtext("(b)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plots 4,5,6: the Residuals vs. all the predictors; should be centered around 0
plot(flood_dat$log10_X24D, flood.rsdl, xlab = bquote(log[10]*"herbicide"~"(ml"%.%m^-2*")"), ylab = "Model residuals",
     pch = 20, cex.lab = 1.15, col = gray(0.3,0.3), las = 1)
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(c)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
plot(flood_dat$cbrt_sprfreq, flood.rsdl, xlab = bquote("sprFreq"^frac(1,3)), ylab = "",
     pch = 20, cex.lab = 1.15, col = gray(0.3,0.3), las = 1)
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(d)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
boxplot(flood.rsdl ~ flood_dat$CY2yrs, xlab = "CY2yrs", ylab = "",
        pch = 20, cex.lab = 1.15, col = gray(0.3,0.3), las = 1)
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(e)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# par(op)
dev.off()

################################################################################

# limited number of CY values for each habitat
# get the unique values of CY and sort the list
cy.list = sort(unique(dat.scl$fCY2yr))
cy.num  = length(cy.list)
cy.chan = sort(unique(chan_dat$fCY2yr))
cy.numc = length(cy.chan)
cy.shal = sort(unique(shal_dat$fCY2yr))
cy.nums = length(cy.shal)
cy.flud = sort(unique(flood_dat$fCY2yr))
cy.numf = length(cy.flud)

# limited number of spray frequency values for each habitat
# get the unique values of sprfreq and sort the list
frq.list = sort(unique(dat.scl$cbrt_sprfreq))
frq.num  = length(frq.list)
frq.chan = sort(unique(chan_dat$cbrt_sprfreq))
frq.numc = length(frq.chan)
frq.shal = sort(unique(shal_dat$cbrt_sprfreq))
frq.nums = length(frq.shal)
frq.flud = sort(unique(flood_dat$cbrt_sprfreq))
frq.numf = length(frq.flud)

#___________________________________________________________________
#   1) Inputs: Build predictor gradient sequences
#___________________________________________________________________

# number of values to be simulated between the 5 to 95 %ile values of each covar
seq_length = 100

# create a sequence of 100 values from 5th to 95 percentile of covar values
log10_X24D_seq = seq(from = quantile(na.omit(dat.scl$log10_X24D), probs = c(0.05)),
                      to = quantile(na.omit(dat.scl$log10_X24D), probs = c(0.95)),
                      length.out = seq_length)

# for all three types of habitat
hab_class = c("Channel", "SlowShallows", "FloodedIsland")

# model for each habitat type
mod_cat = list(chan_mod, shal_mod, flood_mod)

# function definition for getting 5th and 95th percentile values
qprobs = function(x){
  quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)
}
# function definition for getting max, min, and median values
qprobs2 = function(x){
  quantile(x, probs = c(0, 0.5, 1), na.rm = TRUE)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   1) Simulate: Model CY = 0, herbicide = 0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# median, min, max of spray frequency per habitat
med.frq    =   dat.scl %>% filter(cbrt_sprfreq != 0) %>% summarize(across(cbrt_sprfreq, median))
med.frq.ch =  chan_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(across(cbrt_sprfreq, median))
med.frq.ss =  shal_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(across(cbrt_sprfreq, median))
med.frq.fi = flood_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(across(cbrt_sprfreq, median))
min.frq.ch =  chan_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(across(cbrt_sprfreq, min))
min.frq.ss =  shal_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(across(cbrt_sprfreq, min))
min.frq.fi = flood_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(across(cbrt_sprfreq, min))
max.frq.ch =  chan_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(across(cbrt_sprfreq, max))
max.frq.ss =  shal_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(across(cbrt_sprfreq, max))
max.frq.fi = flood_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(across(cbrt_sprfreq, max))
q40.frq.ch =  chan_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(quantile(cbrt_sprfreq, probs = 0.3))
q60.frq.ch =  chan_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(quantile(cbrt_sprfreq, probs = 0.7))
q40.frq.ss =  shal_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(quantile(cbrt_sprfreq, probs = 0.3))
q60.frq.ss =  shal_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(quantile(cbrt_sprfreq, probs = 0.7))
q40.frq.fi = flood_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(quantile(cbrt_sprfreq, probs = 0.3))
q60.frq.fi = flood_dat %>% filter(cbrt_sprfreq != 0) %>% summarize(quantile(cbrt_sprfreq, probs = 0.7))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   2) Data constrain: Identify limits of herbicide application per CY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# store 2,4-D values for 0%ile, 50%ile and 100%ile values for all 3 habitats for each CY (data type: list)
 chan_X24D_qprobs2 = aggregate( chan_dat$log10_X24D ~  chan_dat$CY2yrs, FUN = qprobs2)
 shal_X24D_qprobs2 = aggregate( shal_dat$log10_X24D ~  shal_dat$CY2yrs, FUN = qprobs2)
flood_X24D_qprobs2 = aggregate(flood_dat$log10_X24D ~ flood_dat$CY2yrs, FUN = qprobs2)

# store 2,4-D values for 5%ile and 95%ile values for all 3 habitats for each CY (data type: list)
 chan_X24D_qprobs = aggregate( chan_dat$log10_X24D ~  chan_dat$CY2yrs, FUN = qprobs)
 shal_X24D_qprobs = aggregate( shal_dat$log10_X24D ~  shal_dat$CY2yrs, FUN = qprobs)
flood_X24D_qprobs = aggregate(flood_dat$log10_X24D ~ flood_dat$CY2yrs, FUN = qprobs)

# store 2,4-D values for 0%ile, 50%ile and 100%ile values for all 3 habitats for each spary freq (data type: list)
 chan_X24D_qprobs2.sfr = aggregate( chan_dat$log10_X24D ~  chan_dat$cbrt_sprfreq, FUN = qprobs2)
 shal_X24D_qprobs2.sfr = aggregate( shal_dat$log10_X24D ~  shal_dat$cbrt_sprfreq, FUN = qprobs2)
flood_X24D_qprobs2.sfr = aggregate(flood_dat$log10_X24D ~ flood_dat$cbrt_sprfreq, FUN = qprobs2)

##########   CHANNEL HABITAT   ############

mod.cur = chan_mod

# for each habitat model and CY, get predictions of FAV pixels
eff = allEffects(mod = mod.cur,  
                 xlevels = list(log10_X24D = log10_X24D_seq),
                 fixed.predictors = list(given.values = c(fCY2yr1 = 0, fCY2yr2 = 0, fCY2yr3 = 0, fCY2yr4 = 0,
                                                          cbrt_sprfreq = frq.chan[1])))

# extract the fitted values, and upper and lower CF bounds
eff_fit = invlogit(eff$`log10_X24D`$fit)
eff_upp = invlogit(eff$`log10_X24D`$upper)
eff_low = invlogit(eff$`log10_X24D`$lower)

# matrix to accept the fit and conf. int. of this specific case
pfav.cy0.ch = array(data = NA, dim = c(3, 100),
                    dimnames = list(c("fit", "upper", "lower"), log10_X24D_seq))
pfav.cy0.ch[1,] = eff_fit
pfav.cy0.ch[2,] = eff_upp
pfav.cy0.ch[3,] = eff_low

# matrix to accept the fit and conf. int. of this specific case
pfav.cy1.ch = array(data = NA, dim = c(3, 100, frq.numc),
                    dimnames = list(c("fit", "upper", "lower"), log10_X24D_seq, frq.chan))

# min and max for CY = 1-2
mn = chan_X24D_qprobs[2,2][1]
mx = chan_X24D_qprobs[2,2][2]

for(i in 2:frq.numc){
  
  # for each habitat model and spray frequency, get predictions of FAV pixels
  eff = allEffects(mod = mod.cur,  
                   xlevels = list(log10_X24D = log10_X24D_seq),
                   fixed.predictors = list(given.values = c(fCY2yr1 = 1, fCY2yr2 = 0, fCY2yr3 = 0, fCY2yr4 = 0,
                                                            cbrt_sprfreq = frq.chan[i])))
  
  # extract the fitted values, and upper and lower CF bounds
  eff_fit = invlogit(eff$`log10_X24D`$fit)
  eff_upp = invlogit(eff$`log10_X24D`$upper)
  eff_low = invlogit(eff$`log10_X24D`$lower)
  
  pfav.cy1.ch[1, ,i] = eff_fit
  pfav.cy1.ch[2, ,i] = eff_upp
  pfav.cy1.ch[3, ,i] = eff_low
  
  # get all subscripts with log10(X24D0 values less than min for that spray frequency
  y = (which(log10_X24D_seq < max(mn, chan_X24D_qprobs2.sfr[i, 2][1])))
  # get all subscripts with log10(X24D0 values more than max for that spray frequency
  z = (which(log10_X24D_seq > min(mx, chan_X24D_qprobs2.sfr[i, 2][3])))
  # declare all data values outside of limits existing in dataset as NA's
  # i.e. only possible combinations of spray frequency and herbicide concentration are retained
  if(length(y) > 0) {
    pfav.cy1.ch[, 1:max(y), i] = NA
  }
  if(length(z) > 0){
    pfav.cy1.ch[, min(z):frq.numc, i] = NA
  }
}

# convert the 3D matrix to a 2D table
pfav.cy1.ch.df = as.data.frame(as.table(pfav.cy1.ch))
# name the columns correctly
colnames(pfav.cy1.ch.df) = c("cf_value", "log10_X24D", "SprFreq", "PrFAV")
# pivot wider to make fit, upper and lower into separate columns
pfav.cy1.ch.pw = pivot_wider(pfav.cy1.ch.df, names_from = "cf_value", values_from = "PrFAV")  %>% filter(!(is.na(fit)))

# matrix to accept the fit and conf. int. of this specific case
pfav.cy4.ch = array(data = NA, dim = c(3, 100, frq.numc),
                    dimnames = list(c("fit", "upper", "lower"), log10_X24D_seq, frq.chan))

# min and max for CY > 6
mn = chan_X24D_qprobs[5,2][1]
mx = chan_X24D_qprobs[5,2][2]

for(i in 2:frq.numc){
  
  # for each habitat model and spray frequency, get predictions of FAV pixels
  eff = allEffects(mod = mod.cur,  
                   xlevels = list(log10_X24D = log10_X24D_seq),
                   fixed.predictors = list(given.values = c(fCY2yr1 = 0, fCY2yr2 = 0, fCY2yr3 = 0, fCY2yr4 = 1,
                                                            cbrt_sprfreq = frq.chan[i])))
  
  # extract the fitted values, and upper and lower CF bounds
  eff_fit = invlogit(eff$`log10_X24D`$fit)
  eff_upp = invlogit(eff$`log10_X24D`$upper)
  eff_low = invlogit(eff$`log10_X24D`$lower)
  
  pfav.cy4.ch[1, ,i] = eff_fit
  pfav.cy4.ch[2, ,i] = eff_upp
  pfav.cy4.ch[3, ,i] = eff_low
  
  # get all subscripts with log10(X24D0 values less than min for that spray frequency
  y = (which(log10_X24D_seq < max(mn, chan_X24D_qprobs2.sfr[i, 2][1])))
  # get all subscripts with log10(X24D0 values more than max for that spray frequency
  z = (which(log10_X24D_seq > min(mx, chan_X24D_qprobs2.sfr[i, 2][3])))
  # declare all data values outside of limits existing in dataset as NA's
  # i.e. only possible combinations of spray frequency and herbicide concentration are retained
  if(length(y) > 0) {
    pfav.cy4.ch[, 1:max(y), i] = NA
  }
  if(length(z) > 0){
    pfav.cy4.ch[, min(z):frq.numc, i] = NA
  }
}

# convert the 3D matrix to a 2D table
pfav.cy4.ch.df = as.data.frame(as.table(pfav.cy4.ch))
# name the columns correctly
colnames(pfav.cy4.ch.df) = c("cf_value", "log10_X24D", "SprFreq", "PrFAV")
# pivot wider to make fit, upper and lower into separate columns
pfav.cy4.ch.pw = pivot_wider(pfav.cy4.ch.df, names_from = "cf_value", values_from = "PrFAV")  %>% filter(!(is.na(fit)))

########## SLOW SHALLOWS HABITAT ############

mod.cur = shal_mod

# for each habitat model and CY, get predictions of FAV pixels
eff = allEffects(mod = mod.cur,  
                 xlevels = list(log10_X24D = log10_X24D_seq),
                 fixed.predictors = list(given.values = c(fCY2yr1 = 0, fCY2yr2 = 0, fCY2yr3 = 0, fCY2yr4 = 0,
                                                          cbrt_sprfreq = 0)))

# extract the fitted values, and upper and lower CF bounds
eff_fit = invlogit(eff$`log10_X24D`$fit)
eff_upp = invlogit(eff$`log10_X24D`$upper)
eff_low = invlogit(eff$`log10_X24D`$lower)

# matrix to accept the fit and conf. int. of this specific case
pfav.cy0.ss = array(data = NA, dim = c(3, 100),
                    dimnames = list(c("fit", "upper", "lower"), log10_X24D_seq))
pfav.cy0.ss[1,] = eff_fit
pfav.cy0.ss[2,] = eff_upp
pfav.cy0.ss[3,] = eff_low

# matrix to accept the fit and conf. int. of this specific case
pfav.cy1.ss = array(data = NA, dim = c(3, 100, frq.nums),
                    dimnames = list(c("fit", "upper", "lower"), log10_X24D_seq, frq.shal))

# min and max for CY = 1-2
mn = shal_X24D_qprobs[2,2][1]
mx = shal_X24D_qprobs[2,2][2]

for(i in 2:frq.nums){
  
  # for each habitat model and spray frequency, get predictions of FAV pixels
  eff = allEffects(mod = mod.cur,  
                   xlevels = list(log10_X24D = log10_X24D_seq),
                   fixed.predictors = list(given.values = c(fCY2yr1 = 1, fCY2yr2 = 0, fCY2yr3 = 0, fCY2yr4 = 0,
                                                            cbrt_sprfreq = frq.shal[i])))
  
  # extract the fitted values, and upper and lower CF bounds
  eff_fit = invlogit(eff$`log10_X24D`$fit)
  eff_upp = invlogit(eff$`log10_X24D`$upper)
  eff_low = invlogit(eff$`log10_X24D`$lower)
  
  pfav.cy1.ss[1, ,i] = eff_fit
  pfav.cy1.ss[2, ,i] = eff_upp
  pfav.cy1.ss[3, ,i] = eff_low
  
  # get all subscripts with log10(X24D0 values less than min for that spray frequency
  y = (which(log10_X24D_seq < max(mn, shal_X24D_qprobs2.sfr[i, 2][1])))
  # get all subscripts with log10(X24D0 values more than max for that spray frequency
  z = (which(log10_X24D_seq > min(mx, shal_X24D_qprobs2.sfr[i, 2][3])))
  # declare all data values outside of limits existing in dataset as NA's
  # i.e. only possible combinations of spray frequency and herbicide concentration are retained
  if(length(y) > 0) {
    pfav.cy1.ss[, 1:max(y), i] = NA
  }
  if(length(z) > 0){
    pfav.cy1.ss[, min(z):frq.nums, i] = NA
  }
}

# convert the 3D matrix to a 2D table
pfav.cy1.ss.df = as.data.frame(as.table(pfav.cy1.ss))
# name the columns correctly
colnames(pfav.cy1.ss.df) = c("cf_value", "log10_X24D", "SprFreq", "PrFAV")
# pivot wider to make fit, upper and lower into separate columns
pfav.cy1.ss.pw = pivot_wider(pfav.cy1.ss.df, names_from = "cf_value", values_from = "PrFAV")  %>% filter(!(is.na(fit)))

# matrix to accept the fit and conf. int. of this specific case
pfav.cy4.ss = array(data = NA, dim = c(3, 100, frq.nums),
                    dimnames = list(c("fit", "upper", "lower"), log10_X24D_seq, frq.shal))

# min and max for CY = 1-2
mn = shal_X24D_qprobs[5,2][1]
mx = shal_X24D_qprobs[5,2][2]

for(i in 2:frq.nums){
  
  # for each habitat model and spray frequency, get predictions of FAV pixels
  eff = allEffects(mod = mod.cur,  
                   xlevels = list(log10_X24D = log10_X24D_seq),
                   fixed.predictors = list(given.values = c(fCY2yr1 = 0, fCY2yr2 = 0, fCY2yr3 = 0, fCY2yr4 = 1,
                                                            cbrt_sprfreq = frq.shal[i])))
  
  # extract the fitted values, and upper and lower CF bounds
  eff_fit = invlogit(eff$`log10_X24D`$fit)
  eff_upp = invlogit(eff$`log10_X24D`$upper)
  eff_low = invlogit(eff$`log10_X24D`$lower)
  
  pfav.cy4.ss[1, ,i] = eff_fit
  pfav.cy4.ss[2, ,i] = eff_upp
  pfav.cy4.ss[3, ,i] = eff_low
  
  # get all subscripts with log10(X24D0 values less than min for that spray frequency
  y = (which(log10_X24D_seq < max(mn, shal_X24D_qprobs2.sfr[i, 2][1])))
  # get all subscripts with log10(X24D0 values more than max for that spray frequency
  z = (which(log10_X24D_seq > min(mx, shal_X24D_qprobs2.sfr[i, 2][3])))
  # declare all data values outside of limits existing in dataset as NA's
  # i.e. only possible combinations of spray frequency and herbicide concentration are retained
  if(length(y) > 0) {
    pfav.cy4.ss[, 1:max(y), i] = NA
  }
  if(length(z) > 0){
    pfav.cy4.ss[, min(z):frq.nums, i] = NA
  }
}

# convert the 3D matrix to a 2D table
pfav.cy4.ss.df = as.data.frame(as.table(pfav.cy4.ss))
# name the columns correctly
colnames(pfav.cy4.ss.df) = c("cf_value", "log10_X24D", "SprFreq", "PrFAV")
# pivot wider to make fit, upper and lower into separate columns
pfav.cy4.ss.pw = pivot_wider(pfav.cy4.ss.df, names_from = "cf_value", values_from = "PrFAV")  %>% filter(!(is.na(fit)))

########## FLOODED ISLAND HABITAT ############

mod.cur = flood_mod

# for each habitat model and CY, get predictions of FAV pixels
eff = allEffects(mod = mod.cur,  
                 xlevels = list(log10_X24D = log10_X24D_seq),
                 fixed.predictors = list(given.values = c(fCY2yr1 = 0, fCY2yr2 = 0, fCY2yr3 = 0, fCY2yr4 = 0,
                                                          cbrt_sprfreq = 0)))

# extract the fitted values, and upper and lower CF bounds
eff_fit = invlogit(eff$`log10_X24D`$fit)
eff_upp = invlogit(eff$`log10_X24D`$upper)
eff_low = invlogit(eff$`log10_X24D`$lower)

# matrix to accept the fit and conf. int. of this specific case
pfav.cy0.fi = array(data = NA, dim = c(3, 100),
                    dimnames = list(c("fit", "upper", "lower"), log10_X24D_seq))
pfav.cy0.fi[1,] = eff_fit
pfav.cy0.fi[2,] = eff_upp
pfav.cy0.fi[3,] = eff_low

# matrix to accept the fit and conf. int. of this specific case
pfav.cy1.fi = array(data = NA, dim = c(3, 100, frq.numf),
                    dimnames = list(c("fit", "upper", "lower"), log10_X24D_seq, frq.flud))

# min and max for CY = 1-2
mn = flood_X24D_qprobs[2,2][1]
mx = flood_X24D_qprobs[2,2][2]

for(i in 2:frq.numf){
  
  # for each habitat model and spray frequency, get predictions of FAV pixels
  eff = allEffects(mod = mod.cur,  
                   xlevels = list(log10_X24D = log10_X24D_seq),
                   fixed.predictors = list(given.values = c(fCY2yr1 = 1, fCY2yr2 = 0, fCY2yr3 = 0, fCY2yr4 = 0,
                                                            cbrt_sprfreq = frq.flud[i])))
  
  # extract the fitted values, and upper and lower CF bounds
  eff_fit = invlogit(eff$`log10_X24D`$fit)
  eff_upp = invlogit(eff$`log10_X24D`$upper)
  eff_low = invlogit(eff$`log10_X24D`$lower)
  
  pfav.cy1.fi[1, ,i] = eff_fit
  pfav.cy1.fi[2, ,i] = eff_upp
  pfav.cy1.fi[3, ,i] = eff_low
  
  # get all subscripts with log10(X24D0 values less than min for that spray frequency
  y = (which(log10_X24D_seq < max(mn, flood_X24D_qprobs2.sfr[i, 2][1])))
  # get all subscripts with log10(X24D0 values more than max for that spray frequency
  z = (which(log10_X24D_seq > min(mx, flood_X24D_qprobs2.sfr[i, 2][3])))
  # declare all data values outside of limits existing in dataset as NA's
  # i.e. only possible combinations of spray frequency and herbicide concentration are retained
  if(length(y) > 0) {
    pfav.cy1.fi[, 1:max(y), i] = NA
  }
  if(length(z) > 0){
    pfav.cy1.fi[, min(z):frq.numf, i] = NA
  }
}

# convert the 3D matrix to a 2D table
pfav.cy1.fi.df = as.data.frame(as.table(pfav.cy1.fi))
# name the columns correctly
colnames(pfav.cy1.fi.df) = c("cf_value", "log10_X24D", "SprFreq", "PrFAV")
# pivot wider to make fit, upper and lower into separate columns
pfav.cy1.fi.pw = pivot_wider(pfav.cy1.fi.df, names_from = "cf_value", values_from = "PrFAV")  %>% filter(!(is.na(fit)))

# matrix to accept the fit and conf. int. of this specific case
pfav.cy4.fi = array(data = NA, dim = c(3, 100, frq.numf),
                    dimnames = list(c("fit", "upper", "lower"), log10_X24D_seq, frq.flud))

# min and max for CY = 1-2
mn = flood_X24D_qprobs[5,2][1]
mx = flood_X24D_qprobs[5,2][2]

for(i in 2:frq.numf){

  # for each habitat model and spray frequency, get predictions of FAV pixels
  eff = allEffects(mod = mod.cur,
                   xlevels = list(log10_X24D = log10_X24D_seq),
                   fixed.predictors = list(given.values = c(fCY2yr1 = 0, fCY2yr2 = 0, fCY2yr3 = 0, fCY2yr4 = 1,
                                                            cbrt_sprfreq = frq.flud[i])))

  # extract the fitted values, and upper and lower CF bounds
  eff_fit = invlogit(eff$`log10_X24D`$fit)
  eff_upp = invlogit(eff$`log10_X24D`$upper)
  eff_low = invlogit(eff$`log10_X24D`$lower)

  pfav.cy4.fi[1, ,i] = eff_fit
  pfav.cy4.fi[2, ,i] = eff_upp
  pfav.cy4.fi[3, ,i] = eff_low

  # get all subscripts with log10(X24D0 values less than min for that spray frequency
  y = (which(log10_X24D_seq < max(mn, flood_X24D_qprobs2.sfr[i, 2][1])))
  # get all subscripts with log10(X24D0 values more than max for that spray frequency
  z = (which(log10_X24D_seq > min(mx, flood_X24D_qprobs2.sfr[i, 2][3])))
  # declare all data values outside of limits existing in dataset as NA's
  # i.e. only possible combinations of spray frequency and herbicide concentration are retained
  if(length(y) > 0) {
    pfav.cy4.fi[, 1:max(y), i] = NA
  }
  if(length(z) > 0){
    pfav.cy4.fi[, min(z):frq.numf, i] = NA
  }
}

# convert the 3D matrix to a 2D table
pfav.cy4.fi.df = as.data.frame(as.table(pfav.cy4.fi))
# name the columns correctly
colnames(pfav.cy4.fi.df) = c("cf_value", "log10_X24D", "SprFreq", "PrFAV")
# pivot wider to make fit, upper and lower into separate columns
pfav.cy4.fi.pw = pivot_wider(pfav.cy4.fi.df, names_from = "cf_value", values_from = "PrFAV")  %>% filter(!(is.na(fit)))
 
# when converting a number from factor to numeric, first convert to character
# otherwise, only the factor value will get converted to numeric
pfav.cy1.ch.pw$log10_X24D = as.numeric(as.character(pfav.cy1.ch.pw$log10_X24D))
pfav.cy4.ch.pw$log10_X24D = as.numeric(as.character(pfav.cy4.ch.pw$log10_X24D))
pfav.cy1.ss.pw$log10_X24D = as.numeric(as.character(pfav.cy1.ss.pw$log10_X24D))
pfav.cy4.ss.pw$log10_X24D = as.numeric(as.character(pfav.cy4.ss.pw$log10_X24D))
pfav.cy1.fi.pw$log10_X24D = as.numeric(as.character(pfav.cy1.fi.pw$log10_X24D))
pfav.cy4.fi.pw$log10_X24D = as.numeric(as.character(pfav.cy4.fi.pw$log10_X24D))

pfav.cy1.ch.pw$SprFreq = as.numeric(as.character(pfav.cy1.ch.pw$SprFreq))
pfav.cy4.ch.pw$SprFreq = as.numeric(as.character(pfav.cy4.ch.pw$SprFreq))
pfav.cy1.ss.pw$SprFreq = as.numeric(as.character(pfav.cy1.ss.pw$SprFreq))
pfav.cy4.ss.pw$SprFreq = as.numeric(as.character(pfav.cy4.ss.pw$SprFreq))
pfav.cy1.fi.pw$SprFreq = as.numeric(as.character(pfav.cy1.fi.pw$SprFreq))
pfav.cy4.fi.pw$SprFreq = as.numeric(as.character(pfav.cy4.fi.pw$SprFreq))

# DIVIDE PLOTS INTO 70-100%ILE & 0-30%ILE OF SPRAY FREQUENCY

pfav.cy1.ch.pw.low = pfav.cy1.ch.pw %>% filter(SprFreq <= q40.frq.ch[[1]])
pfav.cy1.ch.pw.hgh = pfav.cy1.ch.pw %>% filter(SprFreq >= q60.frq.ch[[1]])
pfav.cy1.ss.pw.low = pfav.cy1.ss.pw %>% filter(SprFreq <= q40.frq.ss[[1]])
pfav.cy1.ss.pw.hgh = pfav.cy1.ss.pw %>% filter(SprFreq >= q60.frq.ss[[1]])
pfav.cy1.fi.pw.low = pfav.cy1.fi.pw %>% filter(SprFreq <= q40.frq.fi[[1]])
pfav.cy1.fi.pw.hgh = pfav.cy1.fi.pw %>% filter(SprFreq >= q60.frq.fi[[1]])
pfav.cy4.ch.pw.low = pfav.cy4.ch.pw %>% filter(SprFreq <= q40.frq.ch[[1]])
pfav.cy4.ch.pw.hgh = pfav.cy4.ch.pw %>% filter(SprFreq >= q60.frq.ch[[1]])
pfav.cy4.ss.pw.low = pfav.cy4.ss.pw %>% filter(SprFreq <= q40.frq.ss[[1]])
pfav.cy4.ss.pw.hgh = pfav.cy4.ss.pw %>% filter(SprFreq >= q60.frq.ss[[1]])
pfav.cy4.fi.pw.low = pfav.cy4.fi.pw %>% filter(SprFreq <= q40.frq.fi[[1]])
pfav.cy4.fi.pw.hgh = pfav.cy4.fi.pw %>% filter(SprFreq >= q60.frq.fi[[1]])

#_______________________________________________________________________________
#     CREATE FIGURE FOR VARYING AMOUNT OF HERBICIDE, 2 CY & SPRAY FREQ
tiff(filename = "figures/FAV_X24D_2cy_sfqTile_0408_noGlph.tiff", width = 8, height = 6, units = "in", res = 300)
#_______________________________________________________________________________

# set outer (oma(bottom, left, top, right)) and inner margins (mar(bottom, left, top, right))
# in terms of lines of text; + 0.1 is a slight offset
par(mfcol = c(2,3), mar = c(1.5,1.5,1,1) + 0.1, oma = c(3.5,4.5,3,0.5))

# Define a gradient of colors
colors <- colorRampPalette(c("blue", "red3"))(25)  # 25 gradient colors

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    CHANNEL HABITAT; 2,4-D; 2004-2008; MEDIAN SPRAY FREQ; ALL CY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot just for the first value of X24D which is 0
plot.cy0 = c(pfav.cy0.ch[1, 1], pfav.cy0.ch[2, 1], pfav.cy0.ch[3, 1])
# plot a point for fitted value for control sites
plot(x = 0, y = plot.cy0[1], pch = 16, ylim = c(0,0.3),
     las = 0.5, ylab = "", xlab = "", main = "", xlim = range(log10_X24D_seq))
mtext("(a)", side = 3, line = -2, adj = 0.95, cex = 1.2)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# draw the perpendicular lines
arrows(x0 = 0, y0 = plot.cy0[2], x1 = 0, y1 = plot.cy0[3], angle = 90, code = 3, length = 0.1, col = "black")
# draw each fitted point with color changing according to spray frequency
points(pfav.cy1.ch.pw.hgh$log10_X24D, pfav.cy1.ch.pw.hgh$fit, pch = 16, cex = 0.5, col = rgb(1, 0.2, 0))
points(pfav.cy1.ch.pw.low$log10_X24D, pfav.cy1.ch.pw.low$fit, pch = 16, cex = 0.5, col = rgb(0, 0.2, 1))
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.cy1.ch.pw.hgh)) {
  segments(pfav.cy1.ch.pw.hgh$log10_X24D[i], pfav.cy1.ch.pw.hgh$upper[i], pfav.cy1.ch.pw.hgh$log10_X24D[i], pfav.cy1.ch.pw.hgh$lower[i], 
           col = rgb(1,0.85,0.7, alpha = 0.4), lty = 1, lwd = 2)
}
for (i in 1:nrow(pfav.cy1.ch.pw.low)) {
  segments(pfav.cy1.ch.pw.low$log10_X24D[i], pfav.cy1.ch.pw.low$upper[i], pfav.cy1.ch.pw.low$log10_X24D[i], pfav.cy1.ch.pw.low$lower[i], 
           col = rgb(0.7,0.85,1, alpha = 0.4), lty = 1, lwd = 2)
}

#_____________________ LABEL ALL AXES AND ALL PLOTS

mtext(text = "Pr(FAV occurrence)", side = 2, line = 1, outer = TRUE, adj = 0.2)
mtext(text = "Pr(FAV occurrence)", side = 2, line = 1, outer = TRUE, adj = 0.8)

mtext(text = "Treated 1-2 years", side = 2, line = 2.5, outer = TRUE, font = 2, adj = 0.8)
mtext(text = "Treated > 6 years", side = 2, line = 2.5, outer = TRUE, font = 2, adj = 0.2)

mtext(text = bquote(log[10]*"herbicide"~"(ml"%.%m^-2*")"), side = 1, line = 1.5, outer = TRUE, adj = 0.05)
mtext(text = bquote(log[10]*"herbicide"~"(ml"%.%m^-2*")"), side = 1, line = 1.5, outer = TRUE, adj = 0.5)
mtext(text = bquote(log[10]*"herbicide"~"(ml"%.%m^-2*")"), side = 1, line = 1.5, outer = TRUE, adj = 0.95)

mtext(text = "Channels",        side = 3, line = 0.5, outer = TRUE, cex = 1.2, font = 2, adj = 0.1)
mtext(text = "Slow Shallows",   side = 3, line = 0.5, outer = TRUE, cex = 1.2, font = 2, adj = 0.5)
mtext(text = "Flooded Islands", side = 3, line = 0.5, outer = TRUE, cex = 1.2, font = 2, adj = 0.9)

#_______________________________________________________________________________

# plot a point for fitted value for control sites
plot(x = 0, y = plot.cy0[1], pch = 16, ylim = c(0,0.3),
     las = 0.5, ylab = "", xlab = "", main = "", xlim = range(log10_X24D_seq))
mtext("(d)", side = 3, line = -2, adj = 0.95, cex = 1.2)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# draw the perpendicular lines
arrows(x0 = 0, y0 = plot.cy0[2], x1 = 0, y1 = plot.cy0[3], angle = 90, code = 3, length = 0.1, col = "black")
# draw each fitted point with color changing according to spray frequency
points(pfav.cy4.ch.pw.hgh$log10_X24D, pfav.cy4.ch.pw.hgh$fit, pch = 16, cex = 0.5, col = rgb(1, 0.2, 0))
points(pfav.cy4.ch.pw.low$log10_X24D, pfav.cy4.ch.pw.low$fit, pch = 16, cex = 0.5, col = rgb(0, 0.2, 1))
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.cy4.ch.pw.hgh)) {
  segments(pfav.cy4.ch.pw.hgh$log10_X24D[i], pfav.cy4.ch.pw.hgh$upper[i], pfav.cy4.ch.pw.hgh$log10_X24D[i], pfav.cy4.ch.pw.hgh$lower[i], 
           col = rgb(1,0.85,0.7, alpha = 0.4), lty = 1, lwd = 2)
}
for (i in 1:nrow(pfav.cy4.ch.pw.low)) {
  segments(pfav.cy4.ch.pw.low$log10_X24D[i], pfav.cy4.ch.pw.low$upper[i], pfav.cy4.ch.pw.low$log10_X24D[i], pfav.cy4.ch.pw.low$lower[i], 
           col = rgb(0.7,0.85,1, alpha = 0.4), lty = 1, lwd = 2)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   SLOW SHALLOWS HABITAT; GLYPHOSATE; 2014-2018; MEDIAN SPRAY FREQ; ALL CY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot just for the first value of X24D which is 0 - so only for glyphosate
plot.cy0 = c(pfav.cy0.ss[1, 1], pfav.cy0.ss[2, 1], pfav.cy0.ss[3, 1])

#_______________________________________________________________________________

# plot a point for fitted value for control sites
plot(x = 0, y = plot.cy0[1], pch = 16, ylim = c(0,0.3),
     las = 0.5, ylab = "", xlab = "", main = "", xlim = range(log10_X24D_seq))
mtext("(b)", side = 3, line = -2, adj = 0.95, cex = 1.2)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# draw the perpendicular lines
arrows(x0 = 0, y0 = plot.cy0[2], x1 = 0, y1 = plot.cy0[3], angle = 90, code = 3, length = 0.1, col = "black")
# draw each fitted point with color changing according to spray frequency
points(pfav.cy1.ss.pw.hgh$log10_X24D, pfav.cy1.ss.pw.hgh$fit, pch = 16, cex = 0.5, col = rgb(1, 0.2, 0))
points(pfav.cy1.ss.pw.low$log10_X24D, pfav.cy1.ss.pw.low$fit, pch = 16, cex = 0.5, col = rgb(0, 0.2, 1))
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.cy1.ss.pw.hgh)) {
  segments(pfav.cy1.ss.pw.hgh$log10_X24D[i], pfav.cy1.ss.pw.hgh$upper[i], pfav.cy1.ss.pw.hgh$log10_X24D[i], pfav.cy1.ss.pw.hgh$lower[i], 
           col = rgb(1,0.85,0.7, alpha = 0.4), lty = 1, lwd = 2)
}
for (i in 1:nrow(pfav.cy1.ss.pw.low)) {
  segments(pfav.cy1.ss.pw.low$log10_X24D[i], pfav.cy1.ss.pw.low$upper[i], pfav.cy1.ss.pw.low$log10_X24D[i], pfav.cy1.ss.pw.low$lower[i], 
           col = rgb(0.7,0.85,1, alpha = 0.4), lty = 1, lwd = 2)
}

#_______________________________________________________________________________

# plot a point for fitted value for control sites
plot(x = 0, y = plot.cy0[1], pch = 16, ylim = c(0,0.3),
     las = 0.5, ylab = "", xlab = "", main = "", xlim = range(log10_X24D_seq))
mtext("(e)", side = 3, line = -2, adj = 0.95, cex = 1.2)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# draw the perpendicular lines
arrows(x0 = 0, y0 = plot.cy0[2], x1 = 0, y1 = plot.cy0[3], angle = 90, code = 3, length = 0.1, col = "black")
# draw each fitted point with color changing according to spray frequency
points(pfav.cy4.ss.pw.hgh$log10_X24D, pfav.cy4.ss.pw.hgh$fit, pch = 16, cex = 0.5, col = rgb(1, 0.2, 0))
points(pfav.cy4.ss.pw.low$log10_X24D, pfav.cy4.ss.pw.low$fit, pch = 16, cex = 0.5, col = rgb(0, 0.2, 1))
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.cy4.ss.pw.hgh)) {
  segments(pfav.cy4.ss.pw.hgh$log10_X24D[i], pfav.cy4.ss.pw.hgh$upper[i], pfav.cy4.ss.pw.hgh$log10_X24D[i], pfav.cy4.ss.pw.hgh$lower[i], 
           col = rgb(1,0.85,0.7, alpha = 0.4), lty = 1, lwd = 2)
}
for (i in 1:nrow(pfav.cy4.ss.pw.low)) {
  segments(pfav.cy4.ss.pw.low$log10_X24D[i], pfav.cy4.ss.pw.low$upper[i], pfav.cy4.ss.pw.low$log10_X24D[i], pfav.cy4.ss.pw.low$lower[i], 
           col = rgb(0.7,0.85,1, alpha = 0.4), lty = 1, lwd = 2)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   FLOODED ISLAND HABITAT; GLYPHOSATE; 2014-2018; ALL SPRAY FREQ; ALL CY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot just for the first value of X24D which is 0 - so only for glyphosate
plot.cy0 = c(pfav.cy0.fi[1, 1], pfav.cy0.fi[2, 1], pfav.cy0.fi[3, 1])
# plot a point for fitted value for control sites
plot(x = 0, y = plot.cy0[1], pch = 16, ylim = c(0,0.3),
     las = 0.5, ylab = "", xlab = "", main = "", xlim = range(log10_X24D_seq))
mtext("(c)", side = 3, line = -2, adj = 0.95, cex = 1.2)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# draw the perpendicular lines
arrows(x0 = 0, y0 = plot.cy0[2], x1 = 0, y1 = plot.cy0[3], angle = 90, code = 3, length = 0.1, col = "black")
# draw each fitted point with color changing according to spray frequency
points(pfav.cy1.fi.pw.hgh$log10_X24D, pfav.cy1.fi.pw.hgh$fit, pch = 16, cex = 0.5, col = rgb(1, 0.2, 0))
points(pfav.cy1.fi.pw.low$log10_X24D, pfav.cy1.fi.pw.low$fit, pch = 16, cex = 0.5, col = rgb(0, 0.2, 1))
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.cy1.fi.pw.hgh)) {
  segments(pfav.cy1.fi.pw.hgh$log10_X24D[i], pfav.cy1.fi.pw.hgh$upper[i], pfav.cy1.fi.pw.hgh$log10_X24D[i], pfav.cy1.fi.pw.hgh$lower[i], 
           col = rgb(1,0.85,0.7, alpha = 0.4), lty = 1, lwd = 2)
}
for (i in 1:nrow(pfav.cy1.fi.pw.low)) {
  segments(pfav.cy1.fi.pw.low$log10_X24D[i], pfav.cy1.fi.pw.low$upper[i], pfav.cy1.fi.pw.low$log10_X24D[i], pfav.cy1.fi.pw.low$lower[i], 
           col = rgb(0.7,0.85,1, alpha = 0.4), lty = 1, lwd = 2)
}

#_______________________________________________________________________________
# 
# plot a point for fitted value for control sites
plot(x = 0, y = plot.cy0[1], pch = 16, ylim = c(0,0.3),
     las = 0.5, ylab = "", xlab = "", main = "", xlim = range(log10_X24D_seq))
mtext("(f)", side = 3, line = -2, adj = 0.95, cex = 1.2)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# draw the perpendicular lines
arrows(x0 = 0, y0 = plot.cy0[2], x1 = 0, y1 = plot.cy0[3], angle = 90, code = 3, length = 0.1, col = "black")
# draw each fitted point with color changing according to spray frequency
points(pfav.cy4.fi.pw.hgh$log10_X24D, pfav.cy4.fi.pw.hgh$fit, pch = 16, cex = 0.5, col = rgb(1, 0.2, 0))
points(pfav.cy4.fi.pw.low$log10_X24D, pfav.cy4.fi.pw.low$fit, pch = 16, cex = 0.5, col = rgb(0, 0.2, 1))
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.cy4.fi.pw.hgh)) {
  segments(pfav.cy4.fi.pw.hgh$log10_X24D[i], pfav.cy4.fi.pw.hgh$upper[i], pfav.cy4.fi.pw.hgh$log10_X24D[i], pfav.cy4.fi.pw.hgh$lower[i],
           col = rgb(1,0.85,0.7, alpha = 0.4), lty = 1, lwd = 2)
}
for (i in 1:nrow(pfav.cy4.fi.pw.low)) {
  segments(pfav.cy4.fi.pw.low$log10_X24D[i], pfav.cy4.fi.pw.low$upper[i], pfav.cy4.fi.pw.low$log10_X24D[i], pfav.cy4.fi.pw.low$lower[i],
           col = rgb(0.7,0.85,1, alpha = 0.4), lty = 1, lwd = 2)
}

#_______________________________________________________________________________

dev.off()

################################################################################
#####                  LEGACY EFFECT MODELS AND FIGURES                    #####
################################################################################

# only select sample reference polygons and polygons treated 1 year ago
ch.dat.lgcy1 =  chan_dat %>% filter(flgspr == 0 | flgspr == 1)
ss.dat.lgcy1 =  shal_dat %>% filter(flgspr == 0 | flgspr == 1)
fi.dat.lgcy1 = flood_dat %>% filter(flgspr == 0 | flgspr == 1)

# only select sample reference polygons and polygons treated 2 years ago
ch.dat.lgcy2 =  chan_dat %>% filter(flgspr == 0 | flgspr == 2)
ss.dat.lgcy2 =  shal_dat %>% filter(flgspr == 0 | flgspr == 2)
fi.dat.lgcy2 = flood_dat %>% filter(flgspr == 0 | flgspr == 2)

# only select sample reference polygons and polygons treated 2 years ago
ch.dat.lgcy3 =  chan_dat %>% filter(flgspr == 0 | flgspr == 3)
ss.dat.lgcy3 =  shal_dat %>% filter(flgspr == 0 | flgspr == 3)
fi.dat.lgcy3 = flood_dat %>% filter(flgspr == 0 | flgspr == 3)

# only select sample reference polygons to test for legacy effects
ch.dat.lgcy =  chan_dat %>% filter(!is.na(flgspr))
ss.dat.lgcy =  shal_dat %>% filter(!is.na(flgspr))
fi.dat.lgcy = flood_dat %>% filter(!is.na(flgspr))

# combine all 3 habitats
allhab.lgcy1 = rbind(ch.dat.lgcy1, ss.dat.lgcy1, fi.dat.lgcy1)
allhab.lgcy2 = rbind(ch.dat.lgcy2, ss.dat.lgcy2, fi.dat.lgcy2)
allhab.lgcy3 = rbind(ch.dat.lgcy3, ss.dat.lgcy3, fi.dat.lgcy3)
allhab.lgcy  = rbind(ch.dat.lgcy,  ss.dat.lgcy,  fi.dat.lgcy)

# models for each habitat, 1 and 2 year legacy, and for all habitats together
chan.mod.spr1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        + fspr1CY2yrs + (1|fYear) + (1|SID), data = ch.dat.lgcy1, family = binomial)
shal.mod.spr1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        + fspr1CY2yrs + (1|fYear) + (1|SID), data = ss.dat.lgcy1, family = binomial)
flud.mod.spr1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        + fspr1CY2yrs + (1|fYear) + (1|SID), data = fi.dat.lgcy1, family = binomial)
allh.mod.spr1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        + fspr1CY2yrs + (1|fYear) + (1|SID), data = allhab.lgcy1, family = binomial)
chan.mod.spr2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        + fflgspr + (1|fYear) + (1|SID), data = ch.dat.lgcy, family = binomial)
shal.mod.spr2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        + fflgspr + (1|fYear) + (1|SID), data = ss.dat.lgcy, family = binomial)
flud.mod.spr2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        + fflgspr + (1|fYear) + (1|SID), data = fi.dat.lgcy, family = binomial)
allh.mod.spr2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        + fflgspr + (1|fYear) + (1|SID), data = allhab.lgcy, family = binomial)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### formatting and writing the statistics figures to be reported for model
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cm.sum = summary(chan.mod.spr1)
sm.sum = summary(shal.mod.spr1)
fm.sum = summary(flud.mod.spr1)
cm.coeff = cm.sum$coefficients
sm.coeff = sm.sum$coefficients
fm.coeff = fm.sum$coefficients
am.coeff = data.frame(rbind(cm.coeff, sm.coeff, fm.coeff))
am.coeff$varname = row.names(am.coeff)
colnames(am.coeff) = c("Estimate", "StdError", "Zvalue", "Pvalue", "VarName")
am.coeff.fm <- am.coeff %>% mutate(
  # Columns 1, 2, 3 as decimal with 4 digits
  Estimate = sprintf("%.3f", Estimate),
  StdError = sprintf("%.3f", StdError),
  Zvalue   = sprintf("%.2f", Zvalue),
  Pvalue   = sprintf("%.4e", Pvalue)
  # ".2e" means 2 digits after the decimal in scientific notation (e.g., 9.88e+06)
) %>% mutate(EstSE = paste0(Estimate, " (", StdError, ")"))
rand.eff = rbind(unlist(cbind(cm.sum$devcomp$dims[1], cm.sum$ngrps, cm.sum$varcor)),
                 unlist(cbind(sm.sum$devcomp$dims[1], sm.sum$ngrps, sm.sum$varcor)),
                 unlist(cbind(fm.sum$devcomp$dims[1], fm.sum$ngrps, fm.sum$varcor)))
colnames(rand.eff) = c("nObs", "nDud", "nSID", "nYr", "seSID", "seYr")
df.reff = data.frame(rand.eff) %>% mutate(
  # Columns 1, 2, 3 as decimal with 4 digits
  nObs  = sprintf("%d", nObs),
  nDud  = sprintf("%d", nDud),
  nSID  = sprintf("%d", nSID),
  nYr   = sprintf("%d", nYr),
  seSID = sprintf("%.3f", seSID),
  seYr  = sprintf("%.3f", seYr))

write_csv(am.coeff.fm, paste0(workdir, "fixeff_allhab_glmm_lgcy_CY2y_0408_noGlyph.csv"))
write_csv(df.reff,     paste0(workdir, "rndeff_allhab_glmm_lgcy_CY2y_0408_noGlyph.csv"))

###~~~~~~~~~~~~~~~~~~~~ FOR LEGACY EFFECT WITH Tflag ~~~~~~~~~~~~~~~~~~~~~~~~~~~

cm.sum = summary(chan.mod.spr2)
sm.sum = summary(shal.mod.spr2)
fm.sum = summary(flud.mod.spr2)
cm.coeff = cm.sum$coefficients
sm.coeff = sm.sum$coefficients
fm.coeff = fm.sum$coefficients
am.coeff = data.frame(rbind(cm.coeff, sm.coeff, fm.coeff))
am.coeff$varname = row.names(am.coeff)
colnames(am.coeff) = c("Estimate", "StdError", "Zvalue", "Pvalue", "VarName")
am.coeff.fm <- am.coeff %>% mutate(
  # Columns 1, 2, 3 as decimal with 4 digits
  Estimate = sprintf("%.3f", Estimate),
  StdError = sprintf("%.3f", StdError),
  Zvalue   = sprintf("%.2f", Zvalue),
  Pvalue   = sprintf("%.4e", Pvalue)
  # ".2e" means 2 digits after the decimal in scientific notation (e.g., 9.88e+06)
) %>% mutate(EstSE = paste0(Estimate, " (", StdError, ")"))
rand.eff = rbind(unlist(cbind(cm.sum$devcomp$dims[1], cm.sum$ngrps, cm.sum$varcor)),
                 unlist(cbind(sm.sum$devcomp$dims[1], sm.sum$ngrps, sm.sum$varcor)),
                 unlist(cbind(fm.sum$devcomp$dims[1], fm.sum$ngrps, fm.sum$varcor)))
colnames(rand.eff) = c("nObs", "nDud", "nSID", "nYr", "seSID", "seYr")
df.reff = data.frame(rand.eff) %>% mutate(
  # Columns 1, 2, 3 as decimal with 4 digits
  nObs  = sprintf("%d", nObs),
  nDud  = sprintf("%d", nDud),
  nSID  = sprintf("%d", nSID),
  nYr   = sprintf("%d", nYr),
  seSID = sprintf("%.3f", seSID),
  seYr  = sprintf("%.3f", seYr))

write_csv(am.coeff.fm, paste0(workdir, "fixeff_allhab_glmm_lgcy_Tflg_0408_noGlyph.csv"))
write_csv(df.reff,     paste0(workdir, "rndeff_allhab_glmm_lgcy_Tflg_0408_noGlyph.csv"))

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get model residuals
chan.rsdl.spr1 = resid(chan.mod.spr1, type = "pearson")
shal.rsdl.spr1 = resid(shal.mod.spr1, type = "pearson")
flud.rsdl.spr1 = resid(flud.mod.spr1, type = "pearson")
allh.rsdl.spr1 = resid(allh.mod.spr1, type = "pearson")
chan.rsdl.spr2 = resid(chan.mod.spr2, type = "pearson")
shal.rsdl.spr2 = resid(shal.mod.spr2, type = "pearson")
flud.rsdl.spr2 = resid(flud.mod.spr2, type = "pearson")
allh.rsdl.spr2 = resid(allh.mod.spr2, type = "pearson")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MODEL VALIDATION PLOTS FOR THE THREE FINAL HABITAT MODELS FOR LEGACY EFFECTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# validate the final model by plotting residuals and save to a tiff
tiff(filename = paste0("figures/", "allhab_glmm_validation_0408_lgcy_CY2yrs.tiff"), width = 8, height = 8, units = "in", res = 200)
op = par(mfrow = c(3, 3), mar = c(4,4,1,1), oma = c(0.5,2,1,0), cex = 1)

# CHANNEL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 1: Residuals vs. Fitted values; should be centered around 0
plot(fitted(chan.mod.spr1), chan.rsdl.spr1, xlab = "Predicted", ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(a)", side = 3, line = -1.5, adj = 0.95, cex = 1.2)
# Plot 2: histogram of the residuals; should be centered around 0
hist(chan.rsdl.spr1, xlab = 'Model Residuals', main = "", breaks = 15)
mtext("(b)", side = 3, line = -1.5, adj = 0.95, cex = 1.2)
# Plots 4,5,6: the Residuals vs. all the predictors; should be centered around 0
boxplot(chan.rsdl.spr1 ~ ch.dat.lgcy1$fspr1CY2yrs, xlab = bquote(CY2yrs[t-1]), ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(c)", side = 3, line = -1.5, adj = 0.95, cex = 1.2)

mtext(text = "Flooded Island", side = 2, line = 0.5, outer = TRUE, cex = 1.2, font = 2, adj = 0.15)
mtext(text = "Slow Shallows",  side = 2, line = 0.5, outer = TRUE, cex = 1.2, font = 2, adj = 0.55)
mtext(text = "Channel",        side = 2, line = 0.5, outer = TRUE, cex = 1.2, font = 2, adj = 0.9)

# SLOW SHALLOWS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 1: Residuals vs. Fitted values; should be centered around 0
plot(fitted(shal.mod.spr1), shal.rsdl.spr1, xlab = "Predicted", ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(d)", side = 3, line = -1.5, adj = 0.95, cex = 1.2)
# Plot 2: histogram of the residuals; should be centered around 0
hist(shal.rsdl.spr1, xlab = 'Model Residuals', main = "", breaks = 15)
mtext("(e)", side = 3, line = -1.5, adj = 0.95, cex = 1.2)
# Plots 4,5,6: the Residuals vs. all the predictors; should be centered around 0
boxplot(shal.rsdl.spr1 ~ ss.dat.lgcy1$fspr1CY2yrs, xlab = bquote(CY2yrs[t-1]), ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(f)", side = 3, line = -1.5, adj = 0.95, cex = 1.2)

# FLOODED ISLANDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 1: Residuals vs. Fitted values; should be centered around 0
plot(fitted(flud.mod.spr1), flud.rsdl.spr1, xlab = "Predicted", ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(g)", side = 3, line = -1.5, adj = 0.95, cex = 1.2)
# Plot 2: histogram of the residuals; should be centered around 0
hist(flud.rsdl.spr1, xlab = 'Model Residuals', main = "", breaks = 15)
mtext("(h)", side = 3, line = -1.5, adj = 0.95, cex = 1.2)
# Plots 4,5,6: the Residuals vs. all the predictors; should be centered around 0
boxplot(flud.rsdl.spr1 ~ fi.dat.lgcy1$fspr1CY2yrs, xlab = bquote(CY2yrs[t-1]), ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(i)", side = 3, line = -1.5, adj = 0.95, cex = 1.2)

dev.off()

#_______________________________________________________________________________

# validate the final model by plotting residuals and save to a tiff
tiff(filename = paste0("figures/", "allhab_glmm_validation_0408_lgcy_tflag.tiff"), width = 8, height = 8, units = "in", res = 200)
op = par(mfrow = c(3, 3), mar = c(4,4,1,1), oma = c(0.5,2,1,0), cex = 1)

# CHANNEL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 1: Residuals vs. Fitted values; should be centered around 0
plot(fitted(chan.mod.spr2), chan.rsdl.spr2, xlab = "Predicted", ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(a)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plot 2: histogram of the residuals; should be centered around 0
hist(chan.rsdl.spr2, xlab = 'Model Residuals', main = "", breaks = 15)
mtext("(b)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plots 4,5,6: the Residuals vs. all the predictors; should be centered around 0
boxplot(chan.rsdl.spr2 ~ ch.dat.lgcy$fflgspr, xlab = bquote(T[flag]), ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(c)", side = 1, line = -1.8, adj = 0.95, cex = 1.2)

mtext(text = "Flooded Island", side = 2, line = 0.5, outer = TRUE, cex = 1.2, font = 2, adj = 0.15)
mtext(text = "Slow Shallows",  side = 2, line = 0.5, outer = TRUE, cex = 1.2, font = 2, adj = 0.55)
mtext(text = "Channel",        side = 2, line = 0.5, outer = TRUE, cex = 1.2, font = 2, adj = 0.9)

# SLOW SHALLOWS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 1: Residuals vs. Fitted values; should be centered around 0
plot(fitted(shal.mod.spr2), shal.rsdl.spr2, xlab = "Predicted", ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(g)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plot 2: histogram of the residuals; should be centered around 0
hist(shal.rsdl.spr2, xlab = 'Model Residuals', main = "", breaks = 15)
mtext("(h)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plots 4,5,6: the Residuals vs. all the predictors; should be centered around 0
boxplot(shal.rsdl.spr2 ~ ss.dat.lgcy$fflgspr, xlab = bquote(T[flag]), ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(i)", side = 1, line = -1.5, adj = 0.95, cex = 1.2)

# FLOODED ISLANDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 1: Residuals vs. Fitted values; should be centered around 0
plot(fitted(flud.mod.spr2), flud.rsdl.spr2, xlab = "Predicted", ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(d)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plot 2: histogram of the residuals; should be centered around 0
hist(flud.rsdl.spr2, xlab = 'Model Residuals', main = "", breaks = 15)
mtext("(e)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)
# Plots 4,5,6: the Residuals vs. all the predictors; should be centered around 0
boxplot(flud.rsdl.spr2 ~ fi.dat.lgcy$fflgspr, xlab = bquote(T[flag]), ylab = "Model Residuals")
abline(h = 0, lwd = 2, col = gray(0.3,0.3))
mtext("(f)", side = 3, line = -1.2, adj = 0.95, cex = 1.2)

dev.off()

#_______________________________________________________________________________

all.spr1 = sort(unique(allhab.lgcy1$fspr1CY2yrs))
all.spr2 = sort(unique(allhab.lgcy$fflgspr))

# matrix to accept the fit and conf. int. of this specific case
pfav.ch.spr1 = array(data = NA, dim = c(length(all.spr1), 5),
                     dimnames = list(all.spr1, c("fit", "upper", "lower", "fspr1CY2yrs", "N")))
pfav.ss.spr1 = array(data = NA, dim = c(length(all.spr1), 5),
                     dimnames = list(all.spr1, c("fit", "upper", "lower", "fspr1CY2yrs", "N")))
pfav.fi.spr1 = array(data = NA, dim = c(length(all.spr1), 5),
                     dimnames = list(all.spr1, c("fit", "upper", "lower", "fspr1CY2yrs", "N")))
pfav.ch.spr2 = array(data = NA, dim = c(length(all.spr2), 5),
                     dimnames = list(all.spr2, c("fit", "upper", "lower", "fflgspr", "N")))
pfav.ss.spr2 = array(data = NA, dim = c(length(all.spr2), 5),
                     dimnames = list(all.spr2, c("fit", "upper", "lower", "fflgspr", "N")))
pfav.fi.spr2 = array(data = NA, dim = c(length(all.spr2), 5),
                     dimnames = list(all.spr2, c("fit", "upper", "lower", "fflgspr", "N")))

# for each habitat model and spray frequency, get predictions of FAV pixels
eff = allEffects(mod = chan.mod.spr1, xlevels = list(fspr1CY2yrs = all.spr1))
# extract the fitted values, and upper and lower CF bounds
pfav.ch.spr1[,1] = invlogit(eff$`fspr1CY2yrs`$fit)
pfav.ch.spr1[,2] = invlogit(eff$`fspr1CY2yrs`$upper)
pfav.ch.spr1[,3] = invlogit(eff$`fspr1CY2yrs`$lower)
pfav.ch.spr1[,4] = all.spr1
num = unlist(table(ch.dat.lgcy1$fspr1CY2yrs))
pfav.ch.spr1[,5] = num

# for each habitat model and spray frequency, get predictions of FAV pixels
eff = allEffects(mod = chan.mod.spr2, xlevels = list(fflgspr = all.spr2))
# extract the fitted values, and upper and lower CF bounds
pfav.ch.spr2[,1] = invlogit(eff$`fflgspr`$fit)
pfav.ch.spr2[,2] = invlogit(eff$`fflgspr`$upper)
pfav.ch.spr2[,3] = invlogit(eff$`fflgspr`$lower)
pfav.ch.spr2[,4] = all.spr2
num = unlist(table(ch.dat.lgcy$fflgspr))
pfav.ch.spr2[,5] = num

# for each habitat model and spray frequency, get predictions of FAV pixels
eff = allEffects(mod = shal.mod.spr1, xlevels = list(fspr1CY2yrs = all.spr1))
# extract the fitted values, and upper and lower CF bounds
pfav.ss.spr1[,1] = invlogit(eff$`fspr1CY2yrs`$fit)
pfav.ss.spr1[,2] = invlogit(eff$`fspr1CY2yrs`$upper)
pfav.ss.spr1[,3] = invlogit(eff$`fspr1CY2yrs`$lower)
pfav.ss.spr1[,4] = all.spr1
num = unlist(table(ss.dat.lgcy1$fspr1CY2yrs))
pfav.ss.spr1[,5] = num

# for each habitat model and spray frequency, get predictions of FAV pixels
eff = allEffects(mod = shal.mod.spr2, xlevels = list(fflgspr = all.spr2))
# extract the fitted values, and upper and lower CF bounds
pfav.ss.spr2[,1] = invlogit(eff$`fflgspr`$fit)
pfav.ss.spr2[,2] = invlogit(eff$`fflgspr`$upper)
pfav.ss.spr2[,3] = invlogit(eff$`fflgspr`$lower)
pfav.ss.spr2[,4] = all.spr2
num = unlist(table(ss.dat.lgcy$fflgspr))
pfav.ss.spr2[,5] = num

# for each habitat model and spray frequency, get predictions of FAV pixels
eff = allEffects(mod = flud.mod.spr1, xlevels = list(fspr1CY2yrs = all.spr1))
# extract the fitted values, and upper and lower CF bounds
pfav.fi.spr1[,1] = invlogit(eff$`fspr1CY2yrs`$fit)
pfav.fi.spr1[,2] = invlogit(eff$`fspr1CY2yrs`$upper)
pfav.fi.spr1[,3] = invlogit(eff$`fspr1CY2yrs`$lower)
pfav.fi.spr1[,4] = all.spr1
num = unlist(table(fi.dat.lgcy1$fspr1CY2yrs))
pfav.fi.spr1[,5] = num

# for each habitat model and spray frequency, get predictions of FAV pixels
eff = allEffects(mod = flud.mod.spr2, xlevels = list(fflgspr = all.spr2))
# extract the fitted values, and upper and lower CF bounds
pfav.fi.spr2[,1] = invlogit(eff$`fflgspr`$fit)
pfav.fi.spr2[,2] = invlogit(eff$`fflgspr`$upper)
pfav.fi.spr2[,3] = invlogit(eff$`fflgspr`$lower)
pfav.fi.spr2[,4] = all.spr2
num = unlist(table(fi.dat.lgcy$fflgspr))
pfav.fi.spr2[,5] = num

#_______________________________________________________________________________
#     CREATE FIGURE FOR LEGACY EFFECTS OF HERBICIDE TREATMENT
tiff(filename = "figures/FAV_X24D_lgcy_allh_0408_CY2yrs.tiff", width = 8, height = 3, units = "in", res = 300)
#_______________________________________________________________________________

# set outer (oma(bottom, left, top, right)) and inner margins (mar(bottom, left, top, right))
# in terms of lines of text; + 0.1 is a slight offset
par(mfcol = c(1,3), mar = c(1.5,1.5,1,1) + 0.1, oma = c(3.5,4.5,3,0.5))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                  CHANNEL HABITAT; GLYPHOSATE; 2014-2018
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot just for the first value of glph which is 0
plot.cy0 = c(pfav.ch.spr1[1, 1], pfav.ch.spr1[1, 2], pfav.ch.spr1[1, 3])
# plot all points, adjust axis limits, color, labeling, etc.
plot(x = pfav.ch.spr1[,4], y = pfav.ch.spr1[,1], pch = 16, ylim = c(0,0.8), xlim = c(0.5, 4.5),
     las = 0.5, ylab = "", xlab = "", main = "", col = rgb(0.8, 0.2, 0), xaxt="n")
# add x-axis labels without the ticks
axis(1, at = 1:4, labels = c('Ref', '1', '2', '3'), tick = FALSE)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.ch.spr1)) {
  segments(pfav.ch.spr1[i, 4], pfav.ch.spr1[i, 3], pfav.ch.spr1[i, 4], pfav.ch.spr1[i, 2],
           col = rgb(1,0.7,0.5, alpha = 0.8), lty = 1, lwd = 2)
  # make the top and bottom tiny line segments of the distribution
  arrows(x0 = pfav.ch.spr1[i, 4], y0 = pfav.ch.spr1[i, 3], x1 = pfav.ch.spr1[i, 4], y1 = pfav.ch.spr1[i, 2], 
         angle = 90, code = 3, length = 0.1, col = rgb(0.8, 0.2, 0))
  # add sample size on top of each line segment
  text(pfav.ch.spr1[i,4], pfav.ch.spr1[i,2], paste("N =", pfav.ch.spr1[i, 5]), pos = 3, offset = 0.7)
}
# draw each fitted point with color changing according to spray frequency
points(pfav.ch.spr1[,4], pfav.ch.spr1[,1], pch = 16, cex = 1, col = rgb(0.8, 0.2, 0))

#_____________________ LABEL ALL AXES AND ALL PLOTS

mtext(text = "Pr(FAV occurrence)", side = 2, line = 1, outer = TRUE, cex = 0.9, adj = 0.35)

mtext(text = bquote(CY2yrs[t-1]), side = 1, line = 1.5, outer = TRUE, cex = 0.9, adj = 0.15)
mtext(text = bquote(CY2yrs[t-1]), side = 1, line = 1.5, outer = TRUE, cex = 0.9, adj = 0.5)
mtext(text = bquote(CY2yrs[t-1]), side = 1, line = 1.5, outer = TRUE, cex = 0.9, adj = 0.85)

mtext(text = "(a) Channels",        side = 3, line = 0, outer = TRUE, cex = 1, font = 2, adj = 0.1)
mtext(text = "(b) Slow Shallows",   side = 3, line = 0, outer = TRUE, cex = 1, font = 2, adj = 0.5)
mtext(text = "(c) Flooded Islands", side = 3, line = 0, outer = TRUE, cex = 1, font = 2, adj = 0.95)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   SLOW SHALLOWS HABITAT; GLYPHOSATE; 2014-2018; MEDIAN SPRAY FREQ; ALL CY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot just for the first value of glph which is 0
plot.cy0 = c(pfav.ss.spr1[1, 1], pfav.ss.spr1[1, 2], pfav.ss.spr1[1, 3])
plot(x = pfav.ss.spr1[,4], y = pfav.ss.spr1[,1], pch = 16, ylim = c(0,0.8), xlim = c(0.5, 4.5),
     las = 0.5, ylab = "", xlab = "", main = "", col = rgb(0.8, 0.2, 0), xaxt="n")
# add x-axis labels without the ticks
axis(1, at = 1:4, labels = c('Ref', '1', '2', '3'), tick = FALSE)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.ss.spr1)) {
  segments(pfav.ss.spr1[i, 4], pfav.ss.spr1[i, 3], pfav.ss.spr1[i, 4], pfav.ss.spr1[i, 2],
           col = rgb(1,0.7,0.5, alpha = 0.8), lty = 1, lwd = 2)
  # make the top and bottom tiny line segments of the distribution
  arrows(x0 = pfav.ss.spr1[i, 4], y0 = pfav.ss.spr1[i, 3], x1 = pfav.ss.spr1[i, 4], y1 = pfav.ss.spr1[i, 2], 
         angle = 90, code = 3, length = 0.1, col = rgb(0.8, 0.2, 0))
  # add sample size on top of each line segment
  text(pfav.ss.spr1[i,4], pfav.ss.spr1[i,2], paste("N =", pfav.ss.spr1[i, 5]), pos = 3, offset = 0.7)
}
# draw each fitted point with color changing according to spray frequency
points(pfav.ss.spr1[,4], pfav.ss.spr1[,1], pch = 16, cex = 1.2, col = rgb(0.8, 0.2, 0))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   FLOODED ISLAND HABITAT; GLYPHOSATE; 2014-2018; ALL SPRAY FREQ; ALL CY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot just for the first value of glph which is 0
plot.cy0 = c(pfav.fi.spr1[1, 1], pfav.fi.spr1[1, 2], pfav.fi.spr1[1, 3])
# plot a point for fitted value for control sites
plot(x = pfav.fi.spr1[,4], y = pfav.fi.spr1[,1], pch = 16, ylim = c(0,0.8), xlim = c(0.5, 4.5),
     las = 0.5, ylab = "", xlab = "", main = "", col = rgb(0.8, 0.2, 0), xaxt="n")
# add x-axis labels without the ticks
axis(1, at = 1:4, labels = c('Ref', '1', '2', '3'), tick = FALSE)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.fi.spr1)) {
  segments(pfav.fi.spr1[i, 4], pfav.fi.spr1[i, 3], pfav.fi.spr1[i, 4], pfav.fi.spr1[i, 2],
           col = rgb(1,0.7,0.5, alpha = 0.8), lty = 1, lwd = 2)
  # make the top and bottom tiny line segments of the distribution
  arrows(x0 = pfav.fi.spr1[i, 4], y0 = pfav.fi.spr1[i, 3], x1 = pfav.fi.spr1[i, 4], y1 = pfav.fi.spr1[i, 2], 
         angle = 90, code = 3, length = 0.1, col = rgb(0.8, 0.2, 0))
  # add sample size on top of each line segment
  text(pfav.fi.spr1[i,4], pfav.fi.spr1[i,2], paste("N =", pfav.fi.spr1[i, 5]), pos = 3, offset = 0.7)
}
# draw each fitted point with color changing according to spray frequency
points(pfav.fi.spr1[,4], pfav.fi.spr1[,1], pch = 16, cex = 1.2, col = rgb(0.8, 0.2, 0))

dev.off()

#_______________________________________________________________________________
#     CREATE FIGURE FOR LEGACY EFFECTS OF HERBICIDE TREATMENT
tiff(filename = "figures/FAV_X24D_lgcy_allh_0408_Tflag.tiff", width = 8, height = 3, units = "in", res = 300)
#_______________________________________________________________________________

# set outer (oma(bottom, left, top, right)) and inner margins (mar(bottom, left, top, right))
# in terms of lines of text; + 0.1 is a slight offset
par(mfcol = c(1,3), mar = c(1.5,1.5,1,1) + 0.1, oma = c(3.5,4.5,3,0.5))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                  CHANNEL HABITAT; 2,4-D; 2004-2008
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot just for the first value of X24D which is 0
plot.cy0 = c(pfav.ch.spr2[1, 1], pfav.ch.spr2[1, 2], pfav.ch.spr2[1, 3])
# plot all points, adjust axis limits, color, labeling, etc.
plot(x = pfav.ch.spr2[,4], y = pfav.ch.spr2[,1], pch = 16, ylim = c(0,0.25), xlim = c(0.5, 4.5),
     las = 0.5, ylab = "", xlab = "", main = "", col = rgb(0.8, 0.2, 0), xaxt="n")
# add x-axis labels without the ticks
axis(1, at = 1:4, labels = c('Ref', '1', '2', '3'), tick = FALSE)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.ch.spr2)) {
  segments(pfav.ch.spr2[i, 4], pfav.ch.spr2[i, 3], pfav.ch.spr2[i, 4], pfav.ch.spr2[i, 2],
           col = rgb(1,0.7,0.5, alpha = 0.8), lty = 1, lwd = 2)
  # make the top and bottom tiny line segments of the distribution
  arrows(x0 = pfav.ch.spr2[i, 4], y0 = pfav.ch.spr2[i, 3], x1 = pfav.ch.spr2[i, 4], y1 = pfav.ch.spr2[i, 2], 
         angle = 90, code = 3, length = 0.1, col = rgb(0.8, 0.2, 0))
  # add sample size on top of each line segment
  text(pfav.ch.spr2[i,4], pfav.ch.spr2[i,2], paste("N =", pfav.ch.spr2[i, 5]), pos = 3, offset = 0.7)
}
# draw each fitted point with color changing according to spray frequency
points(pfav.ch.spr2[,4], pfav.ch.spr2[,1], pch = 16, cex = 1, col = rgb(0.8, 0.2, 0))

#_____________________ LABEL ALL AXES AND ALL PLOTS

mtext(text = "Pr(FAV occurrence)", side = 2, line = 1, outer = TRUE, cex = 0.9, adj = 0.35)

mtext(text = bquote(T[flag]), side = 1, line = 1.5, outer = TRUE, cex = 0.9, adj = 0.15)
mtext(text = bquote(T[flag]), side = 1, line = 1.5, outer = TRUE, cex = 0.9, adj = 0.5)
mtext(text = bquote(T[flag]), side = 1, line = 1.5, outer = TRUE, cex = 0.9, adj = 0.85)

mtext(text = "(a) Channels",        side = 3, line = 0, outer = TRUE, cex = 1, font = 2, adj = 0.1)
mtext(text = "(b) Slow Shallows",   side = 3, line = 0, outer = TRUE, cex = 1, font = 2, adj = 0.5)
mtext(text = "(c) Flooded Islands", side = 3, line = 0, outer = TRUE, cex = 1, font = 2, adj = 0.95)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   SLOW SHALLOWS HABITAT; GLYPHOSATE; 2014-2018; MEDIAN SPRAY FREQ; ALL CY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot just for the first value of X24D which is 0
plot.cy0 = c(pfav.ss.spr2[1, 1], pfav.ss.spr2[1, 2], pfav.ss.spr2[1, 3])
plot(x = pfav.ss.spr2[,4], y = pfav.ss.spr2[,1], pch = 16, ylim = c(0,0.25), xlim = c(0.5, 4.5),
     las = 0.5, ylab = "", xlab = "", main = "", col = rgb(0.8, 0.2, 0), xaxt="n")
# add x-axis labels without the ticks
axis(1, at = 1:4, labels = c('Ref', '1', '2', '3'), tick = FALSE)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.ss.spr2)) {
  segments(pfav.ss.spr2[i, 4], pfav.ss.spr2[i, 3], pfav.ss.spr2[i, 4], pfav.ss.spr2[i, 2],
           col = rgb(1,0.7,0.5, alpha = 0.8), lty = 1, lwd = 2)
  # make the top and bottom tiny line segments of the distribution
  arrows(x0 = pfav.ss.spr2[i, 4], y0 = pfav.ss.spr2[i, 3], x1 = pfav.ss.spr2[i, 4], y1 = pfav.ss.spr2[i, 2], 
         angle = 90, code = 3, length = 0.1, col = rgb(0.8, 0.2, 0))
  # add sample size on top of each line segment
  text(pfav.ss.spr2[i,4], pfav.ss.spr2[i,2], paste("N =", pfav.ss.spr2[i, 5]), pos = 3, offset = 0.7)
}
# draw each fitted point with color changing according to spray frequency
points(pfav.ss.spr2[,4], pfav.ss.spr2[,1], pch = 16, cex = 1.2, col = rgb(0.8, 0.2, 0))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   FLOODED ISLAND HABITAT; GLYPHOSATE; 2014-2018; ALL SPRAY FREQ; ALL CY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot just for the first value of X24D which is 0
plot.cy0 = c(pfav.fi.spr2[1, 1], pfav.fi.spr2[1, 2], pfav.fi.spr2[1, 3])
# plot a point for fitted value for control sites
plot(x = pfav.fi.spr2[,4], y = pfav.fi.spr2[,1], pch = 16, ylim = c(0,0.25), xlim = c(0.5, 4.5),
     las = 0.5, ylab = "", xlab = "", main = "", col = rgb(0.8, 0.2, 0), xaxt="n")
# add x-axis labels without the ticks
axis(1, at = 1:4, labels = c('Ref', '1', '2', '3'), tick = FALSE)
# extend the lines for the upper and lower confidence interval
abline(h = plot.cy0[2], lwd = 1, lty = 5, col = "gray")
abline(h = plot.cy0[3], lwd = 1, lty = 5, col = "gray")
# Add vertical lines to show the confidence interval around each fitted value
for (i in 1:nrow(pfav.fi.spr2)) {
  segments(pfav.fi.spr2[i, 4], pfav.fi.spr2[i, 3], pfav.fi.spr2[i, 4], pfav.fi.spr2[i, 2],
           col = rgb(1,0.7,0.5, alpha = 0.8), lty = 1, lwd = 2)
  # make the top and bottom tiny line segments of the distribution
  arrows(x0 = pfav.fi.spr2[i, 4], y0 = pfav.fi.spr2[i, 3], x1 = pfav.fi.spr2[i, 4], y1 = pfav.fi.spr2[i, 2], 
         angle = 90, code = 3, length = 0.1, col = rgb(0.8, 0.2, 0))
  # add sample size on top of each line segment
  text(pfav.fi.spr2[i,4], pfav.fi.spr2[i,2], paste("N =", pfav.fi.spr2[i, 5]), pos = 3, offset = 0.7)
}
# draw each fitted point with color changing according to spray frequency
points(pfav.fi.spr2[,4], pfav.fi.spr2[,1], pch = 16, cex = 1.2, col = rgb(0.8, 0.2, 0))

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            FIT THE INTERACTION MODEL FOR FLOODED ISLANDS ONLY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

flood_dat = subset(dat.scl, dat.scl$AggHabitat == "FloodedIsland")

# predict based on each individually - glyphosate, 2,4-D, CY and spray frequency
# random effect of year and site
flood_mod.1.1   = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                          log10_X24D +
                          (1|fYear) + (1|SID), data = flood_dat, family = binomial)
flood_mod.1.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_glph +
                        (1|fYear) + (1|SID), data = flood_dat, family = binomial)
flood_mod.1.3  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         cbrt_sprfreq +
                         (1|fYear) + (1|SID), data = flood_dat, family = binomial)
flood_mod.1.4  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         CY2yrs +
                         (1|fYear) + (1|SID), data = flood_dat, family = binomial)

BIC(flood_mod.1.1, flood_mod.1.2, flood_mod.1.3, flood_mod.1.4)
# Best model: flood_mod.1.1 - first var chosen: log10_X24D

flood_mod.2.1   = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                          log10_X24D + log10_glph +
                          (1|fYear) + (1|SID), data = flood_dat, family = binomial)
flood_mod.2.2  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         log10_X24D + CY2yrs +
                         (1|fYear) + (1|SID), data = flood_dat, family = binomial)
flood_mod.2.3 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + cbrt_sprfreq +
                        (1|fYear) + (1|SID), data = flood_dat, family = binomial)
BIC(flood_mod.1.1, flood_mod.2.1, flood_mod.2.2, flood_mod.2.3)
# Best model: flood_mod.2.3 - second var chosen: cbrt_sprfreq

flood_mod.3.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + cbrt_sprfreq + log10_glph +
                        (1|fYear) + (1|SID), data = flood_dat, family = binomial)
flood_mod.3.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + cbrt_sprfreq + CY2yrs +
                        (1|fYear) + (1|SID), data = flood_dat, family = binomial)
BIC(flood_mod.2.3, flood_mod.3.1, flood_mod.3.2)
# Best model: flood_mod.3.2 - third var chosen: CY2yrs

flood_mod.4.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + cbrt_sprfreq + CY2yrs + log10_glph +
                        (1|fYear) + (1|SID), data = flood_dat, family = binomial)
BIC(flood_mod.3.2, flood_mod.4.1)
# Best model: flood_mod.4.1 - fourth var chosen: log10_glph

# interaction - yes or no
flood_mod_all = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~ log10_X24D + log10_glph +
                        log10_X24D * log10_glph + cbrt_sprfreq + CY2yrs +
                        (1|fYear) + (1|SID), data = flood_dat, family = binomial)

BIC(flood_mod.3.2, flood_mod.4.1, flood_mod_all)
# best model is with everything

#===================================================================================#
# MODEL TESTING FOR LEGACY EFFECTS  #  ARE IMPACTS DETECTABLE 1, 2 OR 3 YRS LATER?  #
#===================================================================================#

ch.dat.lgcy1 =  chan_dat %>% filter(CY == 0) %>% mutate(flgspr = ifelse(spr2cy != 0, 2, ifelse(spr1cy != 0, 1, 0))) %>% filter(spr2cy == 0)
ss.dat.lgcy1 =  shal_dat %>% filter(CY == 0) %>% mutate(flgspr = ifelse(spr2cy != 0, 2, ifelse(spr1cy != 0, 1, 0))) %>% filter(spr2cy == 0)
fi.dat.lgcy1 = flood_dat %>% filter(CY == 0) %>% mutate(flgspr = ifelse(spr2cy != 0, 2, ifelse(spr1cy != 0, 1, 0))) %>% filter(spr2cy == 0)

ch.dat.lgcy2 =  chan_dat %>% filter(CY == 0) %>% mutate(flgspr = ifelse(spr2cy != 0, 2, ifelse(spr1cy != 0, 1, 0))) %>% filter(spr1cy == 0)
ss.dat.lgcy2 =  shal_dat %>% filter(CY == 0) %>% mutate(flgspr = ifelse(spr2cy != 0, 2, ifelse(spr1cy != 0, 1, 0))) %>% filter(spr1cy == 0)
fi.dat.lgcy2 = flood_dat %>% filter(CY == 0) %>% mutate(flgspr = ifelse(spr2cy != 0, 2, ifelse(spr1cy != 0, 1, 0))) %>% filter(spr1cy == 0)

table(ch.dat.lgcy1$flgspr)
table(ss.dat.lgcy1$flgspr)
table(fi.dat.lgcy1$flgspr)

table(ch.dat.lgcy2$flgspr)
table(ss.dat.lgcy2$flgspr)
table(fi.dat.lgcy2$flgspr)

chan_mod.1.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        + flgspr + (1|fYear) + (1|SID), data = ch.dat.lgcy1, family = binomial)
chan_mod.1.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        + spr1cy + (1|fYear) + (1|SID), data = ch.dat.lgcy1, family = binomial)
chan_mod.1.3 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + spr1CY7 + (1|fYear) + (1|SID), data = ch.dat.lgcy1, family = binomial)

BIC(chan_mod.1.1, chan_mod.1.2, chan_mod.1.3)
summary(chan_mod.1.2)

shal_mod.1.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + flgspr + (1|fYear) + (1|SID), data = ss.dat.lgcy1, family = binomial)
shal_mod.1.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + spr1cy + (1|fYear) + (1|SID), data = ss.dat.lgcy1, family = binomial)
shal_mod.1.3 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + spr1CY7 + (1|fYear) + (1|SID), data = ss.dat.lgcy1, family = binomial)

BIC(shal_mod.1.1, shal_mod.1.2, shal_mod.1.3)
summary(shal_mod.1.2)

flud_mod.1.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + flgspr + (1|fYear) + (1|SID), data = fi.dat.lgcy1, family = binomial)
flud_mod.1.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + spr1cy + (1|fYear) + (1|SID), data = fi.dat.lgcy1, family = binomial)
flud_mod.1.3 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + spr1CY7 + (1|fYear) + (1|SID), data = fi.dat.lgcy1, family = binomial)

BIC(flud_mod.1.1, flud_mod.1.2, flud_mod.1.3)
summary(flud_mod.1.1)
summary(flud_mod.1.2)

plot(chan_mod.1.2)
plot(shal_mod.1.2)
plot(flud_mod.1.2)

chan_mod.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + flgspr + (1|fYear) + (1|SID), data = ch.dat.lgcy2, family = binomial)
chan_mod.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + spr2cy + (1|fYear) + (1|SID), data = ch.dat.lgcy2, family = binomial)
chan_mod.3 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                      + spr2CY7 + (1|fYear) + (1|SID), data = ch.dat.lgcy2, family = binomial)

BIC(chan_mod.1, chan_mod.2, chan_mod.3)
summary(chan_mod.1)
summary(chan_mod.2)

shal_mod.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + flgspr + (1|fYear) + (1|SID), data = ss.dat.lgcy2, family = binomial)
shal_mod.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + spr2cy + (1|fYear) + (1|SID), data = ss.dat.lgcy2, family = binomial)
shal_mod.3 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                      + spr2CY7 + (1|fYear) + (1|SID), data = ss.dat.lgcy2, family = binomial)

BIC(shal_mod.1, shal_mod.2, shal_mod.3)
summary(shal_mod.1)
summary(shal_mod.2)

flud_mod.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + flgspr + (1|fYear) + (1|SID), data = fi.dat.lgcy2, family = binomial)
flud_mod.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       + spr2cy + (1|fYear) + (1|SID), data = fi.dat.lgcy2, family = binomial)
flud_mod.3 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                      + spr2CY7 + (1|fYear) + (1|SID), data = fi.dat.lgcy2, family = binomial)

BIC(flud_mod.1, flud_mod.2, flud_mod.3)
summary(flud_mod.1)
summary(flud_mod.2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       FIGURE OUT IF SPRAY FREQUENCY IS JUST A REFLECTION OF PATCH SIZE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# dataset: chan_dat, shal_dat, flood_dat

flood_dat_sbst = flood_dat[(flood_dat$Year != 2004 & flood_dat$Year != 2014),]

# linear regression between patch size and spray frequency
# R2 really low so no direct relationship
lmmod = lm(log10_pszfav ~ cbrt_sprfreq, data = dat.scl)

# predict based on each individually - glyphosate, 2,4-D, CY and spray frequency
# random effect of year and site
flood_mod.1.5  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         log10_pszsav +
                         (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.1.6  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         log10_pszfav +
                         (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.1.7   = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                          log10_X24D +
                          (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.1.8 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_glph +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.1.9  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         cbrt_sprfreq +
                         (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.1.10  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         CY2yrs +
                         (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)

BIC(flood_mod.1.5, flood_mod.1.6, flood_mod.1.7, flood_mod.1.8, flood_mod.1.9, flood_mod.1.10)
# Best model was: flood_mod.1.1 - first var chosen: log10_X24D
# Best model  is: flood_mod.1.7 - first var chosen: log10_X24D

flood_mod.2.4   = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                          log10_X24D + log10_glph +
                          (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.2.5  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         log10_X24D + CY2yrs +
                         (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.2.6   = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                          log10_X24D + cbrt_sprfreq +
                          (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.2.7  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         log10_X24D + log10_pszsav +
                         (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.2.8 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)

# Best model was: flood_mod.2.3 - second var chosen: cbrt_sprfreq
# Best model  is: flood_mod.2.8 - second var chosen: log10_pszfav
BIC(flood_mod.1.7, flood_mod.2.4, flood_mod.2.5, flood_mod.2.6, flood_mod.2.7, flood_mod.2.8)

flood_mod.3.0 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav + log10_X24D*log10_pszfav +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.3.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav + log10_glph +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.3.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav + CY2yrs +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.3.3 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav + cbrt_sprfreq +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.3.4 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav + log10_pszsav +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
# Best model was: flood_mod.3.2 - third var chosen: CY2yrs
# Best model  is: flood_mod.3.3 - third var chosen: cbrt_sprfreq
BIC(flood_mod.2.8, flood_mod.3.0, flood_mod.3.1, flood_mod.3.2, flood_mod.3.3, flood_mod.3.4)

flood_mod.4.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav + cbrt_sprfreq + CY2yrs +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.4.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav + cbrt_sprfreq + log10_glph +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.4.3 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav + cbrt_sprfreq + log10_pszsav +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.4.4 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav + cbrt_sprfreq + log10_X24D*log10_pszfav +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.4.5 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav + cbrt_sprfreq + log10_pszfav*cbrt_sprfreq +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
flood_mod.4.6 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav + cbrt_sprfreq + log10_X24D*cbrt_sprfreq +
                        (1|fYear) + (1|SID), data = flood_dat_sbst, family = binomial)
# Best model was: flood_mod.4.1 - fourth var chosen: log10_glph
# Best model  is: flood_mod.4.3 - fourth var chosen: log10_pszsav
BIC(flood_mod.3.3, flood_mod.4.1, flood_mod.4.2, flood_mod.4.3, flood_mod.4.4, flood_mod.4.5, flood_mod.4.6)

#-----------------------------------------------------------------------------------------------------------

shal_dat_sbst = shal_dat[(shal_dat$Year != 2004 & shal_dat$Year != 2014),]

# linear regression between patch size and spray frequency
# R2 really low so no direct relationship
lmmod = lm(log10_pszfav ~ cbrt_sprfreq, data = shal_dat)

# predict based on each individually - glyphosate, 2,4-D, CY and spray frequency
# random effect of year and site
shal_mod.1.5  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         log10_pszsav +
                         (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.1.6  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         log10_pszfav +
                         (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.1.7   = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                          log10_X24D +
                          (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.1.8 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_glph +
                        (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.1.9  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         cbrt_sprfreq +
                         (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.1.10  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                          CY2yrs +
                          (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)

BIC(shal_mod.1.5, shal_mod.1.6, shal_mod.1.7, shal_mod.1.8, shal_mod.1.9, shal_mod.1.10)
# Best model was: shal_mod.1.1 - first var chosen: log10_X24D
# Best model  is: shal_mod.1.7 - first var chosen: log10_X24D

shal_mod.2.4   = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                          log10_X24D + log10_glph +
                          (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.2.5  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         log10_X24D + CY2yrs +
                         (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.2.6   = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                          log10_X24D + cbrt_sprfreq +
                          (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.2.7  = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                         log10_X24D + log10_pszsav +
                         (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.2.8 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + log10_pszfav +
                        (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)

# Best model was: shal_mod.2.3 - second var chosen: cbrt_sprfreq
# Best model  is: shal_mod.2.5 - second var chosen: CY2yrs
BIC(shal_mod.1.7, shal_mod.2.4, shal_mod.2.5, shal_mod.2.6, shal_mod.2.7, shal_mod.2.8)

shal_mod.3.0 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + CY2yrs + log10_X24D*CY2yrs +
                        (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.3.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + CY2yrs + log10_glph +
                        (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.3.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + CY2yrs + log10_pszfav +
                        (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.3.3 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + CY2yrs + cbrt_sprfreq +
                        (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.3.4 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                        log10_X24D + CY2yrs + log10_pszsav +
                        (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
# Best model was: shal_mod.3.2 - third var chosen: CY2yrs
# Best model  is: shal_mod.3.4 - third var chosen: log10_pszsav
BIC(shal_mod.2.8, shal_mod.3.0, shal_mod.3.1, shal_mod.3.2, shal_mod.3.3, shal_mod.3.4)

shal_mod.4.1 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       log10_X24D + CY2yrs + log10_pszsav + log10_glph +
                       (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.4.2 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       log10_X24D + CY2yrs + log10_pszsav + cbrt_sprfreq +
                       (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.4.3 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       log10_X24D + CY2yrs + log10_pszsav + log10_pszfav +
                       (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.4.4 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       log10_X24D + CY2yrs + log10_pszsav + log10_X24D*CY2yrs +
                       (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.4.5 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       log10_X24D + CY2yrs + log10_pszsav + CY2yrs*log10_pszsav +
                       (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
shal_mod.4.6 = glmer(cbind(FAV_pixel_count, non_FAV_pixel_count) ~
                       log10_X24D + CY2yrs + log10_pszsav + log10_pszsav*log10_X24D +
                       (1|fYear) + (1|SID), data = shal_dat_sbst, family = binomial)
# Best model was: shal_mod.4.1 - fourth var chosen: log10_glph
# Best model  is: shal_mod.4.2 - fourth var chosen: cbrt_sprfreq
BIC(shal_mod.3.4, shal_mod.4.1, shal_mod.4.2, shal_mod.4.3, shal_mod.4.4, shal_mod.4.5, shal_mod.4.6)


