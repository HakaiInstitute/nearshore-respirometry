###############################################################################
# Exploratory Nearshore Respirometry R Script #
# Code adapted from Nate S, UFloridasee package `respR`                       #
# Created by AO and OP on Oct 11, 2020 #
# Automate Respirometry-Related Data. https://github.com/januarharianto/respr,#

rm(list = ls())

#install.packages("devtools")
#devtools::install_github("januarharianto/respR")
#library(respR)
lapply(c("tidyr", "plyr", "dplyr", "ggplot2", "magrittr", "respR",
         "lubridate", "knitr", "tidyverse", "reshape2", "RcppRoll", "readr"), 
       library, character.only = T)
library(tidyverse)

setwd("~/Repository/nearshore-respirometry")



#Framework

#1. open csv for sample1 and corresponding bg csv

#2. apply function1 : slices df in two 
# based on first 30 minutes, dark treatment (500 rows from the start) 
# and the last 30 minutes, light treatment (500 rows from the end)

#3. apply function2 : inspects all 4 chunks for errors
#(sample1_dark, sample1_light, bg_dark, bg_light)

#4. apply function3 : calculates rates for both treatment by 
# a. calculating bg rate for dark, 
# b. calculating bg rate for light,
# c. calculating sample rate for dark,
# d. calculates sample rate for light, 
# e. subtracts sample rate by bg rate for dark,
# f. subtracts sample rate by bg rate for light

#5. apply function4 : converts final rates for light and dark based on metadata volumes




##1. 
metadata <- read_csv("metadata_OP.csv")
S01_M3_002 <- read_csv("01_M3_002.csv") 
S01_C <- read_csv("01_C.csv")

##2. 
S01_M3_002.rd <- S01_M3_002 %>%
  rename(oxygen = Oxygen, datetime = Date, time = "Delta T [min]")
#renames variables - sample data
S01_C.rd <- S01_C %>%
  rename(oxygen = Oxygen, datetime = Date, time = "Delta T [min]")
#renames variables - bg data

S01_M3_002_dark <- subset(S01_M3_002.rd, time >5 & time <30)
S01_M3_002_light <- subset(S01_M3_002.rd, time >35)
#separates dark and light treatment - sample data
S01_C_dark <- subset(S01_C.rd, time >5 & time <30)
S01_C_light <- subset(S01_C.rd, time >35)
#separates dark and light treatment - bg data


##3.
insp_S01_M3_002_dark <- inspect(S01_M3_002_dark, time = 3, oxygen = 7)
insp_S01_M3_002_dark
insp_S01_M3_002_light <- inspect(S01_M3_002_light, time = 3, oxygen = 7)
insp_S01_M3_002_light


insp_S01_C_dark <- inspect(S01_C_dark, time = 3, oxygen = 7)
insp_S01_C_dark
insp_S01_C_light <- inspect(S01_C_light, time = 3, oxygen = 7)
insp_S01_C_light

##4.
rate_S01_M3_002_dark <- auto_rate(insp_S01_M3_002_dark)
print(rate_S01_M3_002_dark)
rate_S01_M3_002_light <- auto_rate(insp_S01_M3_002_light)
print(rate_S01_M3_002_light)


rate_S01_C_dark <- calc_rate.bg(insp_S01_C_dark)
rate_S01_C_light <- calc_rate.bg(insp_S01_C_light)


adj_rate_S01_M3_002_dark <- adjust_rate(rate_S01_M3_002_dark, rate_S01_C_dark)
print(adj_rate_S01_M3_002_dark)
adj_rate_S01_M3_002_light <- adjust_rate(rate_S01_M3_002_light, rate_S01_C_light)
print(adj_rate_S01_M3_002_light)

##5. 
V01_M3_002 <- metadata %>%
  filter(sampleID == "01_M3_002")

dark <- convert_rate(adj_rate_S01_M3_002_dark, 
             o2.unit = "mg/L", 
             time.unit = "min", 
             output.unit = "mg/h/g", 
             volume = (V01_M3_002$volume_chamber - V01_M3_002$volume_kelp),
             mass = V01_M3_002$wet_weight_g) 
print(dark)


light <- convert_rate(adj_rate_S01_M3_002_light, 
                     o2.unit = "mg/L", 
                     time.unit = "min", 
                     output.unit = "mg/h/g", 
                     volume = (V01_M3_002$volume_chamber - V01_M3_002$volume_kelp), 
                     mass = V01_M3_002$wet_weight_g) 
print(light)

