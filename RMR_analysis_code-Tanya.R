

## Oct 15 2020
# Tanya Prinzing Shark Respirometry 
# --- Resting metabolic rate calculation script  -- #




# Load libraries
library("readr")
library("tidyverse")
library("rollRegres")
library("lubridate")
library("respR")


###########################################################################################################################
####################################### ----- FUNCTIONS FOR CALCULATING RMR  ------ ########################################
##############################################################################################################################


# respR -- FUNCTION to tidy raw FIBOX 4 df and convert %O2 to O2mgL, rename columns, 
#Filter columns to keep, add column of time in a more managable format, make cols into "double" numbers where needed
F4_RMR_trial_cleanup <- function(df, lag, end, S) { # Enter the raw data frame, the amount of lag time in minutes, and the salinity in mbar
  
  # Columns to keep
  keep_columns <- c("Date", "Time", "Value", "Temp", "patm") # Create object of columns to keep
  
  # Filter columns to keep, rename "Value" to "O2air", change Temp, patm and O2air into decimal values
  F4_trial_df <- dplyr::select(df, all_of(keep_columns)) %>% 
    mutate(Time = as.character(Time), patm = as.double(patm), 
           Value = as.double(Value), Temp = as.double(Temp)) %>%
    rename(O2air = Value) %>%
    slice(1:n()-1) %>%
    mutate(TimeMinutes = seq(0, by=5/60, length.out=nrow(.))) %>% # Make minutes column for plotting
    filter(TimeMinutes >= lag,
           TimeMinutes <= end) %>% # Filter out lag time when no shark was in the chamber
    mutate(TimeMinutes = seq(0, by=5/60, length.out=nrow(.)))
  
  # Convert %DO to mg/L oxygen with convert DO function
  converted_df <- F4_trial_df %>%
    mutate(O2mgL = ((respR::convert_DO(x=100, from = "%", to = "mg/L", S = S, 
                                       t = Temp, P = (patm/1000))[[2]]) * (O2air/100)))
  
  return(converted_df)
  
  
}


###### ----- respR -- FUNCTION to tidy and convert FIBOX 3 MMR data ------ ##############
F3_RMR_trial_cleanup <- function(df, lag, end, S, patm) {
  # Tidy up the raw data set
  MO2_calculation_df <- df %>% 
    rename(Time = `time/hh:mm:ss`, LogTime = `logtime/min`,
           O2air = `oxygen/% airsatur.`,Temp = `temp/∞C`, Amp = amp, Phase = `phase/∞`)  %>% 
    dplyr::select(-X9, -ErrorMessage) %>%
    mutate(TimeMinutes = as.double(LogTime), O2air = as.double(O2air), Temp = as.double(Temp)) %>%
    filter(TimeMinutes >= lag, TimeMinutes <= end) %>%
    mutate(TimeMinutes = seq(0, by=5/60, length.out=nrow(.)))
  
  
  # Convert %DO to mg/L oxygen with convert DO function
  converted_df <- MO2_calculation_df %>%
    mutate(O2mgL = ((respR::convert_DO(x=100, from = "%", to = "mg/L", S = S, 
                                       t = Temp, P = (patm/1000))[[2]]) * (O2air/100)))
  
  return(converted_df)
  
}


#---- Function to plot basic RESTING MR O2mgL/TimeMinutes for individual sharks
basic_RMR_plot <- function(df, hs, scaleyby, scalexby) {
  
  plot <- ggplot(data = df, aes(x = TimeMinutes, y = O2mgL)) +
    geom_point() + 
    ggtitle(hs)+
    scale_y_continuous(breaks = seq(0, 12, scaleyby)) + 
    scale_x_continuous(breaks = seq(0, 2000, scalexby)) + theme_bw()
  
  return(plot)
}


# --- FUNCTION TO CALCULATE RMR ---- #
# Rolling regression with the width of the measurement window as the width of the regressions
HS_RMR_mgL_hour <- function(df, flush, close, wait, measure, meanofn, vf, vr, hs, rsqrd, 
                            low_n = 2, high_n = 10) {
  
  # Total time of one measurement cycle
  ClosePeriod <- (flush * 12L) + (close *12L)
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts <- (close*12L) - (wait * 12L)
  # Measurement period
  windowwidth <- (measure * 12L)
  mean_loop_vec <- seq(low_n, high_n, 1)
  MeanOfN <- c()
  MeanSlope <- c()
  
  # --- Rolling regression with the width of the measurement window as the width of the regressions
  roll_df <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df, width = windowwidth, 
                                       do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df <- roll_df %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df <- RMR_logtime_df %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod)) %>%
    slice_tail(n = closedEsts) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)%>%
    arrange(desc(coefs.TimeMinutes))  # Arrange slopes so lowest RMR slope is first
  
  # Take mean of lowest n slope values, loop to get a sample of what that mean is using different numbers of slopes
  for(i in mean_loop_vec) {
    
  rmr_loop_df <- rmr_rr_close_df %>%
                 head(i) # Take lowest n estimates
  
  
  # -- Calculate Standard Deviation and use it to remove outliers, if there are any -- #
  # Slopes column only
  slopes3 <- rmr_loop_df$coefs.TimeMinutes
  # Meam of estimates
  mean_unfiltered_slope <- mean(slopes3)
  # Lower and upper SD of the slope estimates
  lower_sd <- mean_unfiltered_slope - (2* (sd(slopes3))) # Define 2*SD bellow mean boundary
  upper_sd <- mean_unfiltered_slope + (2* (sd(slopes3))) # Define 2*SD above mean boundary
  # Filer out any outliers from the lowest slopes
  filtered_slopes_av3 <- rmr_loop_df %>%
    filter(coefs.TimeMinutes > lower_sd,
           coefs.TimeMinutes < upper_sd) 
  
  MeanSlopeVal <- mean(filtered_slopes_av3$coefs.TimeMinutes)
  
  # Mean low slope of the slopes kept from above
  MeanSlope <- c(MeanSlope, MeanSlopeVal)
  # Number of slopes used to calulate the mean
  MeanOfN <- c(MeanOfN, i)
  
  }
  
  HSID <- rep(hs, length(MeanSlope)) # The specimen ID will be coded as a column on the output data frame here
  # Calculate MO2 
  RMR_kghr_df <- data.frame(HSID) %>%
    mutate(MassKG = vf,
           RMR = (((vr - vf) * MeanSlope) / vf) * -60,# Negative b/c want a positive RMR estimate
           MeanOfN = MeanOfN, 
           MeanSlope = MeanSlope,
           ChamVol = vr)
          
  # Return the slope estimate for that shark
  return(RMR_kghr_df)
  
}





# --- FUNCTION TO CALCULATE RMR ** TWO DATA FRAMES ---- #
# Rolling regression with the width of the measurement window as the width of the regressions
HS_RMR_2df_mgL_hour <- function(df1, df2, df3, flush1, close1, flush2, close2,
                                wait, measure, meanofn, vf, vr, hs, rsqrd, low_n = 2, high_n = 10) {
  
  # Total time of one measurement cycle
  ClosePeriod1 <- (flush1 * 12L) + (close1 *12L)
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts1 <- (close1*12L) - (wait * 12L)
  # Total time of one measurement cycle
  ClosePeriod2 <- (flush2 * 12L) + (close2 *12L)
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts2 <- (close2*12L) - (wait * 12L)
  # Measurement period
  windowwidth <- (measure * 12L)
  mean_loop_vec <- seq(low_n, high_n, 1)
  MeanOfN <- c()
  MeanSlope <- c()
  
  # - ONE -- Rolling regression with the width of the measurement window as the width of the regressions
  roll_df1 <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df1, width = windowwidth, 
                                        do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df1 <- roll_df1 %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df1) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df1 <- RMR_logtime_df1 %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod1)) %>%
    slice_tail(n = closedEsts1) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)
  
  # -- TWO - Rolling regression with the width of the measurement window as the width of the regressions
  roll_df2 <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df2, width = windowwidth, 
                                        do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df2 <- roll_df2 %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df2) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df2 <- RMR_logtime_df2 %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod2)) %>%
    slice_tail(n = closedEsts2) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)
  
  
  # Join all the slope estimates from all the rolling regression outputs
  RR_out_slopes_all_df <- rbind(rmr_rr_close_df1, rmr_rr_close_df2) 
  RR_out_slopes_arranged_df <- RR_out_slopes_all_df %>%
    arrange(desc(coefs.TimeMinutes))  # Arrange slopes so lowest RMR slope is first
  
  ################### --------------------- ####################
  
  # Take mean of lowest n slope values, loop to get a sample of what that mean is using different numbers of slopes
  for(i in mean_loop_vec) {
    
    rmr_loop_df <- RR_out_slopes_arranged_df %>%
      head(i) # Take lowest n estimates
    
    
    # -- Calculate Standard Deviation and use it to remove outliers, if there are any -- #
    # Slopes column only
    slopes3 <- rmr_loop_df$coefs.TimeMinutes
    # Meam of estimates
    mean_unfiltered_slope <- mean(slopes3)
    # Lower and upper SD of the slope estimates
    lower_sd <- mean_unfiltered_slope - (2* (sd(slopes3))) # Define 2*SD bellow mean boundary
    upper_sd <- mean_unfiltered_slope + (2* (sd(slopes3))) # Define 2*SD above mean boundary
    # Filer out any outliers from the lowest slopes
    filtered_slopes_av3 <- rmr_loop_df %>%
      filter(coefs.TimeMinutes > lower_sd,
             coefs.TimeMinutes < upper_sd) 
    
    MeanSlopeVal <- mean(filtered_slopes_av3$coefs.TimeMinutes)
    
    # Mean low slope of the slopes kept from above
    MeanSlope <- c(MeanSlope, MeanSlopeVal)
    # Number of slopes used to calulate the mean
    MeanOfN <- c(MeanOfN, i)
    
  }
  
  HSID <- rep(hs, length(MeanSlope)) # The specimen ID will be coded as a column on the output data frame here
  # Calculate MO2 
  RMR_kghr_df <- data.frame(HSID) %>%
    mutate(MassKG = vf,
           RMR = (((vr - vf) * MeanSlope) / vf) * -60,# Negative b/c want a positive RMR estimate
           MeanOfN = MeanOfN, 
           MeanSlope = MeanSlope,
           ChamVol = vr)
  
  # Return the slope estimate for that shark
  return(RMR_kghr_df)
  
}



# --- FUNCTION TO CALCULATE RMR ** THREE DATA FRAMES ---- #
# Rolling regression with the width of the measurement window as the width of the regressions
HS_RMR_3df_mgL_hour <- function(df1, df2, df3, flush1, close1, flush2, close2, flush3, close3 ,
                                wait, measure, meanofn, vf, vr, hs, rsqrd, 
                                low_n = 2, high_n = 10) {
  
  # Total time of one measurement cycle
  ClosePeriod1 <- (flush1 * 12L) + (close1 *12L)
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts1 <- (close1*12L) - (wait * 12L)
  # Total time of one measurement cycle
  ClosePeriod2 <- (flush2 * 12L) + (close2 *12L)
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts2 <- (close2*12L) - (wait * 12L)
  # Total time of one measurement cycle
  ClosePeriod3 <- (flush3 * 12L) + (close3 *12L)
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts3 <- (close3*12L) - (wait * 12L)
  # Measurement period
  windowwidth <- (measure * 12L)
  mean_loop_vec <- seq(low_n, high_n, 1)
  MeanOfN <- c()
  MeanSlope <- c()
  
  # - ONE -- Rolling regression with the width of the measurement window as the width of the regressions
  roll_df1 <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df1, width = windowwidth, 
                                       do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df1 <- roll_df1 %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df1) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df1 <- RMR_logtime_df1 %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod1)) %>%
    slice_tail(n = closedEsts1) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)
  
  # -- TWO - Rolling regression with the width of the measurement window as the width of the regressions
  roll_df2 <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df2, width = windowwidth, 
                                       do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df2 <- roll_df2 %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df2) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df2 <- RMR_logtime_df2 %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod2)) %>%
    slice_tail(n = closedEsts2) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)
  
  # -- THREE - Rolling regression with the width of the measurement window as the width of the regressions
  roll_df3 <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df3, width = windowwidth, 
                                        do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df3 <- roll_df3 %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df3) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df3 <- RMR_logtime_df3 %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod3)) %>%
    slice_tail(n = closedEsts3) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)
  
  
  # Join all the slope estimates from all the rolling regression outputs
  RR_out_slopes_all_df <- rbind(rmr_rr_close_df1, rmr_rr_close_df2, rmr_rr_close_df3) 
  RR_out_slopes_arranged_df <- RR_out_slopes_all_df %>%
                        arrange(desc(coefs.TimeMinutes))  # Arrange slopes so lowest RMR slope is first

  ################### --------------------- ####################
  
  # Take mean of lowest n slope values, loop to get a sample of what that mean is using different numbers of slopes
  for(i in mean_loop_vec) {
    
    rmr_loop_df <- RR_out_slopes_arranged_df %>%
      head(i) # Take lowest n estimates
    
    
    # -- Calculate Standard Deviation and use it to remove outliers, if there are any -- #
    # Slopes column only
    slopes3 <- rmr_loop_df$coefs.TimeMinutes
    # Meam of estimates
    mean_unfiltered_slope <- mean(slopes3)
    # Lower and upper SD of the slope estimates
    lower_sd <- mean_unfiltered_slope - (2* (sd(slopes3))) # Define 2*SD bellow mean boundary
    upper_sd <- mean_unfiltered_slope + (2* (sd(slopes3))) # Define 2*SD above mean boundary
    # Filer out any outliers from the lowest slopes
    filtered_slopes_av3 <- rmr_loop_df %>%
      filter(coefs.TimeMinutes > lower_sd,
             coefs.TimeMinutes < upper_sd) 
    
    MeanSlopeVal <- mean(filtered_slopes_av3$coefs.TimeMinutes)
    
    # Mean low slope of the slopes kept from above
    MeanSlope <- c(MeanSlope, MeanSlopeVal)
    # Number of slopes used to calulate the mean
    MeanOfN <- c(MeanOfN, i)
    
  }
  
  HSID <- rep(hs, length(MeanSlope)) # The specimen ID will be coded as a column on the output data frame here
  # Calculate MO2 
  RMR_kghr_df <- data.frame(HSID) %>%
    mutate(MassKG = vf,
           RMR = (((vr - vf) * MeanSlope) / vf) * -60,# Negative b/c want a positive RMR estimate
           MeanOfN = MeanOfN, 
           MeanSlope = MeanSlope,
           ChamVol = vr)
  
  # Return the slope estimate for that shark
  return(RMR_kghr_df)
  
}



##############################################################################################################################
##############################################################################################################################




#### --- Pre-sets for all sharks -- #####

# Number of slopes to take a mean of for RMR estimate
SlopesForMean <- 5
# Wait period (lag period) in minutes
LagPeriod <- 3
# R-Squared for filtering slopes
HS_RSqrd <- 0.85
# Measurement window length
HSmeasure <- 7
# Lowest and highest number of measurement periods to take mean RMR from
HSlow_n = 2
HShigh_n = 10







########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 29 ----- ##

# Chamber volume
vrHS29 <- 52.5
# Fish mass (mass = volume)
vfHS29 <- 3.72
# FLush cycle length 
FlushHS29 <- 20
# Close cycle length 
CloseHS29 <- 10

# Read in intermittent flow resprometry data data file 
raw_HS29_RMR_df <- read_csv2("RawHornSharkCSV/RMR_HS29_1202_1.csv", skip = 1, # The first row of the raw CSV is junk so this skips it
                             col_types = cols(Value = col_character(), patm = col_character(), # Make these columns into characters so they can actually be read in. Convert to number in next step
                                              Temp = col_character())) 

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS29_converted_df <- F4_RMR_trial_cleanup(df = raw_HS29_RMR_df, lag = 32, end = 992, S = HS_S) 

# Plot the converted df
basic_RMR_plot(df = HS29_converted_df, hs = "HS29", scaleyby = 0.1, scalexby = (FlushHS29 +CloseHS29))

# Calculate RMR - mean of n lowest slope estimates
HS29_RMR_est_df <- HS_RMR_mgL_hour(df = HS29_converted_df, flush = FlushHS29, close = CloseHS29, wait = LagPeriod, 
                                   measure = HSmeasure, meanofn = SlopesForMean, vf = vfHS29, vr = vrHS29, 
                                   hs = "HS29", rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n)









########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 30 ----- ##

# Air Pressure - patm
HS30pressure <- 1018
# Chamber volume
vrHS30 <- 30.2
# Fish mass (mass = volume)
vfHS30 <- 1.7
# FLush cycle length 
FlushHS30 <- 20
# Close cycle length
CloseHS30 <- 10

# Read in intermittent flow resprometry data data file 
raw_HS30_RMR_df <- read_csv2("RawHornSharkCSV/RMR_191202_HS30_intermit1.csv", skip = 57, # First 57 rows are junk, skip
                             col_types = cols(`oxygen/% airsatur.` = col_character(),  
                                              `temp/∞C` = col_character()))

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS30_converted_df <- F3_RMR_trial_cleanup(df = raw_HS30_RMR_df, lag = 0, end = 990, patm =  HS30pressure, S = HS_S) 

# Plot the converted df
basic_RMR_plot(df = HS30_converted_df, hs = "HS30", scaleyby = 0.1, scalexby = (FlushHS30 +CloseHS30))

# Calculate RMR - mean of n lowest slope estimates
HS30_RMR_est_df <- HS_RMR_mgL_hour(df = HS30_converted_df, flush = FlushHS30, close = CloseHS30, wait = LagPeriod, 
                                   measure = HSmeasure, meanofn = SlopesForMean, vf = vfHS30, 
                                   vr = vrHS30, hs = "HS30", rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n)



