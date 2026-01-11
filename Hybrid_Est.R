# This OpenSource R program consists of four functions that comprise the 
# Hybrid_Est.R package that estimates CART bias trees and bias-adjusted 
# impacts under multisite hybrid impact designs. It is described in the package 
# documentation pdf file. It was created by Peter Z. Schochet, Ph.D. 
# (pzschochet@gmail.com) and was posted to github in January 2026.

library(stringr)
library(listr)
library(dplyr)
library(data.tree)
library(ivreg)
library(survey)
library(lmtest)
library(sandwich)
library(clubSandwich)

#####
# This is the bias_tree function that has several sub functions
#####

bias_tree <- function(rct_dat_df, qed_dat_df, rct_qed_site, yvar, xvars_cart, 
                      site_id, t_c, got_treat, ipw_wgt, cace_itt_est, out_cart,
                      minsplit_site = 6, minsplit_tot = 20, 
                      max_depth = 30, holdout = 1, holdout_seed = 42) {
  
  # Initialize err text and counter
  err_count <<- 0
  err_num <- seq(1:100)
  err_txt <- character(100)
  err_df  <- data.frame(err_num,err_txt)
  
  cace_itt     <<- cace_itt_est   # Makes this global
  hold_seed    <<- holdout_seed
  hold_out     <<- holdout
  
  # Create variables to keep for analysis for rct and qed datasets
  # These steps deparse the inputs to make them character strings
  yv     <- deparse(substitute(yvar))
  #yv    <- gsub('[(\")(\\)]', '', yv_temp) # don't need this anymore
  sitev  <- deparse(substitute(site_id))
  tc     <- deparse(substitute(t_c))
  partic <- deparse(substitute(got_treat))
  ipw_wt <- deparse(substitute(ipw_wgt))
  rct_qed_ind <- deparse(substitute(rct_qed_site)) 
  out_txt     <- deparse(substitute(out_cart))   
  
  if (out_txt == "") {
    out_txt <- c('cart_results.txt')
  }
  
  # Define sink function to output results and also writes to the console
  #sink(file = out_txt, split = TRUE)
  
  # Check if names exist
  # Create rct and qed datasets
  if (yv %in% names(rct_dat_df)) {
    rct_dat_df$yv <- rct_dat_df[,yv]
  } else {
    rct_dat_df$yv <- NA
  }
  
  if (sitev %in% names(rct_dat_df)) {
    rct_dat_df$sitev <- rct_dat_df[,sitev]
  } else {
    rct_dat_df$sitev <- NA
  }
  
  if (tc %in% names(rct_dat_df)) {
    rct_dat_df$tc <- rct_dat_df[,tc]
  } else {
    rct_dat_df$tc <- NA
  }
  
  if (partic %in% names(rct_dat_df)) {
    rct_dat_df$partic <- rct_dat_df[,partic]
    bad_particz <- 0
  } else {
    rct_dat_df$partic <- NA
    bad_particz <- 1
  }
  
  # For ITT analysis, set participation variable to 1 for RCT Ts
  # and 0 for RCT Cs
  if (cace_itt == 0) {
    rct_dat_df$partic <- ifelse(rct_dat_df$tc == 1, 1, 0) 
  } else if ((cace_itt == 1) & (bad_particz == 0)) {
    rct_dat_df$partic <- rct_dat_df[,partic]
  } else {
    rct_dat_df$partic <- NA
  }
  
  rct_dat_df$ipw_wt <- 1   
  
  # Process xvars which is a vector with + signs
  xva <- deparse(substitute(xvars_cart))
  xvb <- as.character(gsub('[+]', '', xva))
  # Split the xvb string by spaces
  xvc <- strsplit(xvb, split = "\\s+")
  xv  <- xvc[[1]] 
  
  for (ix in 1:length(xv)) {
    if (xv[ix] %in% names(rct_dat_df)) {
      junk <- 1
    } else {
      xxv <- as.character(xv[ix])
      rct_dat_df[,xxv] <- NA
    }
  }
  
  rcty  <- rct_dat_df[,c("yv","sitev","tc","partic","ipw_wt")]
  rctx  <- data.frame(rct_dat_df[,xv])
  
  rct_dat <- data.frame(rcty,rctx)
  rct_dat$rct_qed_ind <- 1
  
  if (yv %in% names(qed_dat_df)) {
    qed_dat_df$yv <- qed_dat_df[,yv]
  } else {
    qed_dat_df$yv <- NA
  }
  
  if (sitev %in% names(qed_dat_df)) {
    qed_dat_df$sitev <- qed_dat_df[,sitev]
  } else {
    qed_dat_df$sitev <- NA
  }
  
  if (tc %in% names(qed_dat_df)) {
    qed_dat_df$tc <- qed_dat_df[,tc]
  } else {
    qed_dat_df$tc <- NA
  }
  
  if (rct_qed_ind %in% names(qed_dat_df)) {
    qed_dat_df$rct_qed_ind <- qed_dat_df[,rct_qed_ind]
    bad_rqi <- 0
  } else {
    qed_dat_df$rct_qed_ind <- NA
    bad_rqi <- 1
  }
  
  if (ipw_wt %in% names(qed_dat_df)) {
    qed_dat_df$ipw_wt <- qed_dat_df[,ipw_wt]
  } else {
    qed_dat_df$ipw_wt <- NA
  }
  
  # Check for a few errors in QED dataset for tc and rct_qed_ind 
  goodq_tc <- 1
  if (any(is.na(qed_dat_df$tc))) {
    goodq_tc <- 0
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing t_c values in qed_dat_df")
  } else {
    badq_tc <- ifelse(((qed_dat_df$tc == 0) | (qed_dat_df$tc == 1)), 0, 1)
    if (any(badq_tc == 1)) {
      goodq_tc <- 0
      err_count <<- err_count+1
      err_df[err_count,2] <- c("t_c values in qed_dat_df must all be 0 or 1")
    }
  }
  
  if (any(is.na(qed_dat_df$rct_qed_ind))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing rct_qed_site values in qed_dat_df")
  } else {
    badq_rct_qed <- ifelse(((qed_dat_df$rct_qed_ind == 0) | (qed_dat_df$rct_qed_ind == 1)), 0, 1)
    if (any(badq_rct_qed == 1)) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("rct_qed_site values in qed_dat_df must all be 0 or 1")
    }
  }
  
  # For QED, set participation variable to 1 for Ts and 0 for Cs
  qed_dat_df$partic <- ifelse(qed_dat_df$tc == 1,1,0)
  
  qedy  <- qed_dat_df[,c("yv","sitev","tc","partic","ipw_wt","rct_qed_ind")]
  
  for (ix in 1:length(xv)) {
    if (xv[ix] %in% names(qed_dat_df)) {
      junk <- 1
    } else {
      xxv <- as.character(xv[ix])
      qed_dat_df[,xxv] <- NA
    }
  }
  
  qedx  <- data.frame(qed_dat_df[,xv])
  
  # Only use the qed sites that are in the rct_qed sample  
  qed_dat_temp <- data.frame(qedy,qedx)
  qed_dat      <- qed_dat_temp[qed_dat_temp$rct_qed_ind == 1,]
  
  ###
  # Check for errors
  ###
  
  # Sitev 
  if (any(is.na(rct_dat$sitev))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing siteid values in rct_dat_df")
  }
  
  if (any(is.na(qed_dat$sitev)) & (bad_rqi == 0)) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing siteid values in qed_dat_df")
  }
  
  # yv and xv
  if (any(is.na(rct_dat$yv))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing yvar values in rct_dat_df")
  }
  
  if (any(is.na(rct_dat[,xv]))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing xvars_cart values in rct_dat_df")
  }
  
  if (any(is.na(qed_dat$yv)) & (bad_rqi == 0)) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing yvar values in qed_dat_df")
  }
  
  if (any(is.na(qed_dat[,xv])) & (bad_rqi == 0)) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing xvars_cart values in qed_dat_df")
  }
  
  # tc
  if (any(is.na(rct_dat$tc))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing t_c values in rct_dat_df")
  } else {
    badr_tc <- ifelse(((rct_dat$tc == 0) | (rct_dat$tc == 1)), 0, 1)
    if (any(badr_tc == 1)) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("t_c values in rct_dat_df must all be 0 or 1")
    }
  }

  # cace_itt and partic
  if ((cace_itt == 0) | (cace_itt == 1)) {
    bad_cace_itt <- 0
  } else {
    bad_cace_itt <- 1
    err_count <<- err_count+1
    err_df[err_count,2] <- c("cace_itt value must be 0 or 1")
  }
  
  if (cace_itt == 1) {
    if (any(is.na(rct_dat$partic))) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing got_treat values in rct_dat_df")
    } else {
      badr_partic <- ifelse(((rct_dat$partic == 0) | (rct_dat$partic == 1)), 0, 1)
      if (any(badr_partic == 1)) {
        err_count <<- err_count+1
        err_df[err_count,2] <- c("got_treat values in rct_dat_df must all be 0 or 1")
      }
    }
  }
  
  # ipw_wt
  if (goodq_tc <- 1) {
    if (any(is.na(qed_dat$ipw_wt)) & (bad_rqi == 0)) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing ipw_wgt values in qed_dat_df")
    } else if (bad_rqi == 0) {
      badq_wt <- ifelse(((qed_dat$ipw_wt < 0)), 1, 0)
      if (any(badq_wt == 1)) {
        err_count <<- err_count+1
        err_df[err_count,2] <- c("ipw_wgt values in qed_dat_df must all be nonnegative")
      }
    }
  }
  
  # holdout, minsplit_site, minsplit_tot, max_depth
  if ((holdout == 0) | (holdout == 1)) {
    bad_hold <- 0
  } else {
    bad_hold <- 1
    err_count <<- err_count+1
    err_df[err_count,2] <- c("holdout value must be 0 or 1")
  }
  
  if (minsplit_site < 4) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("minsplit_site value must be at least 4")
  }
  
  if (minsplit_tot < 10) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("minsplit_tot value must be at least 10")
  }
  
  if (max_depth < 2) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("max_depth value must be at least 2")
  }
  
  # Return function if there are any errors
  
  if (err_count >= 1) {
    err_dfg <- err_df[1:err_count,2]
    cat("\n")
    print("ERRORS TO FIX IN BIAS_TREE FUNCTION")
    print(err_dfg)
    return(err_dfg)
  }
  
  # Define datasets for RCT control nonparticipants, RCT treatment nonparticipants,
  # and QED comparisons
  qed_datc <- qed_dat[qed_dat$tc == 0,]
  rct_datc <- rct_dat[rct_dat$tc == 0,]
  rct_datt <- rct_dat[rct_dat$tc == 1,]
  rct_datc_np <- rct_datc[rct_datc$partic == 0,]
  rct_datt_p  <- rct_datt[rct_datt$partic == 1,]
  
  # Check for missing site IDs and check that the same ones are on
  # the rct and qed datasets. Do separately for the cace and itt analyses
  if (cace_itt == 0) {
    
    site_qed <- data.frame(table(qed_datc$sitev,useNA= "always"))
    site_qed <- site_qed %>%
      rename(site_var = Var1,
             qed_freq = Freq)
    site_qed <- site_qed[,c("site_var","qed_freq")]
    
    site_rct <- data.frame(table(rct_datc$sitev,useNA= "always"))
    site_rct <- site_rct %>%
      rename(site_var = Var1,
             rct_freq = Freq)
    site_rct <- site_rct[,c("site_var","rct_freq")]
    
    sitev_m <- merge(site_rct, site_qed, by = "site_var",
                     all = TRUE)
    
    # Replace NA with 0 in site frequencies
    sitev_m$rct_freq <- ifelse(is.na(sitev_m$rct_freq),0,sitev_m$rct_freq)
    sitev_m$qed_freq <- ifelse(is.na(sitev_m$qed_freq),0,sitev_m$qed_freq)
    
    sitev_mz1 <- sitev_m
    
    cat("\n")
    pr1 <- "TABLE 1. BIAS_TREE FUNCTION: RCT AND QED SITE SAMPLE SIZES IN THE"
    pr2 <- "RCT-QED SITES FOR THE FULL SAMPLE ITT ANALYSIS"
    cat(pr1,"\n",pr2,"\n") 
    print(sitev_m)
    
    # Check for missing site IDs
    miss_site_var <- sitev_m[is.na(sitev_m$site_var),]
    
    msite_ind <- 0
    if ((miss_site_var$rct_freq > 0) | (miss_site_var$qed_freq > 0)) {
      msite_ind <- 1
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing Site IDs")
    }
    
    # Check for mismatched RCT and QED site IDs
    nonmiss_site_var <- sitev_m[!is.na(sitev_m$site_var),]
    
    rct_zero_freq <- any(nonmiss_site_var$rct_freq == 0)
    qed_zero_freq <- any(nonmiss_site_var$qed_freq == 0)
    
    if ((rct_zero_freq == TRUE) | (qed_zero_freq == TRUE)) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Some RCT and QED Site IDs in the RCT-QED sites do not match")
      #print(err_df[1:err_count,],right='F')
    }
    
    # Check if there are enough RCT and QED sample members in each site and in total
    # using the minsplit_site and minsplit_tot inputs
    rct_site_freq <- any(nonmiss_site_var$rct_freq < minsplit_site)
    qed_site_freq <- any(nonmiss_site_var$qed_freq < minsplit_site)
    
    rct_site_tot  <- sum(nonmiss_site_var$rct_freq)
    qed_site_tot  <- sum(nonmiss_site_var$qed_freq)
    
    if (msite_ind == 0) {
      if ((rct_site_freq == TRUE) | (qed_site_freq == TRUE)) {
        err_count <<- err_count+1
        err_df[err_count,2] <- paste(c("Some RCT-QED sites have fewer RCT control or QED comparison sample members than specified in the minsplit_site input of"),minsplit_site)
      }
      
      if ((rct_site_tot < minsplit_tot) | (qed_site_tot < minsplit_tot)) {
        err_count <<- err_count+1
        err_df[err_count,2] <- paste(c("Total RCT control or QED comparison sample size is less than the minsplit_tot input of"), minsplit_tot)
      }
    }
    
  } else if (cace_itt == 1) {
    
    site_qedc <- data.frame(table(qed_datc$sitev,useNA= "always"))
    site_qedc <- site_qedc %>%
      rename(site_var = Var1,
             qed_freqc = Freq)
    site_qedc <- site_qedc[,c("site_var","qed_freqc")]
    #site_qedc <- site_qedc[,c("site_var","rct_freq")]
    
    site_rctc <- data.frame(table(rct_datc$sitev,useNA= "always"))
    site_rctc <- site_rctc %>%
      rename(site_var = Var1,
             rct_freqc = Freq)
    site_rctc <- site_rctc[,c("site_var","rct_freqc")]
    #site_rctc <- site_rctc[,c("site_var","rct_freq")]
    
    site_rctt <- data.frame(table(rct_datt$sitev,useNA= "always"))
    site_rctt <- site_rctt %>%
      rename(site_var = Var1,
             rct_freqt = Freq)
    site_rctt <- site_rctt[,c("site_var","rct_freqt")]
    
    site_rctc_np <- data.frame(table(rct_datc_np$sitev,useNA= "always"))
    site_rctc_np <- site_rctc_np %>%
      rename(site_var = Var1,
             rct_freqc_np = Freq)
    site_rctc_np <- site_rctc_np[,c("site_var","rct_freqc_np")]
    
    site_rctt_p <- data.frame(table(rct_datt_p$sitev,useNA= "always"))
    site_rctt_p <- site_rctt_p %>%
      rename(site_var = Var1,
             rct_freqt_p = Freq)
    site_rctt_p <- site_rctt_p[,c("site_var","rct_freqt_p")]
    
    sitev_m1 <- merge(site_rctc, site_rctt, , by = "site_var", all = TRUE)
    sitev_m2 <- merge(sitev_m1, site_qedc, by = "site_var", all = TRUE)
    sitev_m3 <- merge(sitev_m2, site_rctc_np, by = "site_var", all = TRUE)
    sitev_m  <- merge(sitev_m3, site_rctt_p, by = "site_var", all = TRUE)
    
    # Replace NA with 0 in site frequencies
    sitev_m$rct_freqc <- ifelse(is.na(sitev_m$rct_freqc),0,sitev_m$rct_freqc)
    sitev_m$rct_freqt <- ifelse(is.na(sitev_m$rct_freqt),0,sitev_m$rct_freqt)
    sitev_m$qed_freqc <- ifelse(is.na(sitev_m$qed_freqc),0,sitev_m$qed_freqc)
    sitev_m$rct_freqc_np <- ifelse(is.na(sitev_m$rct_freqc_np),0,
                                   sitev_m$rct_freqc_np)
    sitev_m$rct_freqt_p  <- ifelse(is.na(sitev_m$rct_freqt_p),0,
                                   sitev_m$rct_freqt_p)
    
    sitev_mz1 <- sitev_m
    
    cat("\n")
    pr1 <- "TABLE 1. BIAS_TREE FUNCTION: RCT CONTROL, RCT TREATMENT, AND QED COMPARISON"
    pr2 <- "SITE SAMPLE SIZES IN THE RCT-QED SITES FOR THE FULL SAMPLE CACE ANALYSIS"
    cat(pr1,"\n",pr2,"\n") 
    
    print(sitev_m[,c("site_var","rct_freqc","rct_freqc_np",
                     "rct_freqt","rct_freqt_p","qed_freqc")])
    fn1 <- "Notes: Figures in order show samples for: (i) the full RCT control group,"
    fn2 <- "(ii) RCT control nonparticipants, (iii) the full RCT treatment group,"
    fn3 <- "(iii) RCT treatment participants, and (iv) the full QED comparison group."
    cat(fn1,"\n",fn2,"\n",fn3,"\n")
    
    # Check for missing site IDs - do not do for rct_np datasets since implied
    # using the larger full datasets
    miss_site_var <- sitev_m[is.na(sitev_m$site_var),]
    
    msite_ind <- 0
    if ((miss_site_var$rct_freqc > 0) | (miss_site_var$rct_freqt > 0) |
        (miss_site_var$qed_freqc > 0)) {
      msite_ind <- 1
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing Site IDs")
    }
    
    # Check for mismatched RCT and QED site IDs
    nonmiss_site_var <- sitev_m[!is.na(sitev_m$site_var),]
    
    rct_zero_freqc <- any(nonmiss_site_var$rct_freqc == 0)
    rct_zero_freqt <- any(nonmiss_site_var$rct_freqt == 0)
    qed_zero_freqc <- any(nonmiss_site_var$qed_freqc == 0)
    rct_zero_freqc_np <- any(nonmiss_site_var$rct_freqc_np == 0)
    rct_zero_freqt_np <- any(nonmiss_site_var$rct_freqt_np == 0)
    
    if ((rct_zero_freqc == TRUE) | (rct_zero_freqt == TRUE) | 
        (qed_zero_freqc == TRUE)) {
      err_count <<- err_count+1
      cat("\n")
      err_df[err_count,2] <- c("Some RCT and QED Site IDs in the RCT-QED sites do not match")
      #print(err_df[1:err_count,],right='F')
    }
    
    # Check if there are enough RCT and QED sample members in each site and in total
    # using the minsplit_site and minsplit_tot inputs
    rct_site_freqc_np <- any(nonmiss_site_var$rct_freqc_np < minsplit_site)
    rct_site_freqt_p  <- any(nonmiss_site_var$rct_freqt_p < minsplit_site)
    qed_site_freqc    <- any(nonmiss_site_var$qed_freqc < minsplit_site)
    
    rct_site_totc_np  <- sum(nonmiss_site_var$rct_freqc_np)
    rct_site_tott_p   <- sum(nonmiss_site_var$rct_freqt_p)
    qed_site_totc     <- sum(nonmiss_site_var$qed_freqc)
    
    if (msite_ind == 0) {
      if ((rct_site_freqc_np == TRUE) | (rct_site_freqt_p == TRUE) | 
          (qed_site_freqc == TRUE)) {
        err_count <<- err_count+1
        err_df[err_count,2] <- paste(c("Some RCT-QED sites have fewer required sample members than specified in the minsplit_site input of"),minsplit_site)
      }
      
      if ((rct_site_totc_np < minsplit_tot) | (rct_site_tott_p < minsplit_tot) | 
          (qed_site_totc < minsplit_tot)) {
        err_count <<- err_count+1
        err_df[err_count,2] <- paste(c("Total required sample size across sites is less than the minsplit_tot input of"), minsplit_tot)
      }
    }
    
    # Check if the participation rate is greater for RCT treatments than controls
    # in all sites
    
    if (err_count == 0) {
    part_ratet <- data.frame(aggregate(partic ~ sitev, data = rct_datt, FUN = mean))
    part_ratec <- data.frame(aggregate(partic ~ sitev, data = rct_datc, FUN = mean))
    
    partic_df <- merge(part_ratet, part_ratec, by = "sitev")
    colnames(partic_df) <- c("site_var", "part_ratet", "part_ratec")
    
    partic_df$compl_rate <- partic_df$part_ratet - partic_df$part_ratec
    
    compl_rate_tot <<- mean(rct_datt$partic) - mean(rct_datc$partic) # makes global
    
    pc_df1 <- partic_df
    pc_df1$part_ratet <- format(round(pc_df1$part_ratet,4), nsmall=4)
    pc_df1$part_ratec <- format(round(pc_df1$part_ratec,4), nsmall=4)
    pc_df1$compl_rate <- format(round(pc_df1$compl_rate,4), nsmall=4)
    
    pc_df1zz <- pc_df1
    
    cat("\n")
    pr1 <- "TABLE 2. BIAS_TREE FUNCTION: RCT TREATMENT AND CONTROL GROUP PARTICIPATION RATES"
    pr2 <- "AND THE COMPLIANCE RATE FOR THE FULL SAMPLE CACE ANALYSIS"
    cat(pr1,"\n",pr2,"\n") 
    print(pc_df1)
    
    any_bad_compl <- any(partic_df$compl_rate <= 0)
    
    if (any_bad_compl == TRUE) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("In some sites, the participation rate is smaller for RCT treatments than RCT controls")
    }
    } # if err_count==0
    
  } #end cace_itt
  
  # Create site ID that goes from 1 to se. The cur_group_id() variable is an
  # internal dplr function that retains the current group counter
  
  rct_datca <- rct_datc %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
    rct_datca <- data.frame(rct_datca)
  
  qed_datca <- qed_datc %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
    qed_datca <- data.frame(qed_datca)
  
  rct_datta <- rct_datt %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
    rct_datta <- data.frame(rct_datta)
  
  rct_datc_npa <- rct_datca[rct_datca$partic == 0,]
  rct_datt_pa  <- rct_datta[rct_datta$partic == 1,]
  
  #print(table(rct_datca$siteID,rct_datca$sitev))
  #print(table(qed_datca$siteID,qed_datca$sitev))
  
  # Create training and holdout samples by site
  
  if (holdout == 1) {
    
    # This is better but is not reproducible with same seed if site order is changed
    #train_rct_datc <- rct_datca %>%
    #  group_by(siteID) %>%
    #  sample_frac(size = 0.7) # Samples 70% of rows from each group for training
    
    #hold_rct_datc <- anti_join(rct_datca, train_rct_datc, 
    #                           by = names(rct_datca))
    
    set.seed(hold_seed)
    
    # Add obs number to sort back to original order later
    rct_datca$persID <- seq_len(nrow(rct_datca))
    rct_datta$persID <- seq_len(nrow(rct_datta))
    qed_datca$persID <- seq_len(nrow(qed_datca))
    
    rct_dat_fullc <- rct_datca
    rct_dat_fullt <- rct_datta
    qed_dat_fullc <- qed_datca
    
    countj <- 1
    for (j in unique(rct_datca$siteID)) {
      
      #print(j)
      
      # RCT controls
      rct_datcaj <- rct_datca[rct_datca$siteID == j,]
      
      train_rct_sizecj <- floor(0.7 * nrow(rct_datcaj))
      # Generate random indices for the training set
      train_rct_indcj  <- sample(seq_len(nrow(rct_datcaj)), size = train_rct_sizecj)
      train_rct_datcj  <- rct_datcaj[train_rct_indcj,]
      hold_rct_datcj   <- rct_datcaj[-train_rct_indcj,]
      
      # QED comparisons
      qed_datcaj <- qed_datca[qed_datca$siteID == j,]
      
      train_qed_sizecj <- floor(0.7 * nrow(qed_datcaj))
      # Generate random indices for the training set
      train_qed_indcj  <- sample(seq_len(nrow(qed_datcaj)), size = train_qed_sizecj)
      train_qed_datcj  <- qed_datcaj[train_qed_indcj,]
      hold_qed_datcj   <- qed_datcaj[-train_qed_indcj,]
      
      # RCT treatments
      rct_dattaj <- rct_datta[rct_datta$siteID == j,]
      
      train_rct_sizetj <- floor(0.7 * nrow(rct_dattaj))
      # Generate random indices for the training set
      train_rct_indtj  <- sample(seq_len(nrow(rct_dattaj)), size = train_rct_sizetj)
      train_rct_dattj  <- rct_dattaj[train_rct_indtj,]
      hold_rct_dattj   <- rct_dattaj[-train_rct_indtj,]
      
      if (countj == 1) {
        train_rct_datc <- train_rct_datcj
        train_qed_datc <- train_qed_datcj
        train_rct_datt <- train_rct_dattj
        
        hold_rct_datc1 <- hold_rct_datcj
        hold_qed_datc1 <- hold_qed_datcj
        hold_rct_datt1 <- hold_rct_dattj
      } else if (countj > 1) {
        train_rct_datc <- rbind(train_rct_datc,train_rct_datcj)
        train_qed_datc <- rbind(train_qed_datc,train_qed_datcj)
        train_rct_datt <- rbind(train_rct_datt,train_rct_dattj)
        
        hold_rct_datc1 <- rbind(hold_rct_datc1,hold_rct_datcj)
        hold_qed_datc1 <- rbind(hold_qed_datc1,hold_qed_datcj)
        hold_rct_datt1 <- rbind(hold_rct_datt1,hold_rct_dattj)
      }
      
      countj <- countj + 1
    }
    
    # Sort back to original order
    
    rct_datc <- train_rct_datc[order(train_rct_datc$persID),]
    qed_datc <- train_qed_datc[order(train_qed_datc$persID),]
    rct_datt <- train_rct_datt[order(train_rct_datt$persID),]
    
    hold_rct_datc <- hold_rct_datc1[order(hold_rct_datc1$persID),]
    hold_qed_datc <- hold_qed_datc1[order(hold_qed_datc1$persID),]
    hold_rct_datt <- hold_rct_datt1[order(hold_rct_datt1$persID),]
    
    rct_datc_np <- rct_datc[rct_datc$partic == 0,]
    rct_datt_p  <- rct_datt[rct_datt$partic == 1,]
    
    # Calculate SE = the number of RCT-QED sites
    se <<- max(rct_datc$siteID)  # Makes this global
    
    if (cace_itt == 0) {
    
    # Need to check again about enough sample sizes by site in training sample
    site_rct <- data.frame(table(rct_datc$sitev))
    site_rct <- site_rct %>%
      rename(site_var = Var1,
             rct_freq = Freq)
    site_rct <- site_rct[,c("site_var","rct_freq")]
    
    site_qed <- data.frame(table(qed_datc$sitev))
    site_qed <- site_qed %>%
      rename(site_var = Var1,
             qed_freq = Freq)
    site_qed <- site_qed[,c("site_var","qed_freq")]
    
    sitev_m <- merge(site_rct, site_qed, by = "site_var",
                     all = TRUE)
    
    # Replace NA with 0 in site frequencies
    sitev_m$rct_freq <- ifelse(is.na(sitev_m$rct_freq),0,sitev_m$rct_freq)
    sitev_m$qed_freq <- ifelse(is.na(sitev_m$qed_freq),0,sitev_m$qed_freq)
    
    sitev_mz2 <- sitev_m
    
    cat("\n")
    pr1 <- "TABLE 1A. BIAS_TREE FUNCTION: RCT AND QED SITE SAMPLE SIZES IN THE" 
    pr2 <- "RCT-QED SITES FOR THE 70 PERCENT TRAINING SAMPLE ITT ANALYSIS"
    cat(pr1,"\n",pr2,"\n")
    print(sitev_m)
    
    # Check if there are enough RCT and QED sample members in each site and in total
    # using the minsplit_site and minsplit_tot inputs
    rct_site_freq <- any(sitev_m$rct_freq < minsplit_site)
    qed_site_freq <- any(sitev_m$qed_freq < minsplit_site)
    
    rct_site_tot  <- sum(sitev_m$rct_freq)
    qed_site_tot  <- sum(sitev_m$qed_freq)
    
    if (msite_ind == 0) {
      if ((rct_site_freq == TRUE) | (qed_site_freq == TRUE)) {
        err_count <<- err_count+1
        err_df[err_count,2] <- paste(c("Some RCT-QED sites have fewer RCT control or QED comparison sample members in the training sample than specified in the minsplit_site input of"),minsplit_site)
      }
      
      if ((rct_site_tot < minsplit_tot) | (qed_site_tot < minsplit_tot)) {
        err_count <<- err_count+1
        err_df[err_count,2] <- paste(c("The total RCT control or QED comparison sample size in the training sample is less than the minsplit_tot input of"), minsplit_tot)
      }
    }
    
    } else if (cace_itt == 1) {
      
      # Need to check again about enough sample sizes by site in training sample
      site_qedc <- data.frame(table(qed_datc$sitev))
      site_qedc <- site_qedc %>%
        rename(site_var = Var1,
               qed_freqc = Freq)
      site_qedc <- site_qedc[,c("site_var","qed_freqc")]
      
      site_rctc <- data.frame(table(rct_datc$sitev))
      site_rctc <- site_rctc %>%
        rename(site_var = Var1,
               rct_freqc = Freq)
      site_rctc <- site_rctc[,c("site_var","rct_freqc")]
      
      site_rctt <- data.frame(table(rct_datt$sitev))
      site_rctt <- site_rctt %>%
        rename(site_var = Var1,
               rct_freqt = Freq)
      site_rctt <- site_rctt[,c("site_var","rct_freqt")]
      
      site_rctc_np <- data.frame(table(rct_datc_np$sitev))
      site_rctc_np <- site_rctc_np %>%
        rename(site_var = Var1,
               rct_freqc_np = Freq)
      site_rctc_np <- site_rctc_np[,c("site_var","rct_freqc_np")]
      
      site_rctt_p <- data.frame(table(rct_datt_p$sitev))
      site_rctt_p <- site_rctt_p %>%
        rename(site_var = Var1,
               rct_freqt_p = Freq)
      site_rctt_p <- site_rctt_p[,c("site_var","rct_freqt_p")]
      
      sitev_m1 <- merge(site_rctc, site_rctt, , by = "site_var", all = TRUE)
      sitev_m2 <- merge(sitev_m1, site_qedc, by = "site_var", all = TRUE)
      sitev_m3 <- merge(sitev_m2, site_rctc_np, by = "site_var", all = TRUE)
      sitev_m  <- merge(sitev_m3, site_rctt_p, by = "site_var", all = TRUE)
      
      # Replace NA with 0 in site frequencies
      sitev_m$rct_freqc <- ifelse(is.na(sitev_m$rct_freqc),0,sitev_m$rct_freqc)
      sitev_m$rct_freqt <- ifelse(is.na(sitev_m$rct_freqt),0,sitev_m$rct_freqt)
      sitev_m$qed_freqc <- ifelse(is.na(sitev_m$qed_freqc),0,sitev_m$qed_freqc)
      sitev_m$rct_freqc_np <- ifelse(is.na(sitev_m$rct_freqc_np),0,
                                     sitev_m$rct_freqc_np)
      sitev_m$rct_freqt_p  <- ifelse(is.na(sitev_m$rct_freqt_p),0,
                                     sitev_m$rct_freqt_p)
      
      sitev_mz2 <- sitev_m
      
      cat("\n")
      pr1 <- "TABLE 1A. BIAS_TREE FUNCTION: RCT AND QED SITE SAMPLE SIZES IN THE" 
      pr2 <- "RCT-QED SITES FOR THE 70 PERCENT TRAINING SAMPLE CACE ANALYSIS"
      cat(pr1,"\n",pr2,"\n")
      print(sitev_m[,c("site_var","rct_freqc","rct_freqc_np",
                       "rct_freqt","rct_freqt_p","qed_freqc")])
      fn1 <- "Notes: Figures in order show samples for: (i) the RCT control group,"
      fn2 <- "(ii) RCT control nonparticipants, (iii) the RCT treatment group,"
      fn3 <- "(iii) RCT treatment participants, and (iv) the QED comparison group."
      cat(fn1,"\n",fn2,"\n",fn3,"\n")
      
      # Check if there are enough RCT and QED sample members in each site and in total
      # using the minsplit_site and minsplit_tot inputs
      rct_site_freqc_np <- any(sitev_m$rct_freqc_np < minsplit_site)
      rct_site_freqt_p  <- any(sitev_m$rct_freqt_p < minsplit_site)
      qed_site_freqc    <- any(sitev_m$qed_freqc < minsplit_site)
      
      rct_site_totc_np  <- sum(sitev_m$rct_freqc_np)
      rct_site_tott_p   <- sum(sitev_m$rct_freqt_p)
      qed_site_totc     <- sum(sitev_m$qed_freqc)
      
      if (msite_ind == 0) {
        if ((rct_site_freqc_np == TRUE) | (rct_site_freqt_p == TRUE) | 
            (qed_site_freqc == TRUE)) {
          err_count <<- err_count+1
          err_df[err_count,2] <- paste(c("Some RCT-QED sites have fewer required sample members in the training sample than specified in the minsplit_site input of"),minsplit_site)
        }
        
        if ((rct_site_totc_np < minsplit_tot) | (rct_site_tott_p < minsplit_tot) | 
            (qed_site_totc < minsplit_tot)) {
          err_count <<- err_count+1
          err_df[err_count,2] <- paste(c("Total required sample size across sites in the training sample is less than the minsplit_tot input of"), minsplit_tot)
        }
      }
      
      # Check if the participation rate is greater for RCT treatments than controls
      # in all sites
      if (err_count == 0) {
      
      part_ratet <- data.frame(aggregate(partic ~ sitev, data = rct_datt, FUN = mean))
      part_ratec <- data.frame(aggregate(partic ~ sitev, data = rct_datc, FUN = mean))
      
      partic_df <- merge(part_ratet, part_ratec, by = "sitev")
      colnames(partic_df) <- c("site_var", "part_ratet", "part_ratec")
      
      partic_df$compl_rate <- partic_df$part_ratet - partic_df$part_ratec
      
      pc_df1 <- partic_df
      pc_df1$part_ratet <- format(round(pc_df1$part_ratet,4), nsmall=4)
      pc_df1$part_ratec <- format(round(pc_df1$part_ratec,4), nsmall=4)
      pc_df1$compl_rate <- format(round(pc_df1$compl_rate,4), nsmall=4)
      
      pc_df1zz2 <- pc_df1
      
      cat("\n")
      pr1 <- "TABLE 2A. BIAS_TREE FUNCTION: RCT TREATMENT AND CONTROL GROUP PARTICIPATION RATES"
      pr2 <- "AND THE COMPLIANCE RATE FOR THE 70 PERCENT TRAINING SAMPLE CACE ANALYSIS"
      cat(pr1,"\n",pr2,"\n")
      print(pc_df1)
      
      any_bad_compl <- any(partic_df$compl_rate <= 0)
      
      if (any_bad_compl == TRUE) {
        err_count <<- err_count+1
        err_df[err_count,2] <- c("Some sites in the training sample where the participation rate is smaller for RCT treatments than controls")
      }
      } # err_count==0
    
    } # end if cace_itt
    
  } else if (holdout == 0) {
    rct_datc <- rct_datca
    qed_datc <- qed_datca
    rct_datt <- rct_datta
    
    rct_dat_fullc <- rct_datc
    rct_dat_fullt <- rct_datt
    qed_dat_fullc <- qed_datc
    
    rct_datc_np <- rct_datc_npa 
    rct_datt_p  <- rct_datt_pa
    
  } # end holdout
  
  # Define dat_noder and dat_nodeq
  dat_noder  <- rct_datc
  dat_nodeq  <- qed_datc
  dat_nodert <- rct_datt
  
  # Return function if there are any errors
  
  if (err_count >= 1) {
    err_dfg <- err_df[1:err_count,2]
    cat("\n")
    print("ERRORS TO FIX IN BIAS_TREE FUNCTION")
    print(err_dfg)
    return(err_dfg)
  }
  
  custom_rpart <- function() {
    
  ####
  # Start tree building analysis
  ####
  
  # Prepare the response variable (target) and predictors (features)
  yr  <- rct_datc$yv  # Target variable (dependent)
  yq  <- qed_datc$yv  # Target variable (dependent)
  yrt <- rct_datt$yv  # Use t suffix for treatments but no suffix for controls 
  
  Xr  <- rct_datc[,xv] # Predictor variables (independent)
  Xq  <- qed_datc[,xv] # Predictor variables (independent)
  Xrt <- rct_datt[,xv]
  
  ssr  <- rct_datc$siteID
  ssq  <- qed_datc$siteID
  ppr  <- rct_datc$partic
  wwq  <- qed_datc$ipw_wt
  
  ssrt <- rct_datt$siteID
  pprt <- rct_datt$partic
  
  # Initialize the tree object
  tree <- list()
  
  # Create a recursive function to build the tree
 
  node_n   <<- -1   #This <<- notation makes node_num global!!
  node_num <<- -1
  
  recursive_split <- function(yr, yq, yrt, Xr, Xq, Xrt, 
                              ssr, ssq, ssrt, ppr, pprt, wwq, 
                              depth = 1, parent_bias = NULL, 
                              leftv = "Root") {
    node_n   <<- node_n + 1
    node_num <<- rbind(node_num,node_n)
    
    if (is.null(parent_bias)) {
      
      # Calculate bias to return
      bias_res  <- calc_bias(yr, yq, yrt, ssr, ssq, ssrt, ppr, pprt, wwq)
      bias      <- bias_res$bias_node
      parent_bias <<- bias  # Make this global
    }
    
    # Create RCT and QED data frames
    rct_datcf <- data.frame(yr,Xr,ssr,ppr)
    qed_datcf <- data.frame(yq,Xq,ssq,wwq)
    rct_dattf <- data.frame(yrt,Xrt,ssrt,pprt)
    
    # Create stopping rules
    stop_rule <- 0
    bad_part  <- 0
    
    # Check whether sample sizes are too small by site and overall
    if (cace_itt == 0) {
      
      site_rct <- data.frame(table(ssr))
      site_rct <- site_rct %>%
        rename(site_var = ssr, rct_freq = Freq)
      site_rct <- site_rct[,c("site_var","rct_freq")]
      
      site_qed <- data.frame(table(ssq))
      site_qed <- site_qed %>%
        rename(site_var = ssq, qed_freq = Freq)
      site_qed <- site_qed[,c("site_var","qed_freq")]
      
      sitev_m <- merge(site_rct, site_qed, by = "site_var", all = TRUE)
      
      # Replace NA with 0 in site frequencies
      sitev_m$rct_freq <- ifelse(is.na(sitev_m$rct_freq),0,sitev_m$rct_freq)
      sitev_m$qed_freq <- ifelse(is.na(sitev_m$qed_freq),0,sitev_m$qed_freq)
      
      rct_site_freq <- any(sitev_m$rct_freq < minsplit_site)
      qed_site_freq <- any(sitev_m$qed_freq < minsplit_site)
      
      rct_site_tot  <- sum(sitev_m$rct_freq)
      qed_site_tot  <- sum(sitev_m$qed_freq)
      
      if ((rct_site_freq == TRUE) | (qed_site_freq == TRUE)) {
        stop_rule <- 1}
      
      if ((rct_site_tot < minsplit_tot) | (qed_site_tot < minsplit_tot)) {
        stop_rule <- 1}
      
      # Unique ys
      if ((length(unique(yr)) == 1) & (length(unique(yq)) == 1)) {
        stop_rule <- 1}
    
    } else if (cace_itt == 1) {
      
      rct_datc_npf <- rct_datcf[rct_datcf$ppr == 0,]
      rct_datt_pf  <- rct_dattf[rct_dattf$pprt == 1,]
      
      site_qedc <- data.frame(table(qed_datcf$ssq))
      site_qedc <- site_qedc %>%
        rename(site_var = Var1,
               qed_freqc = Freq)
      site_qedc <- site_qedc[,c("site_var","qed_freqc")]
      
      site_rctc_np <- data.frame(table(rct_datc_npf$ssr))
      site_rctc_np <- site_rctc_np %>%
        rename(site_var = Var1,
               rct_freqc_np = Freq)
      site_rctc_np <- site_rctc_np[,c("site_var","rct_freqc_np")]
      
      site_rctt_p <- data.frame(table(rct_datt_pf$ssrt))
      site_rctt_p <- site_rctt_p %>%
        rename(site_var = Var1,
               rct_freqt_p = Freq)
      site_rctt_p <- site_rctt_p[,c("site_var","rct_freqt_p")]
      
      sitev_m1 <- merge(site_qedc, site_rctc_np, , by = "site_var", all = TRUE)
      sitev_m  <- merge(sitev_m1, site_rctt_p, by = "site_var", all = TRUE)
      
      # Replace NA with 0 in site frequencies
      sitev_m$qed_freqc <- ifelse(is.na(sitev_m$qed_freqc),0,sitev_m$qed_freqc)
      sitev_m$rct_freqc_np <- ifelse(is.na(sitev_m$rct_freqc_np),0,
                                     sitev_m$rct_freqc_np)
      sitev_m$rct_freqt_p  <- ifelse(is.na(sitev_m$rct_freqt_p),0,
                                     sitev_m$rct_freqt_p)
      
      # Check if there are enough RCT and QED sample members in each site and in total
      # using the minsplit_site and minsplit_tot inputs
      rct_site_freqc_np <- any(sitev_m$rct_freqc_np < minsplit_site)
      rct_site_freqt_p  <- any(sitev_m$rct_freqt_p < minsplit_site)
      qed_site_freqc    <- any(sitev_m$qed_freqc < minsplit_site)
      
      rct_site_totc_np  <- sum(sitev_m$rct_freqc_np)
      rct_site_tott_p   <- sum(sitev_m$rct_freqt_p)
      qed_site_totc     <- sum(sitev_m$qed_freqc)
      
      if ((rct_site_freqc_np == TRUE) | (rct_site_freqt_p == TRUE) | 
          (qed_site_freqc == TRUE)) {
        stop_rule <- 1 }
      
      if ((rct_site_totc_np < minsplit_tot) | (rct_site_tott_p < minsplit_tot) | 
          (qed_site_totc < minsplit_tot)) {
        stop_rule <- 1 }
      
      # Check if the participation rate is greater for RCT treatments than controls
      # in all sites
      
      if (stop_rule == 0) {
      
      part_ratet <- data.frame(aggregate(pprt ~ ssrt, data = rct_dattf, FUN = mean))
      part_ratec <- data.frame(aggregate(ppr  ~ ssr,  data = rct_datcf, FUN = mean))
      
      compl_rate <- part_ratet$pprt - part_ratec$ppr
      
      any_bad_compl <- any(compl_rate <= 0)
      
      if (any_bad_compl == TRUE) {
        bad_part  <- 1
        stop_rule <- 1 }
      
      } # if stop_rule==0
      
      # Unique ys
      if ((length(unique(rct_datc_npf$yr)) == 1) & (length(unique(yq)) == 1)) {
        stop_rule <- 1}
      
    } # end if cace_itt
      
    # Depth rule
    if (depth > max_depth) {
      stop_rule <- 1}
    
    # If stopping conditions are met (e.g., no further splits), 
    # return a leaf node
    
    if (stop_rule == 1) {
      
      # Calculate bias to return
      bias_res <- calc_bias(yr, yq, yrt, ssr, ssq, ssrt, ppr, pprt, wwq)
      bias <- bias_res$bias_node
      return(list(type = "leaf", prediction = bias, nsampr = length(yr), 
                  nsampq = length(yq)))
    } # if stop_rule == 1
    
    # Calculate the best split 
    best_split <- find_best_split(yr, yq, yrt, Xr, Xq, Xrt, 
                                  ssr, ssq, ssrt, ppr, pprt, wwq)
    
    #bbb <- 1
    
    if (is.null(best_split)) {
      
      # Calculate bias to return
      bias_res <- calc_bias(yr, yq, yrt, ssr, ssq, ssrt, ppr, pprt, wwq)
      bias <- bias_res$bias_node
      return(list(type = "leaf", prediction = bias, nsampr = length(yr), 
                  nsampq = length(yq)))
    } # if no best split
    
    # Recursively split the data
    l_leftr  <- length(best_split$left_indicesr)
    l_rightr <- length(best_split$right_indicesr)
    l_leftq  <- length(best_split$left_indicesq)
    l_rightq <- length(best_split$right_indicesq)
    
    l_leftrt  <- length(best_split$left_indicesrt)
    l_rightrt <- length(best_split$right_indicesrt)
    
    #cat("\n")
    if (holdout == 0) {
      if (cace_itt == 0) {
        if (node_n == 0) {
          cat("\n")
          pr1 <- "BUILDING THE BIAS TREE: THE NUMBER OF RCT CONTROLS AND QED COMPARISONS"
          pr2 <- "SEQUENTIALLY ASSIGNED TO THE LEFT AND RIGHT TREE NODES FOR THE ITT ANALYSIS"
          cat(pr1,"\n",pr2,"\n")
        }
        print(paste(l_leftr, l_leftq, l_rightr, l_rightq))
      } else if (cace_itt == 1) {
        if (node_n == 0) {
          cat("\n")
          pr1 <- "BUILDING THE BIAS TREE: THE NUMBER OF RCT CONTROLS, QED COMPARISONS, AND RCT"
          pr2 <- "TREATMENTS ASSIGNED TO THE LEFT AND RIGHT TREE NODES FOR THE CACE ANALYSIS"
          cat(pr1,"\n",pr2,"\n")
        }
        print(paste(l_leftr, l_leftq, l_leftrt, l_rightr, l_rightq, l_rightrt))
      }
    } else if (holdout == 1) {
      if (cace_itt == 0) {
        if (node_n == 0) {
          cat("\n")
          pr1 <- "BUILDING THE BIAS TREE: THE NUMBER OF RCT CONTROLS AND QED COMPARISONS"
          pr2 <- "IN THE 70 PERCENT TRAINING SAMPLE SEQUENTIALLY ASSIGNED TO THE"
          pr3 <- "LEFT AND RIGHT TREE NODES FOR THE ITT ANALYSIS"
          cat(pr1,"\n",pr2,"\n",pr3,"\n")
        }
        print(paste(l_leftr, l_leftq, l_rightr, l_rightq))
      } else if (cace_itt == 1) {
        if (node_n == 0) {
          cat("\n")
          pr1 <- "BUILDING THE BIAS TREE: THE NUMBER OF RCT CONTROLS, QED COMPARISONS,"
          pr2 <- "AND RCT TREATMENTS IN THE 70 PERCENT TRAINING SAMPLE SEQUENTIALLY"
          pr3 <- "ASSIGNED TO THE LEFT AND RIGHT TREE NODES FOR THE CACE ANALYSIS"
          cat(pr1,"\n",pr2,"\n",pr3,"\n")
        }
        print(paste(l_leftr, l_leftq, l_leftrt, l_rightr, l_rightq, l_rightrt))
      } # end cace_itt
    } # end holdout
    
    left_datar    <- Xr[best_split$left_indicesr, , drop = FALSE]
    right_datar   <- Xr[best_split$right_indicesr, , drop = FALSE]
    left_targetr  <- yr[best_split$left_indicesr]
    right_targetr <- yr[best_split$right_indicesr]
    left_ssr      <- ssr[best_split$left_indicesr]
    right_ssr     <- ssr[best_split$right_indicesr]
    left_ppr      <- ppr[best_split$left_indicesr]
    right_ppr     <- ppr[best_split$right_indicesr]
    
    left_dataq    <- Xq[best_split$left_indicesq, , drop = FALSE]
    right_dataq   <- Xq[best_split$right_indicesq, , drop = FALSE]
    left_targetq  <- yq[best_split$left_indicesq]
    right_targetq <- yq[best_split$right_indicesq]
    left_ssq      <- ssq[best_split$left_indicesq]
    right_ssq     <- ssq[best_split$right_indicesq]
    left_wwq      <- wwq[best_split$left_indicesq]
    right_wwq     <- wwq[best_split$right_indicesq]
    
    left_datart    <- Xrt[best_split$left_indicesrt, , drop = FALSE]
    right_datart   <- Xrt[best_split$right_indicesrt, , drop = FALSE]
    left_targetrt  <- yrt[best_split$left_indicesrt]
    right_targetrt <- yrt[best_split$right_indicesrt]
    left_ssrt      <- ssrt[best_split$left_indicesrt]
    right_ssrt     <- ssrt[best_split$right_indicesrt]
    left_pprt      <- pprt[best_split$left_indicesrt]
    right_pprt     <- pprt[best_split$right_indicesrt]
    
    # Recursively create subtrees for the left and right nodes
    
    # First create left and right biases to input for parent_bias
    left_bias_res  <- calc_bias(left_targetr, left_targetq, left_targetrt, 
                                left_ssr, left_ssq, left_ssrt, 
                                left_ppr, left_pprt, left_wwq)
    left_bias <- left_bias_res$bias_node
    
    right_bias_res <- calc_bias(right_targetr, right_targetq, right_targetrt, 
                                right_ssr, right_ssq, right_ssrt, 
                                right_ppr, right_pprt, right_wwq)
    right_bias <- right_bias_res$bias_node
    
    left_tree  <- recursive_split(left_targetr, left_targetq, left_targetrt, 
                                  left_datar, left_dataq, left_datart, 
                                  left_ssr, left_ssq, left_ssrt, 
                                  left_ppr, left_pprt, left_wwq,
                                  depth + 1, 
                                  parent_bias <<- left_bias,
                                  leftv = "Left")
    right_tree <- recursive_split(right_targetr, right_targetq, right_targetrt,
                                  right_datar, right_dataq, right_datart,
                                  right_ssr, right_ssq, right_ssrt, 
                                  right_ppr, right_pprt, right_wwq,
                                  depth + 1, 
                                  parent_bias <<- right_bias, 
                                  leftv = "Right")
    
    # Return the current node (split)
    
    return(list(type = "node", 
                split_var = best_split$split_var,
                split_value = best_split$split_value,
                left = left_tree,
                right = right_tree)) 
    
  } # end recursive split function
  
  # Function to calculate the best split using MSE for regression
  
  find_best_split <- function(yr, yq, yrt, Xr, Xq, Xrt,
                              ssr, ssq, ssrt, ppr, pprt, wwq) {
    
    rct_datcff <- data.frame(yr,Xr,ssr,ppr)
    qed_datcff <- data.frame(yq,Xq,ssq,wwq)
    rct_dattff <- data.frame(yrt,Xrt,ssrt,pprt)
    rct_datc_npff <- rct_datcff[rct_datcff$ppr == 0,]
    rct_datt_pff  <- rct_dattff[rct_dattff$pprt == 1,]
    
    # Define Xs for RCT control nonparticpants and RCT treatment participants.
    # This is tricky since Xr and Xrt are matrices so pull off the right columns
    ncolrc <- ncol(rct_datcff)
    ncolrt <- ncol(rct_dattff)
    
    Xr_np <- rct_datc_npff[,2:(ncolrc-2)]
    Xrt_p <- rct_datt_pff[,2:(ncolrt-2)]
    
    best_split <- NULL
    best_score <- -10000 #Inf
    
    # Loop through all possible splits and find the one with the best score
    
    # First create Xrq that stacks Xr and Xq
    Xrq <- rbind(Xr,Xq)
    
      for (var in names(Xrq)) {
        
        # Calculate midpoints of values like rpart
        valuesa <- unique(Xrq[[var]])
        valuesb <- valuesa[order(valuesa)]
        values  <- valuesb[-length(valuesb)] + (diff(valuesb)/2)
        
        for (value in values) {
    
        # Split data based on the variable and its value
        left_indicesra  <- which(Xr[[var]] <= value)
        right_indicesra <- which(Xr[[var]] > value)
        
        left_indicesqa  <- which(Xq[[var]] <= value)
        right_indicesqa <- which(Xq[[var]] > value)
        
        left_indicesrta  <- which(Xrt[[var]] <= value)
        right_indicesrta <- which(Xrt[[var]] > value)
        
        left_indicesr_npa  <- which(Xr_np[[var]] <= value)
        right_indicesr_npa <- which(Xr_np[[var]] > value)
        
        left_indicesrt_pa  <- which(Xrt_p[[var]] <= value)
        right_indicesrt_pa <- which(Xrt_p[[var]] > value)
        
        l_leftr   <- length(left_indicesra)
        l_rightr  <- length(right_indicesra)
        l_leftq   <- length(left_indicesqa)
        l_rightq  <- length(right_indicesqa)
        l_leftrt  <- length(left_indicesrta)
        l_rightrt <- length(right_indicesrta)
        
        l_leftr_np  <- length(left_indicesr_npa)
        l_rightr_np <- length(right_indicesr_npa)
        l_leftrt_p  <- length(left_indicesrt_pa)
        l_rightrt_p <- length(right_indicesrt_pa)
        
        #print(l_leftq)
        #print(l_rightq)
        #print(l_leftr_np)
        #print(l_rightr_np)
        #print(l_leftrt_p)
        #print(l_rightrt_p)
        
        stop_rule <- 0
        no_ind    <- 0
        bad_part  <- 0
        
        if (cace_itt == 0) {
          if ((l_leftq == 0) | (l_rightq == 0) | 
              (l_leftr == 0) | (l_rightr == 0) ) {
            stop_rule <- 1
            no_ind    <- 1
          }
        } else if (cace_itt == 1) {
          if ((l_leftq == 0) | (l_rightq == 0) |
              (l_leftr_np == 0) | (l_rightr_np == 0) |
              (l_leftrt_p == 0) | (l_rightrt_p == 0)) {
            stop_rule <- 1
            no_ind    <- 1
          }
        } # end cace_itt
        
        # Create minbucket stop rules
        if (no_ind == 0) {
        
        minbucket_site <- round(minsplit_site/3)
        minbucket_tot  <- round(minsplit_tot/3)
        
        if (cace_itt == 0) {
          
          # See if sample sizes are too small by site and overall using minbucket
          site_rctl <- data.frame(table(ssr[left_indicesra]))
          site_rctl <- site_rctl %>%
            rename(site_var = Var1, rct_freql = Freq)
          site_rctl <- site_rctl[,c("site_var","rct_freql")]
          
          site_rctr <- data.frame(table(ssr[right_indicesra]))
          site_rctr <- site_rctr %>%
            rename(site_var = Var1, rct_freqr = Freq)
          site_rctr <- site_rctr[,c("site_var","rct_freqr")]
          
          site_qedl <- data.frame(table(ssq[left_indicesqa]))
          site_qedl <- site_qedl %>%
            rename(site_var = Var1, qed_freql = Freq)
          site_qedl <- site_qedl[,c("site_var","qed_freql")]
          
          site_qedr <- data.frame(table(ssq[right_indicesqa]))
          site_qedr <- site_qedr %>%
            rename(site_var = Var1, qed_freqr = Freq)
          site_qedr <- site_qedr[,c("site_var","qed_freqr")]
          
          sitev_ma <- merge(site_rctl, site_rctr, by = "site_var", 
                            all = TRUE)
          sitev_mb <- merge(sitev_ma, site_qedl, by = "site_var", 
                            all = TRUE)
          sitev_m  <- merge(sitev_mb, site_qedr, by = "site_var", 
                            all = TRUE)
          
          # Replace NA with 0 in site frequencies
          sitev_m$rct_freql <- ifelse(is.na(sitev_m$rct_freql),0,
                                      sitev_m$rct_freql)
          sitev_m$rct_freqr <- ifelse(is.na(sitev_m$rct_freqr),0,
                                      sitev_m$rct_freqr)
          sitev_m$qed_freql <- ifelse(is.na(sitev_m$qed_freql),0,
                                      sitev_m$qed_freql)
          sitev_m$qed_freqr <- ifelse(is.na(sitev_m$qed_freqr),0,
                                      sitev_m$qed_freqr)
          
          rct_site_freql <- any(sitev_m$rct_freql < minbucket_site)
          rct_site_freqr <- any(sitev_m$rct_freqr < minbucket_site)
          qed_site_freql <- any(sitev_m$qed_freql < minbucket_site)
          qed_site_freqr <- any(sitev_m$qed_freqr < minbucket_site)
          
          rct_site_totl  <- sum(sitev_m$rct_freql)
          rct_site_totr  <- sum(sitev_m$rct_freqr)
          qed_site_totl  <- sum(sitev_m$qed_freql)
          qed_site_totr  <- sum(sitev_m$qed_freqr)
          
          if ((rct_site_freql == TRUE) | (rct_site_freqr == TRUE) |
              (qed_site_freql == TRUE) | (qed_site_freqr == TRUE)) {
            stop_rule <- 1
          }
          
          if ((rct_site_totl < minbucket_tot) | (rct_site_totr < minbucket_tot) |
              (qed_site_totl < minbucket_tot) | (qed_site_totr < minbucket_tot)) {
            stop_rule <- 1
          }
          
        } else if (cace_itt == 1) {
          
          # See if sample sizes are too small by site and overall using minbucket
          site_qedcl <- data.frame(table(qed_datcff$ssq[left_indicesqa]))
          site_qedcl <- site_qedcl %>%
            rename(site_var = Var1, qed_freqcl = Freq)
          site_qedcl <- site_qedcl[,c("site_var","qed_freqcl")]
          
          site_qedcr <- data.frame(table(qed_datcff$ssq[right_indicesqa]))
          site_qedcr <- site_qedcr %>%
            rename(site_var = Var1, qed_freqcr = Freq)
          site_qedcr <- site_qedcr[,c("site_var","qed_freqcr")]
          
          site_rctc_npl <- data.frame(table(rct_datc_npff$ssr[left_indicesr_npa]))
          site_rctc_npl <- site_rctc_npl %>%
            rename(site_var = Var1,
                   rct_freqc_npl = Freq)
          site_rctc_npl <- site_rctc_npl[,c("site_var","rct_freqc_npl")]
          
          site_rctc_npr <- data.frame(table(rct_datc_npff$ssr[right_indicesr_npa]))
          site_rctc_npr <- site_rctc_npr %>%
            rename(site_var = Var1,
                   rct_freqc_npr = Freq)
          site_rctc_npr <- site_rctc_npr[,c("site_var","rct_freqc_npr")]
          
          site_rctt_pl <- data.frame(table(rct_datt_pff$ssrt[left_indicesrt_pa]))
          site_rctt_pl <- site_rctt_pl %>%
            rename(site_var = Var1,
                   rct_freqt_pl = Freq)
          site_rctt_pl <- site_rctt_pl[,c("site_var","rct_freqt_pl")]
          
          site_rctt_pr <- data.frame(table(rct_datt_pff$ssrt[right_indicesrt_pa]))
          site_rctt_pr <- site_rctt_pr %>%
            rename(site_var = Var1,
                   rct_freqt_pr = Freq)
          site_rctt_pr <- site_rctt_pr[,c("site_var","rct_freqt_pr")]
          
          sitev_m1 <- merge(site_qedcl, site_qedcr, by = "site_var", all = TRUE)
          sitev_m2 <- merge(sitev_m1, site_rctc_npl, by = "site_var", all = TRUE)
          sitev_m3 <- merge(sitev_m2, site_rctc_npr, by = "site_var", all = TRUE)
          sitev_m4 <- merge(sitev_m3, site_rctt_pl, by = "site_var", all = TRUE)
          sitev_m  <- merge(sitev_m4, site_rctt_pr, by = "site_var", all = TRUE)
          
          # Replace NA with 0 in site frequencies
          sitev_m$qed_freqcl <- ifelse(is.na(sitev_m$qed_freqcl),0,
                                       sitev_m$qed_freqcl)
          sitev_m$qed_freqcr <- ifelse(is.na(sitev_m$qed_freqcr),0,
                                       sitev_m$qed_freqcr)
          sitev_m$rct_freqc_npl <- ifelse(is.na(sitev_m$rct_freqc_npl),0,
                                          sitev_m$rct_freqc_npl)
          sitev_m$rct_freqc_npr <- ifelse(is.na(sitev_m$rct_freqc_npr),0,
                                          sitev_m$rct_freqc_npr)
          sitev_m$rct_freqt_pl  <- ifelse(is.na(sitev_m$rct_freqt_pl),0,
                                          sitev_m$rct_freqt_pl)
          sitev_m$rct_freqt_pr  <- ifelse(is.na(sitev_m$rct_freqt_pr),0,
                                          sitev_m$rct_freqt_pr)
          
          # Check if there are enough RCT and QED sample members in each site and in total
          # using the minbucket_site and minbucket_tot inputs
          qed_site_freqcl    <- any(sitev_m$qed_freqcl < minbucket_site)
          qed_site_freqcr    <- any(sitev_m$qed_freqcr < minbucket_site)
          rct_site_freqc_npl <- any(sitev_m$rct_freqc_npl < minbucket_site)
          rct_site_freqc_npr <- any(sitev_m$rct_freqc_npr < minbucket_site)
          rct_site_freqt_pl  <- any(sitev_m$rct_freqt_pl < minbucket_site)
          rct_site_freqt_pr  <- any(sitev_m$rct_freqt_pr < minbucket_site)
          
          qed_site_totcl     <- sum(sitev_m$qed_freqcl)
          qed_site_totcr     <- sum(sitev_m$qed_freqcr)
          rct_site_totc_npl  <- sum(sitev_m$rct_freqc_npl)
          rct_site_totc_npr  <- sum(sitev_m$rct_freqc_npr)
          rct_site_tott_pl   <- sum(sitev_m$rct_freqt_pl)
          rct_site_tott_pr   <- sum(sitev_m$rct_freqt_pr)
          
          if ((qed_site_freqcl == TRUE) | (qed_site_freqcr == TRUE) |
              (rct_site_freqc_npl == TRUE) | (rct_site_freqc_npr == TRUE) | 
              (rct_site_freqt_pl == TRUE) | (rct_site_freqt_pr == TRUE)) {
            stop_rule <- 1 
            }
          
          if ((qed_site_totcl < minbucket_tot) | (qed_site_totcr < minbucket_tot) |
              (rct_site_totc_npl < minbucket_tot) | 
              (rct_site_totc_npr < minbucket_tot) |
              (rct_site_tott_pl < minbucket_tot) | 
              (rct_site_tott_pr < minbucket_tot)) {
            stop_rule <- 1 
          }
          
          # Check if the participation rate is greater for RCT treatments than 
          # controls in all sites in left and right nodes if stop_rule is 0
          
          if (stop_rule == 0) {
          
          rct_trtl <- rct_dattff[left_indicesrta,]
          rct_conl <- rct_datcff[left_indicesra,]
          
          part_ratetl <- data.frame(aggregate(pprt ~ ssrt, data = rct_trtl, 
                                             FUN = mean))
          part_ratecl <- data.frame(aggregate(ppr  ~ ssr,  data = rct_conl, 
                                             FUN = mean))
          
          rct_trtr <- rct_dattff[right_indicesrta,]
          rct_conr <- rct_datcff[right_indicesra,]
          
          part_ratetr <- data.frame(aggregate(pprt ~ ssrt, data = rct_trtr, 
                                              FUN = mean))
          part_ratecr <- data.frame(aggregate(ppr  ~ ssr,  data = rct_conr, 
                                              FUN = mean))
          
          compl_ratel <- part_ratetl$pprt - part_ratecl$ppr
          compl_rater <- part_ratetr$pprt - part_ratecr$ppr
          
          any_bad_compll <- any(compl_ratel <= 0)
          any_bad_complr <- any(compl_rater <= 0)
          
          if ((any_bad_compll == TRUE) | (any_bad_complr == TRUE)) {
            bad_part  <- 1
            stop_rule <- 1 }
          
          #print(paste(stop_rule,any_bad_compll,any_bad_complr))
          
          } # end if stop_rule == 0
          
        } # if cace_itt
        
        } # end NO_IND == 0 
        
        if (stop_rule == 1) {
          next
        }
        
        # Calculate score = Between SSE as in (16) in article
        
        biasl_res <- calc_bias(yr[left_indicesra], yq[left_indicesqa], 
                               yrt[left_indicesrta],
                               ssr[left_indicesra], ssq[left_indicesqa],
                               ssrt[left_indicesrta],
                               ppr[left_indicesra], pprt[left_indicesrta],
                               wwq[left_indicesqa])
        biasl  <- biasl_res$bias_node
        wgtl   <- biasl_res$sw
        
        biasr_res <- calc_bias(yr[right_indicesra], yq[right_indicesqa],
                               yrt[right_indicesrta],
                               ssr[right_indicesra], ssq[right_indicesqa],
                               ssrt[right_indicesrta],
                               ppr[right_indicesra], pprt[right_indicesrta], 
                               wwq[right_indicesqa])
        biasr  <- biasr_res$bias_node
        wgtr   <- biasr_res$sw
        
        w1 <- wgtl/(wgtl+wgtr)
        w2 <- wgtr/(wgtl+wgtr)
        
        score <- (w1*((biasl - parent_bias)^2)) +
          (w2*((biasr - parent_bias)^2))
        
        #print(paste("FIND_BEST_SPLIT CACE_ITT",cace_itt,w1,biasl,w2,biasr))
        
        # If the score is the best (max) so far, store it
        if (score > best_score) {
          best_score <- score
          
          best_split <- list(split_var = var, 
                             split_value = value,
                             left_indicesr = left_indicesra, 
                             right_indicesr = right_indicesra, 
                             left_indicesq = left_indicesqa,
                             right_indicesq = right_indicesqa,
                             left_indicesrt = left_indicesrta, 
                             right_indicesrt = right_indicesrta)
        }
      }  # for value
    } # for var
    
    #print(best_split)

    return(best_split)
    
  } # end of find_best_fit 
  
  # Start the tree construction
  tree$root <- recursive_split(yr, yq, yrt, Xr, Xq, Xrt,
                               ssr, ssq, ssrt, ppr, pprt, wwq)
  
  # Return the constructed tree
  return(tree)
  
} # End custom_rpart function

# Run custom_rpart function
rt <- custom_rpart()

#print(rt)

# Create a recursive function to process rt that contains tree results 
# in a complicated nested list

tree_results <<- -1  # <<- means it is global!
recursive_nested_processor <- function(my_list) {
  for (item in my_list) {
    if (is.list(item)) {
      # If it's a sublist, recurse
      recursive_nested_processor(item) 
    } else {
      # If it's a non-list element, process it
      #print(paste("Processing element:", item))
      #print(item)
      tree_results <<- rbind(tree_results,item)
    }
  }
}

#print("TREE_RESULTS")
#print(tree_results)

#print(paste("RT LIST:",rt))

# Run recursive nested function
recursive_nested_processor(rt)

# Get Right Left tree movement patterns using listr and stringr libraries 

# Flatten rt list and save names in list
flat_rt    <- list_flatten(rt, max_depth = -1)
left_right <- names(flat_rt)

# Only keep the left_right_type names and then
# remove "_type" from the names
left_right1 <- left_right[str_detect(left_right,"type") == TRUE]
left_right1 <- str_replace(left_right1,"_type","")
left_right1 <- str_replace(left_right1,"root_","")

left_right1 <- str_replace_all(left_right1,"right","R")
left_right1 <- str_replace_all(left_right1,"left","L")
left_right  <- str_replace_all(left_right1,"_","")

left_rightp <- format(left_right, justify = "left") # Left justified for printing

nleft_right <- length(left_right)

# Only continue if there is more than one split
if (nleft_right >= 3) {
  
  # Get vector of node_types from tree_results and create data frames
  # to merge on split_var and split_val and predictions later
  node_type    <- tree_results[tree_results == "node" | tree_results == "leaf"]
  node_type_df <- data.frame(node_type)
  
  # Create ID variables for sorting later
  node_type_df$idnum <- seq_len(nrow(node_type_df))
  
  # Create node and leaf data frames
  node_df <- node_type_df[node_type_df$node_type == "node",]
  leaf_df <- node_type_df[node_type_df$node_type == "leaf",]
  
  # Store split vars and values from tree_results for non-leaf nodes
  split_var <- ifelse(tree_results == 'node',lead(tree_results),"")
  split_val <- ifelse(tree_results == 'node',lead(tree_results,n=2),"")
  
  split_var <- split_var[split_var != ""]
  split_val <- as.numeric(split_val[split_val != ""])
  
  split_valp <- format(round(split_val,4), nsmall=4) # For printing
  
  # Add split vars and values to node_df and initialize predictions and nsamp 
  node_df$split_var   <- split_var
  node_df$split_val   <- split_val
  node_df$split_valp  <- split_valp
  node_df$prediction  <- ""
  node_df$predictionp <- ""
  node_df$nsampr      <- ""
  node_df$nsampq      <- ""
  
  # Get predictions and sample sizes for leaf nodes
  prediction <- ifelse(tree_results == 'leaf',lead(tree_results),"")
  nsampr     <- ifelse(tree_results == 'leaf',lead(tree_results,n=2),"")
  nsampq     <- ifelse(tree_results == 'leaf',lead(tree_results,n=3),"")
  
  nsampr     <- nsampr[nsampr != ""]
  nsampq     <- nsampq[nsampq != ""]
  prediction <- as.numeric(prediction[prediction != ""])
  
  predictionp <- format(round(prediction,4), nsmall=4) # For printing
  
  leaf_df$split_var   <- ""
  leaf_df$split_val   <- ""
  leaf_df$split_valp  <- ""
  leaf_df$prediction  <- prediction
  leaf_df$predictionp <- predictionp
  leaf_df$nsampr      <- nsampr
  leaf_df$nsampq      <- nsampq
  
  # Stack node and leaf data frames
  nodeleaf_df <- rbind(node_df,leaf_df)
  
  # Sort back to original order
  treesum_df <- nodeleaf_df[order(nodeleaf_df$idnum), ]
  
  # Get correct vector of node numbers from node_num that removes first -1 value
  node_num <- node_num[node_num != -1]
  
  # Add left_right and node_num to tree data frame of results and 
  # remove id which isn't needed anymore
  treesum_df$left_right  <- left_right
  treesum_df$left_rightp <- left_rightp
  treesum_df$node_num    <- node_num
  
  treesum_df$idnum <- NULL
  
  # Number of rows in tree
  n_treesum <- nrow(treesum_df)
  
  # Create if statements
  
  # Number of characters (length) of left-right entry
  treesum_df$nlr <- nchar(trimws(treesum_df$left_right))
  
  # Initialize if conditions
  treesum_df$ifcondr  <- c("")
  treesum_df$ifcondq  <- c("")
  treesum_df$ifcondrt <- c("")
  
  # Store if variables and values for the root - that is, for the first split
  var_ifr  <- paste0("dat_noder$",treesum_df$split_var[1])
  var_ifq  <- paste0("dat_nodeq$",treesum_df$split_var[1])
  var_ifrt <- paste0("dat_nodert$",treesum_df$split_var[1])
  val_if   <- treesum_df$split_val[1] 
  
  # Calculate if conditions for the L and R left-right entries and find their
  # row numbers for later
  treesum_df$ifcondr[treesum_df$left_right == 'L'] <- 
    paste("(",var_ifr,"<",val_if,")",sep=" ")
  
  treesum_df$ifcondr[treesum_df$left_right == 'R'] <- 
    paste("(",var_ifr,">=",val_if,")",sep=" ")
  
  treesum_df$ifcondq[treesum_df$left_right == 'L'] <- 
    paste("(",var_ifq,"<",val_if,")",sep=" ")
  
  treesum_df$ifcondq[treesum_df$left_right == 'R'] <- 
    paste("(",var_ifq,">=",val_if,")",sep=" ")
  
  treesum_df$ifcondrt[treesum_df$left_right == 'L'] <- 
    paste("(",var_ifrt,"<",val_if,")",sep=" ")
  
  treesum_df$ifcondrt[treesum_df$left_right == 'R'] <- 
    paste("(",var_ifrt,">=",val_if,")",sep=" ")
  
  lrow <- treesum_df$node_num[treesum_df$left_right == 'L']+1
  rrow <- treesum_df$node_num[treesum_df$left_right == 'R']+1
  
  # Loop over left-right entries to create if conditions, separately for
  # initial L and R splits
  
  # First for left initial splits
  startl <- lrow+1
  endl   <- rrow-1
  if (startl <= endl) {
    for (i in startl:endl) {
      
      # Find the substring that takes off the last element of left_right
      lri   <- treesum_df$left_right[i]
      nlri  <- treesum_df$nlr[i]
      lri_nolast  <- substring(lri,1,(nlri-1))
      
      # Find the last element of left_right
      lri_last  <- substring(lri,nlri,nlri)
      
      #print(paste(lri,lri_nolast,lri_last))
      
      # Do over the left_right entries until it matches lri_nolast
      for (j in lrow:(i-1)) {
        
        lrj       <- treesum_df$left_right[j]
        type_lrj  <- treesum_df$node_type[j]
        var_lrj   <- treesum_df$split_var[j]
        val_lrj   <- treesum_df$split_val[j]  
        if_lrjr   <- treesum_df$ifcondr[j]
        if_lrjq   <- treesum_df$ifcondq[j]
        if_lrjrt  <- treesum_df$ifcondrt[j]
        
        if (lrj == lri_nolast) {
          
          # Create new split vars and values 
          var_ifr  <- paste0("dat_noder$",var_lrj) 
          var_ifq  <- paste0("dat_nodeq$",var_lrj)
          var_ifrt <- paste0("dat_nodert$",var_lrj) 
          val_if   <- val_lrj 
          
          # Calculate if conditions for the L and R last entries 
          if (lri_last == 'L') {
            treesum_df$ifcondr[i] <- 
              paste(if_lrjr,"&","(",var_ifr,"<",val_if,")",sep=" ")
            
            treesum_df$ifcondq[i] <- 
              paste(if_lrjq,"&","(",var_ifq,"<",val_if,")",sep=" ")
            
            treesum_df$ifcondrt[i] <- 
              paste(if_lrjrt,"&","(",var_ifrt,"<",val_if,")",sep=" ")
            
          } else if (lri_last == 'R') {
            treesum_df$ifcondr[i] <- 
              paste(if_lrjr,"&","(",var_ifr,">=",val_if,")",sep=" ")
            
            treesum_df$ifcondq[i] <- 
              paste(if_lrjq,"&","(",var_ifq,">=",val_if,")",sep=" ")
            
            treesum_df$ifcondrt[i] <- 
              paste(if_lrjrt,"&","(",var_ifrt,">=",val_if,")",sep=" ")
          } 
          
          break # exits for j loop
          
        } #if lrj=slri
      } # for j
    } # for i 
  } # if (startl <=endl)
  
  # Now repeat for right initial splits
  startr <- rrow+1
  endr   <- n_treesum
  if (startr <= endr) {
    for (i in startr:endr) {
      
      # Find the substring that takes off the last element of left_right
      lri   <- treesum_df$left_right[i]
      nlri  <- treesum_df$nlr[i]
      lri_nolast  <- substring(lri,1,(nlri-1))
      
      # Find the last element of left_right
      lri_last  <- substring(lri,nlri,nlri)
      
      # Do over the left_right entries until it matches lri_nolast
      
      for (j in rrow:(i-1)) {
        
        lrj       <- treesum_df$left_right[j]
        type_lrj  <- treesum_df$node_type[j]
        var_lrj   <- treesum_df$split_var[j]
        val_lrj   <- treesum_df$split_val[j]  
        if_lrjr   <- treesum_df$ifcondr[j]
        if_lrjq   <- treesum_df$ifcondq[j]
        if_lrjrt  <- treesum_df$ifcondrt[j]
        
        if (lrj == lri_nolast) {
          
          # Create new split vars and values 
          var_ifr  <- paste0("dat_noder$",var_lrj) 
          var_ifq  <- paste0("dat_nodeq$",var_lrj)
          var_ifrt <- paste0("dat_nodert$",var_lrj) 
          val_if   <- val_lrj 
          
          # Calculate if conditions for the L and R last entries 
          if (lri_last == 'L') {
            treesum_df$ifcondr[i] <- 
              paste(if_lrjr,"&","(",var_ifr,"<",val_if,")",sep=" ")
            
            treesum_df$ifcondq[i] <- 
              paste(if_lrjq,"&","(",var_ifq,"<",val_if,")",sep=" ")
            
            treesum_df$ifcondrt[i] <- 
              paste(if_lrjrt,"&","(",var_ifrt,"<",val_if,")",sep=" ")
            
          } else if (lri_last == 'R') {
            treesum_df$ifcondr[i] <- 
              paste(if_lrjr,"&","(",var_ifr,">=",val_if,")",sep=" ")
            
            treesum_df$ifcondq[i] <- 
              paste(if_lrjq,"&","(",var_ifq,">=",val_if,")",sep=" ")
            
            treesum_df$ifcondrt[i] <- 
              paste(if_lrjrt,"&","(",var_ifrt,">=",val_if,")",sep=" ")
          } 
          
          break # exits for j loop
          
        } #if lrj=slri
      } # for j
    } # for i 
  } # if (startr <= endr)
  
  # Create data frames for printing
  treesum_df$ifcondr  <- format(treesum_df$ifcondr, justify = "left") # to print
  treesum_df$ifcondq  <- format(treesum_df$ifcondq, justify = "left") # to print
  treesum_df$ifcondrt <- format(treesum_df$ifcondrt, justify = "left") 
  
  #treep_df  <- treesum_df[,c("left_rightp","node_num","node_type","split_var",
  #                           "split_valp","predictionp","nsampr",
  #                           "nsampq")]
  
  #treep1_df <- treesum_df[,c("left_rightp","ifcondr","ifcondq")]
  
  #print(treep_df)
  #print(treep1_df)
  
  # Parse if conditions, which means translate if statements into r code and  
  # then use eval to run that r code
  
  treesum_df$node_bias  <- rep(0,n_treesum)
  treesum_df$node_sw    <- rep(0,n_treesum)
  treesum_df$node_sampr <- rep(0,n_treesum)
  treesum_df$node_sampq <- rep(0,n_treesum)
  treesum_df$node_bss   <- rep(0,n_treesum)
  
  for (i in 1:n_treesum) {
    
    if (i == 1) {
      
      yr  <- dat_noder$yv  # Target variable (dependent)
      yq  <- dat_nodeq$yv  # Target variable (dependent)
      yrt <- dat_nodert$yv  # Target variable (dependent)
      Xr  <- dat_noder[,xv] # Predictor variables (independent)
      Xq  <- dat_nodeq[,xv] # Predictor variables (independent)
      Xrt <- dat_nodert[,xv] # Predictor variables (independent)
      ssr <- dat_noder$siteID
      ssq <- dat_nodeq$siteID
      ssrt <- dat_nodert$siteID
      ppr  <- dat_noder$partic
      pprt <- dat_nodert$partic
      wwq  <- dat_nodeq$ipw_wt
      
      bias_res  <- calc_bias(yr, yq, yrt, ssr, ssq, ssrt, ppr, pprt, wwq)
      bias      <- bias_res$bias_node
      sumwt     <- bias_res$sw
      
      treesum_df$node_bias[1]  <- bias   #mean(dat_node[,yvar])
      treesum_df$node_sw[1]    <- sumwt  
      treesum_df$node_sampr[1] <- NROW(yr)
      treesum_df$node_sampq[1] <- NROW(yq)
      treesum_df$node_bss[1]   <- sumwt*(bias^2)
      #treesum_df$node_rss[1]  <- var(dat_node[,yvar])*(NROW(dat_node[,yvar])-1)
      
    } else if (i > 1) {
      
      # Pull off if conditions for node
      ifnoder  <- treesum_df$ifcondr[i] 
      ifnodeq  <- treesum_df$ifcondq[i]
      ifnodert <- treesum_df$ifcondrt[i]
      
      # Define samp_in_node which subsets the observations in dat_node 
      # who are in node i
      
      ifparse_yr   <- paste0("samp_yr <- ", "yr", "[", ifnoder, "]")
      ifparse_yq   <- paste0("samp_yq <- ", "yq", "[", ifnodeq, "]")
      ifparse_yrt  <- paste0("samp_yrt <- ", "yrt", "[", ifnodert, "]")
      ifparse_ssr  <- paste0("samp_ssr <- ", "ssr", "[", ifnoder, "]")
      ifparse_ssq  <- paste0("samp_ssq <- ", "ssq", "[", ifnodeq, "]")
      ifparse_ssrt <- paste0("samp_ssrt <- ", "ssrt", "[", ifnodert, "]")
      ifparse_ppr  <- paste0("samp_ppr <- ", "ppr", "[", ifnoder, "]")
      ifparse_pprt <- paste0("samp_pprt <- ", "pprt", "[", ifnodert, "]")
      ifparse_wwq  <- paste0("samp_wwq <- ", "wwq", "[", ifnodeq, "]")
      
      # Parse the ifparse text and then eval it
      parse_code_yr <- parse(text=ifparse_yr)
      eval(parse_code_yr)
      
      parse_code_yq <- parse(text=ifparse_yq)
      eval(parse_code_yq)
      
      parse_code_yrt <- parse(text=ifparse_yrt)
      eval(parse_code_yrt)
      
      parse_code_ssr <- parse(text=ifparse_ssr)
      eval(parse_code_ssr)
      
      parse_code_ssq <- parse(text=ifparse_ssq)
      eval(parse_code_ssq)
      
      parse_code_ssrt <- parse(text=ifparse_ssrt)
      eval(parse_code_ssrt)
      
      parse_code_ppr <- parse(text=ifparse_ppr)
      eval(parse_code_ppr)
      
      parse_code_pprt <- parse(text=ifparse_pprt)
      eval(parse_code_pprt)
      
      parse_code_wwq <- parse(text=ifparse_wwq)
      eval(parse_code_wwq)
      
      # Calculate means, nsamp, and bias in node using the samp_in_node sample
      
      bias_res  <- calc_bias(samp_yr, samp_yq, samp_yrt,
                             samp_ssr, samp_ssq, samp_ssrt, 
                             samp_ppr, samp_pprt,samp_wwq)
      bias      <- bias_res$bias_node
      sumwt     <- bias_res$sw
      
      treesum_df$node_bias[i]  <- bias #mean(sampy_in_node)
      treesum_df$node_sw[i]    <- sumwt
      treesum_df$node_sampr[i] <- NROW(samp_yr)
      treesum_df$node_sampq[i] <- NROW(samp_yq)
      treesum_df$node_bss[i]   <- sumwt*(bias^2)
      #treesum_df$node_rss[i]  <- var(sampy_in_node)*(NROW(sampy_in_node)-1)
      
    } # if i>1
  } # for i
  
  treesum_df$pred_biasp <- treesum_df$predictionp
  treesum_df$node_biasp <- format(round(treesum_df$node_bias,4), nsmall=4)
  treesum_df$node_bssp  <- format(round(treesum_df$node_bss,2), nsmall=2) 
  treesum_df$node_swp   <- format(round(treesum_df$node_sw,0)) 
  
  # Prune tree using cost complexity, weakest link pruning. The idea is to
  # calculate ALPHA = (SSE(CURRENT TREE) - SSE(TREE WITH PRUNED NODE)))/(|T|-1),
  # Where |T| = number of nodes in current tree, SSE = between ssqs measured as the 
  # sum of the bias squared, and ALPHA = Complexity parameter. We do this for each node in the current tree 
  # and find the node with the smallest CP and then chop that off and
  # We keep on going until we hit the root node. Idea is Loss function is
  # C = SSE - ALPHA*|T| and this method finds the subtree with the maximum C
  # OLD for regular trees: C = SSE + ALPHA*|T| and find the subtree 
  # with the minimum C. The above is for regular CART - in our bias context, 
  # we instead use the between sum of squares, BSS = node_sw*(node_bias^2) which 
  # turns out to be the same thing as above for pruning, but TSSY is not observed so
  # we instead estimate it by specifying the R2 and noting that 
  # R2 ~ avg(node_sw)*BSS/TSS so we can solve for TSS. We use the avg comparison group
  # sample size for node_sw
  # This doesn't work since the estimated tss_bias increases with sample size 
  # so not easy to select a cutoff alpha value. So we rely on holdout sample
  
  # Initialize current BSS for the full tree
  bss_curr <- sum(treesum_df$node_bss[treesum_df$node_type == "leaf"]) 
  
  # First print unpruned tree
  
  treesum_df$node_n  <- seq_len(nrow(treesum_df))
  treesum_df$node_id <- format(treesum_df$node_n, justify = "left")
  
  # Assign better varable names
  df_temp <- treesum_df
  
  df_temp$left_right     <- df_temp$left_rightp
  df_temp$pred_bias      <- df_temp$node_biasp
  df_temp$bss            <- df_temp$node_bssp
  df_temp$eff_samp       <- df_temp$node_swp
  df_temp$sampr          <- df_temp$node_sampr
  df_temp$sampq          <- df_temp$node_sampq
  df_temp$split_val      <- df_temp$split_valp
  
  # Write to output file
  # Define sink function
  sink(file = out_txt)  #, append = TRUE)
  
  if (cace_itt == 0) {
    cat("\n")
    pr1 <- "TABLE 1. BIAS_TREE FUNCTION: RCT AND QED SITE SAMPLE SIZES IN"
    pr2 <- "THE RCT-QED SITES FOR THE FULL SAMPLE ITT ANALYSIS"
    cat(pr1,"\n",pr2,"\n") 
    print(sitev_mz1)
    
    if (holdout == 1) {
      cat("\n")
      pr1 <- "TABLE 1A. BIAS_TREE FUNCTION: RCT AND QED SITE SAMPLE SIZES IN THE" 
      pr2 <- "RCT-QED SITES FOR THE 70 PERCENT TRAINING SAMPLE ITT ANALYSIS"
      cat(pr1,"\n",pr2,"\n")
      print(sitev_mz2)
    }
  }
  
  if (cace_itt == 1) {
    
    cat("\n")
    pr1 <- "TABLE 1. BIAS_TREE FUNCTION: RCT CONTROL, RCT TREATMENT, AND QED COMPARISON"
    pr2 <- "SITE SAMPLE SIZES IN THE RCT-QED SITES FOR THE FULL SAMPLE CACE ANALYSIS"
    cat(pr1,"\n",pr2,"\n")
    print(sitev_mz1[,c("site_var","rct_freqc","rct_freqc_np",
                       "rct_freqt","rct_freqt_p","qed_freqc")])
    fn1 <- "Notes: Figures in order show samples for: (i) the full RCT control group,"
    fn2 <- "(ii) RCT control nonparticipants, (iii) the full RCT treatment group,"
    fn3 <- "(iii) RCT treatment participants, and (iv) the full QED comparison group."
    cat(fn1,"\n",fn2,"\n",fn3,"\n")
    
    cat("\n")
    pr1 <- "TABLE 2. BIAS_TREE FUNCTION: RCT TREATMENT AND CONTROL GROUP PARTICIPATION"
    pr2 <- "RATES AND THE COMPLIANCE RATE FOR THE FULL SAMPLE CACE ANALYSIS"
    cat(pr1,"\n",pr2,"\n") 
    print(pc_df1zz)
    
    if (holdout == 1) {
      cat("\n")
      pr1 <- "TABLE 1A. BIAS_TREE FUNCTION: RCT AND QED SITE SAMPLE SIZES IN THE" 
      pr2 <- "RCT-QED SITES FOR THE 70 PERCENT TRAINING SAMPLE CACE ANALYSIS"
      cat(pr1,"\n",pr2,"\n")
      print(sitev_mz2[,c("site_var","rct_freqc","rct_freqc_np",
                         "rct_freqt","rct_freqt_p","qed_freqc")])
      fn1 <- "Notes: Figures in order show samples for: (i) the RCT control group,"
      fn2 <- "(ii) RCT control nonparticipants, (iii) the RCT treatment group,"
      fn3 <- "(iii) RCT treatment participants, and (iv) the QED comparison group."
      cat(fn1,"\n",fn2,"\n",fn3,"\n")
      
      cat("\n")
      pr1 <- "TABLE 2A. BIAS_TREE FUNCTION: RCT TREATMENT AND CONTROL GROUP PARTICIPATION RATES"
      pr2 <- "AND THE COMPLIANCE RATE FOR THE 70 PERCENT TRAINING SAMPLE CACE ANALYSIS"
      cat(pr1,"\n",pr2,"\n")
      print(pc_df1zz2)
    }
  }
  
  sink()
  sink(file = out_txt, split = TRUE, append = TRUE)
  
  #cat("\n")
  if (cace_itt == 0) {
    cat("\n")
    cat("TABLE 2. BIAS_TREE FUNCTION: INITIAL UNPRUNED TREE FOR THE ITT ANALYSIS","\n")
  } else if (cace_itt == 1) {
    cat("\n")
    cat("TABLE 3. BIAS_TREE FUNCTION: INITIAL UNPRUNED TREE FOR THE CACE ANALYSIS","\n")
  }
  
  print(df_temp[,c("node_id","left_right",
                   "node_type","pred_bias","sampr","sampq","eff_samp", 
                   "split_var", "split_val")],row.names = FALSE,
        right='F')
  
  strp1 <- "Notes: The left_right variable displays the movement along the binary" 
  strp2 <- "tree branches where L signifies a left move if the split variable is "
  strp3 <- "less than the split value at the parent node, and R signifies a right "
  strp4 <- "move otherwise. The pred_bias variable is the node bias, node_type "
  strp5 <- "signifies whether the node is an internal tree node or a terminal leaf, "
  strp6 <- "sampr and sampq are RCT control and QED comparison sample sizes, and "
  strp7 <- "eff_samp is the total effective sample size accounting for design"
  strp8 <- "effects due to inverse probability weighting."
  cat(strp1,"\n",strp2,"\n",strp3,"\n", strp4,"\n", strp5,"\n", 
      strp6,"\n",strp7,"\n",strp8,"\n") 
  
  sink()
  
  # Initialize current treesum dataset
  ctreesum_df <- treesum_df 
  
  # Create current tree_sum_leaf_df with leaf nodes only
  ctreesum_leaf_df <- ctreesum_df[ctreesum_df$node_type == 'leaf',]
  nctree_sum_leaf  <- nrow(ctreesum_leaf_df)
  
  # Recursively create pruned trees where we find the minimum cp (alpha)
  # and save the associated tree until nleaves <= 2
  
  # Initialize variables and dfs for storing results
  cw  <- 0
  nleaves <- nctree_sum_leaf
  
  cp_sum <- data.frame(matrix(nrow = 1000, ncol = 5))
  colnames(cp_sum) <- c("subtree","cp", "node_chop", "nleaves","nsplits")
  
  cp_sum$subtree[1]      <- 1
  cp_sum$cp[1]           <- 0
  cp_sum$node_chop[1]    <- c("None") 
  cp_sum$nleaves[1]      <- nleaves
  cp_sum$nsplits[1]      <- nleaves - 1
  
  ctreesum_df_temp         <- ctreesum_df
  ctreesum_df_temp$subtree <- 1 
  
  subtrees_df <- ctreesum_df_temp
  
  # This is the WHILE LOOP until the number of leaves in the pruned tree 
  # is at least 3 (root L R)!!
  
  while (nleaves >= 3) {
    
    # Increase cw counter  
    cw <- cw + 1  
    
    #for (cw in 1:30) {
    
    min_alpha  <- Inf  
    min_leaves <- Inf  
    
    lri_nolast_prev <- c("JunkZZZ")
    
    # Number of characters (length) of left-right entry in current trees
    ctreesum_df$nlr      <- nchar(trimws(ctreesum_df$left_right))
    ctreesum_leaf_df$nlr <- nchar(trimws(ctreesum_leaf_df$left_right))
    
    # Do over all leaf nodes to find weakest link in current tree
    
    ii <- 0
    for (i in 1:nleaves) {                  
      
      # Initialize tsl_prune_df which contains leaves only and ts_prune
      # which contains leaves and nodes
      tsl_prune_df <- ctreesum_leaf_df
      ts_prune_df  <- ctreesum_df
  
      # Find the length of left_right for leaf node i
      nlri  <- tsl_prune_df$nlr[i]
      
      # Only continue if nlri > 1
      if (nlri > 1) {
        
        # Find the left_right string of tsl_prune_df that excludes 
        # the last L or R entry of leaf i
        lri        <- tsl_prune_df$left_right[i]
        lri_nolast <- substring(lri,1,(nlri-1))
        
        # If duplicate strings (e.g., LLR and LLL) then only do the first one
        if (lri_nolast == lri_nolast_prev) {
          next
        }
        
        lri_nolast_prev <- lri_nolast
        
        ii <- ii+1
        
        # Find the if condition and node_bias for lri_nolast using ctreesum_df
        ifcondr_lri_nolast    <- ctreesum_df$ifcondr[ctreesum_df$left_right ==
                                                     lri_nolast]
        ifcondq_lri_nolast    <- ctreesum_df$ifcondq[ctreesum_df$left_right ==
                                                       lri_nolast]
        ifcondrt_lri_nolast   <- ctreesum_df$ifcondrt[ctreesum_df$left_right ==
                                                         lri_nolast]
        node_bias_lri_nolast  <- ctreesum_df$node_bias[ctreesum_df$left_right == 
                                                        lri_nolast]
        node_sampr_lri_nolast <- ctreesum_df$node_sampr[ctreesum_df$left_right == 
                                                        lri_nolast]
        node_sampq_lri_nolast <- ctreesum_df$node_sampq[ctreesum_df$left_right == 
                                                        lri_nolast]
        node_bss_lri_nolast   <- ctreesum_df$node_bss[ctreesum_df$left_right == 
                                                       lri_nolast]
        
        # Find all leaves whose first elements equal lri_nolast
        # and set their left_right equal to lri_last,
        # their if condition to ifcondr_lri_nolast, and their bias to
        # node_bias_lri_nolast
        
        first_lr_chars <- substring(tsl_prune_df$left_right,1,(nlri-1)) 
        
        tsl_prune_df$left_right <- ifelse(first_lr_chars == lri_nolast, 
                                          lri_nolast, 
                                          tsl_prune_df$left_right)
        tsl_prune_df$ifcondr    <- ifelse(first_lr_chars == lri_nolast, 
                                      ifcondr_lri_nolast, 
                                      tsl_prune_df$ifcondr)
        tsl_prune_df$ifcondq    <- ifelse(first_lr_chars == lri_nolast, 
                                          ifcondq_lri_nolast, 
                                          tsl_prune_df$ifcondq)
        tsl_prune_df$ifcondrt   <- ifelse(first_lr_chars == lri_nolast, 
                                          ifcondrt_lri_nolast, 
                                          tsl_prune_df$ifcondrt)
        
        
        tsl_prune_df$node_bias  <- ifelse(first_lr_chars == lri_nolast,
                                         node_bias_lri_nolast,
                                         tsl_prune_df$node_bias)
        tsl_prune_df$node_sampr <- ifelse(first_lr_chars == lri_nolast,
                                         node_sampr_lri_nolast,
                                         tsl_prune_df$node_sampr)
        tsl_prune_df$node_sampq <- ifelse(first_lr_chars == lri_nolast,
                                          node_sampq_lri_nolast,
                                          tsl_prune_df$node_sampq)
        tsl_prune_df$node_bss   <- ifelse(first_lr_chars == lri_nolast,
                                         node_bss_lri_nolast,
                                         tsl_prune_df$node_bss)
        
        # Now do this for the data frame with leaves and nodes
        first_lr_all_chars <- substring(ts_prune_df$left_right,1,(nlri-1)) 
        
        ts_prune_df$left_right <- ifelse(first_lr_all_chars == lri_nolast, 
                                         lri_nolast, 
                                         ts_prune_df$left_right)
        ts_prune_df$ifcondr    <- ifelse(first_lr_all_chars == lri_nolast, 
                                     ifcondr_lri_nolast, 
                                     ts_prune_df$ifcondr)
        ts_prune_df$ifcondq    <- ifelse(first_lr_all_chars == lri_nolast, 
                                         ifcondq_lri_nolast, 
                                         ts_prune_df$ifcondq)
        ts_prune_df$ifcondrt   <- ifelse(first_lr_all_chars == lri_nolast, 
                                         ifcondrt_lri_nolast, 
                                         ts_prune_df$ifcondrt)
        ts_prune_df$node_bias  <- ifelse(first_lr_all_chars == lri_nolast,
                                        node_bias_lri_nolast,
                                        ts_prune_df$node_bias)
        ts_prune_df$node_sampr <- ifelse(first_lr_all_chars == lri_nolast,
                                        node_sampr_lri_nolast,
                                        ts_prune_df$node_sampr)
        ts_prune_df$node_sampq <- ifelse(first_lr_all_chars == lri_nolast,
                                        node_sampq_lri_nolast,
                                        ts_prune_df$node_sampq)
        ts_prune_df$node_bss   <- ifelse(first_lr_all_chars == lri_nolast,
                                        node_bss_lri_nolast,
                                        ts_prune_df$node_bss)
        
        # Delete affected nodes from ts_prune_df
        del_node <- ifelse(((first_lr_all_chars == lri_nolast) & 
                              (ts_prune_df$node_type == 'node')),1,0) 
        ts_prune_df <- ts_prune_df[del_node == 0,]
        
        # Find unique rows in pruned trees to deal with duplicate left_right
        tsl_prune_df <- tsl_prune_df[!duplicated(tsl_prune_df$left_right),]
        ts_prune_df  <- ts_prune_df[!duplicated(ts_prune_df$left_right),]
        
        # Calculate BSS (Old RSS) for pruned tree using tsl_prune_df
        
        # Number of rows in new pruned tree with leaves
        n_tsl_prune <- nrow(tsl_prune_df)
        
        # Calculate pruned between (old residual) sum of squares
        bss_prune <- sum(tsl_prune_df$node_bss)
        
        # Calculate alpha and save it in alpha_df
        alpha <- (bss_curr - bss_prune) / (nleaves - n_tsl_prune)
        #alpha <- alpha/tss_bias
        
        numb_leaves <- nrow(tsl_prune_df)
        
        #print(paste(cw,lri,alpha,bss_prune,bss_curr))
        #print(paste(alpha,min_alpha,numb_leaves,min_leaves))
        
        #if (is.na(alpha)) {alpha <- 1000000}
        
        # Update the minimum alpha
        if ((alpha < min_alpha) | ((alpha == min_alpha) & 
                                   (numb_leaves < min_leaves))) {
          min_alpha  <- alpha
          node_chop  <- lri_nolast
          min_bss    <- bss_prune
          min_leaves <- numb_leaves
          tsl_prune_best_df <- tsl_prune_df
          ts_prune_best_df  <- ts_prune_df
        }
        
      } # if nlri > 1
      
    } # for i over leaves in tree
    
    ts_prune_best_df$subtree <- cw+1
    
    subtrees_df <- rbind(subtrees_df,ts_prune_best_df)
    
    # Reset nleaves, bss (old rss), and current pruned trees
    nleaves  <- nrow(tsl_prune_best_df)
    bss_curr <- min_bss
    
    ctreesum_leaf_df <- tsl_prune_best_df
    ctreesum_df      <- ts_prune_best_df
    
    cp_sum$subtree[cw+1]   <- (cw+1)
    cp_sum$cp[cw+1]        <- min_alpha
    cp_sum$node_chop[cw+1] <- node_chop
    cp_sum$nleaves[cw+1]   <- min_leaves
    cp_sum$nsplits[cw+1]   <- min_leaves - 1
    
  } # end While
  
  # Deal with last tree with only root, L, and R only
  # First define ts_prune_best_df and cp_sum when nleft_right == 3
  
  if (nleft_right == 3) {
    
    cp_sum <- data.frame(matrix(nrow = 2, ncol = 5))
    colnames(cp_sum) <- c("subtree","cp", "node_chop", "nleaves","nsplits")
    
    cp_sum$subtree[1]      <- 1
    cp_sum$cp[1]           <- 0
    cp_sum$node_chop[1]    <- c("None")
    cp_sum$nleaves[1]      <- 2
    cp_sum$nsplits[1]      <- 1
    
    subtrees_df <- treesum_df
    subtrees_df$subtree <- 1
    
    ts_prune_best_df  <- subtrees_df
    tsl_prune_best_df <- ts_prune_best_df[ts_prune_best_df$node_type == 'leaf',]
  }
  
  # Now check if the best root, L, R tree needs to be pruned to the root only
  if (nrow(ts_prune_best_df) == 3) {
    
    bss_prune <- ts_prune_best_df$node_bss[ts_prune_best_df$left_right == 'root']
    bss_curr  <- sum(tsl_prune_best_df$node_bss)
    alpha     <- (bss_curr - bss_prune) / (2-1)
    #alpha     <- alpha/tss_bias
    
    numb_leaves <- 1
    
    ncp <- nrow(cp_sum[!is.na(cp_sum$cp),])
    
    cp_sum$subtree[ncp+1]   <- ncp+1
    cp_sum$cp[ncp+1]        <- alpha
    cp_sum$node_chop[ncp+1] <- c("None")
    cp_sum$nleaves[ncp+1]   <- numb_leaves
    cp_sum$nsplits[ncp+1]   <- numb_leaves - 1
    
    root_info <- subtrees_df[(subtrees_df$subtree == ncp & 
                                subtrees_df$left_right == "root"),]
    
    root_info$subtree <- ncp+1
    root_info$split_var <- ""
    root_info$split_val <- ""
    
    subtrees_df <- rbind(subtrees_df,root_info)
    
    # Calculate nlr in subtrees_df
    subtrees_df$nlr      <- nchar(trimws(subtrees_df$left_right))
    
  } # if last tree root
  
  # Assign better varable names to subtrees_df for output 
  colnames(subtrees_df)[colnames(subtrees_df) == "node_bias"] <- "pred_bias"
  colnames(subtrees_df)[colnames(subtrees_df) == "node_sw"]   <- "eff_samp"
  colnames(subtrees_df)[colnames(subtrees_df) == "node_sampr"]  <- "sampr"
  colnames(subtrees_df)[colnames(subtrees_df) == "node_sampq"]  <- "sampq"
  colnames(subtrees_df)[colnames(subtrees_df) == "subtree"] <- "subtree_num"
  
  # Calculate best tree using holdout sample using holdout_analysis function
  # and recalculate biases using the full sample if holdout == 1
  
  if (holdout == 1) {
    
    subtrees_full_df <- subtrees_full(rct_dat_fullc, qed_dat_fullc,
                                      rct_dat_fullt, treeinfo_df=subtrees_df)
    subtrees_df <- subtrees_full_df
    
    best_subtree_hold <- holdout_analysis(hold_rct_datc, hold_qed_datc, 
                                          hold_rct_datt,
                                          treeinfo_df = subtrees_df)
    
    best_subtree_num <- best_subtree_hold$best_subtree
    #best_subtree_sqerr <- best_subtree_hold$min_sq_diff
    
    # Print cp information
    
    # Remove NA values in cp_sum and also remove only the root node since
    # the alpha value does not make sense for this since 
    # total sum of squares is unobserved
    
    cp_sum  <- cp_sum[!is.na(cp_sum$cp),]
    cp_suma <- cp_sum[-length(cp_sum$subtree),]
    #cp_suma <- cp_sum
    
    cp_suma$pruned_node <- cp_suma$node_chop
    
    ncp_trees <- nrow(cp_suma)
    
    cp_suma$holdout_tree <- character(ncp_trees)
    cp_suma$holdout_tree[best_subtree_num] <- c("*")
    cp_suma$subtree_num <- cp_suma$subtree
    
    sink(file = out_txt, split = TRUE, append = TRUE)
    
    cat("\n")
    if (cace_itt == 0) {
      pr1 <- "TABLE 3. BIAS_TREE FUNCTION: CP VALUES, SIZE OF NESTED SUBTREES, AND"
      pr2 <- "PRUNED NODES FROM WEAKEST LINK PRUNING"
    } else if (cace_itt == 1) {
      pr1 <- "TABLE 4. BIAS_TREE FUNCTION: CP VALUES, SIZE OF NESTED SUBTREES, AND"
      pr2 <- "PRUNED NODES FROM WEAKEST LINK PRUNING"
    }
    cat(pr1,"\n",pr2,"\n") 
    #print("BIAS_TREE FUNCTION: CP VALUES, SIZE OF NESTED SUBTREES, AND PRUNED NODES FROM WEAKEST LINK PRUNING",
    #      row.names = FALSE)
    print(cp_suma[,c("subtree_num","cp","nleaves","nsplits","pruned_node",
                     "holdout_tree")], row.names = FALSE, right='F')
    
    strp1a <- "Notes: The pruned nodes are cumulative moving from one nested" 
    strp2a <- "subtree to the next. The holdout_tree variable displays the selected"
    strp3a <- "subtree based on the 30 percent holdout sample that minimizes the "
    strp4a <- "squared difference between the actual bias and predicted bias from"
    strp5a <- "each subtree. The cp (complexity parameter) values are not scaled."
    
    cat(strp1a,"\n",strp2a,"\n",strp3a,"\n",strp4a,"\n",strp5a,"\n")
    
    sink()
    
  } else if (holdout == 0) {
    
    #subtrees_df$pred_bias_full <- subtrees_df$pred_bias
    #subtrees_df$sampr_full     <- subtrees_df$sampr
    #subtrees_df$sampq_full     <- subtrees_df$sampq
    
    cp_sum  <- cp_sum[!is.na(cp_sum$cp),]
    cp_suma <- cp_sum[-length(cp_sum$subtree),]
    
    cp_suma$pruned_node <- cp_suma$node_chop
    
    ncp_trees <- nrow(cp_suma)
    cp_suma$holdout_tree <- character(ncp_trees)
    
    cp_suma$subtree_num <- cp_suma$subtree
    
    sink(file = out_txt, split = TRUE, append = TRUE)
    
    cat("\n")
    if (cace_itt == 0) {
      pr1 <- "TABLE 3. BIAS_TREE FUNCTION: CP VALUES, SIZE OF NESTED SUBTREES, AND"
      pr2 <- "PRUNED NODES FROM WEAKEST LINK PRUNING"
    } else if (cace_itt == 1) {
      pr1 <- "TABLE 4. BIAS_TREE FUNCTION: CP VALUES, SIZE OF NESTED SUBTREES, AND"
      pr2 <- "PRUNED NODES FROM WEAKEST LINK PRUNING"
    }
    #pr1 <- "BIAS_TREE FUNCTION: CP VALUES, SIZE OF NESTED SUBTREES, AND PRUNED NODES"
    #pr2 <- "FROM WEAKEST LINK PRUNING"
    cat(pr1,"\n",pr2,"\n") 
    #print("BIAS_TREE FUNCTION: CP VALUES, SIZE OF NESTED SUBTREES, AND PRUNED NODES FROM WEAKEST LINK PRUNING",
    #      row.names = FALSE)
    print(cp_suma[,c("subtree_num","cp","nleaves","nsplits","pruned_node")], 
          row.names = FALSE, right='F')
    
    strp1a <- "Notes: The pruned nodes are cumulative moving from one nested subtree" 
    strp2a <- "to the next. The cp (complexity parameter) values are not scaled."
    
    cat(strp1a,"\n",strp2a,"\n")
    
    sink()
  }
  
} else if (nleft_right == 1) {
  
  # This occurs if the initial tree has only the root which is a leaf
  
  # Create treesum_df from tree_results
  treesum_df <- data.frame(left_right)
  treesum_df$left_rightp <- format(treesum_df$left_right,justify = "left")
  
  treesum_df$node_type <- "leaf"
  
  treesum_df$split_var  <- NA
  treesum_df$split_val  <- NA
  treesum_df$split_valp <- NA
  
  tree_results <- tree_results[tree_results != -1]
  
  prediction <- ifelse(tree_results == 'leaf',lead(tree_results),"")
  nsampr     <- ifelse(tree_results == 'leaf',lead(tree_results,n=2),"")
  nsampq     <- ifelse(tree_results == 'leaf',lead(tree_results,n=3),"")
  
  treesum_df$node_bias   <- as.numeric(prediction[1])
  treesum_df$node_sampr  <- nsampr[1]
  treesum_df$node_sampq  <- nsampq[1]
  treesum_df$node_biasp  <- format(round(treesum_df$node_bias,4), 
                                   nsmall=4) 
  
  treesum_df$node_n  <- seq_len(nrow(treesum_df))
  treesum_df$node_id <- format(treesum_df$node_n, justify = "left")
  
  # Assign better variable names
  df_temp <- treesum_df
  
  df_temp$left_right     <- df_temp$left_rightp
  df_temp$pred_bias      <- df_temp$node_biasp
  df_temp$sampr          <- df_temp$node_sampr
  df_temp$sampq          <- df_temp$node_sampq
  df_temp$split_val      <- df_temp$split_valp
  
  sink(file = out_txt, split = TRUE, append = TRUE)
  
  cat("\n")
  print("BIAS_TREE FUNCTION: INITIAL UNPRUNED TREE",row.names = FALSE)
  print(df_temp[,c("node_id","left_right",
                   "node_type","pred_bias","sampr","sampq", 
                   "split_var", "split_val")],row.names = FALSE,
        right='F')
  sink()
  
  # Create subtree_df and ts_prune_best_df
  subtrees_df           <- treesum_df
  subtrees_df$subtree   <- 1
  subtrees_df$eff_samp  <- NA
  subtrees_df$ifcondr   <- NA
  subtrees_df$ifcondq   <- NA
  
  colnames(subtrees_df)[colnames(subtrees_df) == "node_bias"] <- "pred_bias"
  colnames(subtrees_df)[colnames(subtrees_df) == "node_sampr"]  <- "sampr"
  colnames(subtrees_df)[colnames(subtrees_df) == "node_sampq"]  <- "sampq"
  colnames(subtrees_df)[colnames(subtrees_df) == "subtree"] <- "subtree_num"
  
  ts_prune_best_df  <- subtrees_df
  tsl_prune_best_df <- ts_prune_best_df[ts_prune_best_df$node_type == 'leaf',]
  
  #print(str(subtrees_df))
  
} # end of nleft_right conditions

# Renumber the nodes
# Create a node sequence within each subtree and name it node_num
subtrees_df$node_num <- with(subtrees_df, ave(subtree_num,
                                              subtree_num, FUN = seq_along))

if (holdout == 0) {
  subtrees_df$pred_bias_full <- subtrees_df$pred_bias
  subtrees_df$sampr_full     <- subtrees_df$sampr
  subtrees_df$sampq_full     <- subtrees_df$sampq
  subtrees_df$eff_samp_full  <- subtrees_df$eff_samp
}

#print(hold_out)
#print(cace_itt)

if (as.numeric(hold_out) == 1) {
  subtrees_df$holdout  <- 1
} else {
  subtrees_df$holdout  <- 0
}

if (as.numeric(cace_itt) == 1) {
  subtrees_df$cace_itt  <- 1
} else {
  subtrees_df$cace_itt  <- 0
}

subtreesg_df <- subtrees_df[,c("subtree_num","node_num","left_right","node_type",
                               "pred_bias", "sampr", "sampq", "eff_samp", "nlr", 
                               "split_var", "split_val","ifcondr", "ifcondq",
                               "pred_bias_full","sampr_full","sampq_full",
                               "eff_samp_full","holdout","cace_itt")]

#sink()

return(subtreesg_df)

}  # End bias_tree function

# Function to calculate biases

calc_bias <- function(yr, yq, yrt, ssr, ssq, ssrt, ppr, pprt, wwq) {
  
  rct_datcfg <- data.frame(yr,ssr,ppr)
  qed_datcfg <- data.frame(yq,ssq,wwq)
  rct_dattfg <- data.frame(yrt,ssrt,pprt)
  
  # Create data sets with RCT control and treatment non-participants 
  rct_datc_npfg <- rct_datcfg[rct_datcfg$ppr == 0,]
  rct_datt_npfg <- rct_dattfg[rct_dattfg$pprt == 0,]
  
  se1 <- max(ssr) 
  
  # Initialize and do over sites
  bias_site <- rep(0,se1)
  wgt_site  <- rep(0,se1)
  badccc <- 0
  
  for (j in 1:se1) {
    
    yrcs  <- rct_datcfg$yr[rct_datcfg$ssr==j]
    
    yqcs  <- qed_datcfg$yq[qed_datcfg$ssq==j]
    wwqs  <- qed_datcfg$wwq[qed_datcfg$ssq==j]
    
    nrcs  <- length(yrcs)
    nqcs  <- length(yqcs)
    
    # QED comparison mean in site using ipw weights
    meanqcs <- sum(wwqs*yqcs) / sum(wwqs)
    
    # QED design effect
    wgt_term  <- sum(wwqs^2)/((sum(wwqs))^2)  
    deffqcs   <- wgt_term*nqcs
    
    if (cace_itt == 0) {
      
      # RCT control mean in site
      meanrcs <- mean(yrcs)
      
      # Bias in site
      bias_site[j] <- meanrcs - meanqcs
      
      # Weights for pooling across sites
      wgt_site[j] <- nrcs + (nqcs/deffqcs)
      
    } else if (cace_itt == 1) {
      
      # Calculate site participation rate and complier proportion
      drts  <- rct_dattfg$pprt[rct_dattfg$ssrt==j]
      drcs  <- rct_datcfg$ppr[rct_datcfg$ssr==j]
      
      part_ratets <- mean(drts)
      part_ratecs <- mean(drcs)
      
      compl_rates <- part_ratets - part_ratecs
      
      # Pull off y for RCT control nonparticipants and calculate mean
      yrc_nps    <- rct_datc_npfg$yr[rct_datc_npfg$ssr==j]
      meanrc_nps <- mean(yrc_nps)
      
      # Pull off y for RCT treatment nonparticipants and calculate mean
      # if there are at least 2 people - if not set the mean to 0
      yrt_nps <- rct_datt_npfg$yrt[rct_datt_npfg$ssrt==j]
      nrt_nps <- length(yrt_nps)
      
      if (nrt_nps >= 2) {
        meanrt_nps <- mean(yrt_nps)
      } else {
        meanrt_nps <- 0
      }
      
      # Calculate the estimated RCT control complier mean in site
      meanrccs1 <- ((1-part_ratecs)*meanrc_nps) - ((1-part_ratets)*meanrt_nps)
      
      # For some reason the screen for complier rates is not full-proof
      # so adding this check and fix to avoid the program bombing
      
      #meanrccs  <-  meanrccs1 / compl_rates
      
      if (compl_rates > 0) {
        meanrccs  <-  meanrccs1 / compl_rates
      } else {
        badccc <- 1
        meanrccs  <-  meanrccs1 / compl_rate_tot
      }
      
      # Bias in site
      bias_site[j] <- meanrccs - meanqcs
      
      if ((bias_site[j] == Inf) | (is.na(bias_site[j]))) {
        #print(paste(compl_rates,length(drts),length(drcs),length(yrcs),compl_rate_tot))
        #print(paste(length(yr),length(yq)))
        junk <- 1
      }
      
      # Weights for pooling across sites
      wgt_site[j] <- (nrcs*compl_rates) + (nqcs/deffqcs)
      
    } # end if cace_itt
        
  } # for sites
  
  # Calculate the pooled bias across sites
  bias_node    <- sum(wgt_site*bias_site)/sum(wgt_site)
  
  # Calculate sum of weights to calculate CART score
  sw <- sum(wgt_site)
  
  bn <- data.frame(bias_node,sw)
  
  return(bn)
  
} # End of Calc Bias function


###
# Function to update subtrees_df to calculate biases using the FULL SAMPLE
# if holdout == 1 
###

subtrees_full <- function(rct_dat_fullc_df, qed_dat_full_df, rct_dat_fullt_df,
                          treeinfo_df=subtrees_df) {
  
  #treeinfol_df <- treeinfo_df[treeinfo_df$node_type == "leaf",]
  treeinfof_df <- data.frame(treeinfo_df)
  
  dat_noder  <- rct_dat_fullc_df
  dat_nodeq  <- qed_dat_full_df
  dat_nodert <- rct_dat_fullt_df
  
  # Prepare data for the calc_bias macro
  
  # Number of subtrees
  nsubtrees <- max(treeinfof_df$subtree_num)
  
  #print(nsubtrees)
  #print(treeinfof_df$subtree_num)
  
  # Run through each subtree
  
  cnt_full <- 1
  for (i in 1:nsubtrees) {
    
    tinfo <- treeinfof_df[treeinfof_df$subtree_num == i,]
    tree_leaves <- nrow(tinfo)
    
    #tz <- tinfo[,c("left_right","node_type","ifcondr")]
    #print(tz)
    
    if (tree_leaves > 1) {
      
      for (k in 1:tree_leaves) {
        
        if (tinfo$left_right[k] != "root") {
          
          #Pull off if conditions for node k
          ifnoder  <- tinfo$ifcondr[k]
          ifnodeq  <- tinfo$ifcondq[k]
          ifnodert <- tinfo$ifcondrt[k]
          
          # Define samp_in_node which subsets the observations in dat_node
          # who lie in node i
          
          ifparser    <- paste0("samp_in_noder <- ", ifnoder)
          parse_coder <- parse(text=ifparser)
          eval(parse_coder)
          
          ifparseq    <- paste0("samp_in_nodeq <- ", ifnodeq)
          parse_codeq <- parse(text=ifparseq)
          eval(parse_codeq)
          
          ifparsert    <- paste0("samp_in_nodert <- ", ifnodert)
          parse_codert <- parse(text=ifparsert)
          eval(parse_codert)
          
          # Define inputs to calc_bias
          yra   <- dat_noder$yv[samp_in_noder]  
          yqa   <- dat_nodeq$yv[samp_in_nodeq] 
          yrta  <- dat_nodert$yv[samp_in_nodert]  
          ssra  <- dat_noder$siteID[samp_in_noder]
          ssqa  <- dat_nodeq$siteID[samp_in_nodeq]
          ssrta <- dat_nodert$siteID[samp_in_nodert]
          ppra  <- dat_noder$partic[samp_in_noder]
          pprta <- dat_nodert$partic[samp_in_nodert]
          wwqa  <- dat_nodeq$ipw_wt[samp_in_nodeq]
          
        } else if (tinfo$left_right[k] == "root") {
          
          yra   <- dat_noder$yv  
          yqa   <- dat_nodeq$yv 
          yrta  <- dat_nodert$yv  
          ssra  <- dat_noder$siteID
          ssqa  <- dat_nodeq$siteID
          ssrta <- dat_nodert$siteID
          ppra  <- dat_noder$partic
          pprta <- dat_nodert$partic
          wwqa  <- dat_nodeq$ipw_wt
        } # end if root
                
        bias_full  <- calc_bias(yra, yqa, yrta, ssra, ssqa, ssrta, 
                                ppra, pprta, wwqa)
        
        treeinfof_df$pred_bias_full[cnt_full] <- bias_full$bias_node
        treeinfof_df$eff_samp_full[cnt_full]  <- bias_full$sw
        treeinfof_df$sampr_full[cnt_full] <- NROW(yra)
        treeinfof_df$sampq_full[cnt_full] <- NROW(yqa)
        
        cnt_full <- cnt_full + 1
        
      } #for k
      
    } else if (tree_leaves == 1) {
      yra   <- dat_noder$yv  
      yqa   <- dat_nodeq$yv 
      yrta  <- dat_nodert$yv  
      ssra  <- dat_noder$siteID
      ssqa  <- dat_nodeq$siteID
      ssrta <- dat_nodert$siteID
      ppra  <- dat_noder$partic
      pprta <- dat_nodert$partic
      wwqa  <- dat_nodeq$ipw_wt
      
      bias_full  <- calc_bias(yra, yqa, yrta, ssra, ssqa, ssrta, 
                              ppra, pprta, wwqa)
      
      treeinfof_df$pred_bias_full[cnt_full] <- bias_full$bias_node
      treeinfof_df$sampr_full[cnt_full] <- NROW(yra)
      treeinfof_df$sampq_full[cnt_full] <- NROW(yqa)
      
      cnt_full <- cnt_full + 1
      
    } # end if tree_leaves
    
  } # for i over subtrees
  
  #ttt <- treeinfof_df[,c("subtree_num","left_right","node_type",
  #                       "pred_bias","pred_bias_full",
  #                       "sampr","sampq","sampr_full","sampq_full")]
  #print("TTT INFO")
  #print(ttt)
  
return(treeinfof_df)
  
} # end subtrees_full function 


#
# Function to Merge Predicted Biases onto Holdout Datasets and find
# the subtree with the smallest squared difference between actual and 
# predicted biases. 

holdout_analysis <- function(rct_datc_dfh, qed_dat_dfh, rct_datt_dfh,
                             treeinfo_df=subtrees_df) {
  
  treeinfol_df <- treeinfo_df[treeinfo_df$node_type == "leaf",]
  
  #print(treeinfol_df)
  
  dat_nodera  <- rct_datc_dfh
  dat_nodeqa  <- qed_dat_dfh
  dat_noderta <- rct_datt_dfh
  
  # First calculate actual bias in holdout sample.  Do this by site to 
  # get the weights needed for the predicted biases later
  
  # Prepare data for the calc_bias macro
  dat_nodera$site_temp <- 1
  dat_nodeqa$site_temp <- 1
  dat_noderta$site_temp <- 1
  
  actual_site_bias   <- rep(0,se)
  site_wtz           <- rep(0,se)
  for (ss in 1:se) {
  
    # Calculate actual bias in holdout sample
    yra   <- dat_nodera$yv[dat_nodera$siteID == ss]  
    yqa   <- dat_nodeqa$yv[dat_nodeqa$siteID == ss] 
    yrta  <- dat_noderta$yv[dat_noderta$siteID == ss]  
    ssra  <- dat_nodera$siteID[dat_nodera$siteID == ss]
    ssqa  <- dat_nodeqa$siteID[dat_nodeqa$siteID == ss]
    ssrta <- dat_noderta$siteID[dat_noderta$siteID == ss]
    ppra  <- dat_nodera$partic[dat_nodera$siteID == ss]
    pprta <- dat_noderta$partic[dat_noderta$siteID == ss]
    wwqa  <- dat_nodeqa$ipw_wt[dat_nodeqa$siteID == ss]
    sssra <- dat_nodera$site_temp[dat_nodera$siteID == ss]
    sssqa <- dat_nodeqa$site_temp[dat_nodeqa$siteID == ss]
    sssrta <- dat_noderta$site_temp[dat_noderta$siteID == ss]
    
    bias_res  <- calc_bias(yra, yqa, yrta, sssra, sssqa, sssrta, ppra, pprta, wwqa)
    actual_site_bias[ss]  <- bias_res$bias_node
    site_wtz[ss]          <- bias_res$sw
    
  } # end ss
  
  # Now calculate predicted biases for the QED sample in each subtree
  
  # Number of subtrees
  nsubtrees <- max(treeinfol_df$subtree_num)
  
  #print(nsubtrees)
  
  # Run through each subtree, merge on predicted bias for the
  # QED holdout sample, and calculate squared bias
  min_sq_diff <- Inf
  
  for (i in 1:nsubtrees) {
    
    # Initialize data sets
    dat_nodeq <- qed_dat_dfh
    dat_nodeq$left_right <- NA
    dat_nodeq$pred_bias  <- NA
    
    tinfo <- treeinfol_df[treeinfol_df$subtree_num == i,]
    tree_leaves <- nrow(tinfo)
    
    if (tree_leaves > 1) {
      
      for (k in 1:tree_leaves) {
        
        # Pull off if conditions for node k
        ifnodeq <- tinfo$ifcondq[k]
        
        # Define samp_in_node which subsets the observations in dat_node
        # who lie in node i
        ifparseq   <- paste0("samp_in_node <- ", ifnodeq)
        parse_codeq <- parse(text=ifparseq)
        eval(parse_codeq)
        
        # Add left_right and its predicted bias to dat_nodeq dataset
        #dat_nodeq$left_right[samp_in_node]    <- tinfo$left_right[k]
        
        dat_nodeq$pred_bias[samp_in_node] <- tinfo$pred_bias[k]

      } # for k
      
    } else if (tree_leaves <= 1) {
      dat_nodeq$pred_bias <- treeinfo_df$pred_bias[1]
    }
    
    # Calculate predicted bias by site for QED sample 
    qed_pred_site_bias <- rep(0,se)
    
    for (ss in 1:se) {
      
      # Calculate predicted bias for QED sample by site
      bias_ss <- dat_nodeq$pred_bias[dat_nodeq$siteID == ss]
      ipw_ss  <- dat_nodeq$ipw_wt[dat_nodeq$siteID == ss]
      
      qed_pred_site_bias[ss] <- sum(ipw_ss*bias_ss)/sum(ipw_ss)
      
    } # end ss
    
    #print(paste("SUBTREE",i))
    #print("Predicted QED Biases by RCT-QED site")
    #print(qed_pred_site_bias)
    
    # Calculate overall predicted and actual biases across sites
    qed_pred_bias <- sum(site_wtz*qed_pred_site_bias)/sum(site_wtz)
    actual_bias   <- sum(site_wtz*actual_site_bias)/sum(site_wtz)
    
    # Calculate minimum (actual - predicted)^2
    sq_diff <- (qed_pred_bias-actual_bias)^2
    
    if (sq_diff <= min_sq_diff) {
      min_sq_diff <- sq_diff
      best_subtree <- i
    }
    
    #print(paste ("PREDICTED AND ACTUAL BIAS AND SQ DIFF:",qed_pred_bias,
    #             actual_bias,sq_diff))
    #print(paste("MIN_SQ_DIFF and BEST_SUBTREE",min_sq_diff,best_subtree))
    
  } # end i over subtrees 
  
  #print(paste("FINAL MIN_SQ_DIFF and BEST_SUBTREE",min_sq_diff,best_subtree))
  
  holdout_info <- data.frame(best_subtree,min_sq_diff)
  
  return(holdout_info)
}

#
# Function to print and graph selected subtree
#

view_selected_tree <- function(subtree_num = NA, treeinfo_df, out_view) {
  
  out_txt1 <- deparse(substitute(out_view))   
  
  if (out_txt1 == "") {
    out_txt1 <- c('view_tree_results.txt')
  }
  
  # Define sink function to output file and also writes to the console
  sink(file = out_txt1, split = TRUE)
  
  blank <- c("")
  
  cat("\n")
  if (is.na(subtree_num)) {
    return(print("ERROR in the view_selected_tree function: No value provided for the required subtree_num input"))
  }
  
  # Check if subtree_num is valid
  sn    <- subtree_num
  sp_df <- treeinfo_df[treeinfo_df$subtree_num == sn,]
  
  if (nrow(sp_df) == 0) {
    return(print("ERROR in the view_selected_tree function: Invalid value provided for the required subtree_num input"))
  }
  
  # Format and print selected tree
  sp_df$pred_bias <-  format(round(sp_df$pred_bias,4),nsmall=4)
  sp_df$eff_samp  <-  format(round(as.numeric(sp_df$eff_samp),0)) 
  
  sp_df$pred_bias_full <-  format(round(sp_df$pred_bias_full,4),nsmall=4)
  sp_df$eff_samp_full  <-  format(round(as.numeric(sp_df$eff_samp_full),0)) 
  
  #print(sp_df$holdout)
  #print(sp_df$cace_itt)
  
  if (any(sp_df$holdout == 1)) {
    hh <- 1
  } else {hh <- 0}
  
  if (any(sp_df$cace_itt == 1)) {
    cace <- 1
  } else {cace <- 0}
  
  #print(paste("HH",hh))
  #print(paste("CACE",cace))
  
  if (hh == 1) {
    
    #cat("\n")
    if (cace == 0) {
      tit <- sprintf("TABLE 1. VIEW_TREE FUNC: INFO ON PRUNED SUBTREE %s FOR THE ITT ANALYSIS \n USING THE TRAINING SAMPLE",sn)
      cat(blank,tit, sep="\n")
    } else if (cace == 1) {
      tit <- sprintf("TABLE 1. VIEW_TREE FUNC: INFO ON PRUNED SUBTREE %s FOR THE CACE ANALYSIS \n USING THE TRAINING SAMPLE",sn)
      cat(blank,tit, sep="\n")
    }
    
    print(sp_df[,c("node_num","left_right","node_type","pred_bias","sampr","sampq",
                   "eff_samp", "split_var", "split_val")],
          row.names = FALSE, right='F')
    
    strp1 <- "Notes: The left_right variable displays the movement along the binary" 
    strp2 <- "tree branches where L signifies a left move if the split variable is "
    strp3 <- "less than the split value at the parent node, and R signifies a right "
    strp4 <- "move otherwise. The pred_bias variable is the node bias, node_type "
    strp5 <- "signifies whether the node is an internal tree node or a terminal leaf, "
    strp6 <- "sampr and sampq are RCT control and QED comparison sample sizes, and "
    strp7 <- "eff_samp is the total effective sample size accounting for design"
    strp8 <- "effects due to inverse probability weighting."
    cat(strp1,"\n",strp2,"\n",strp3,"\n", strp4,"\n", strp5,"\n", 
        strp6,"\n",strp7,"\n",strp8,"\n") 
    
    cat("\n")
    if (cace == 0) {
      tit <- sprintf("TABLE 1A. VIEW_TREE FUNC: INFO ON PRUNED SUBTREE %s FOR THE ITT ANALYSIS \n USING THE FULL SAMPLE",sn)
      cat(blank,tit, sep="\n")
    } else if (cace == 1) {
      tit <- sprintf("TABLE 1A. VIEW_TREE FUNC: INFO ON PRUNED SUBTREE %s FOR THE CACE ANALYSIS \n USING THE FULL SAMPLE",sn)
      cat(blank,tit, sep="\n")
    }
    
    print(sp_df[,c("node_num","left_right","node_type","pred_bias_full",
                   "sampr_full","sampq_full", "eff_samp_full")],
          row.names = FALSE, right='F')
    
  } else if (hh == 0) {
    
    cat("\n")
    if (cace == 0) {
      tit <- sprintf("TABLE 1. VIEW_TREE FUNC: INFO ON PRUNED SUBTREE %s FOR THE ITT ANALYSIS",sn)
      cat(blank,tit, sep="\n")
      
    } else if (cace ==1) {
      tit <- sprintf("TABLE 1. VIEW_TREE FUNC: INFO ON PRUNED SUBTREE %s FOR THE CACE ANALYSIS",sn)
      cat(blank,tit, sep="\n")
    }
    
    print(sp_df[,c("node_num","left_right","node_type","pred_bias","sampr","sampq",
                   "eff_samp", "split_var", "split_val")],
          row.names = FALSE, right='F')
    
    strp1 <- "Notes: The left_right variable displays the movement along the binary" 
    strp2 <- "tree branches where L signifies a left move if the split variable is "
    strp3 <- "less than the split value at the parent node, and R signifies a right "
    strp4 <- "move otherwise. The pred_bias variable is the node bias, node_type "
    strp5 <- "signifies whether the node is an internal tree node or a terminal leaf, "
    strp6 <- "sampr and sampq are RCT control and QED comparison sample sizes, and "
    strp7 <- "eff_samp is the total effective sample size accounting for design"
    strp8 <- "effects due to inverse probability weighting."
    cat(strp1,"\n",strp2,"\n",strp3,"\n", strp4,"\n", strp5,"\n", 
        strp6,"\n",strp7,"\n",strp8,"\n")
    
  }
  
  #
  # Visualize tree
  #
  
  sp_df$pred_bias <- sp_df$pred_bias_full
  sp_df$sampr     <- sp_df$sampr_full
  sp_df$sampq     <- sp_df$sampq_full
  sp_df$eff_samp  <- sp_df$eff_samp_full
  
  ntree_nodesz  <- nrow(sp_df)
  
  if (ntree_nodesz > 3) {
    
    # Initialize tree structure at Root for data.tree function and set the 
    # split vars and values for the L and R left-right entries equal to the 
    # Root node values and also initialize the plot labels
    
    tree_p <- Node$new("Root")  
    sp_df$plot_lab <- "None"
    
    root_svar  <- sp_df$split_var[sp_df$left_right == "root"]
    root_sval  <- sp_df$split_val[sp_df$left_right == "root"]
    
    rootlab1 <- paste(root_svar,"<",root_sval,sep=" ")
    sp_df$plot_lab[sp_df$left_right == "L"] <- rootlab1
    
    rootlab2 <- paste(root_svar,">=",root_sval,sep=" ")
    sp_df$plot_lab[sp_df$left_right == "R"] <- rootlab2
    
    sp_df$obs_num <- seq_len(nrow(sp_df))
    
    lrow <- sp_df$obs_num[sp_df$left_right == 'L']
    rrow <- sp_df$obs_num[sp_df$left_right == 'R']
    
    # Define the tree structure for the L and R entries
    
    tree_temp <- paste0("tree_p",lrow-1)
    assign(tree_temp, tree_p$AddChild(rootlab1)) 
    
    tree_temp <- paste0("tree_p",rrow-1)
    assign(tree_temp, tree_p$AddChild(rootlab2)) 
    
    # Loop over left-right entries to create the tree structure and plot labels
    for (i in 1:ntree_nodesz) {
      
      i1 <- i-1
      
      # Find the substring that takes off the last element of left_right
      lri <- sp_df$left_right[i]
      
      if ((lri != "root") & (lri != "L") & (lri != "R")) {
        
        nlri        <- sp_df$nlr[i]
        lri_nolast  <- substring(lri,1,(nlri-1))
        
        # Find the last element of left_right
        lri_last  <- substring(lri,nlri,nlri)
        
        # Find the split vars and vals for the string without that last element
        temp_svar  <- sp_df$split_var[sp_df$left_right == lri_nolast]
        temp_sval  <- sp_df$split_val[sp_df$left_right == lri_nolast]
        
        #print(paste(lri,nlri, lri_nolast,lri_last,temp_svar,temp_sval))
        
        # Do over the left_right entries until it matches lri_nolast
        for (j in 1:ntree_nodesz) {
          
          lrj       <- sp_df$left_right[j]
          lrj_obs   <- sp_df$obs_num[j]
          lrj_obs1  <- lrj_obs-1
          
          if (lrj == lri_nolast) {
            
            # Create plot label
            if (lri_last == 'L') {
              temp_lab <- paste(temp_svar,"<",temp_sval,sep=" ")
            } else if (lri_last == 'R') {
              temp_lab <- paste(temp_svar,">=",temp_sval,sep=" ")
            }
            
            sp_df$plot_lab[sp_df$left_right == lri] <- temp_lab
            
            # Create tree structure parent and child nodes
            
            tree_text <- paste0("tree_p",i1, "<-","tree_p",lrj_obs1,"$AddChild(temp_lab)")
            parse_code <- parse(text=tree_text)
            eval(parse_code)
            
            #print(parse_code)
            
          } # for j
        } # for lri not root, L, or R
      } # for i
    }
    
    # Add on node predictions and sample sizes
    
    tree_p$pred_bias <- sp_df$pred_bias[1]
    tree_p$sampr     <- sp_df$sampr[1]
    tree_p$sampq     <- sp_df$sampq[1]
    tree_p$eff_samp  <- sp_df$eff_samp[1]
    
    df_list <- list()
    for (i in 1:(ntree_nodesz-1)) {
      
      df_name      <- get(paste0("tree_p",i))  # Get creates a data frame
      df_list[[i]] <- df_name
    }
    
    for (i in seq_along(df_list)) {
      
      current_df <- df_list[[i]]
      
      type_temp  <- sp_df$node_type[i+1]
      current_df$node_type  <- ifelse(type_temp == 'leaf',type_temp,"")
      current_df$pred_bias  <- sp_df$pred_bias[i+1]
      current_df$sampr      <- sp_df$sampr[i+1]
      current_df$sampq      <- sp_df$sampq[i+1]
      current_df$eff_samp   <- sp_df$eff_samp[i+1]
      
      df_list[[i]] <- current_df
    }
    
    cat("\n")
    if (cace == 0) {
      cat("\n")
      tit <- sprintf("TABLE 2. VIEW_TREE FUNC: VISUALIZATION OF PRUNED SUBTREE %s FOR THE ITT ANALYSIS \n USING THE FULL SAMPLE",sn)
      cat(blank,tit, sep="\n")
    } else if (cace == 1) {
      tit <- sprintf("TABLE 2. VIEW_TREE FUNC: VISUALIZATION OF PRUNED SUBTREE %s FOR THE CACE ANALYSIS \n USING THE FULL SAMPLE",sn)
      cat(blank,tit, sep="\n")
    }
    
    print(tree_p,"node_type","pred_bias","sampr","sampq","eff_samp",
          row.names = FALSE)
    
  } else if (ntree_nodesz <= 3) {
    cat("\n")
    print("Plots are not provided for trees with fewer than 3 leaves")
  }
  
  #return(v_tree)
  
  sink()

} # end view_selected_tree function

###
# Function to calculate predicted values for the best tree
###

calc_pred_bias <- function(pred_dat_df, subtree_num = NA, treeinfo_df, 
                           out_pred_bias) {
  
  out_txtp <- deparse(substitute(out_pred_bias))   
  
  if (out_txtp == "") {
    out_txtp <- c('pred_bias_results.txt')
  }
  
  # Define sink function to output file and also writes to the console
  sink(file = out_txtp, split = TRUE)
  
  input_df <- deparse(substitute(pred_dat_df))
  
  blank <- c("")
  
  cat("\n")
  if (is.na(subtree_num)) {
    return(print("ERROR in the calc_pred_bias function: No value provided for the required subtree_num input"))
  }
  
  # Check if subtree_num is valid
  sn    <- subtree_num
  sp_df <- treeinfo_df[treeinfo_df$subtree_num == sn,]
  
  if (nrow(sp_df) == 0) {
    return(print("ERROR in the calc_pred_bias function: Invalid value provided for the required subtree_num input"))
  }
  
  treeinfol_df <- sp_df[sp_df$node_type == "leaf",]
  
  # Now calculate predicted biases for the input data in the selected subtree
  
  # Number of rows in new pruned tree with leaves
  tree_leaves <- nrow(treeinfol_df)
  
  # Initialize the data set 
  
  dat_nodeq <- pred_dat_df
  
  dat_nodeq$pred_bias <- dat_nodeq$pred_bias_full
  
  # Parse if statements for each leaf to merge right_left and 
  # node mean onto dat_pred
  
  if (tree_leaves > 1) {
    
    for (k in 1:tree_leaves) {
      
      # Pull off if condition for node k
      ifnode <- treeinfol_df$ifcondq[k] 
      
      # Define samp_in_node which subsets the observations in dat_node
      # who lie in node i
      
      ifparse   <- paste0("samp_in_node <- ", ifnode)
      parse_code <- parse(text=ifparse)
      eval(parse_code)
      
      # Add left_right and its mean to dat_pred dataset
      
      dat_nodeq$left_right[samp_in_node] <- treeinfol_df$left_right[k]
      dat_nodeq$pred_bias[samp_in_node]  <- treeinfol_df$pred_bias[k]
      
    } # for k
    
  } else if (tree_leaves <= 1) {
    dat_nodeq$pred_bias <- treeinfo_df$pred_bias[1]
  }
  
  freq_pred <- table(dat_nodeq$pred_bias,useNA = "always")
  
  freq_pred_df <- data.frame(freq_pred)
  freq_pred_df$Predicted_bias <- freq_pred_df$Var1
  
  #cat("\n")
  tit <- sprintf("TABLE 1. CALC_PRED_BIAS FUNC: FREQS OF FULL SAMPLE PREDICTED BIASES \n USING PRUNED SUBTREE %s AND INPUT DATA %s",sn,
                 input_df)
  cat(blank,tit, sep="\n")
  print(freq_pred_df[,c("Predicted_bias","Freq")])
  
  sink()
  
  return(dat_nodeq)
  
} # end calc_pred_bias function


###
# Hybrid Impacts Function
###

hybrid_impacts <- function(rct_dat_df, qed_dat_df, rct_qed_site, yvar, 
                           xvars_rct_adj = 0, site_id, t_c, got_treat, 
                           ipw_wgt, cace_itt_est, out_impacts, inv_var_agg_wgt = 1,  
                           xvars_qed_adj = 0, cluster_id = 0) {
  
  out_txt2 <- deparse(substitute(out_impacts))   
  
  if (out_txt2 == "") {
    out_txt2 <- c('impact_results.txt')
  }
  
  # Define sink function to output file and also writes to the console
  sink(file = out_txt2, split = TRUE)
  #cat("Testing sink")
  sink()
  
  # Initialize err text and counter
  err_count <<- 0
  err_num <- seq(1:100)
  err_txt <- character(100)
  err_df  <- data.frame(err_num,err_txt)
  
  cace_itt     <<- cace_itt_est   # Makes this global
  
  # Create variables to keep for analysis for rct and qed datasets
  # These steps deparse the inputs to make them character strings
  yv     <- deparse(substitute(yvar))
  #yv    <- gsub('[(\")(\\)]', '', yv_temp) # don't need this anymore
  sitev  <- deparse(substitute(site_id))
  tc     <- deparse(substitute(t_c))
  partic <- deparse(substitute(got_treat))
  ipw_wt <- deparse(substitute(ipw_wgt))
  rct_qed_ind <- deparse(substitute(rct_qed_site))
  clusv  <- deparse(substitute(cluster_id))
  
  # Check if names exist
  # Create rct and qed datasets
  if (yv %in% names(rct_dat_df)) {
    rct_dat_df$yv <- rct_dat_df[,yv]
  } else {
    rct_dat_df$yv <- NA
  }
  
  if (sitev %in% names(rct_dat_df)) {
    rct_dat_df$sitev <- rct_dat_df[,sitev]
  } else {
    rct_dat_df$sitev <- NA
  }
  
  if (tc %in% names(rct_dat_df)) {
    rct_dat_df$tc <- rct_dat_df[,tc]
  } else {
    rct_dat_df$tc <- NA
  }
  
  if (partic %in% names(rct_dat_df)) {
    rct_dat_df$partic <- rct_dat_df[,partic]
    bad_particz <- 0
  } else {
    rct_dat_df$partic <- NA
    bad_particz <- 1
  }
  
  if ("pred_bias" %in% names(rct_dat_df)) {
    rct_dat_df$pred_bias   <- rct_dat_df[,c("pred_bias")]
  } else {
    rct_dat_df$pred_bias <- NA
  }
  
  if ("left_right" %in% names(rct_dat_df)) {
    rct_dat_df$left_right  <- rct_dat_df[,c("left_right")]
  } else {
    rct_dat_df$left_right <- NA
  }
  
  if (clusv == 0) {
    clus_dat <- 0 
    rct_dat_df$clusv <- NA
  } else {
    clus_dat <- 1
    if (clusv %in% names(rct_dat_df)) {
      rct_dat_df$clusv <- rct_dat_df[,clusv]
    } else {
      rct_dat_df$clusv <- NA
    }
  }
  
  # For ITT analysis, set participation variable to 1 for RCT Ts
  # and 0 for RCT Cs
  if (cace_itt == 0) {
    rct_dat_df$partic <- ifelse(rct_dat_df$tc == 1, 1, 0) 
  } else if ((cace_itt == 1) & (bad_particz == 0)) {
    rct_dat_df$partic <- rct_dat_df[,partic]
  } else {
    rct_dat_df$partic <- NA
  }
  
  rct_dat_df$ipw_wt <- 1   
  
  rcty  <- rct_dat_df[,c("yv","sitev","tc","partic","ipw_wt","pred_bias",
                         "left_right","clusv")]
  
  # Process xvars which is a vector with + signs
  xva <- deparse(substitute(xvars_rct_adj))
  xvb <- as.character(gsub('[+]', '', xva))
  # Split the xvb string by spaces
  xvc <- strsplit(xvb, split = "\\s+")
  xv  <- xvc[[1]] 
  
  if (xv[1] == 0) {
    adjx <- 0
    nadjx <- 0
  } else {
    adjx <- 1
    nadjx <- length(xv)
  }
  
  # If adjx is 1 then check if input x names are on the data frame 
  if (adjx == 1) {
    for (ix in 1:nadjx) {
      if (xv[ix] %in% names(rct_dat_df)) {
        junk <- 1
      } else {
        xxv <- as.character(xv[ix])
        rct_dat_df[,xxv] <- NA
      }
    }
    rctx  <- data.frame(rct_dat_df[,xv])
    rct_dat <- data.frame(rcty,rctx)
    
  } else if (adjx == 0) {
    rct_dat <- data.frame(rcty)
  }
  
  rct_dat$rct_qed_ind <- 1
  
  # Now for qed
  
  if (yv %in% names(qed_dat_df)) {
    qed_dat_df$yv <- qed_dat_df[,yv]
  } else {
    qed_dat_df$yv <- NA
  }
  
  if (sitev %in% names(qed_dat_df)) {
    qed_dat_df$sitev <- qed_dat_df[,sitev]
  } else {
    qed_dat_df$sitev <- NA
  }
  
  if (tc %in% names(qed_dat_df)) {
    qed_dat_df$tc <- qed_dat_df[,tc]
  } else {
    qed_dat_df$tc <- NA
  }
  
  if (rct_qed_ind %in% names(qed_dat_df)) {
    qed_dat_df$rct_qed_ind <- qed_dat_df[,rct_qed_ind]
    bad_rqi <- 0
  } else {
    qed_dat_df$rct_qed_ind <- NA
    bad_rqi <- 1
  }
  
  if (ipw_wt %in% names(qed_dat_df)) {
    qed_dat_df$ipw_wt <- qed_dat_df[,ipw_wt]
  } else {
    qed_dat_df$ipw_wt <- NA
  }
  
  if ("pred_bias" %in% names(qed_dat_df)) {
    qed_dat_df$pred_bias   <- qed_dat_df[,c("pred_bias")]
  } else {
    qed_dat_df$pred_bias <- NA
  }
  
  if ("left_right" %in% names(qed_dat_df)) {
    qed_dat_df$left_right  <- qed_dat_df[,c("left_right")]
  } else {
    qed_dat_df$left_right <- NA
  }
  
  if (clus_dat == 0) {
    qed_dat_df$clusv <- NA
  } else {
    if (clusv %in% names(qed_dat_df)) {
      qed_dat_df$clusv <- qed_dat_df[,clusv]
    } else {
      qed_dat_df$clusv <- NA
    }
  }
  
  # Check for a few errors in QED dataset for tc and rct_qed_ind 
  goodq_tc <- 1
  if (any(is.na(qed_dat_df$tc))) {
    goodq_tc <- 0
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing t_c values in qed_dat_df")
  } else {
    badq_tc <- ifelse(((qed_dat_df$tc == 0) | (qed_dat_df$tc == 1)), 0, 1)
    if (any(badq_tc == 1)) {
      goodq_tc <- 0
      err_count <<- err_count+1
      err_df[err_count,2] <- c("t_c values in qed_dat_df must all be 0 or 1")
    }
  }
  
  if (any(is.na(qed_dat_df$rct_qed_ind))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing rct_qed_site values in qed_dat_df")
  } else {
    badq_rct_qed <- ifelse(((qed_dat_df$rct_qed_ind == 0) | (qed_dat_df$rct_qed_ind == 1)), 0, 1)
    if (any(badq_rct_qed == 1)) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("rct_qed_site values in qed_dat_df must all be 0 or 1")
    }
  }
  
  # For QED, set participation variable to 1 for Ts and 0 for Cs
  qed_dat_df$partic <- ifelse(qed_dat_df$tc == 1,1,0)
  
  qedy  <- qed_dat_df[,c("yv","sitev","tc","partic","ipw_wt","rct_qed_ind",
                         "pred_bias","left_right","clusv")]
  
  # Process xvars which is a vector with + signs
  xvaq <- deparse(substitute(xvars_qed_adj))
  xvbq <- as.character(gsub('[+]', '', xvaq))
  # Split the xvb string by spaces
  xvcq <- strsplit(xvbq, split = "\\s+")
  xvq  <- xvcq[[1]] 
  
  if (xvq[1] == 0) {
    adjxq <- 0
    nadjxq <- 0
  } else {
    adjxq <- 1
    nadjxq <- length(xvq)
  }
  
  # If adjxq is 1 then check if input x names are on the data frame 
  if (adjxq == 1) {
    for (ix in 1:nadjxq) {
      if (xvq[ix] %in% names(qed_dat_df)) {
        junk <- 1
      } else {
        xxvq <- as.character(xvq[ix])
        qed_dat_df[,xxvq] <- NA
      }
    }
    qedx  <- data.frame(qed_dat_df[,xvq])
    qed_dat <- data.frame(qedy,qedx)
    
  } else if (adjxq == 0) {
    qed_dat <- data.frame(qedy)
  }
  
  # Create qed sites that are in the rct_qed sample and those in the qed-only sample  
  #qed_dat_temp <- data.frame(qedy,qedx)
  #qed_dat  <- data.frame(qedy)
  qed_datr <- qed_dat[qed_dat$rct_qed_ind == 1,]
  qed_datq <- qed_dat[qed_dat$rct_qed_ind == 0,]
  
  ###
  # Check for errors
  ###
  
  # Sitev 
  if (any(is.na(rct_dat$sitev))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing site_id values in rct_dat_df")
  }
  
  if (any(is.na(qed_dat$sitev))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing site_id values in qed_dat_df")
  }
  
  # yv and xv
  if (any(is.na(rct_dat$yv))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing yvar values in rct_dat_df")
  }
  
  if (adjx == 1) {
    if (any(is.na(rct_dat[,xv]))) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing xvars_rct_adj values in rct_dat_df")
    }
  }
  
  if (any(is.na(qed_dat$yv))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing yvar values in qed_dat_df")
  }
  
  if (adjxq == 1) {
    if (any(is.na(qed_dat[,xvq]))) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing xvars_qed_adj values in qed_dat_df")
    }
  }
  
  # tc
  if (any(is.na(rct_dat$tc))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing t_c values in rct_dat_df")
  } else {
    badr_tc <- ifelse(((rct_dat$tc == 0) | (rct_dat$tc == 1)), 0, 1)
    if (any(badr_tc == 1)) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("t_c values in rct_dat_df must all be 0 or 1")
    }
  }
  
  # pred-bias and left_right
  if (any(is.na(qed_dat$pred_bias))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing pred_bias values in qed_dat_df from the calc_pred_bias function")
  }
  
  if (any(is.na(rct_dat$pred_bias))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing pred_bias values in rct_dat_df from the calc_pred_bias function")
  }
  
  if (any(is.na(qed_dat$left_right))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing left_right values in qed_dat_df from the calc_pred_bias function")
  }
  
  if (any(is.na(rct_dat$left_right))) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing left_right values in rct_dat_df from the calc_pred_bias function")
  }
  
  # clusv
  if (clus_dat == 1) {
    if (any(is.na(rct_dat$clusv))) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing cluster_id values in rct_dat_df")
    }
    
    if (any(is.na(qed_dat$clusv))) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing cluster_id values in qed_dat_df")
    }
  }
  
  # cace_itt and partic
  if ((cace_itt == 0) | (cace_itt == 1)) {
    bad_cace_itt <- 0
  } else {
    bad_cace_itt <- 1
    err_count <<- err_count+1
    err_df[err_count,2] <- c("cace_itt value must be 0 or 1")
  }
  
  if (cace_itt == 1) {
    if (any(is.na(rct_dat$partic))) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing got_treat values in rct_dat_df")
    } else {
      badr_partic <- ifelse(((rct_dat$partic == 0) | (rct_dat$partic == 1)), 0, 1)
      if (any(badr_partic == 1)) {
        err_count <<- err_count+1
        err_df[err_count,2] <- c("got_treat values in rct_dat_df must all be 0 or 1")
      }
    }
  }
  
  # ipw_wt
  if (goodq_tc <- 1) {
    if (any(is.na(qed_dat$ipw_wt))) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing ipw_wgt values in qed_dat_df")
    } else {
      badq_wt <- ifelse(((qed_dat$ipw_wt < 0)), 1, 0)
      if (any(badq_wt == 1)) {
        err_count <<- err_count+1
        err_df[err_count,2] <- c("ipw_wgt values in qed_dat_df must all be nonnegative")
      }
    }
  }
  
  # inv_var_agg_wgt
  if ((inv_var_agg_wgt == 0) | (inv_var_agg_wgt == 1)) {
    bad_inv_wgt <- 0
  } else {
    bad_inv_wgt <- 1
    err_count <<- err_count+1
    err_df[err_count,2] <- c("inv_var_agg_wgt value must be 0 or 1")
  }
  
  # Return function if there are any errors
  
  if (err_count >= 1) {
    err_dfg <- err_df[1:err_count,2]
    cat("\n")
    print("ERRORS TO FIX IN HYBRID_IMPACTS FUNCTION")
    print(err_dfg)
    return(err_dfg)
  }
  
  # Define datasets for RCT control nonparticipants, RCT treatment nonparticipants,
  # and QED comparisons
  qed_datc  <- qed_dat[qed_dat$tc == 0,]
  qed_datt  <- qed_dat[qed_dat$tc == 1,]
  qed_datrc <- qed_datr[qed_datr$tc == 0,]
  qed_datrt <- qed_datr[qed_datr$tc == 1,]
  qed_datqc <- qed_datq[qed_datq$tc == 0,]
  qed_datqt <- qed_datq[qed_datq$tc == 1,]
  
  rct_datc <- rct_dat[rct_dat$tc == 0,]
  rct_datt <- rct_dat[rct_dat$tc == 1,]
  rct_datc_np <- rct_datc[rct_datc$partic == 0,]
  rct_datt_p  <- rct_datt[rct_datt$partic == 1,]
  
  # Check for missing site IDs and check that the same ones are on
  # the rct and qed datasets. Do separately for cace and itt analyses
  
  #sink(file = out_txt2, split = TRUE)
  
  if (cace_itt == 0) {
    
    site_qedrc <- data.frame(table(qed_datrc$sitev,useNA= "always"))
    site_qedrc <- site_qedrc %>%
      rename(site_var = Var1,
             qed_freqrc = Freq)
    site_qedrc <- site_qedrc[,c("site_var","qed_freqrc")]
    
    site_qedrt <- data.frame(table(qed_datrt$sitev,useNA= "always"))
    site_qedrt <- site_qedrt %>%
      rename(site_var = Var1,
             qed_freqrt = Freq)
    site_qedrt <- site_qedrt[,c("site_var","qed_freqrt")]
    
    site_rctc <- data.frame(table(rct_datc$sitev,useNA= "always"))
    site_rctc <- site_rctc %>%
      rename(site_var = Var1,
             rct_freqc = Freq)
    site_rctc <- site_rctc[,c("site_var","rct_freqc")]
    
    site_rctt <- data.frame(table(rct_datt$sitev,useNA= "always"))
    site_rctt <- site_rctt %>%
      rename(site_var = Var1,
             rct_freqt = Freq)
    site_rctt <- site_rctt[,c("site_var","rct_freqt")]
    
    sitev_m1 <- merge(site_rctt, site_rctc, by = "site_var", all = TRUE)
    sitev_m2 <- merge(sitev_m1, site_qedrt, by = "site_var", all = TRUE)
    sitev_m  <- merge(sitev_m2, site_qedrc, by = "site_var", all = TRUE)
    
    # Replace NA with 0 in site frequencies
    sitev_m$rct_freqt  <- ifelse(is.na(sitev_m$rct_freqt),0,sitev_m$rct_freqt)
    sitev_m$rct_freqc  <- ifelse(is.na(sitev_m$rct_freqc),0,sitev_m$rct_freqc)
    sitev_m$qed_freqrt <- ifelse(is.na(sitev_m$qed_freqrt),0,sitev_m$qed_freqrt)
    sitev_m$qed_freqrc <- ifelse(is.na(sitev_m$qed_freqrc),0,sitev_m$qed_freqrc)
    
    sitev_mr <- sitev_m[!is.na(sitev_m$site_var),]
    
    sink(file = out_txt2, split = TRUE, append = TRUE)
    
    cat("\n")
    pr1 <- "TABLE 1. HYBRID_IMPACT FUNCTION: RCT AND QED SITE SAMPLE SIZES IN THE"
    pr2 <- "RCT-QED SITES FOR THE FULL SAMPLE ITT ANALYSIS"
    cat(pr1,"\n",pr2,"\n") 
    print(sitev_m)
    
    sink()
    
    # Check for missing site IDs
    miss_site_var <- sitev_m[is.na(sitev_m$site_var),]
    
    if ((miss_site_var$rct_freqt > 0) | (miss_site_var$rct_freqc > 0) | 
        (miss_site_var$qed_freqrt > 0) | (miss_site_var$qed_freqrc > 0)) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing Site IDs in the RCT-QED sites")
    }
    
    # Check for mismatched RCT and QED site IDs
    nonmiss_site_var <- sitev_m[!is.na(sitev_m$site_var),]
    
    rct_zero_freqt  <- any(nonmiss_site_var$rct_freqt == 0)
    rct_zero_freqc  <- any(nonmiss_site_var$rct_freqc == 0)
    qed_zero_freqrt <- any(nonmiss_site_var$qed_freqrt == 0)
    qed_zero_freqrc <- any(nonmiss_site_var$qed_freqrc == 0)
    
    msite_ind <- 0
    if ((rct_zero_freqt == TRUE) | (rct_zero_freqc == TRUE) | 
        (qed_zero_freqrt == TRUE) | (qed_zero_freqrc == TRUE)) {
      msite_ind <- 1
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Some RCT and QED Site IDs in the RCT-QED sites do not match")
    }
    
  } else if (cace_itt == 1) {
    
    site_qedrc <- data.frame(table(qed_datrc$sitev,useNA= "always"))
    site_qedrc <- site_qedrc %>%
      rename(site_var = Var1,
             qed_freqrc = Freq)
    site_qedrc <- site_qedrc[,c("site_var","qed_freqrc")]
    
    site_qedrt <- data.frame(table(qed_datrt$sitev,useNA= "always"))
    site_qedrt <- site_qedrt %>%
      rename(site_var = Var1,
             qed_freqrt = Freq)
    site_qedrt <- site_qedrt[,c("site_var","qed_freqrt")]
    
    site_rctc <- data.frame(table(rct_datc$sitev,useNA= "always"))
    site_rctc <- site_rctc %>%
      rename(site_var = Var1,
             rct_freqc = Freq)
    site_rctc <- site_rctc[,c("site_var","rct_freqc")]
    
    site_rctt <- data.frame(table(rct_datt$sitev,useNA= "always"))
    site_rctt <- site_rctt %>%
      rename(site_var = Var1,
             rct_freqt = Freq)
    site_rctt <- site_rctt[,c("site_var","rct_freqt")]
    
    sitev_m1 <- merge(site_rctt, site_rctc, by = "site_var", all = TRUE)
    sitev_m2 <- merge(sitev_m1, site_qedrt, by = "site_var", all = TRUE)
    sitev_m  <- merge(sitev_m2, site_qedrc, by = "site_var", all = TRUE)
    
    # Replace NA with 0 in site frequencies
    sitev_m$rct_freqt  <- ifelse(is.na(sitev_m$rct_freqt),0,sitev_m$rct_freqt)
    sitev_m$rct_freqc  <- ifelse(is.na(sitev_m$rct_freqc),0,sitev_m$rct_freqc)
    sitev_m$qed_freqrt <- ifelse(is.na(sitev_m$qed_freqrt),0,sitev_m$qed_freqrt)
    sitev_m$qed_freqrc <- ifelse(is.na(sitev_m$qed_freqrc),0,sitev_m$qed_freqrc)
    
    sitev_mr <- sitev_m[!is.na(sitev_m$site_var),]
    
    sink(file = out_txt2, split = TRUE, append = TRUE)
    
    cat("\n")
    pr1 <- "TABLE 1. HYBRID_IMPACT FUNCTION: RCT AND QED SITE SAMPLE SIZES IN THE"
    pr2 <- " RCT-QED SITES FOR THE FULL SAMPLE CACE ANALYSIS"
    cat(pr1,"\n",pr2,"\n") 
    print(sitev_m)
    
    sink()
    
    # Check for missing site IDs
    miss_site_var <- sitev_m[is.na(sitev_m$site_var),]
    
    if ((miss_site_var$rct_freqt > 0) | (miss_site_var$rct_freqc > 0) | 
        (miss_site_var$qed_freqrt > 0) | (miss_site_var$qed_freqrc > 0)) {
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Missing Site IDs in the RCT-QED sites")
    }
    
    # Check for mismatched RCT and QED site IDs
    nonmiss_site_var <- sitev_m[!is.na(sitev_m$site_var),]
    
    rct_zero_freqt  <- any(nonmiss_site_var$rct_freqt == 0)
    rct_zero_freqc  <- any(nonmiss_site_var$rct_freqc == 0)
    qed_zero_freqrt <- any(nonmiss_site_var$qed_freqrt == 0)
    qed_zero_freqrc <- any(nonmiss_site_var$qed_freqrc == 0)
    
    msite_ind <- 0
    if ((rct_zero_freqt == TRUE) | (rct_zero_freqc == TRUE) | 
        (qed_zero_freqrt == TRUE) | (qed_zero_freqrc == TRUE)) {
      msite_ind <- 1
      err_count <<- err_count+1
      err_df[err_count,2] <- c("Some RCT and QED Site IDs in the RCT-QED sites do not match")
    }
    
    # Check if there are enough degrees of freedom for RCT sample members 
    # across all sites to run regressions
    nonmiss_site_var$rct_freqtc <- nonmiss_site_var$rct_freqt + 
      nonmiss_site_var$rct_freqc
    
    nonmiss_var_rct <- sum(nonmiss_site_var$rct_freqtc)
    
    rct_freq <- 0
    if (nonmiss_var_rct <= (nadjx+2)) {rct_freq <- 1}
    
    if (msite_ind == 0) {
      if (rct_freq == 1) {
        err_count <<- err_count+1
        err_df[err_count,2] <- paste(c("The RCT-QED sites have fewer RCT members than available degrees of freedom"))
      }
    }
    
    # Check if the participation rate is greater for RCT treatments than controls
    # in all sites
    
    if (err_count == 0) {
      part_ratet <- data.frame(aggregate(partic ~ sitev, data = rct_datt, FUN = mean))
      part_ratec <- data.frame(aggregate(partic ~ sitev, data = rct_datc, FUN = mean))
      
      partic_df <- merge(part_ratet, part_ratec, by = "sitev")
      colnames(partic_df) <- c("site_var", "part_ratet", "part_ratec")
      
      partic_df$compl_rate <- partic_df$part_ratet - partic_df$part_ratec
      
      pc_df1 <- partic_df
      pc_df1$part_ratet <- format(round(pc_df1$part_ratet,4), nsmall=4)
      pc_df1$part_ratec <- format(round(pc_df1$part_ratec,4), nsmall=4)
      pc_df1$compl_rate <- format(round(pc_df1$compl_rate,4), nsmall=4)
      
      sink(file = out_txt2, split = TRUE, append = TRUE)
      
      cat("\n")
      pr1 <- "TABLE 1A. HYBRID_IMPACT FUNCTION: RCT TREATMENT AND CONTROL GROUP PARTICIPATION"
      pr2 <- "RATES AND THE COMPLIANCE RATE FOR THE FULL SAMPLE CACE ANALYSIS"
      cat(pr1,"\n",pr2,"\n") 
      print(pc_df1)
      
      sink()
      
      any_bad_compl <- any(partic_df$compl_rate <= 0)
      
      if (any_bad_compl == TRUE) {
        err_count <<- err_count+1
        err_df[err_count,2] <- c("In some sites, the participation rate is smaller for RCT treatments than RCT controls")
      }
    } # if err_count==0
    
  } #end cace_itt
  
  # Do this for the qed sample in the QED-only sites
  
  site_qedqc <- data.frame(table(qed_datqc$sitev,useNA= "always"))
  site_qedqc <- site_qedqc %>%
    rename(site_var = Var1,
           qed_freqqc = Freq)
  site_qedqc <- site_qedqc[,c("site_var","qed_freqqc")]
  
  site_qedqt <- data.frame(table(qed_datqt$sitev,useNA= "always"))
  site_qedqt <- site_qedqt %>%
    rename(site_var = Var1,
           qed_freqqt = Freq)
  site_qedqt <- site_qedqt[,c("site_var","qed_freqqt")]
  
  sitev_m  <- merge(site_qedqt, site_qedqc, by = "site_var", all = TRUE)
  
  # Replace NA with 0 in site frequencies
  sitev_m$qed_freqqt <- ifelse(is.na(sitev_m$qed_freqqt),0,sitev_m$qed_freqqt)
  sitev_m$qed_freqqc <- ifelse(is.na(sitev_m$qed_freqqc),0,sitev_m$qed_freqqc)
  
  sitev_mq <- sitev_m[!is.na(sitev_m$site_var),]
  
  sink(file = out_txt2, split = TRUE, append = TRUE)
  
  cat("\n")
  pr1 <- "TABLE 2. HYBRID_IMPACT FUNCTION: QED SITE SAMPLE SIZES IN THE QED-ONLY SITES"
  pr2 <- "FOR THE FULL SAMPLE ANALYSIS"
  cat(pr1,"\n",pr2,"\n") 
  print(sitev_m)
  
  sink()
  
  # Check for missing site IDs
  miss_site_var <- sitev_m[is.na(sitev_m$site_var),]
  
  if ((miss_site_var$qed_freqqt > 0) | (miss_site_var$qed_freqqc > 0)) {
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Missing Site IDs in the QED-only sites")
  }
  
  # Check for mismatched QED and comparison site IDs
  nonmiss_site_var <- sitev_m[!is.na(sitev_m$site_var),]
  
  qed_zero_freqqt <- any(nonmiss_site_var$qed_freqqt == 0)
  qed_zero_freqqc <- any(nonmiss_site_var$qed_freqqc == 0)
  
  msite_indq <- 0
  if ((qed_zero_freqqt == TRUE) | (qed_zero_freqqc == TRUE)) {
    msite_indq <- 1
    err_count <<- err_count+1
    err_df[err_count,2] <- c("Some QED Site IDs for treatments and comparisons in the QED-only sites do not match")
  }
  
  # Check for degrees of freedom
  nonmiss_site_var$qed_freqqtc <- nonmiss_site_var$qed_freqqt + 
    nonmiss_site_var$qed_freqqc
 
  nonmiss_var_qed <- sum(nonmiss_site_var$qed_freqqtc)
 
  qed_freq <- 0
  if (nonmiss_var_qed <= (nadjxq+2)) {qed_freq <- 1}
  
  if (msite_indq == 0) {
    if (qed_freq == 1) {
      err_count <<- err_count+1
      err_df[err_count,2] <- paste(c("The RCT-QED sites have fewer RCT members than available degrees of freedom"))
    }
  }
  
  # Create site ID that goes from 1 to # of sites. The cur_group_id() variable is an
  # internal dplr function that retains the current group counter
  
  rct_dat <- rct_dat %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
    rct_dat <- data.frame(rct_dat)
  
  rct_datc <- rct_datc %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
    rct_datc <- data.frame(rct_datc)
  
  rct_datt <- rct_datt %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
    rct_datt <- data.frame(rct_datt)
  
  #qed_dat <- qed_dat %>%
  #  group_by(sitev) %>%
  #  mutate(siteID = cur_group_id()) %>%
  #  ungroup()
  #qed_dat <- data.frame(qed_dat)
    
  qed_datr <- qed_datr %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
  qed_datr <- data.frame(qed_datr)
    
  qed_datrc <- qed_datrc %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
    qed_datrc <- data.frame(qed_datrc)
  
  qed_datrt <- qed_datrt %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
    qed_datrt <- data.frame(qed_datrt)
  
  qed_datq <- qed_datq %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
    qed_datq <- data.frame(qed_datq)
  
  qed_datqc <- qed_datqc %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
    qed_datqc <- data.frame(qed_datqc)
  
  qed_datqt <- qed_datqt %>%
    group_by(sitev) %>%
    mutate(siteID = cur_group_id()) %>%
    ungroup()
    qed_datqt <- data.frame(qed_datqt)
    
  # Create cluster variable that differs across sites
  # since some might code clusters 1,.. within each site
    
  if (clus_dat == 1) {
    
    rdatclus <- rct_dat[,c("sitev","clusv")]
    qdatclus <- qed_dat[,c("sitev","clusv")]
    
    rqdatclus <- rbind(rdatclus,qdatclus)
    
    rqclus <- data.frame(ftable(rqdatclus$sitev,rqdatclus$clusv))
    rqclus <- rqclus[rqclus$Freq > 0,]
    
    colnames(rqclus)[which(names(rqclus) == "Var1")] <- "sitev"
    colnames(rqclus)[which(names(rqclus) == "Var2")] <- "clusv"
    rqclus$clusvID <- 1:nrow(rqclus)
    rqclus <- rqclus[,c("sitev","clusv","clusvID")]
    
    rct_dat <- merge(rct_dat, rqclus, by = c("sitev","clusv"), all.x = TRUE)
    rct_dat$clusv_orig <- rct_dat$clusv
    rct_dat$clusv      <- rct_dat$clusvID
    #print("RCT_DAT")
    #print(table(rct_dat$clusv,useNA = "always"))
    
    rct_datc <- merge(rct_datc, rqclus, by = c("sitev","clusv"), all.x = TRUE)
    rct_datc$clusv_orig <- rct_datc$clusv
    rct_datc$clusv      <- rct_datc$clusvID
    #print(rct_datc[,c("clusv","clusvID","clusv_orig")])
    #print("RCT_DATC")
    #print(table(rct_datc$clusv,useNA = "always"))
    
    rct_datt <- merge(rct_datt, rqclus, by = c("sitev","clusv"), all.x = TRUE)
    rct_datt$clusv_orig <- rct_datt$clusv
    rct_datt$clusv      <- rct_datt$clusvID
    #print("RCT_DATT")
    #print(table(rct_datt$clusv,useNA = "always"))
    
    qed_dat <- merge(qed_dat, rqclus, by = c("sitev","clusv"), all.x = TRUE)
    qed_dat$clusv_orig <- qed_dat$clusv
    qed_dat$clusv      <- qed_dat$clusvID
    #print("QED_DAT")
    #print(table(qed_dat$clusv,useNA = "always"))
    
    qed_datr <- merge(qed_datr, rqclus, by = c("sitev","clusv"), all.x = TRUE)
    qed_datr$clusv_orig <- qed_datr$clusv
    qed_datr$clusv      <- qed_datr$clusvID
    #print("QED_DATR")
    #print(table(qed_datr$clusv,useNA = "always"))
    
    qed_datrc <- merge(qed_datrc, rqclus, by = c("sitev","clusv"), all.x = TRUE)
    qed_datrc$clusv_orig <- qed_datrc$clusv
    qed_datrc$clusv      <- qed_datrc$clusvID
    #print("QED_DATRC")
    #print(table(qed_datrc$clusv,useNA = "always"))
    
    qed_datrt <- merge(qed_datrt, rqclus, by = c("sitev","clusv"), all.x = TRUE)
    qed_datrt$clusv_orig <- qed_datrt$clusv
    qed_datrt$clusv      <- qed_datrt$clusvID
    #print("QED_DATRT")
    #print(table(qed_datrt$clusv,useNA = "always"))
    
    qed_datq <- merge(qed_datq, rqclus, by = c("sitev","clusv"), all.x = TRUE)
    qed_datq$clusv_orig <- qed_datq$clusv
    qed_datq$clusv      <- qed_datq$clusvID
    #print("QED_DATQ")
    #print(table(qed_datq$clusv,useNA = "always"))
    
    qed_datqc <- merge(qed_datqc, rqclus, by = c("sitev","clusv"), all.x = TRUE)
    qed_datqc$clusv_orig <- qed_datqc$clusv
    qed_datqc$clusv      <- qed_datqc$clusvID
    #print("QED_DATQC")
    #print(table(qed_datqc$clusv,useNA = "always"))
    
    qed_datqt <- merge(qed_datqt, rqclus, by = c("sitev","clusv"), all.x = TRUE)
    qed_datqt$clusv_orig <- qed_datqt$clusv
    qed_datqt$clusv      <- qed_datqt$clusvID
    #print("QED_DATQT")
    #print(table(qed_datqt$clusv,useNA = "always"))
    
  }
    
  rct_datc_np <- rct_datc[rct_datc$partic == 0,]
  rct_datt_p  <- rct_datt[rct_datt$partic == 1,]
  rct_datt_np <- rct_datt[rct_datt$partic == 0,]
  
  # Return function if there are any errors
  
  if (err_count >= 1) {
    err_dfg <- err_df[1:err_count,2]
    cat("\n")
    print("ERRORS TO FIX IN HYBRID_IMPACTS FUNCTION")
    print(err_dfg)
    return(err_dfg)
  }
  
  se_site <- max(rct_datt$siteID)
  sq_site <- max(qed_datqt$siteID)
  
  # Run RCT IV Analysis for CACE analysis or regression analysis for ITT analysis
  
  rct_imp <- rep(0,se_site)
  rct_se  <- rep(0,se_site)
  rct_ne  <- rep(0,se_site)
  rct_net <- rep(0,se_site)
  rct_nec <- rep(0,se_site)
  rct_p10 <- rep(0,se_site)
  
  meancr  <- rep(0,se_site)
  rct_s2t <- rep(0,se_site)
  
  # Run interacted model
  
  # Create site indicators
  rct_dat$siteIDf <- as.factor(rct_dat$siteID)
  
  # Create regression formulas
  
  if (cace_itt == 1) {
    if (se_site > 1) {
      if (adjx == 1) {
        formula_iv <- paste("yv", " ~ ", "siteIDf", "+", "partic:siteIDf", "+", 
                            xva, "+", "0", "|", "siteIDf", "+", "tc:siteIDf", "+",
                            xva, "+", "0")
      } else if (adjx == 0) {
        formula_iv <- paste("yv", " ~ ", "siteIDf", "+", "partic:siteIDf", "+", 
                            "0", "|", "siteIDf", "+", "tc:siteIDf", "+", "0")
      } # end adjx
    } else if (se_site == 1) {
      if (adjx == 1) {
        formula_iv <- paste("yv", " ~ ", "partic", "+", xva, "|", "tc", "+", xva)
      } else if (adjx == 0) {
        formula_iv <- paste("yv", " ~ ", "partic", "|", "tc")
      }
    }
  } else if (cace_itt == 0) {
    if (se_site > 1) {
      if (adjx == 1) {
        formula_lm <- paste("yv", " ~ ", "siteIDf", "+", "tc:siteIDf", "+", 
                            xva, "+", "0")
      } else if (adjx == 0) {
        formula_lm <- paste("yv", " ~ ", "siteIDf", "+", "tc:siteIDf", "+", "0")
      } # end adjx
    } else if (se_site == 1) {
      if (adjx == 1) {
        formula_lm <- paste("yv", " ~ ", "tc", "+", xva)
      } else if (adjx == 0) {
        formula_lm <- paste("yv", " ~ ", "tc")
      }
    }
  }
  
  if (cace_itt == 1) {
    rct_agg <- ivreg(formula_iv, data = rct_dat)
  } else if (cace_itt == 0) {
    rct_agg <- lm(formula_lm, data = rct_dat)
  }
  
  # Calculate robust standard errs and do it again if there is clustering
  # to calculate design effects
  # NEED TO TEST
  
  robust_vcov <- vcovHC(rct_agg, type = "HC1")
  #robust_res  <- coeftest(rct_agg, vcov = vcovHC(rct_agg, type = 'HC1'))
  
  # Extract the square root of the diagonal elements to get the standard errs
  rct_se  <- sqrt(diag(robust_vcov))
  rct_imp <- rct_agg$coefficients
  
  if (clus_dat == 1) {
    robust_vcov_clus <- vcovCR(rct_agg, cluster = rct_dat$clusv, type = "CR1")
    rct_se_clus  <- sqrt(diag(robust_vcov_clus))
  }
  
  # Pull off the impacts and ses for interaction terms
  nimpr    <- NROW(rct_imp)
  
  if (se_site > 1) {
    rct_imp  <- rct_imp[(se_site+nadjx+1):nimpr]
    rct_se   <- rct_se[(se_site+nadjx+1):nimpr]
    if (clus_dat == 1) {
      rct_se_clus <- rct_se_clus[(se_site+nadjx+1):nimpr]
    }
  } else if (se_site == 1) {
    rct_imp  <- rct_imp[2]
    rct_se   <- rct_se[2]
    if (clus_dat == 1) {
      rct_se_clus <- rct_se_clus[2]
    }
  }
  
  # Calculate control group mean and treatment group standard deviation
  # to facilitate the computations as the biases complicate the QED
  # control group std 
  for (i in 1:se_site) {
    
    rct_dats   <- rct_dat[rct_dat$siteID==i,]
    rct_ne[i]  <- nrow(rct_dats)
    rct_net[i] <- nrow(rct_dats[rct_dats$tc == 1,])
    rct_nec[i] <- nrow(rct_dats[rct_dats$tc == 0,])
    
    # Calculate treatment group s2 - assumed the same for CACE and ITT
    rct_s2t[i] <- var(rct_dats$yv[rct_dats$tc == 1])
    
    if (cace_itt == 1) {
      
      # Use unadjusted complier rate to calculate control complier mean
      rct_d1 <- mean(rct_dats$partic[rct_dats$tc == 1])
      rct_d0 <- mean(rct_dats$partic[rct_dats$tc == 0])
      
      rct_p10[i] <- rct_d1 - rct_d0
      
      # Calculate complier mean in site
      mcnp <- mean(rct_dats$yv[((rct_dats$tc == 0) & (rct_dats$partic == 0))])
      
      if (rct_d1 < 1) {
        mtnp <- mean(rct_dats$yv[(rct_dats$tc == 1) & (rct_dats$partic == 0)])
      } else if (rct_d1 == 1) {
        mtnp <- 0
      }
      
      meancr[i] <- (((1-rct_d0)*mcnp) - ((1-rct_d1)*mtnp)) / (rct_d1 - rct_d0)
      
    } else if (cace_itt == 0) {
      
      rct_p10[i] <- 1
      meancr[i]  <- mean(rct_dats$yv[rct_dats$tc == 0])
    }
    
  } # End for loop for IV or White regressions by site
  
  # Calculate the RCT weighted mean and variance across sites
  # Where the weights are the estimated number of compliers for IV or individuals 
  # for ITT in the site
  
  rct_s2t_pool <- sum(rct_net*rct_s2t)/sum(rct_net)
  rct_sdt_pool <- rct_s2t_pool^.5
  
  rct_wgt_site <- rct_ne*rct_p10
  rct_wgt_sum  <- sum(rct_wgt_site)
  
  rct_impw     <- rct_imp*rct_wgt_site
  rct_impw_sum <- sum(rct_impw)
  rct_imp_pool <- rct_impw_sum/rct_wgt_sum
  
  rct_meancw     <- meancr*rct_wgt_site
  rct_meancw_sum <- sum(rct_meancw)
  rct_meanc_pool <- rct_meancw_sum/rct_wgt_sum
  
  rct_var  <- rct_se^2
  rct_w2   <- rct_wgt_site^2
  rct_varw <- rct_var*rct_w2
  rct_varw_sum <- sum(rct_varw)
  
  rct_var_pool <- rct_varw_sum/(rct_wgt_sum^2)
  rct_se_pool  <- rct_var_pool^.5
  
  if (clus_dat == 1) {
    
    rct_var_clus  <- rct_se_clus^2
    rct_varw_clus <- rct_var_clus*rct_w2
    rct_varw_sum_clus <- sum(rct_varw_clus)
    
    rct_var_pool_clus <- rct_varw_sum_clus/(rct_wgt_sum^2)
    rct_se_pool_clus  <- rct_var_pool_clus^.5
    
    rct_deff_clus <- rct_var_pool_clus / rct_var_pool
    
    rct_var_pool <- rct_var_pool_clus
    rct_se_pool  <- rct_se_pool_clus
    rct_se       <- rct_se_clus
    rct_var      <- rct_var_clus
  }
  
  # Significance testing
  
  if (clus_dat == 1) {
    ncluse <- NROW(table(rct_dat$clusv))
    dofr   <- ncluse - 2*se_site - nadjx
  } else if (clus_dat == 0) {
    dofr <- rct_wgt_sum - 2*se_site - nadjx
  }
  
  if (dofr <= 0) {
    dofr <- 1
  }
  
  t_statr  <- abs(rct_imp_pool)/rct_se_pool
  p_valuer <- 2*pt(t_statr, df = dofr,  lower.tail=FALSE)
  
  # Print RCT impact results
  rct_wgt_sitez    <- format(rct_wgt_site,digits=2, nsmall=2)
  rct_wgt_sumz     <- format(rct_wgt_sum,digits=2, nsmall=2)
  rct_impz         <- format(rct_imp,digits=4, nsmall=4)
  rct_sez          <- format(rct_se,digits=4, nsmall=4)
  rct_meanc_poolz  <- format(rct_meanc_pool,digits=4, nsmall=4)
  rct_imp_poolz    <- format(rct_imp_pool,digits=4, nsmall=4)
  rct_se_poolz     <- format(rct_se_pool,digits=4, nsmall=4)
  rct_sdt_poolz    <- format(rct_sdt_pool,digits=4, nsmall=4)
  meancrz          <- format(meancr,digits=4, nsmall=4)
  t_statrz         <- format(t_statr,digits=2, nsmall=2)
  p_valuerz        <- format(p_valuer,digits=3, nsmall=3)
  
  ###
  # CONDUCT QED ANALYSIS FOR QED-ONLY SITES: USES SURVEY SVY_DESIGN and SVYGLM
  ###
  
  # Create site indicators
  qed_datq$siteIDf <- as.factor(qed_datq$siteID)
  
  # Create regression formulas
  if (sq_site > 1) {
    if (adjxq == 1) {
      formula_qed <- paste("yv", " ~ ", "siteIDf", "+", "tc:siteIDf", "+",
                           xvaq, "+", "0")
    } else if (adjxq == 0) {
      formula_qed <- paste("yv", " ~ ", "siteIDf", "+", "tc:siteIDf", "+", "0")
    }
  } else if (sq_site == 1) {
    if (adjx == 1) {
      formula_qed <- paste("yv", " ~ ", "tc", "+", xvaq)
    } else if (adjx == 0) {
      formula_qed <- paste("yv", " ~ ", "tc")
    }
  }
  
  # Define design depending on whether clustering or not
  dsgn_spec    <- svydesign(id = ~1, weights = ~ipw_wt, data = qed_datq)
  if (clus_dat == 1) {
    dsgn_spec_clus <- svydesign(id = ~clusv, weights = ~ipw_wt, data = qed_datq)
  }
  
  surv_reg <- svyglm(formula_qed, design = dsgn_spec)
  s_out    <- summary(surv_reg)
  qed_coef <- s_out$coefficients
  
  if (clus_dat == 1) {
    surv_reg_clus <- svyglm(formula_qed, design = dsgn_spec_clus)
    s_out_clus    <- summary(surv_reg_clus)
    qed_coef_clus <- s_out_clus$coefficients
  }
  
  if (sq_site > 1) {
    nimpq   <- nrow(qed_coef)
    qed_imp <- qed_coef[(sq_site+nadjxq+1):nimpq,c("Estimate")]
    qed_se  <- qed_coef[(sq_site+nadjxq+1):nimpq,c("Std. Error")]
    if (clus_dat == 1) {
      qed_se_clus <- qed_coef_clus[(sq_site+nadjxq+1):nimpq,c("Std. Error")]
    }
  } else if (sq_site == 1) {
    qed_imp <- qed_coef[1,c("Estimate")]
    qed_se  <- qed_coef[1,c("Std. Error")]
    if (clus_dat == 1) {
      qed_se_clus  <- qed_coef_clus[1,c("Std. Error")]
    }
  }
  
  # Run over each QED-only site to get design effects and eff samp sizes
  # Some of this is old to get impacts and std. errs separately by site
  # Yields very similar results to svyglm
  meanqt      <- rep(0,sq_site)
  s2yqt       <- rep(0,sq_site)
  seqt        <- rep(0,sq_site)
  
  meanqc      <- rep(0,sq_site)
  s2yqc       <- rep(0,sq_site)
  seqc        <- rep(0,sq_site)
  deffqc      <- rep(0,sq_site)
  
  qed_eff_samp <- rep(0,sq_site)
  qed_nq      <- rep(0,sq_site)
  qed_nqt     <- rep(0,sq_site)   
  qed_nqc     <- rep(0,sq_site)
  
  qed_s2t     <- rep(0,sq_site)
  
  for (i in 1:sq_site) {
    
    yqti <- qed_datqt$yv[qed_datqt$siteID==i]
    yqci <- qed_datqc$yv[qed_datqc$siteID==i]
    ipwi <- qed_datqc$ipw_wt[qed_datqc$siteID==i]
    
    # Calculate means and standard errors for Ts
    meanqt[i]  <- mean(yqti)
    qed_nqt[i] <- length(yqti)
    s2yqt[i]   <- sum((yqti-meanqt[i])^2)/(qed_nqt[i]-1)
    seqt[i]    <- (s2yqt[i] / qed_nqt[i])^.5
    
    # Calculate treatment group s2 - assumed the same for CACE and ITT
    qed_s2t[i] <- var(yqti)
    
    # Calculate means and standard errors for Comparisons
    # Weighted means for Comparisons
    meanqc[i] <- sum(ipwi*yqci) / sum(ipwi)
    
    # Weighted variances for Comparisons
    qed_nqc[i] <- length(yqci)
    vv         <- sum(ipwi*((yqci - meanqc[i])^2)) / sum(ipwi)
    s2yqc[i]   <- qed_nqc[i]*vv/(qed_nqc[i]-1)
    
    # Standard error and design effects for Comparisons
    vv1 <- sum((ipwi^2)*((yqci - meanqc[i])^2)) / ((sum(ipwi))^2)
    vv1 <- qed_nqc[i]*vv1/(qed_nqc[i]-1)
    seqc[i]   <- vv1^.5
    
    wwq        <- sum(ipwi^2)/((sum(ipwi))^2)  
    deffqc[i]  <- wwq*qed_nqc[i]
    
    # Effective sample size needed for pooling
    qed_eff_samp[i] <- (qed_nqc[i]/deffqc[i]) + qed_nqt[i]
    qed_nq[i]   <- qed_nqt[i] + qed_nqc[i]
    
  } # End loop over QED-only sites 
  
  # Calculate pooled std, mean, and se
  qed_s2t_pool <- sum(qed_nqt*qed_s2t)/sum(qed_nqt)
  qed_sdt_pool <- qed_s2t_pool^.5
  
  qed_var <- qed_se^2  
  qed_imp_pool <- sum(qed_eff_samp*qed_imp)/ sum(qed_eff_samp)
  qed_var_pool <- sum((qed_eff_samp^2)*qed_var) / ((sum(qed_eff_samp))^2)
  qed_se_pool  <- qed_var_pool^.5
  
  qed_meanc_pool  <- sum(qed_eff_samp*meanqc)/ sum(qed_eff_samp)
  
  if (clus_dat == 1) {
    qed_var_clus <- qed_se_clus^2  
    qed_var_pool_clus <- sum((qed_eff_samp^2)*qed_var_clus) / ((sum(qed_eff_samp))^2)
    qed_se_pool_clus  <- qed_var_pool_clus^.5
    
    qed_deff_clus <- qed_var_pool_clus / qed_var_pool
    
    qed_se  <- qed_se_clus
    qed_var <- qed_var_clus
    qed_var_pool <- qed_var_pool_clus
    qed_se_pool  <- qed_se_pool_clus
  }
  
  if (clus_dat == 1) {
    nclusq <- NROW(table(qed_datq$clusv))
    dofq   <- nclusq - 2*sq_site - nadjxq
  } else if (clus_dat == 0) {
    dofq <- sum(qed_nq) - 2*sq_site - nadjxq
  }
  
  if (dofq <= 0) {
    dofq <- 1
  }
  
  t_statq  <- abs(qed_imp_pool)/qed_se_pool
  p_valueq <- 2*pt(t_statq, df = dofq,  lower.tail=FALSE)
  
  # Print QED impact results
  qed_impz  <- format(qed_imp,digits=4, nsmall=4)
  qed_sez   <- format(qed_se,digits=4, nsmall=4)
  eff_sampz <- format(qed_eff_samp,digits=2, nsmall=2)
  meanqtz   <- format(meanqt,digits=4, nsmall=4)
  meanqcz   <- format(meanqc,digits=4, nsmall=4)
  seqtz     <- format(seqt,digits=4, nsmall=4)
  seqcz     <- format(seqc,digits=4, nsmall=4)
  qed_imp_poolz   <- format(qed_imp_pool,digits=4, nsmall=4)
  qed_meanc_poolz <- format(qed_meanc_pool,digits=4, nsmall=4)
  qed_se_poolz    <- format(qed_se_pool,digits=4, nsmall=4)
  qed_sdt_poolz   <- format(qed_sdt_pool,digits=4, nsmall=4)
  t_statqz        <- format(t_statq,digits=2, nsmall=2)
  p_valueqz       <- format(p_valueq,digits=3, nsmall=3)
  
  ###
  # CONDUCT QED ANALYSIS FOR RCT-QED SITES: USES SURVEY SVY_DESIGN and SVYGLM
  ###
  
  # Create site indicators
  
  qed_datr$siteIDf <- as.factor(qed_datr$siteID)
  
  # Create regression formulas
  if (se_site > 1) {
    if (adjxq == 1) {
      formula_qedr <- paste("yv", " ~ ", "siteIDf", "+", "tc:siteIDf", "+",
                            xvaq, "+", "0")
    } else if (adjxq == 0) {
      formula_qedr <- paste("yv", " ~ ", "siteIDf", "+", "tc:siteIDf", "+", "0")
    }
  } else if (se_site == 1) {
    if (adjx == 1) {
      formula_qedr <- paste("yv", " ~ ", "tc", "+", xvaq)
    } else if (adjx == 0) {
      formula_qedr <- paste("yv", " ~ ", "tc")
    }
  }
  
  # Define design depending on whether clustering or not
  dsgn_specr    <- svydesign(id = ~1, weights = ~ipw_wt, data = qed_datr)
  if (clus_dat == 1) {
    dsgn_spec_clusr <- svydesign(id = ~clusv, weights = ~ipw_wt, data = qed_datr)
  }
  
  surv_regr <- svyglm(formula_qedr, design = dsgn_specr)
  s_outr    <- summary(surv_regr)
  qed_coefr <- s_outr$coefficients
  
  if (clus_dat == 1) {
    surv_reg_clusr <- svyglm(formula_qedr, design = dsgn_spec_clusr)
    s_out_clusr    <- summary(surv_reg_clusr)
    qed_coef_clusr <- s_out_clusr$coefficients
  }
  
  if (se_site > 1) {
    nimpqr   <- nrow(qed_coefr)
    qed_impr <- qed_coefr[(se_site+nadjxq+1):nimpqr,c("Estimate")]
    qed_ser  <- qed_coefr[(se_site+nadjxq+1):nimpqr,c("Std. Error")]
    if (clus_dat == 1) {
      qed_se_clusr <- qed_coef_clusr[(se_site+nadjxq+1):nimpqr,c("Std. Error")]
    }
  } else if (se_site == 1) {
    qed_impr <- qed_coefr[1,c("Estimate")]
    qed_ser  <- qed_coefr[1,c("Std. Error")]
    if (clus_dat == 1) {
      qed_se_clusr  <- qed_coef_clusr[1,c("Std. Error")]
    }
  }
  
  # Run over each RCT-QED site to get design effects and eff samp sizes
  # Some of this is old to get impacts and std. errs separately by site
  # Yields very similar results to svyglm
  meanqtr      <- rep(0,se_site)
  s2yqtr       <- rep(0,se_site)
  seqtr        <- rep(0,se_site)
  
  meanqcr      <- rep(0,se_site)
  s2yqcr       <- rep(0,se_site)
  seqcr        <- rep(0,se_site)
  deffqcr      <- rep(0,se_site)
  
  qed_eff_sampr <- rep(0,se_site)
  qed_nqr      <- rep(0,se_site)
  qed_nqtr     <- rep(0,se_site)   
  qed_nqcr     <- rep(0,se_site)
  
  qed_s2tr     <- rep(0,se_site)
  
  for (i in 1:se_site) {
    
    yqtir <- qed_datrt$yv[qed_datrt$siteID==i]
    yqcir <- qed_datrc$yv[qed_datrc$siteID==i]
    ipwir <- qed_datrc$ipw_wt[qed_datrc$siteID==i]
    
    # Calculate means and standard errors for Ts
    meanqtr[i]  <- mean(yqtir)
    qed_nqtr[i] <- length(yqtir)
    s2yqtr[i]   <- sum((yqtir-meanqtr[i])^2)/(qed_nqtr[i]-1)
    seqtr[i]    <- (s2yqtr[i] / qed_nqtr[i])^.5
    
    # Calculate treatment group s2 - assumed the same for CACE and ITT
    qed_s2tr[i] <- var(yqtir)
    
    # Calculate means and standard errors for Comparisons
    # Weighted means for Comparisons
    meanqcr[i] <- sum(ipwir*yqcir) / sum(ipwir)
    
    # Weighted variances for Comparisons
    qed_nqcr[i] <- length(yqcir)
    vvr         <- sum(ipwir*((yqcir - meanqcr[i])^2)) / sum(ipwir)
    s2yqcr[i]   <- qed_nqcr[i]*vvr/(qed_nqcr[i]-1)
    
    # Standard error and design effects for Comparisons
    vv1r <- sum((ipwir^2)*((yqcir - meanqcr[i])^2)) / ((sum(ipwir))^2)
    vv1r <- qed_nqcr[i]*vv1r/(qed_nqcr[i]-1)
    seqcr[i]   <- vv1r^.5
    
    wwqr        <- sum(ipwir^2)/((sum(ipwir))^2)  
    deffqcr[i]  <- wwqr*qed_nqcr[i]
    
    # Effective sample size needed for pooling
    qed_eff_sampr[i] <- (qed_nqcr[i]/deffqcr[i]) + qed_nqtr[i]
    qed_nqr[i]   <- qed_nqtr[i] + qed_nqcr[i]
    
  } # End loop over RCT-QED sites 
  
  # Calculate pooled std, mean, and se
  qed_s2t_poolr <- sum(qed_nqtr*qed_s2tr)/sum(qed_nqtr)
  qed_sdt_poolr <- qed_s2t_poolr^.5
  
  qed_varr <- qed_ser^2  
  qed_imp_poolr <- sum(qed_eff_sampr*qed_impr)/ sum(qed_eff_sampr)
  qed_var_poolr <- sum((qed_eff_sampr^2)*qed_varr) / ((sum(qed_eff_sampr))^2)
  qed_se_poolr  <- qed_var_poolr^.5
  
  qed_meanc_poolr  <- sum(qed_eff_sampr*meanqcr)/ sum(qed_eff_sampr)
  
  if (clus_dat == 1) {
    qed_var_clusr <- qed_se_clusr^2  
    qed_var_pool_clusr <- sum((qed_eff_sampr^2)*qed_var_clusr) / 
      ((sum(qed_eff_sampr))^2)
    qed_se_pool_clusr  <- qed_var_pool_clusr^.5
    
    qed_deff_clusr <- qed_var_pool_clusr / qed_var_poolr
    
    qed_ser  <- qed_se_clusr
    qed_varr <- qed_var_clusr
    qed_var_poolr <- qed_var_pool_clusr
    qed_se_poolr  <- qed_se_pool_clusr
  }
  
  if (clus_dat == 1) {
    nclusqr <- NROW(table(qed_datr$clusv))
    dofqr   <- nclusqr - 2*se_site - nadjxq
  } else if (clus_dat == 0) {
    dofqr <- sum(qed_nqr) - 2*se_site - nadjxq
  }
  
  if (dofqr <= 0) {
    dofqr <- 1
  }
  
  t_statqr  <- abs(qed_imp_poolr)/qed_se_poolr
  p_valueqr <- 2*pt(t_statqr, df = dofqr,  lower.tail=FALSE)
  
  # Print QED impact results
  qed_impzr  <- format(qed_impr,digits=4, nsmall=4)
  qed_sezr   <- format(qed_ser,digits=4, nsmall=4)
  eff_sampzr <- format(qed_eff_sampr,digits=2, nsmall=2)
  meanqtzr   <- format(meanqtr,digits=4, nsmall=4)
  meanqczr   <- format(meanqcr,digits=4, nsmall=4)
  seqtzr     <- format(seqtr,digits=4, nsmall=4)
  seqczr     <- format(seqcr,digits=4, nsmall=4)
  qed_imp_poolzr   <- format(qed_imp_poolr,digits=4, nsmall=4)
  qed_meanc_poolzr <- format(qed_meanc_poolr,digits=4, nsmall=4)
  qed_se_poolzr    <- format(qed_se_poolr,digits=4, nsmall=4)
  qed_sdt_poolzr   <- format(qed_sdt_poolr,digits=4, nsmall=4)
  t_statqzr        <- format(t_statqr,digits=2, nsmall=2)
  p_valueqzr       <- format(p_valueqr,digits=3, nsmall=3)
  
  ###
  # CONDUCT QED ANALYSIS FOR ALL QED SITES
  ###
  
  # Create qed dataset and new site indicator
  
  qed_datr$siteIDa <- qed_datr$siteID
  qed_datq$siteIDa <- qed_datq$siteID + se_site
  
  qed_datqa     <- rbind(qed_datr,qed_datq)
  
  # Create site indicators for regression models
  qed_datqa$siteIDf <- as.factor(qed_datqa$siteIDa)
  
  # Create regression formulas
  seq_site <- se_site + sq_site
  
  if (seq_site > 1) {
    if (adjxq == 1) {
      formula_qeda <- paste("yv", " ~ ", "siteIDf", "+", "tc:siteIDf", "+",
                           xvaq, "+", "0")
    } else if (adjxq == 0) {
      formula_qeda <- paste("yv", " ~ ", "siteIDf", "+", "tc:siteIDf", "+", "0")
    }
  } else if (seq_site == 1) {
    if (adjx == 1) {
      formula_qeda <- paste("yv", " ~ ", "tc", "+", xvaq)
    } else if (adjx == 0) {
      formula_qeda <- paste("yv", " ~ ", "tc")
    }
  }
  
  # Define design depending on whether clustering or not
  dsgn_speca    <- svydesign(id = ~1, weights = ~ipw_wt, data = qed_datqa)
  if (clus_dat == 1) {
    dsgn_spec_clusa <- svydesign(id = ~clusv, weights = ~ipw_wt, data = qed_datqa)
  }
  
  surv_rega <- svyglm(formula_qeda, design = dsgn_speca)
  s_outa    <- summary(surv_rega)
  qed_coefa <- s_outa$coefficients
  
  if (clus_dat == 1) {
    surv_reg_clusa <- svyglm(formula_qeda, design = dsgn_spec_clusa)
    s_out_clusa    <- summary(surv_reg_clusa)
    qed_coef_clusa <- s_out_clusa$coefficients
  }
  
  if (seq_site > 1) {
    nimpqa   <- nrow(qed_coefa)
    qed_impa <- qed_coefa[(seq_site+nadjxq+1):nimpqa,c("Estimate")]
    qed_sea  <- qed_coefa[(seq_site+nadjxq+1):nimpqa,c("Std. Error")]
    if (clus_dat == 1) {
      qed_se_clusa <- qed_coef_clusa[(seq_site+nadjxq+1):nimpqa,c("Std. Error")]
    }
  } else if (seq_site == 1) {
    qed_impa <- qed_coefa[1,c("Estimate")]
    qed_sea  <- qed_coefa[1,c("Std. Error")]
    if (clus_dat == 1) {
      qed_se_clusa  <- qed_coef_clusa[1,c("Std. Error")]
    }
  }
  
  # Run over each QED site to get design effects and eff samp sizes
  # Some of this is old to get impacts and std. errs separately by site
  # Yields very similar results to svyglm
  meanqta      <- rep(0,seq_site)
  s2yqta       <- rep(0,seq_site)
  seqta        <- rep(0,seq_site)
  
  meanqca      <- rep(0,seq_site)
  s2yqca       <- rep(0,seq_site)
  seqca        <- rep(0,seq_site)
  deffqca      <- rep(0,seq_site)
  
  qed_eff_sampa <- rep(0,seq_site)
  qed_nqa      <- rep(0,seq_site)
  qed_nqta     <- rep(0,seq_site)   
  qed_nqca     <- rep(0,seq_site)
  
  qed_s2ta     <- rep(0,seq_site)
  
  for (i in 1:seq_site) {
    
    yqtia <- qed_datqa$yv[qed_datqa$siteIDa==i & qed_datqa$tc == 1]
    yqcia <- qed_datqa$yv[qed_datqa$siteIDa==i & qed_datqa$tc == 0]
    ipwia <- qed_datqa$ipw_wt[qed_datqa$siteIDa==i & qed_datqa$tc == 0]
    
    # Calculate means and standard errors for Ts
    meanqta[i]  <- mean(yqtia)
    qed_nqta[i] <- length(yqtia)
    s2yqta[i]   <- sum((yqtia-meanqta[i])^2)/(qed_nqta[i]-1)
    seqta[i]    <- (s2yqta[i] / qed_nqta[i])^.5
    
    # Calculate treatment group s2 - assumed the same for CACE and ITT
    qed_s2ta[i] <- var(yqtia)
    
    # Calculate means and standard errors for Comparisons
    # Weighted means for Comparisons
    meanqca[i] <- sum(ipwia*yqcia) / sum(ipwia)
    
    # Weighted variances for Comparisons
    qed_nqca[i] <- length(yqcia)
    vva         <- sum(ipwia*((yqcia - meanqca[i])^2)) / sum(ipwia)
    s2yqca[i]   <- qed_nqca[i]*vva/(qed_nqca[i]-1)
    
    # Standard error and design effects for Comparisons
    vv1a <- sum((ipwia^2)*((yqcia - meanqca[i])^2)) / ((sum(ipwia))^2)
    vv1a <- qed_nqca[i]*vv1a/(qed_nqca[i]-1)
    seqca[i]   <- vv1a^.5
    
    wwqa        <- sum(ipwia^2)/((sum(ipwia))^2)  
    deffqca[i]  <- wwqa*qed_nqca[i]
    
    # Effective sample size needed for pooling
    qed_eff_sampa[i] <- (qed_nqca[i]/deffqca[i]) + qed_nqta[i]
    qed_nqa[i]   <- qed_nqta[i] + qed_nqca[i]
    
  } # End loop over QED sites 
  
  # Calculate pooled std, mean, and se
  qed_s2t_poola <- sum(qed_nqta*qed_s2ta)/sum(qed_nqta)
  qed_sdt_poola <- qed_s2t_poola^.5
  
  qed_vara <- qed_sea^2  
  qed_imp_poola <- sum(qed_eff_sampa*qed_impa)/ sum(qed_eff_sampa)
  qed_var_poola <- sum((qed_eff_sampa^2)*qed_vara) / ((sum(qed_eff_sampa))^2)
  qed_se_poola  <- qed_var_poola^.5
  
  qed_meanc_poola  <- sum(qed_eff_sampa*meanqca)/ sum(qed_eff_sampa)
  
  if (clus_dat == 1) {
    qed_var_clusa <- qed_se_clusa^2  
    qed_var_pool_clusa <- sum((qed_eff_sampa^2)*qed_var_clusa) / 
      ((sum(qed_eff_sampa))^2)
    qed_se_pool_clusa  <- qed_var_pool_clusa^.5
    
    qed_deff_clusa <- qed_var_pool_clusa / qed_var_poola
    
    qed_sea  <- qed_se_clusa
    qed_vara <- qed_var_clusa
    qed_var_poola <- qed_var_pool_clusa
    qed_se_poola  <- qed_se_pool_clusa
  }
  
  if (clus_dat == 1) {
    nclusqa <- NROW(table(qed_datqa$clusv))
    dofqa   <- nclusqa - 2*seq_site - nadjxq
  } else if (clus_dat == 0) {
    dofqa <- sum(qed_nqa) - 2*seq_site - nadjxq
  }
  
  if (dofqa <= 0) {
    dofqa <- 1
  }
  
  t_statqa  <- abs(qed_imp_poola)/qed_se_poola
  p_valueqa <- 2*pt(t_statqa, df = dofqa,  lower.tail=FALSE)
  
  # Print QED impact results
  qed_impza  <- format(qed_impa,digits=4, nsmall=4)
  qed_seza   <- format(qed_sea,digits=4, nsmall=4)
  eff_sampza <- format(qed_eff_sampa,digits=2, nsmall=2)
  meanqtza   <- format(meanqta,digits=4, nsmall=4)
  meanqcza   <- format(meanqca,digits=4, nsmall=4)
  seqtza     <- format(seqta,digits=4, nsmall=4)
  seqcza     <- format(seqca,digits=4, nsmall=4)
  qed_imp_poolza   <- format(qed_imp_poola,digits=4, nsmall=4)
  qed_meanc_poolza <- format(qed_meanc_poola,digits=4, nsmall=4)
  qed_se_poolza    <- format(qed_se_poola,digits=4, nsmall=4)
  qed_sdt_poolza   <- format(qed_sdt_poola,digits=4, nsmall=4)
  t_statqza        <- format(t_statqa,digits=2, nsmall=2)
  p_valueqza       <- format(p_valueqa,digits=3, nsmall=3)
  
  ###
  # CALCULATE BIAS VARIANCE
  ###
  
  # Create bias_grp number that goes from 1 to number of bias tree leaves
  leaves_df <- data.frame(table(rct_dat$left_right))
  n_leaves   <- nrow(leaves_df)
  
  leaves_df$left_right <- leaves_df$Var1
  leaves_df$bias_grp   <- seq_len(nrow(leaves_df))
  leaves_df            <- leaves_df[,c("left_right","bias_grp")]
  
  # Merge leaves onto rct_dat and qed_datrc datasets to 
  # calculate biases and their variances in the RCT-QED sites:
  # Use r suffix for specific tree nodes
  rct_datt    <- merge(rct_datt, leaves_df, by = "left_right", all.x = TRUE)
  rct_datc    <- merge(rct_datc, leaves_df, by = "left_right", all.x = TRUE)
  rct_datc_np <- merge(rct_datc_np, leaves_df, by = "left_right", all.x = TRUE)
  rct_datt_p  <- merge(rct_datt_p, leaves_df, by = "left_right", all.x = TRUE)
  rct_datt_np <- merge(rct_datt_np, leaves_df, by = "left_right", all.x = TRUE)
  qed_datrc   <- merge(qed_datrc, leaves_df, by = "left_right", all.x = TRUE)
  
  bias_node      <- rep(0,n_leaves) 
  var_bias_node  <- rep(0,n_leaves)
  mwrsu          <- rep(0,n_leaves) 
  
  for (i in 1:n_leaves) {
    
    meanccr  <- rep(0,se_site)
    seccr    <- rep(0,se_site)
    dbartr   <- rep(0,se_site)
    dbarcr   <- rep(0,se_site)
    dbarr    <- rep(0,se_site)
    meantnpr <- rep(0,se_site)  
    meancnpr <- rep(0,se_site)
    
    meanqcr  <- rep(0,se_site)
    s2yqcr   <- rep(0,se_site)
    seqcr    <- rep(0,se_site)
    deffqc1r <- rep(0,se_site)
    deffqcr  <- rep(0,se_site)
    
    nyqcr    <- rep(0,se_site)
    nycr     <- rep(0,se_site)
    nytr     <- rep(0,se_site)
    nycnpr   <- rep(0,se_site)
    nytnpr   <- rep(0,se_site)
    nytpr    <- rep(0,se_site)
    
    wgtr     <- rep(0,se_site)
    
    site_bias <- rep(0,se_site)
    
    # Loop over RCT-QED sites
    for (j in 1:se_site) {
      
      # Read in RCT and QED data in node and site
      
      yqcr <- qed_datrc[qed_datrc$bias_grp==i & qed_datrc$siteID==j,c("yv")]
      ipwr <- qed_datrc[qed_datrc$bias_grp==i & qed_datrc$siteID==j,c("ipw_wt")]
      
      ycr   <- rct_datc[rct_datc$bias_grp==i & rct_datc$siteID==j,c("yv","partic")]
      ytr   <- rct_datt[rct_datt$bias_grp==i & rct_datt$siteID==j,c("yv","partic")]
      ycnpr <- rct_datc_np[rct_datc_np$bias_grp==i & rct_datc_np$siteID==j,c("yv")]
      ytnpr <- rct_datt_np[rct_datt_np$bias_grp==i & rct_datt_np$siteID==j,c("yv")]
      ytpr  <- rct_datt_p[rct_datt_p$bias_grp==i & rct_datt_p$siteID==j,c("yv")]
      
      nyqcr[j]  <- length(yqcr)    # QED comparisons
      nycr[j]   <- length(ycr$yv)  # RCT Cs
      nytr[j]   <- length(ytr$yv)  # RCT Ts 
      nycnpr[j] <- length(ycnpr)   # RCT C NPs
      nytnpr[j] <- length(ytnpr)   # RCT T NPs
      nytpr[j]  <- length(ytpr)    # RCT T Ps
      
      # Conduct CACE analysis only if there are enough RCT C NPs and T Ps 
      # and QED comparisons in the node for variance calculations
      
      conduct_anal <- 0
      if (cace_itt == 1) {
        if ((nyqcr[j] > 1) & (nycnpr[j] > 1) & (nytpr[j] > 0)) {
          conduct_anal <- 1
        } 
      } 
      
      if (cace_itt == 0) {
        if ((nycr[j] > 1) & (nyqcr[j] > 1)) {
          conduct_anal <- 1
        }
      }
      
      if (conduct_anal == 1) {
        
        # Calculate RCT control mean and variance if cace_itt == 0 
        if (cace_itt == 0) {
          meanccr[j] <- mean(ycr$yv)
          varccr     <- var(ycr$yv)/nycr[j]
          
          if (varccr > 0) {
            seccr[j] <- varccr^.5
          } else {
            seccr[j] <- NA
          }
        } else if (cace_itt == 1) {
          
          # ***
          # *** Calculate RCT complier mean and standard errs if cace_itt == 1
          # ***
          
          # Calculate RCT participation rates for Ts and Cs: dbart and dbarc
          # Note that we allow for crossovers
          dbartr[j] <- mean(ytr$partic)
          dbarcr[j] <- mean(ycr$partic)
          
          dbarr[j]  <- dbartr[j] - dbarcr[j]
          
          # Only proceed if the compliance rate is positive
          if (dbarr[j] > 0) {
            
            # Calculate mean for RCT control nonparticipants
            meancnpr[j] <- mean(ycnpr)
            
            # Calculate mean for RCT treatment nonparticipants if there are 
            # at least 2 people - if not set the mean to 0
            if (nytnpr[j] >= 2) {
              meantnpr[j] <- mean(ytnpr) 
            } else {
              meantnpr[j] <- 0 
            }
            
            # Calculate RCT control complier mean and variance 
            meanccr[j] <- ( ((1-dbarcr[j])*meancnpr[j]) 
                            - ((1-dbartr[j])*meantnpr[j]) )/dbarr[j]
            
            varcnpr  <- var(ycnpr)/nycnpr[j]
            
            if (nytnpr[j] >= 2) {
              vartnpr <- var(ytnpr)/nytnpr[j] 
            } else {
              vartnpr <- 0 
            }
            
            varccr <- ( (((1-dbarcr[j])^2)*varcnpr) 
                        + (((1-dbartr[j])^2)*vartnpr) )/(dbarr[j]^2)
            
            if (varccr > 0) {
              seccr[j] <- varccr^.5
            } else {
              seccr[j] <- NA
            }
            
          } else if (dbarr[j] <= 0) {
            
            meanccr[j] <- NA
            seccr[j]   <- NA
            
          } # end if dbarr > 0
          
        } # end if cace_itt == 1
        
        # ***
        # QED comparison group means and standard errs in node  
        # ***
        
        # Weighted means for comparisons in node
        meanqcr[j] <- sum(ipwr*yqcr) / sum(ipwr)
        
        # Weighted s2 variances for QED Comparisons in node
        vv        <- sum(ipwr*((yqcr - meanqcr[j])^2)) / sum(ipwr)
        s2yqcr[j] <- nyqcr[j]*vv/(nyqcr[j]-1)
        
        # Standard error and design effects for Comparisons in node
        vv1 <- sum((ipwr^2)*((yqcr - meanqcr[j])^2)) / ((sum(ipwr))^2)
        vv1 <- nyqcr[j]*vv1/(nyqcr[j]-1)
        seqcr[j]   <- vv1^.5
        
        wwqr        <- sum(ipwr^2)/((sum(ipwr))^2)  
        deffqcr[j]  <- wwqr*nyqcr[j]
        
        # wgtrs needed for pooling
        wgtr[j]     <- nycr[j]*dbarr[j] + (nyqcr[j]/deffqcr[j]) 
        
      } else if (conduct_anal == 0) {
        
        meanccr[j]  <- NA
        seccr[j]    <- NA
        dbartr[j]   <- NA
        dbarcr[j]   <- NA
        meantnpr[j] <- NA  
        meancnpr[j] <- NA
        
        meanqcr[j]  <- NA
        s2yqcr[j]   <- NA
        seqcr[j]    <- NA
        deffqcr[j]  <- NA
        wgtr[j]     <- NA
        
        nyqcr[j]    <- NA
        nycr[j]     <- NA
        nytr[j]     <- NA 
        nycnpr[j]   <- NA
        nytnpr[j]   <- NA
        nytpr[j]    <- NA
        
      }  # End conduct_anal
      
    } # End j site loop
    
    # Calculate biases and their variances in leaf across the RCT sites
    bias_rs     <- meanccr - meanqcr
    var_bias_rs <- (seccr^2) + (seqcr^2) 
    
    # Calculate weighted biases and variances using wgtr weight
    # but only if there are some biases that are not NA
    nbiasg <- sum(!is.na(var_bias_rs))
    
    if (nbiasg > 0) {
      swr <- sum(wgtr,na.rm = TRUE)
      bias_node[i]     <- sum(wgtr*bias_rs,na.rm = TRUE)/swr 
      var_bias_node[i] <- sum((wgtr^2)*var_bias_rs,na.rm = TRUE)/(swr^2) 
      mwrsu[i] <- mean(wgtr,na.rm = TRUE)
    } else {
      bias_node[i]     <- NA 
      var_bias_node[i] <- NA
      mwrsu[i]         <- NA
    }
    
    # Print bias stuff out
    meanccrz    <- format(meanccr,digits=4, nsmall=4)
    seccrz      <- format(seccr,digits=4, nsmall=4)
    dbarrz      <- format(dbarr,digits=4, nsmall=4)
    dbartrz     <- format(dbartr,digits=4, nsmall=4)
    dbarcrz     <- format(dbarcr,digits=4, nsmall=4)
    wgtrz       <- format(wgtr,digits=4, nsmall=4)
    
    meanqcrz    <- format(meanqcr,digits=4, nsmall=4)
    seqcrz      <- format(seqcr,digits=4, nsmall=4)
    deffqcrz    <- format(deffqcr,digits=4, nsmall=4)
    bias_rsz    <- format(bias_rs,digits=4, nsmall=4)
    var_bias_rsz  <- format(var_bias_rs,digits=4, nsmall=4)
    
  }  # End i grp loop
  
  bias_nodez      <- format(bias_node,digits=4, nsmall=4)
  var_bias_nodez  <- format(var_bias_node,digits=4, nsmall=4)
  se_bias_node    <- var_bias_node^.5
  se_bias_nodez   <- format(se_bias_node,digits=4, nsmall=4)
  
  # Merge var_bias_node to leaves and qed_datqc and qt datasets
  leaves_df$bias_node <- bias_node
  leaves_df$var_bias_node <- var_bias_node
  
  qed_datqc <- merge(qed_datqc, leaves_df, by = "left_right", all.x = TRUE)
  qed_datqt <- merge(qed_datqt, leaves_df, by = "left_right", all.x = TRUE)
  
  # Calculate bias-adjusted QED comparison mean at individual level
  qed_datqc$yq_adj <- qed_datqc$yv - qed_datqc$pred_bias #qed_datqc$bias_node  
  
  # Calculate mean bias and variance at the site level using IPW weights
  # in QED-only sites using QED comparisons
  # This uses (D.7) in the paper
  
  mean_biasq <- 0
  var_biasq  <- 0
  deff_biasq <- 0
  var_u_node  <- rep(0,n_leaves)
  var_biasqu <- 0
  
  for (r in 1:n_leaves) {
    
    vvrk    <- 0
    swqrku <- 0
    
    for (k in 1:sq_site) {
      
      # Pull off C and T outcomes in bias tree leaf r and QED-only site k
      anyrkc <- qed_datqc[qed_datqc$bias_grp==r & qed_datqc$siteID==k,c("yv")]
      anyrkt <- qed_datqt[qed_datqt$bias_grp==r & qed_datqt$siteID==k,c("yv")]
      
      if (length(anyrkc) > 0) {
        
        # Pull off the C datasets in rk and site k 
        qbrk <- qed_datqc[qed_datqc$bias_grp==r & qed_datqc$siteID==k,]
        qbk  <- qed_datqc[qed_datqc$siteID==k,]
        
        # wgtk is the wQEDk weight in the paper
        wgtk <- qed_eff_samp[k]
        
        # sumwbr is sum of wgtk*I*wstar / sum(wstar) over Cs in rk  
        sumwbr <- sum(wgtk*qbrk$ipw_wt)/sum(qbk$ipw_wt)
        
        # mean_biasq is calculating bhatQED in (D.6) in the paper (pooled bias) 
        mean_biasq <- mean_biasq + sumwbr*bias_node[r]
        
        # vvrk sums over k needed to get the variance of bhatQED 
        vvrk <- vvrk + sumwbr
        
        # This is for Var(u) for forecasts
        # Calculate design effect needed to define weight in denominator of Var(ur)
        # in (D.15) in the paper
        n0qrku <- length(anyrkc)
        wwqrku <- sum(qbrk$ipw_wt^2)/((sum(qbrk$ipw_wt))^2)
        deffqrku  <- wwqrku*n0qrku
        
        n1qrku <- length(anyrkt)
        
        if ((n0qrku > 0) & (n1qrku > 0)) {
          effqrku <- n1qrku + (n0qrku/deffqrku)
          swqrku <- swqrku + effqrku
        } 
        
      }  # anyrkc
      
    } # end k QED-only sites
    
    # Add in Var(ur) to var_bias_node[r] using (D.15) in paper
    if (swqrku > 0) {
      
      var_u_node[r] <- se*mwrsu[r]*var_bias_node[r]/swqrku
    }
    
    # Calculate variances in (D.7) with and without the forecast error 
    var_biasq  <- var_biasq + (vvrk^2)*var_bias_node[r]
    deff_biasq <- deff_biasq + (vvrk^2)
    
    var_biasqu <- var_biasqu + (vvrk^2)*(var_bias_node[r] + var_u_node[r])
    
  } # end r nodes
  
  # Adjust denominators of pooled biases and variances by sum of site weights
  bias_poolq <- mean_biasq / sum(qed_eff_samp)
  
  bias_var_poolq  <- var_biasq / ((sum(qed_eff_samp))^2)
  bias_var_poolqu <- var_biasqu / ((sum(qed_eff_samp))^2)
  
  # If clustered design, multiply this bias variance by the average of the 
  # design effects due to clustering for the RCT and QED-only estimates
  if (clus_dat == 1) {
    avg_deff_clus  <- (rct_deff_clus + qed_deff_clus)/2
    bias_var_poolq  <- bias_var_poolq*avg_deff_clus
    bias_var_poolqu <- bias_var_poolqu*avg_deff_clus
  }
  
  bias_se_poolq   <- bias_var_poolq^.5
  bias_se_poolqu  <- bias_var_poolqu^.5
  
  ###
  # Pool RCT, QED, and Bias findings
  ###
  
  tot_qed_imp   <- qed_imp_pool - bias_poolq
  tot_qed_meanc <- qed_meanc_pool + bias_poolq
  
  tot_qed_var  <- qed_var_pool + bias_var_poolq
  tot_qed_varu <- qed_var_pool + bias_var_poolqu
  
  tot_qed_se   <- tot_qed_var^.5
  tot_qed_seu  <- tot_qed_varu^.5
  
  tot_t_statq  <- abs(tot_qed_imp)/tot_qed_se
  tot_p_valueq <- 2*pt(tot_t_statq, df = dofq,  lower.tail=FALSE)
  
  tot_t_statqu  <- abs(tot_qed_imp)/tot_qed_seu
  tot_p_valuequ <- 2*pt(tot_t_statqu, df = dofq,  lower.tail=FALSE)
  
  # Weight RCT and Adjusted QED estimates
  # Calculate weight assigned to RCT and QED components
  
  # std
  www_sdt <- sum(rct_net)/(sum(rct_net)+sum(qed_nqt))
  s2t_rct_adjqed <- www_sdt*rct_s2t_pool + (1-www_sdt)*qed_s2t_pool
  sdt_rct_adjqed <- s2t_rct_adjqed^.5
  
  if (inv_var_agg_wgt == 1) {
    www_rct  <- (tot_qed_var)/(tot_qed_var + rct_var_pool)
    www_rctu <- (tot_qed_varu)/(tot_qed_varu + rct_var_pool)
    www_rctn <- (qed_var_pool)/(qed_var_pool + rct_var_pool)
  } else if (inv_var_agg_wgt == 0) {
    www_rct  <- sum(rct_ne) / (sum(rct_ne) + sum(qed_eff_samp))
    www_rctu <- www_rct
    www_rctn <- sum(rct_ne) / (sum(rct_ne) + sum(qed_eff_samp))
  } 
  
  # Calculate combined estimates across rct and adj QED-only sites
  meanc_rct_adjqed <- www_rct*rct_meanc_pool + (1-www_rct)*tot_qed_meanc
  imp_rct_adjqed   <- www_rct*rct_imp_pool + (1-www_rct)*tot_qed_imp
  var_rct_adjqed   <- (www_rct^2)*rct_var_pool + ((1-www_rct)^2)*tot_qed_var
  se_rct_adjqed    <- var_rct_adjqed^.5
  
  var_rct_adjqedu  <- (www_rctu^2)*rct_var_pool + ((1-www_rctu)^2)*tot_qed_varu
  se_rct_adjqedu   <- var_rct_adjqedu^.5
  
  # Combined without bias adjusment
  meanc_rct_noadjqed <- www_rct*rct_meanc_pool + (1-www_rct)*qed_meanc_pool
  imp_rct_noadjqed   <- www_rctn*rct_imp_pool + (1-www_rctn)*qed_imp_pool
  var_rct_noadjqed   <- (www_rctn^2)*rct_var_pool + ((1-www_rctn)^2)*qed_var_pool
  se_rct_noadjqed    <- var_rct_noadjqed^.5
  
  # Hypothesis testing for combined RCT and adjusted QED-only impact estimates
  if (clus_dat == 1) {
    dof_tot  <- ncluse + nclusq - 2*se_site - 2*sq_site - nadjx - nadjxq
  } else if (clus_dat == 0) {
    dof_tot <- rct_wgt_sum + sum(qed_nq)- 2*se_site - 2*sq_site - 
      nadjx - nadjxq
  }
  
  if (dof_tot <= 0) {
    dof_tot <- 1
  }
  
  t_stat_tot   <- abs(imp_rct_adjqed)/se_rct_adjqed
  p_value_tot  <- 2*pt(t_stat_tot, df = dof_tot,  lower.tail=FALSE)
  t_stat_totu  <- abs(imp_rct_adjqed)/se_rct_adjqedu
  p_value_totu <- 2*pt(t_stat_totu, df = dof_tot,  lower.tail=FALSE)
  
  t_stat_totn   <- abs(imp_rct_noadjqed)/se_rct_noadjqed
  p_value_totn  <- 2*pt(t_stat_totn, df = dof_tot,  lower.tail=FALSE)
  
  control_meanz    <- format(meanc_rct_adjqed,digits=4, nsmall=4)
  agg_impactz      <- format(imp_rct_adjqed,digits=4, nsmall=4)
  senou_impactz    <- format(se_rct_adjqed,digits=4, nsmall=4)
  se_impactz       <- format(se_rct_adjqedu,digits=4, nsmall=4)
  t_statnouz       <- format(t_stat_tot,digits=2, nsmall=2)
  p_valuenouz      <- format(p_value_tot,digits=3, nsmall=3)
  t_statz          <- format(t_stat_totu,digits=2, nsmall=2)
  p_valuez         <- format(p_value_totu,digits=3, nsmall=3)
  sd_outcomez      <- format(sdt_rct_adjqed,digits=4, nsmall=4)
  dofz             <- round(dof_tot)
  
  meanc_rct_noadjqedz <- format(meanc_rct_noadjqed,digits=4, nsmall=4)
  imp_rct_noadjqedz <- format(imp_rct_noadjqed,digits=4, nsmall=4)
  se_rct_noadjqedz  <- format(se_rct_noadjqed,digits=4, nsmall=4)
  t_stat_totnz      <- format(t_stat_totn,digits=2, nsmall=2)
  p_value_totnz     <- format(p_value_totn,digits=3, nsmall=3)
  
  tot_qed_meancz <- format(tot_qed_meanc,digits=4, nsmall=4)
  tot_qed_impz   <- format(tot_qed_imp,digits=4, nsmall=4)
  tot_qed_seuz   <- format(tot_qed_seu,digits=4, nsmall=4)
  tot_t_statqz   <- format(tot_t_statqu,digits=2, nsmall=2)
  tot_p_valueqz  <- format(tot_p_valuequ,digits=3, nsmall=3)
  
  ###
  # Print stuff out
  ###
  
  # Aggregated
  imp_res <- data.frame(control_meanz,
                        agg_impactz,
                        se_impactz,
                        sd_outcomez,
                        t_statz, 
                        p_valuez)
  
  imp_res$ctrl_mean    <- imp_res$control_meanz
  imp_res$mean_bias    <- NA
  imp_res$impact       <- imp_res$agg_impactz
  imp_res$std_err      <- imp_res$se_impactz
  imp_res$sd_out       <- imp_res$sd_outcomez
  imp_res$t_stat       <- imp_res$t_statz
  imp_res$p_value      <- imp_res$p_valuez
  
  imp_res$Design       <- c("Pooled RCT & Bias-Adj QED Imps")
  
  # RCT pooled results
  imp_res_rct <- data.frame(rct_meanc_poolz, rct_imp_poolz,rct_se_poolz,
                            rct_sdt_poolz,t_statrz, p_valuerz)
  
  imp_res_rct$ctrl_mean    <- imp_res_rct$rct_meanc_poolz
  imp_res_rct$mean_bias    <- NA
  imp_res_rct$impact       <- imp_res_rct$rct_imp_poolz
  imp_res_rct$std_err      <- imp_res_rct$rct_se_poolz
  imp_res_rct$sd_out       <- imp_res_rct$rct_sdt_poolz
  imp_res_rct$t_stat       <- imp_res_rct$t_statrz
  imp_res_rct$p_value      <- imp_res_rct$p_valuerz
  
  imp_res_rct$Design       <- c("RCT Imp in RCT-QED Sites")
  
  # QED in QED-only sites Unadjusted
  imp_res_qed <- data.frame(qed_meanc_poolz, qed_imp_poolz, qed_se_poolz,
                            qed_sdt_poolz, t_statqz, p_valueqz)
  
  imp_res_qed$ctrl_mean    <- imp_res_qed$qed_meanc_poolz
  imp_res_qed$mean_bias    <- NA
  imp_res_qed$impact       <- imp_res_qed$qed_imp_poolz
  imp_res_qed$std_err      <- imp_res_qed$qed_se_poolz
  imp_res_qed$sd_out       <- imp_res_qed$qed_sdt_poolz
  imp_res_qed$t_stat       <- imp_res_qed$t_statqz
  imp_res_qed$p_value      <- imp_res_qed$p_valueqz
  
  imp_res_qed$Design       <- c("Unadj QED Imp in QED-Only Sites")
  
  # QED in QED-only sites adjusted
  imp_res_qed_adj <- data.frame(tot_qed_meancz, tot_qed_impz, tot_qed_seuz,
                                qed_sdt_poolz, tot_t_statqz, tot_p_valueqz)
  
  imp_res_qed_adj$ctrl_mean    <- imp_res_qed_adj$tot_qed_meancz
  imp_res_qed_adj$mean_bias    <- NA
  imp_res_qed_adj$impact       <- imp_res_qed_adj$tot_qed_impz
  imp_res_qed_adj$std_err      <- imp_res_qed_adj$tot_qed_seuz
  imp_res_qed_adj$sd_out       <- imp_res_qed_adj$qed_sdt_poolz
  imp_res_qed_adj$t_stat       <- imp_res_qed_adj$tot_t_statqz
  imp_res_qed_adj$p_value      <- imp_res_qed_adj$tot_p_valueqz
  
  imp_res_qed_adj$Design       <- c("Adj QED Imp in QED-Only Sites")
  
  # QED in RCT-QED sites Unadjusted
  imp_res_qedr <- data.frame(qed_meanc_poolzr, qed_imp_poolzr, qed_se_poolzr,
                            qed_sdt_poolzr, t_statqzr, p_valueqzr)
  
  imp_res_qedr$ctrl_mean    <- imp_res_qedr$qed_meanc_poolzr
  imp_res_qedr$mean_bias    <- NA
  imp_res_qedr$impact       <- imp_res_qedr$qed_imp_poolzr
  imp_res_qedr$std_err      <- imp_res_qedr$qed_se_poolzr
  imp_res_qedr$sd_out       <- imp_res_qedr$qed_sdt_poolzr
  imp_res_qedr$t_stat       <- imp_res_qedr$t_statqzr
  imp_res_qedr$p_value      <- imp_res_qedr$p_valueqzr
  
  imp_res_qedr$Design       <- c("Unadj QED Imp in RCT-QED Sites")
  
  # QED in All QED sites Unadjusted
  imp_res_qeda <- data.frame(qed_meanc_poolza, qed_imp_poolza, qed_se_poolza,
                             qed_sdt_poolza, t_statqza, p_valueqza)
  
  imp_res_qeda$ctrl_mean    <- imp_res_qeda$qed_meanc_poolza
  imp_res_qeda$mean_bias    <- NA
  imp_res_qeda$impact       <- imp_res_qeda$qed_imp_poolza
  imp_res_qeda$std_err      <- imp_res_qeda$qed_se_poolza
  imp_res_qeda$sd_out       <- imp_res_qeda$qed_sdt_poolza
  imp_res_qeda$t_stat       <- imp_res_qeda$t_statqza
  imp_res_qeda$p_value      <- imp_res_qeda$p_valueqza
  
  imp_res_qeda$Design       <- c("Unadj QED Imp in All QED Sites")
  
  # Bias
  if (clus_dat == 1) {
    dof_bias  <- ncluse - 2*se_site
  } else if (clus_dat == 0) {
    dof_bias <- rct_wgt_sum - 2*se_site
  }
  
  if (dof_bias <= 0) {
    dof_bias <- 1
  }
  
  dof_bias     <- rct_wgt_sum - 2*se_site 
  t_stat_bias  <- abs(bias_poolq)/bias_se_poolqu
  p_value_bias <- 2*pt(t_stat_bias, df = dof_bias, lower.tail=FALSE)
  
  bias_poolqz      <- format(bias_poolq,digits=4, nsmall=4)
  bias_se_poolquz  <- format(bias_se_poolqu,digits=4, nsmall=4)
  t_stat_biasz     <- format(t_stat_bias,digits=2, nsmall=2)
  p_value_biasz    <- format(p_value_bias,digits=3, nsmall=3)
  dofz_bias        <- round(dof_bias)
  
  res_bias <- data.frame(bias_poolqz,bias_se_poolquz,t_stat_biasz,p_value_biasz)
  
  res_bias$ctrl_mean    <- NA
  res_bias$mean_bias    <- res_bias$bias_poolqz
  res_bias$impact       <- res_bias$bias_poolqz
  res_bias$std_err      <- res_bias$bias_se_poolquz
  res_bias$sd_out       <- NA
  res_bias$t_stat       <- res_bias$t_stat_biasz
  res_bias$p_value      <- res_bias$p_value_biasz
  
  res_bias$Design        <- c("QED Bias in QED-only sites")
  
  res_biasf <- res_bias[,c("Design","ctrl_mean","impact","std_err",
                           "sd_out", "t_stat", "p_value")]
  
  imp_resf  <- imp_res[,c("Design","ctrl_mean","impact","std_err",
                          "sd_out", "t_stat", "p_value")]
  
  imp_res_rctf  <- imp_res_rct[,c("Design","ctrl_mean","impact",
                                  "sd_out", "std_err","t_stat", "p_value")]
  
  imp_res_qedf  <- imp_res_qed[,c("Design","ctrl_mean","impact",
                                  "sd_out", "std_err","t_stat", "p_value")]
  
  imp_res_qed_adjf  <- imp_res_qed_adj[,c("Design","ctrl_mean","impact",
                                          "sd_out", "std_err","t_stat", "p_value")]
  
  imp_res_qedrf  <- imp_res_qedr[,c("Design","ctrl_mean","impact",
                                    "sd_out", "std_err","t_stat", "p_value")]
  
  imp_res_qedaf  <- imp_res_qeda[,c("Design","ctrl_mean","impact",
                                    "sd_out", "std_err","t_stat", "p_value")]
  
  imp_resg <- rbind(imp_resf,imp_res_rctf)
  imp_resg <- rbind(imp_resg,imp_res_qedf)
  imp_resg <- rbind(imp_resg,imp_res_qed_adjf)
  imp_resg <- rbind(imp_resg,res_biasf)
  imp_resg <- rbind(imp_resg,imp_res_qedrf)
  imp_resg <- rbind(imp_resg,imp_res_qedaf)
  
  sink(file = out_txt2, split = TRUE, append = TRUE)
  
  cat("\n")
  if (cace_itt == 1) {
    pr1 <- "TABLE 3. HYBRID_IMPACT FUNCTION: OVERALL CACE IMPACT FINDINGS, BY DESIGN" 
    cat(pr1,"\n")
  } else if (cace_itt == 0) {
    pr1 <- "TABLE 3. HYBRID_IMPACT FUNCTION: OVERALL ITT IMPACT FINDINGS, BY DESIGN" 
    cat(pr1,"\n")
  }
  
  cat("\n")
  print(imp_resg[1:5,c("Design","ctrl_mean","impact",
                       "std_err","t_stat", "p_value")],row.names = FALSE,
        right='F')
  fn1 <- "Notes: The ctrl_mean column is the mean of the outcome variable for the control"
  fn2 <- "and/or comparison group. The impact and std_err columns are estimated impacts and"  
  fn3 <- "standard errors in nominal units, where Table 5 in the output file shows the" 
  fn4 <- "standard deviation of the outcome to convert the impacts into effect size units."  
  fn5 <- "The remaining columns are t-statistics and p-values. The impacts may not align across"
  fn6 <- "designs due to differences in the samples and their relative weights for pooling"
  fn7 <- "the site estimates."
  
  cat("\n")
  cat(fn1,"\n",fn2,"\n",fn3,"\n",fn4,"\n",fn5,"\n",fn6,"\n",fn7,"\n")
  
  sink()
  
  sink(file = out_txt2, append = TRUE)
  
  cat("\n")
  pr1 <- "ADDITIONAL HYBRID_IMPACT FUNCTION RESULTS" 
  cat(pr1,"\n")
  blank <- c("")
  
  cat("\n")
  if (cace_itt == 1) {
    pr1 <- "TABLE 4. FULL SET OF CACE IMPACT FINDINGS, BY DESIGN" 
    cat(pr1,"\n")
  } else if (cace_itt == 0) {
    pr1 <- "TABLE 4. FULL SET OF ITT IMPACT FINDINGS, BY DESIGN" 
    cat(pr1,"\n")
  }
  
  cat("\n")
  print(imp_resg[,c("Design","ctrl_mean","impact",
                    "std_err","t_stat", "p_value")],row.names = FALSE,
        right='F')
  fn1 <- "Notes: The ctrl_mean column is the mean of the outcome variable for the control"
  fn2 <- "and/or comparison group. The impact and std_err columns are estimated impacts and"  
  fn3 <- "standard errors in nominal units, where Table 5 in the output file shows the" 
  fn4 <- "standard deviation of the outcome to convert the impacts into effect size units."  
  fn5 <- "The remaining columns are t-statistics and p-values. The impacts may not align across"
  fn6 <- "designs due to differences in the samples and their relative weights for pooling"
  fn7 <- "the site estimates."
  
  cat("\n")
  cat(fn1,"\n",fn2,"\n",fn3,"\n",fn4,"\n",fn5,"\n",fn6,"\n",fn7,"\n")
  
  cat("\n")
  pr1 <- "TABLE 5. STANDARD DEVIATION OF THE OUTCOME, BY DESIGN" 
  cat(pr1,"\n")
  
  cat("\n")
  print(imp_resg[c(1, 2, 3, 6, 7),c("Design","sd_out")],row.names = FALSE,
        right='F')
  fn1 <- "Notes: The sd_out column shows the standard deviation of the outcome to convert the"
  fn2 <- "impacts in Tables 3 and 4 into effect size units. It is assumed the CACE and ITT" 
  fn3 <- "standard deviations are the same, which are obtained using the treatment"
  fn4 <- "group to facilitate the calculations."
  
  cat("\n")
  cat(fn1,"\n",fn2,"\n",fn3,"\n",fn4,"\n")
  
  cat("\n")
  cat("TABLE 6. RCT FINDINGS IN RCT-QED SITES")
  s_str <- paste0(meancrz, collapse = " ")
  tit <- sprintf("RCT CONTROL MEANS BY SITE: %s",s_str)
  cat(blank,tit, sep="\n")
  
  s_str <- paste0(rct_impz, collapse = " ")
  tit <- sprintf("RCT IMPACTS BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(rct_sez, collapse = " ")
  tit <- sprintf("RCT STD ERRS BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(rct_nte, collapse = " ")
  tit <- sprintf("RCT TREATMENT SAMPLES BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(rct_nce, collapse = " ")
  tit <- sprintf("RCT CONTROL SAMPLES BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  cat("\n")
  cat("TABLE 7. QED FINDINGS IN QED-ONLY SITES")
  s_str <- paste0(meanqcz, collapse = " ")
  tit <- sprintf("QED COMPARISON MEANS BY SITE: %s",s_str)
  cat(blank,tit, sep="\n")
  
  s_str <- paste0(qed_impz, collapse = " ")
  tit <- sprintf("QED IMPACTS BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(qed_sez, collapse = " ")
  tit <- sprintf("QED STD ERRS BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(qed_nqt, collapse = " ")
  tit <- sprintf("QED TREATMENT SAMPLES BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(qed_nqc, collapse = " ")
  tit <- sprintf("QED COMPARISON SAMPLES BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(eff_sampz, collapse = " ")
  tit <- sprintf("QED TOTAL EFFECTIVE SAMPLE BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  cat("\n")
  cat("TABLE 8. QED BIAS FINDINGS IN RCT-QED SITES")
  s_str <- paste0(bias_nodez, collapse = " ")
  tit <- sprintf("QED BIASES FROM THE CART BIAS TREE LEAVES: %s",s_str)
  cat(blank,tit, sep="\n")
  
  s_str <- paste0(se_bias_nodez, collapse = " ")
  tit <- sprintf("STD ERRS OF QED BIASES FROM THE CART TREE: %s",s_str)
  cat(tit, sep="\n")
  
  cat("\n")
  cat("TABLE 9. QED FINDINGS IN RCT-QED SITES")
  s_str <- paste0(meanqczr, collapse = " ")
  tit <- sprintf("QED COMPARISON MEANS BY SITE: %s",s_str)
  cat(blank,tit, sep="\n")
  
  s_str <- paste0(qed_impzr, collapse = " ")
  tit <- sprintf("QED IMPACTS BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(qed_sezr, collapse = " ")
  tit <- sprintf("QED STD ERRS BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(qed_nqtr, collapse = " ")
  tit <- sprintf("QED TREATMENT SAMPLES BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(qed_nqcr, collapse = " ")
  tit <- sprintf("QED COMPARISON SAMPLES BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(eff_sampzr, collapse = " ")
  tit <- sprintf("QED TOTAL EFFECTIVE SAMPLE BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  cat("\n")
  cat("TABLE 10. QED FINDINGS IN ALL QED SITES")
  s_str <- paste0(meanqcza, collapse = " ")
  tit <- sprintf("QED COMPARISON MEANS BY SITE: %s",s_str)
  cat(blank,tit, sep="\n")
  
  s_str <- paste0(qed_impza, collapse = " ")
  tit <- sprintf("QED IMPACTS BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(qed_seza, collapse = " ")
  tit <- sprintf("QED STD ERRS BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(qed_nqta, collapse = " ")
  tit <- sprintf("QED TREATMENT SAMPLES BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(qed_nqca, collapse = " ")
  tit <- sprintf("QED COMPARISON SAMPLES BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  s_str <- paste0(eff_sampza, collapse = " ")
  tit <- sprintf("QED TOTAL EFFECTIVE SAMPLE BY SITE: %s",s_str)
  cat(tit, sep="\n")
  
  sink()
  
} # end of HYBRID_IMPACT function

