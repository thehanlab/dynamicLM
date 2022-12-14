#########################################################################
#### BUILDING SUPER DATASET ####
#########################################################################

# -----------------------------------------------------------------------
# packyears2: helper function to calculate packyears, only used if fup
#             occurs between IPLC and prediction time point
# -----------------------------------------------------------------------

pack_years <- function(IPLC_py, smkstatus, smkstatus_fup, cigday, time_IPLC_to_fup, pred_time){
  time_fup_to_pred = pred_time - time_IPLC_to_fup

  ## if fup info is missing use BL info
  cigday_fup = cigday
  smkstatus_fup = if_else(is.na(smkstatus_fup), smkstatus, smkstatus_fup)

  ## smkstatus=2,smkstatus_fup=2 => py = IPLC_py
  ## smkstatus=2,smkstatus_fup=3 => py = IPLC_py + (time_fup to pred point)*cigday_fup/20
  ## smkstatus=3,smkstatus_fup=2 => py = IPLC_py + (time IPLC to fup)*cigday/20
  ## smkstatus=3,smkstatus_fup=3 => py = IPLC_py + (time IPLC to fup)*cigday/20 + (time_fup to pred point)*cigday_fup/20

  py = if_else(smkstatus==2,
                if_else(smkstatus_fup==2,
                        IPLC_py,
                        IPLC_py + time_fup_to_pred*cigday_fup/20),
                #smkstatus==3
                if_else(smkstatus_fup==2,
                        IPLC_py + time_IPLC_to_fup*cigday/20,
                        IPLC_py + time_IPLC_to_fup*cigday/20 + time_fup_to_pred*cigday_fup/20)
                )
  return(py)
}


# -----------------------------------------------------------------------
# update: updates covariates in LM data frame to be LM dependent
# -----------------------------------------------------------------------
update_df <- function(LMDAT){
  agevar="age_ix"
  pyvar ="packyears2"
  quitvar = "quityears2"
  statvar="smkstatus2"
  startage = 55
  stopage = 80
  py.thred = 30

  # TODO: lots of NA handling to think about

  if (class(LMDAT)=="LMdataframe") { DAT <- LMDAT$LMdata }
  else { DAT <- LMDAT }

  DAT <- DAT %>%
    mutate(
      age_ix = age_ix + LM,

        # diff_fup_IPLC = time fup->IPLC
        # so -diff_fup_IPLC gives time IPLC(Time)->fup (>0 is relevant)
        # and -diff_fup_IPLC-LM gives time from pred time point -> fup
      fup_af_IPLC_bef_LM = if_else(is.na(diff_fup_IPLC),          # no fup
                                   F, if_else(-diff_fup_IPLC>0 & -diff_fup_IPLC-LM<0,  # fup after IPLC fup before LM
                                              T,F)),

      smkstatus_changed = if_else(fup_af_IPLC_bef_LM,
                                  if_else(is.na(smkstatus_fup),F,
                                          if_else(smkstatus2!=smkstatus_fup, T, F)), F),
      smkstatus2 = if_else(fup_af_IPLC_bef_LM,
                           if_else(is.na(smkstatus_fup),
                                   smkstatus2,
                                   smkstatus_fup),
                           smkstatus2
      ),

      # stick to cigday2 - do not update
      # cigday2 =    if_else(fup_af_IPLC_bef_LM,
      #                      if_else(is.na(cigday_fup),
      #                              cigday2,
      #                              cigday_fup),
      #                      cigday2
      # ),

      ## smkstatus2 has been updated

      quityears2 = if_else(fup_af_IPLC_bef_LM,

                          ## Need to update
                          if_else(is.na(smkstatus2),
                                  NA_real_, ## assume current smoker
                                  if_else(smkstatus2 == 2, ## former smoker
                                          ## now need to check if they become a former smoker or have always been
                                          if_else(smkstatus_changed,
                                                  diff_fup_IPLC+LM, ## quityears are from fup time to pred time
                                                  if_else(is.na(quityears2), LM, quityears2+LM)), ## otherwise prev amount + time to pred
                                          NA_real_)), ## current smoker

                          ## only using info from baseline
                          if_else(is.na(smkstatus2),
                                  NA_real_, # current smoker
                                  if_else(smkstatus2 == 2, # former smoker
                                          if_else(is.na(quityears2),# no BL quit info
                                                  LM,               # assume at least LM
                                                  quityears2+LM),   # BL info+LM
                                          NA_real_)) # current smoker
      ),

      packyears2 = if_else(fup_af_IPLC_bef_LM,
                           pack_years(packyears2, smkstatus2, smkstatus_fup, cigday2, -diff_fup_IPLC, LM),
                           if_else(smkstatus2==2, # TODO: add an if_else for NA
                                   packyears2,
                                   packyears2 + LM * (cigday2/20))
      )

    ) %>% select(-fup_af_IPLC_bef_LM,-smkstatus_changed)

  DAT[["USPSTF2"]] = myEligibility2(agevar, pyvar, quitvar, statvar, DAT, startage, stopage, py.thred)
  DAT[["USPSTF.stage"]] = DAT$USPSTF2 * DAT$stage2.ix

  TESTDAT = FALSE
  if (TESTDAT) {
    if (missing(add_vars)){
      stop("If wanting test data, add_vars that need to be removed must be specified.")
    }
    DAT <- DAT %>% select(-all_of(add_vars),-LM)
  }

  if (class(LMDAT)=="LMdataframe") { LMDAT$LMdata <- DAT }
  else { LMDAT <-  DAT }


  return(LMDAT)
}

#############################################################################
## CREATING NEW DATA
## This is the same at cutLM but without censoring time at prediction window
#############################################################################

cutLMnewdata <- function (data, outcome, LM, horizon, covs, format = c("wide",
                                                       "long"), id, rtime, right = TRUE)
{
  format <- match.arg(format)
  if (format == "wide") {
    LMdat <- data
    if (!is.null(covs$varying))
      LMdat[[covs$varying]] <- 1 - as.numeric(LMdat[[covs$varying]] >
                                                 LM)
  }
  else {
    if (missing(id))
      stop("argument 'id' should be specified for long format data")
    if (missing(rtime))
      stop("argument 'rtime' should be specified for long format data")
    ord <- order(data[[id]], data[[rtime]])
    data <- data[ord, ]
    ids <- unique(data[[id]])
    n <- length(ids)
    LMdat <- data[which(!duplicated(data[[id]])), ]
    for (i in 1:n) {
      wh <- which(data[[id]] == ids[i])
      di <- data[wh, ]
      idx <- cut(LM, c(data[[rtime]][wh], Inf), right = right,
                 labels = FALSE)
      if (!is.na(idx))
        LMdat[i, ] <- di[idx, ]
      else {
        LMdat[i, ] <- di[1, ]
        LMdat[[covs$varying]][i] <- NA
        LMdat[[rtime]][i] <- NA
      }
    }
  }
  LMdat <- LMdat[LMdat[[outcome$time]] > LM, ]
  if (format == "long")
    LMdat <- LMdat[!is.na(LMdat[[id]]), ]
  LMdat[outcome$status] <- LMdat[[outcome$status]] * as.numeric(LMdat[[outcome$time]] <=
                                                                    horizon)
  # LMdat[outcome$time] <- pmin(as.vector(LMdat[[outcome$time]]),
  #                              horizon)
  LMdat$LM <- LM
  if (format == "long")
    cols <- match(c(id, outcome$time, outcome$status, covs$fixed,
                    covs$varying, rtime, "LM"), names(LMdat))
  else cols <- match(c(outcome$time, outcome$status, covs$fixed,
                       covs$varying, "LM"), names(LMdat))
  return(LMdat[, cols])
}



#########################################################################
#### PRINT FUNCTIONS FOR COEFFICIENTS OF MODELS  ####
#########################################################################

CI <- function(mu, se, n_dec){
  paste("(",round(mu-1.96*se,n_dec),", ",round(mu+1.96*se,n_dec),")",sep="")
}


# -----------------------------------------------------------------------
# For CSC model with all three coxph models using the same covars
# -----------------------------------------------------------------------
get_coef_info <- function(fm, n_dec, CI=TRUE, pVal=FALSE){
  bet <- lapply(1:3,function(i) fm$models[[i]]$coefficients)

  # Get coefficients and desired quantities
  if(CI){
    se <- lapply(1:3,function(i) summary(fm$models[[i]])$coefficients[,4])
    coef_info <- data.frame(
      bet1=bet[[1]],bet2=bet[[2]],bet3=bet[[3]],
      CI1 = CI(bet1,se[[1]],n_dec),
      CI2 = CI(bet2,se[[2]],n_dec),
      CI3 = CI(bet3,se[[3]],n_dec)
    ) %>%
      transmute(
        # `Coef 1`=round(bet1, n_dec),
        `HR 1`=exp(bet[[1]]),
        `95% CI 1`=CI1,
        # `Coef 2`=round(bet2, n_dec),
        `HR 2`=exp(bet[[2]]),
        `95% CI 2`=CI2,
        # `Coef 3`=round(bet3, n_dec),
        `HR 3`=exp(bet[[3]]),
        `95% CI 3`=CI3
      )
    coef_info <- coef_info[ order(row.names(coef_info)), ]

  } else if (pVal) {
    pval <- lapply(1:3,function(i) summary(fm$models[[i]])$coefficients[,6])
    coef_info <- data.frame(
      # `Coef 1`=bet[[1]],
      `HR 1`=exp(bet[[1]]), `p-value 1`=pval[[1]],
      # `Coef 2`=bet[[2]],
      `HR 2`=exp(bet[[2]]), `p-value 2`=pval[[2]],
      # `Coef 3`=bet[[3]],
      `HR 3`=exp(bet[[3]]), `p-value 3`=pval[[3]]
    )
    coef_info <- round(coef_info[ order(row.names(coef_info)), ],n_dec)
  }

  # # Get AIC and BIC
  # aic <- lapply(1:3,function(i) round(AIC(fm$models[[i]]),n_dec))
  # bic <- lapply(1:3,function(i) round(BIC(fm$models[[i]]),n_dec))
  # model_info <- as.data.frame(matrix(
  #   c("","",aic[[1]],bic[[1]],"","",aic[[2]],bic[[2]],"","",aic[[3]],bic[[3]]),
  #   ncol = 6,
  #   dimnames = list(c("AIC:","BIC:"),
  #                   colnames(coef_info))))
  # coef_info <- rbind(coef_info, model_info)

  # Combine and format
  coef_info <- coef_info %>%
    addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7")) %>%
    htmlTable(align = "rrrr",
              cgroup=c("Cause 1: SPLC","Cause 2: IPLC death","Cause 3: Other-cause death"),
              n.cgroup=c(2,2,2),
              # tspanner = c("","Model information criteria:"),
              n.tspanner = c(nrow(coef_info)-2,2)) #,
  return(coef_info)
}


# -----------------------------------------------------------------------
# For a single coxph model (one of the cause-specific functions)
# -----------------------------------------------------------------------
get_coef_info2 <- function(fm, cause, n_dec, CI=TRUE, pVal=FALSE){
  coef_info <- data.frame(
    Coef = fm$coefficients,
    pval = summary(fm)$coefficients[,6]
  )
  coef_info <- round(coef_info[ order(row.names(coef_info)), ],n_dec)

  # Get AIC and BIC
  model_info <- as.data.frame(matrix(
    c("","",round(AIC(fm),n_dec),round(BIC(fm),n_dec)),
    ncol = 2,
    dimnames = list(c("AIC:","BIC:"),
                    colnames(coef_info))))

  # Combine and format
  coef_info <- rbind(coef_info, model_info)
  coef_info <- coef_info %>%
    addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7")) %>%
    htmlTable(align = "rr",
              tspanner = c("","Model information criteria:"),
              n.tspanner = c(nrow(coef_info)-2,2),
              caption=paste("Cause:",cause))
  return(coef_info)
}

######################################
#### USPSTF ELIGIBILITY FUNCTION  ####
######################################
### 7/14/2020 update the function so it can handle missing data ########
### By Summer and Eunji

myEligibility2 = function(agevar, pyvar, quitvar, statvar, mydata, startage, stopage, py.thred){

  ## statvar (smoking staus should be coded as 1, 2, 3 for never, ever, and current smokers


  vars = c(agevar, pyvar, quitvar, statvar)
  dat = mydata[,vars]
  #> dat[1:20,]
  #     age_ix packyears2 quityears2 smkstatus2
  #1  75.13425       31.8    0.00000          3
  #2  78.97260       32.5   28.17534          2
  #3  77.80274        0.0         NA          1
  #4  89.39726         NA   21.00000          2
  #5  76.80000       31.8    0.50000          2
  #6  75.47123       19.8   13.00000          2
  #7  80.21644       45.3   18.00000          2
  #8  68.29863       19.8    0.00000          3
  #9  87.55616       45.3    0.00000          3
  #10 91.97808         NA         NA       <NA>
  #11 52.45205       19.8    0.00000          3
  #12 63.12877       55.0    0.00000          3
  #13 74.30411       16.4    0.00000          3
  #14 71.29589       52.3    0.00000          3
  #15 69.38082       12.0   29.50411          2
  #16 87.72603       45.3   13.00000          2
  #17 88.55890        0.0         NA          1
  #18 87.72329       27.5    0.00000          3
  #19 78.55890       27.5    0.00000          3
  #20 86.39178       60.3   22.42466          2

  ### OK, define ineligibility for each criterion ###

  ## ineligible due to age
  ix.age = (dat[,agevar] < startage | dat[,agevar] > stopage) #& is.na(dat[,agevar])==F

  sum(ix.age)


  ### ineligible due to packyears ##


  ix.py = 	(dat[,pyvar] < py.thred) #& is.na(dat[,pyvar])==F


  #### ineligible due to quityears ###

  ## make sure this should be former smokers only criterion ##

  ix.quityears =  (dat[,statvar]==2 ) & (dat[,quitvar] > 15 )

  #ix.quityears =  (dat[,statvar]==2 & is.na(dat[,statvar])==F) & (dat[,quitvar] > 15 & is.na(dat[,quitvar]) == F)

  tmp=cbind(ix.age=ix.age, ix.py=ix.py, ix.quityears=ix.quityears)
  #> (1*tmp)[1:20,]
  #      ix.age ix.py ix.quityears
  # [1,]      0     0            0
  # [2,]      0     0            1
  # [3,]      0     1            0
  # [4,]      1    NA            1
  # [5,]      0     0            0
  # [6,]      0     1            0
  # [7,]      1     0            1
  # [8,]      0     1            0
  # [9,]      1     0            0
  #[10,]      1    NA           NA
  #[11,]      1     1            0
  #[12,]      0     0            0
  #[13,]      0     1            0
  #[14,]      0     0            0
  #[15,]      0     1            1
  #[16,]      1     0            0
  #[17,]      1     1            0
  #[18,]      1     1            0
  #[19,]      0     1            0
  #[20,]      1     0            1



  ix.eligible = rep(NA, nrow(dat))


  ### (1) those who don't have NA ######

  mysum = rowSums(tmp)
  #> mysum[1:10]
  # [1]  0  1  1 NA  0  1  2  1  1 NA

  ### if no criaterion is violated, then eligible


  ix.eligible[is.na(mysum)==F & mysum ==0]=1



  ### if at least one violated, not eligible
  ix.eligible[is.na(mysum)==F & mysum >0 ]=0



  ## (2) those who have at least one 1 should be ineligible


  ix.eligible[is.na(mysum)==T & rowSums(tmp,na.rm=T)>0]=0

  #rowSums(tmp,na.rm=T)[1:20]
  #[1] 0 1 1 2 0 1 2 1 1 1 2 0 1 0 2 1 2 2 1 2
  #> length(rowSums(tmp,na.rm=T))
  #[1] 7293
  #> length(rowSums(tmp))
  #[1] 7293
  ## even though the rest of them dont violate, we cannot say this person is eligible


  ## only current smokers or never smokers with missing quityears

  table(ix.eligible,useNA="always")
  #ix.eligible
  #   0    1 <NA>
  #5100 2058  135



  ix.eligible*1


}#	myEligibility2
