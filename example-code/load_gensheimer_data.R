library(feather)
library(rhdf5)
library(dplyr)

RemoveDups <- function(df, cols) {
  inds = sample(1:nrow(df))
  df   = df[inds, ]

  dups = duplicated(df[cols])
  df   = df[!dups, ]
  inds = inds[!dups]

  df[sort(inds, index=T)$ix, ]
}

# prefix="Documents/Stanford/RA/RA-reg-landmark/initial-articles/gensheimer/prognosis-model/jnci_paper_2019/"
get_simulated_data <- function(){
  data_dir = "jnci_paper_2019/"
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cox_probs <- c(0,.16,.5,.84,1) # as suggested in https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-13-33

  visits <- read_feather(paste(data_dir,'visits.feather',sep=''))
  visits$age <- rnorm(nrow(visits),mean=60,sd=10)
  visits$sex <- sample(0:1,size=nrow(visits),replace=TRUE)
  text <- h5read(paste(data_dir,'text.h5',sep=''),'visits_tfidf_lasso',compoundAsDataFrame=FALSE)
  text <- t(text)
  n_text <- ncol(text)
  labsvitals <- h5read(paste(data_dir,'labsvitals.h5',sep=''),'visits_labs',compoundAsDataFrame=FALSE)
  labsvitals <- t(labsvitals)
  n_labsvitals <- ncol(labsvitals)
  diag_proc_medi <- h5read(paste(data_dir,'diag_proc_medi.h5',sep=''),'diag_proc_medi',compoundAsDataFrame=FALSE)
  diag_proc_medi <- t(diag_proc_medi)
  n_diag_proc_medi <- ncol(diag_proc_medi)
  data <- cbind(visits$age, visits$sex, text, labsvitals, diag_proc_medi)
  text <- 0
  labsvitals <- 0
  diag_proc_medi <- 0
  gc()

  #exclude children and visits with no f/u
  idx <- visits$age>=18 & visits$days_to_last_contact_or_death>0
  data <- data[idx,]
  visits <- visits[idx,]

  visits$visit_date <- as.Date(visits$visit_date)
  first_visit <- visits %>% group_by(patient_id) %>% summarise(first_visit=min(visit_date))
  visits <- inner_join(visits,first_visit,by='patient_id')
  visits$time1 <- as.numeric(visits$visit_date-visits$first_visit) #time1 = time since landmark time t0 (first visit after metastatic cancer diagnosis)
  visits_train <- visits[visits$set==0,]
  visits_test <- visits[visits$set==2,]
  dataMean <- apply(data[visits$set==0,], 2, mean)
  dataStd <- apply(data[visits$set==0,], 2, sd)
  data <- (data - do.call('rbind',rep(list(dataMean),dim(data)[1]))) / do.call('rbind',rep(list(dataStd),dim(data)[1]))

  data_df <- cbind(visits$patient_id, as.Date(visits$first_visit), as.Date(visits$visit_date), as.Date(visits$date_last_contact_or_death), visits$dead, data.frame(data))
  x_cols <- c("age", "sex", paste0("text", 1:n_text), paste0("labsvitals", 1:n_labsvitals), paste0("diag_proc_medi", 1:n_diag_proc_medi))
  colnames(data_df) <- c("id", "first_visit", "visit_date", "date_last_contact_or_death", "status", x_cols)
  data_df$time = data_df$date_last_contact_or_death - data_df$first_visit
  data_df$visit = data_df$visit_date - data_df$first_visit
  data_df = RemoveDups(data_df, c("id", "visit"))

  return(list(df=data_df, x_cols=x_cols))
}
