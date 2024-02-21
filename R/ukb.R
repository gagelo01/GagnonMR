#' Leverage hospital admission record for a subset of three level ICD10
#'
#' @param code_inclusion The code that you which to include if you include I10, it will include I101 and so on.
#' @param code_exclusion The code that you wish to exclude. If participants have this code before inclusion code they will be tagges as toexclude
#' @param fieldid_date the fieldid for the date in ukb either 41280, or 41282
#' @param fieldid_code the fieldid for the code in ukb either 41270, or 41272
#' @param nrows the number of row. When testing I suggest only including 1e4 row as it speeds the process.
#' @param ukbpath the path to the ukb phenotype data where c("f.41280.0", "f.41270.0") are present
#'
#' @return
#' @export

ukb_format_HAR_OPCS <- function(code_inclusion, code_exclusion,
                                fieldid_date = "f.41280.0", #41282
                                fieldid_code = "f.41270.0", #41272
                                nrows = Inf,
                                ukbpath="/home/couchr02/Mendel_UKB/Source/Phenotype/February_2023_Update/ukb671338.tab") {

  columnname <- fread(ukbpath, nrows = 2)
  columnname <- names(columnname)
  tolocate<-paste(unique(c(code_inclusion, code_exclusion)), collapse = "|")
  fieldid <- c(date= fieldid_date, code =fieldid_code) #date, code
  data_disease <- fread(input = ukbpath,
                        select = c(1, which(grepl(paste(fieldid, collapse = "|"), columnname))),
                        nrows = nrows)


  data_disease <- data_disease[2:.N,]

  measure.vars<-colnames(data_disease)[grepl(fieldid_date, colnames(data_disease))]
  data_disease[, (measure.vars):=lapply(.SD,as.IDate), .SDcols = measure.vars]

  index<-data_disease[,apply(.SD, 1, function(x) any(grepl(tolocate, x))), .SDcols = colnames(data_disease)[grepl(fieldid_code, colnames(data_disease))]]
  k <- data_disease[index, ]
  if(nrow(k)==0){return(data.table(f.eid = integer(), exclude_date = Date(), outcome_date = Date(), toexclude = logical()))}
  list_res<-vector(mode = "list", length = length(fieldid))
  for(i in 1:length(fieldid)) {
    measure.vars<-colnames(data_disease)[grepl(fieldid[i], colnames(data_disease))]
    long <- melt(k, id.vars = "f.eid", measure.vars = measure.vars)
    long[,c("todump1", "data_field", "todump2", "index") := data.table::tstrsplit(variable, split = ".", fixed = TRUE)]
    long[,c("todump1","todump2"):=NULL]
    long[,data_field := names(fieldid)[i]]
    wide<- dcast(long, f.eid + index ~ data_field, value.var = "value")
    list_res[[i]]<- wide
  }

  wide <- merge(list_res[[1]], list_res[[2]], by = c("f.eid", "index"))
  rm(list="list_res")
  wide <- wide[grepl(tolocate, code), ]
  wide <- dcast(wide, f.eid ~ code, value.var = "date")
  my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
  my.lowerthan <- function(k1, k2) {
    totransformNA<- is.na(k1)&is.na(k2)
    k1[is.na(k1)]<- Inf
    k2[is.na(k2)]<- Inf
    res <- k1 < k2
    res[totransformNA]<-NA
    return(res)
  }

  if(!is.null(code_exclusion)) {
    col_exclude <- colnames(wide)[grepl(paste(code_exclusion, collapse = "|"), colnames(wide))]
    wide[, exclude_date := apply(.SD, 1, function(x) my.min(x)),.SDcols = col_exclude]} else
    {col_exclude<-NULL; wide[, exclude_date := NA]}
  therest <-  setdiff(colnames(wide), c(col_exclude, "f.eid", "exclude_date"))
  wide[, outcome_date := apply(.SD, 1, function(x) my.min(x)),.SDcols = therest]
  wide[, toexclude := my.lowerthan(exclude_date, outcome_date) ]
  return(wide)
}



#' Title
#'
#' @param fieldid_inclusion_code the ukb fieldid of the inclusion code
#' @param fieldid_exclusion_code the ukb fieldid of the exclusion code can be either in the form  c("f.131270.0.0", "f.131272.0.0") or  "f.131270.0.0|f.131272.0.0"
#' @param outcome_name the name of the outcome for example "aortic_stenosis"
#' @param start_followup either "date_birth" or "date_assessment"
#' @param nrow the number of row. When testing I suggest only including 1e4 row as it speeds the process.
#' @param ukbpath the path to the ukb phenotype data where c("f.41280.0", "f.41270.0") are present
#'
#' @return
#' @export

ukb_disease_to_coxph <- function(fieldid_inclusion_code, fieldid_exclusion_code, outcome_name,
start_followup = c("date_birth", "date_assessment"),
nrows=Inf,
ukbpath="/home/couchr02/Mendel_UKB/Source/Phenotype/February_2023_Update/ukb671338.tab") {

columnname <- fread(ukbpath, nrows = 2)
columnname <- names(columnname)
fieldid <- c("f.34.0.0", "f.52.0.0", "f.53.0.0", "f.40000.0.0",#year of birth, date of birth, and date of attending assessment center, #date of death
             fieldid_inclusion_code, fieldid_exclusion_code)

data_disease <- fread(input = ukbpath,
                      select = c(1, which(grepl(paste(fieldid, collapse = "|"), columnname))),
                      nrows = nrows)
data_disease <- data_disease[2:.N,]

data_disease[, date_birth := ymd(paste0(f.34.0.0, "-", f.52.0.0, "-01")) %>% as.IDate(.)]
data_disease[, age_enrollment := f.53.0.0 - date_birth]
data_disease[, outcome_date := apply(.SD, 1, function(x) min(x, na.rm = TRUE)), .SDcols = c(fieldid_inclusion_code)]
data_disease[, outcome_date := outcome_date %>% as.IDate(.)]

date_end_followup <- data_disease[, max(outcome_date, na.rm = TRUE)]
date_end_followup <- round_date(date_end_followup, unit = "month")+30 #Here is the date the follow up ends.

data_disease[,censoring_date := dplyr::if_else(is.na(f.40000.0.0), ymd(date_end_followup) %>% as.IDate(.), f.40000.0.0) %>%
               dplyr::if_else(!is.na(outcome_date), outcome_date, .)]
data_disease[, outcome_censored := ifelse(is.na(outcome_date), 0, 1)]

if(start_followup=="date_birth"){data_disease[, outcome_time := as.IDate(censoring_date) - as.IDate(date_birth)]}
if(start_followup=="date_assessement"){data_disease[, outcome_time := as.IDate(censoring_date) - as.IDate(f.53.0.0)]}

data_disease[,toinclude:=1]
if(!is.null(fieldid_exclusion_code)) {
data_disease[,exclusion_disease_before := apply(.SD, 1, function(x) any(outcome_date>x)),.SDcols = c(fieldid_exclusion_code)]
data_disease[exclusion_disease_before==TRUE,toinclude := 0] #y revenir
} else {data_disease[,exclusion_disease_before := NA]}
data_disease[outcome_time<1,toinclude := 0] #y revenir
data_disease <- data_disease[,.(f.eid, outcome_censored, outcome_time, toinclude, exclusion_disease_before)]
colnom<-c("outcome_censored", "outcome_time", "toinclude", "exclusion_disease_before")
setnames(data_disease, colnom, paste0(outcome_name, "_", colnom))
return(data_disease) }


