

# Data cleaning

#This file is for downloading and cleaning the data required for estimating
# multimorbidity incidence and prevalence by IMD quintile.

#Our patient sample was derived from May 2020 Aurum database

#Inputs:
#1. Aurum lookup files (July 2020)
#2. 1,00,140 patients from May 2020 Aurum database; using: patient, practice,
# observation, & linked IMD files
#3. Codelists derived from CALIBER lists (available here:
# https://github.com/annalhead/CPRD_multimorbidity_codelists)

#######################################################################

#### SET UP ####

#######################################################################

library(data.table)
#library(foreign)
#library(stringr)
#library(ggplot2)
library(fst) #for quick writing of files
data.table::setDTthreads(10) #this is so that don't use all the processors

# Shortcuts for the directories
data_dir_CPRD <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Data May 2020/", x)


data_dir_lookup <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Dictionaries etc/", x)

data_dir_DL <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Disease_lists/", x)



#Chris' package for data management
library(Rcpp)
sourceCpp(data_dir_lookup("shift_bypid.cpp"),
          cacheDir = "/mnt/alhead/.cache")
shift_bypid <- #for quick shifting of cols by id
  function(x, lag, id, replace = NA) {
    if (typeof(x) == "integer") {
      return(shift_bypidInt(x, lag, replace, id))
    } else if (typeof(x) == "logical") {
      return(shift_bypidBool(x, lag, replace, id))
    } else if (typeof(x) == "double") {
      return(shift_bypidNum(x, lag, replace, id))
    } else
      stop("class of x not supported")
  }



#######################################################################

#### LOAD LOOKUP FILES - AURUM DOCUMENTATION ####

#######################################################################


# Aurum dictionary
aurumdict <-
  fread(
    file = data_dir_lookup("202007_EMISMedicalDictionary.txt"),
    header = TRUE,
    sep = "\t" ,
    colClasses = 'character'
  ) #reads everything in as characters
setnames(aurumdict, tolower(names(aurumdict))) #making all col names lower case

# Patient type
patienttype <- fread(
  file = data_dir_lookup("PatientType.txt"),
  header = TRUE,
  sep = "\t" ,
  colClasses = 'character'
)

# Unit type for numerical values
numunit <- fread(
  file = data_dir_lookup("NumUnit.txt"),
  header = TRUE,
  sep = "\t" ,
  colClasses = 'character'
)

# Observation type
obstype <- fread(
  file = data_dir_lookup("ObsType.txt"),
  header = TRUE,
  sep = "\t",
  colClasses = 'character'
)

#######################################################################

#### LOAD MEDCODEIDS FOR CONDITIONS OF INTEREST ####

#######################################################################


# Loading the list of conditions of interest (incl disease system and gives
# number to each disease) & the codelists
#read in the summary list of diseases of interest
diseasesum <- fread(data_dir_DL("DiseaseSumm.csv"))

fullcodelistmedid <-
  read_fst(data_dir_DL("fullcodelistmedidJAN2021.fst"),
           as.data.table = TRUE)
fullcodelistmedid[, disease_num := as.numeric(disease_num)]

codelistsinglemedid <-
  read_fst(data_dir_DL("codelistsinglemedidJAN2021.fst"),
           as.data.table = TRUE)
codelistsinglemedid[, disease_num := as.numeric(disease_num)]

codelistmultimedid <-
  read_fst(data_dir_DL("codelistmultimedidJAN2021.fst"),
           as.data.table = TRUE)
codelistmultimedid[, disease_num := as.numeric(disease_num)]

# Medcodeids for test value results
testcodes <- fread(
  file = data_dir_DL("testcodes.csv"),
  header = TRUE,
  sep = ",",
  colClasses = 'character'
)

#######################################################################

####  LOAD PRACTICE FILES ####

#pracid: Practice identifier

#region: Region - Value to indicate where in the UK the practice is based.
#The region denotes the Strategic Health Authority for practices within England,
#and the country i.e. Wales, Scotland, or Northern Ireland for the rest

#lcd: Last Collection Date - Date of the last collection for the practice

#uts: Up To Standard Date - Date at which the practice data is deemed to be of
#research quality.
#Derived using a CPRD algorithm that primarily looks at practice death
#recording and gaps in the data

#######################################################################

if (file.exists(data_dir_CPRD("practice.fst"))) {
  practice <-
    read_fst(data_dir_CPRD("practice.fst"), as.data.table = TRUE)
} else {
  require(tidyverse)
  require(data.table)
  
  #Puts together a list of all the files matching the pattern
  pracfiles <-
    list.files(
      data_dir_CPRD("Patient_Practice_Staff_Referral"),
      pattern = ".*Practice_001.txt[.]zip",
      recur = T,
      full = T
    )
  
  
  practice <- data.table()
  for (fn in pracfiles[1:2]) {
    subtab <-
      readr::read_tsv(fn, col_type = cols(.default = col_character())) %>%
      as.data.table() %>%
      mutate(pracid = as.numeric(pracid)) %>%
      mutate(lcd = as.IDate(lcd, format = "%d/%m/%Y"))
    practice <- rbind(practice, subtab)
    rm(subtab)
  }
  rm(fn)
  
  #Checking for missing values
  practice[is.na(region),] # prac 20957 has no region.
  
  #dropping uts col as empty - not yet available in Aurum
  practice[, uts := NULL]
  
  #only need 1 record per practice (there are 2 because 2 extractions)
  setkey(practice, pracid, lcd)
  practice <-
    practice[practice[, .I[1], by = c("pracid", "lcd")]$V1]
  
  write_fst(practice, data_dir_CPRD("practice.fst"), 100)
}

#######################################################################

####  LOAD PATIENT & IMD FILES, COMBINE WITH PRACID INFO  ####
# & ADD IN BLACK ETHNICITY (FOR CKD) FROM OBSERVATION FILES

#patid: Patient Identifier. Last 3 digits are practice identifiers
#pracid: practice id

#gender: 1 = male (recode to "M"), 2 = female ("F" ), 3 = indeterminate ("I")
#NB: in practice this is sex

#yob: Birth Year

#regstartdate: The date that the patient registered with
#the CPRD contributing practice.
#Most recent date the patient is recorded as having registered at the practice.
#If a patient deregistered for a period of time and returned,
#the return date would be recorded.

#patienttypeid: The category that the patient has been assigned
#to e.g. private, regular, temporary. Lookup: PatientType.txt.

#emis_ddate: Date of death as recorded in the EMIS software. Researchers are
#advised to treat the emis_ddate with caution and consider using the
#cprd_ddate variable below.

#regenddate: Date the patient's registration at the practice ended.
#This may represent a transfer-out date or death date.
#Transferred out period is the time between a patient transferring out and
#re-registering at the same practice. 0 is continuous registration.

#cprd_ddate: Estimated date of death of patient- derived using a CPRD algorithm


#######################################################################


if (file.exists(data_dir_CPRD("patient.fst"))) {
  #patient <- read_fst(data_dir_CPRD("patient.fst"), as.data.table = TRUE)
  patient <- read_fst(
    data_dir_CPRD("patient.fst"),
    columns = c(
      "patid",
      "gender",
      "imd",
      "black",
      "yob" ,
      "dob",
      "reg1yr",
      "censordate",
      "censorreason",
      "region"
    ),
    as.data.table = TRUE
  )
} else {
  require(tidyverse)
  require(data.table)
  
  #Puts together a list of all the files matching the pattern
  patfiles <-
    list.files(
      data_dir_CPRD("Patient_Practice_Staff_Referral"),
      pattern = ".*Patient_001.txt[.]zip",
      recur = T,
      full = T
    )
  
  patient <- data.table()
  for (fn in patfiles[1:2]) {
    subtab <-
      readr::read_tsv(fn, col_type = cols(.default = col_character())) %>%
      as.data.table()
    subtab[, `:=` (
      pracid = as.numeric(pracid),
      yob = as.numeric(yob),
      #correcting date formats:
      emis_ddate = as.IDate(emis_ddate, format = "%d/%m/%Y"),
      regstartdate = as.IDate(regstartdate, format = "%d/%m/%Y"),
      regenddate = as.IDate(regenddate, format = "%d/%m/%Y"),
      cprd_ddate = as.IDate(cprd_ddate, format = "%d/%m/%Y"),
      gender = factor(
        gender,
        levels = c("1", "2", "3"),
        labels = c("M", "F", "I")
      )
    )][, c("mob", "usualgpstaffid") := NULL]
    patient <- rbind(patient, subtab)
    rm(subtab)
  }
  rm(fn, patfiles)
  
  
  #patient acceptability flag & patient type
  table(patient$acceptable) # everybody 'acceptable'
  table(patient$patienttypeid) # everybody type 3: "regular"
  #view(patienttype)
  rm(patienttype)
  
  #including only patients with
  #acceptable flags and regular patienttype
  patient <- patient[acceptable == 1 & patienttypeid == 3,]
  
  acceptpatno <- patient[, length(patid)] #total no. participants
  
  
  #Looking at gender
  patient[, .(.N, "%" = .N / acceptpatno), by = gender] #tablulating gender
  
  #Removing cols no longer need
  #Creating 3 new variables: censordate, dob, reg1yr
  patient[, `:=` (
    patienttypeid = NULL,
    acceptable = NULL,
    #creating new end date of reg: this is regenddate for those where recorded:
    censordate = regenddate,
    #making everybody born on the 1st July:
    dob = as.IDate(paste0(yob, "0701"), format = "%Y%m%d"),
    #creating variable for date that have 1 year of registration:
    reg1yr = regstartdate + 365
  )]
  
  #ignoring emis deathdate as per CPRD recommendations
  patient[!is.na(cprd_ddate) &
            (cprd_ddate != regenddate | is.na(regenddate)),
          .N] # 66,114: cprd_ddate recorded but different to regenddate
  patient[cprd_ddate > regenddate,
          .N] #431 cprd_ddate after regenddate
  patient[cprd_ddate - regenddate > 30,
          .N] #257 cprd_ddate more than 1 month after regenddate
  patient[regenddate > cprd_ddate,
          .N] #65,683 cprd_ddate before regenddate
  patient[regenddate - cprd_ddate  > 30,
          .N] #11,216 cprd_ddate more than 1 month before regenddate
  patient[regenddate - cprd_ddate  > 365,
          .N] #654 cprd_ddate more than 1 month before regenddate
  
  #if cprd_ddate recorded and different to regenddate, using the earliest.
  #Although there are some observations recorded within the time between
  #cprd_ddate and regendate, these will be ignored in the analysis
  patient[!is.na(cprd_ddate) &
            (cprd_ddate != regenddate | is.na(regenddate)),
          censordate := ifelse(cprd_ddate > regenddate, regenddate, cprd_ddate)]
  
  
  
  #creating censor reason 1 = death, 0 = other
  patient[, `:=` (censorreason = factor(ifelse(is.na(cprd_ddate), 0, 1)))]
  
  
  
  #Adding in info from practice file
  patient[practice, on = 'pracid', `:=` (region = i.region, lcd = i.lcd)]
  
  #for people with no regenddate, changing censordate
  #to the practice last collection date
  patient[, `:=` (#emis_ddate = NULL, #keep this for now....
    censordate = as.IDate(ifelse(is.na(censordate), lcd, censordate),
                          format = "%d/%m/%Y"))]
  
  patient[, lcd := NULL] #dropping lcd as no longer needed
  
  #adding in IMD data - 1 least deprived; 5 most deprived
  IMD <- fread(
    file = data_dir_CPRD(
      "IMD Oct2020/Aurum_linked/patient_imd2015_19_173_request2.txt"
    ),
    header = TRUE,
    sep = "\t",
    colClasses = 'character'
  )
  patient[IMD, on = 'patid', imd := i.imd2015_5]
  table(patient$imd) #1,345 have no imd recorded
  patient[, .(.N, "%" = .N / acceptpatno), keyby = imd] #tablulating imd
  
  rm(acceptpatno, IMD)
  write_fst(patient, data_dir_CPRD("patient.fst"), 100)
  
  
  
  
  #Extracting codes for black ethnicity and adding to patient,
  #this is for calculating egfr
  black <- fread(
    file = data_dir_lookup("ckdepi_black.txt"),
    header = TRUE,
    sep = "\t" ,
    colClasses = 'character'
  ) #reads everything in as characters
  black_mc <- black$medcodeid
  
  require(tidyverse)
  require(data.table)
  
  #Puts together a list of all the files matching the pattern
  obsfiles <- list.files(
    data_dir_CPRD("Observation") ,
    pattern = ".*txt[.]zip",
    recur = T,
    full = T
  )
  
  obstabblack <- data.table()
  for (fn in obsfiles[1:111]) {
    subtab <-
      readr::read_tsv(fn, col_type = cols(.default = col_character())) %>%
      select(patid, medcodeid) %>%
      filter(medcodeid %in% black_mc) %>% as.data.table()
    obstabblack <- rbind(obstabblack, subtab)
    rm(subtab)
  }
  rm(fn)
  
  patient[, black := ifelse(patid %in% obstabblack$patid, 1, 0)][
    , black := factor(black)]
  rm(obstabblack, black_mc)
  write_fst(patient, data_dir_CPRD("patient.fst"), 100)
}

#Total number of patient years
patient[censordate - reg1yr >= 0, sum(censordate - reg1yr) / 365.24]





#######################################################################

####  LOAD OBSERVATION FILE DATA FOR MEDCODEIDS OF INTEREST  ####

#patid: Patient identifier

#consid: Consultation identifier

#pracid: CPRD Practice identifier

#obsid: Observation identifier

#obsdate: Event date

#enterdate: Entered date

#staffid: Staff identifier

#parentobsid: Parent observation identifier. Observation identifier (obsid)
#that is the parent to the observation. This enables grouping of multiple
#observations, such as systolic and diastolic blood pressure, or
#blood test results.

#medcodeid: Medical code. CPRD unique code for the medical term selected by
#the GP. Lookup: Medical dictionary

#value: Value. Measurement or test value

#numunitid:Numeric unit identifier. Unit of measurement. Lookup: NumUnit.txt

#obstypeid: Observation type identifier. Type of observation (allergy,
#family history, observation). Lookup: ObsType.txt

#numrangelow: Numeric range low. Value representing the low boundary of the
#normal range for this measurement

#numrangehigh: Numeric range high. Value representing the high boundary of the
#normal range for this measurement

#probobsid: Problem observation identifier. Observation identifier (obsid) of
#any problem that an observation is associated with. An example of this might
#be an overarching condition that the current observation is associated with
#such as 'wheezing' with the problem observation identifier that links to an
#observation of 'asthma', that in turn contains information in the
#problem table. Link Observation table

#######################################################################

if (file.exists(data_dir_CPRD("obstaball_raw.fst"))) {
  obstaball <-
    read_fst(data_dir_CPRD("obstaball_raw.fst"), as.data.table = TRUE)
} else {
  require(tidyverse)
  require(data.table)
  
  #Puts together a list of all the files matching the pattern
  obsfiles <- list.files(
    data_dir_CPRD("Observation"),
    pattern = ".*txt[.]zip",
    recur = T,
    full = T
  )
  
  # Medcodeids of interest
  obscodes <-
    unique(c(fullcodelistmedid$medcodeid, testcodes$medcodeid))
  
  obstaball <- data.table()
  for (fn in obsfiles[1:111]) {
    subtab <-
      readr::read_tsv(fn, col_type = cols(.default = col_character())) %>%
      select(
        patid,
        obsid,
        obsdate,
        enterdate,
        medcodeid,
        value,
        numunitid,
        obstypeid,
        numrangelow,
        numrangehigh
      ) %>%
      filter(medcodeid %in% obscodes) %>%
      as.data.table() %>%
      mutate(enterdate = as.IDate(enterdate, format = "%d/%m/%Y")) %>%
      mutate(obsdate = as.IDate(obsdate, format = "%d/%m/%Y"))
    obstaball <- rbind(obstaball, subtab)
    rm(subtab)
  }
  rm(fn, obscodes, obsfiles)
  
  obstaball[, `:=` (
    value = as.numeric(value),
    numrangelow = as.numeric(numrangelow),
    numrangehigh = as.numeric(numrangehigh)
  )]
  
  #looking at enterdate and obsdate
  obstaball[is.na(obsdate), .N] #44,999 (0.01% of all obs)
  obstaball[is.na(enterdate), .N] #0
  obstaball[obsdate > enterdate, .N] #48,8431
  head(obstaball[obsdate > enterdate,
                 .(obsdate, enterdate)],
       100) #lots of these have enterdates before 1900
  obstaball[obsdate > enterdate &
              year(enterdate) > 1900,
            .N] #7,771 obsdates are after sensible enterdates
  obstaball[year(enterdate) > 1900 &
              obsdate - enterdate > 30,
            .N] #5,813 are more than a month apart 
  obstaball[year(enterdate) > 1900 & obsdate - enterdate > 365,
            .N] #3,184 are more than a year apart
  
  obstaball[enterdate > obsdate,
            .N] #12,486,515 (27.2%) were entered after obsdate.
  obstaball[enterdate - obsdate > 30,
            .N] #4,286,698 (9.3%) less problematic
  obstaball[enterdate - obsdate > 365,
            .N] #3,166,846 (6.9%)
  
  #Creating an eventdate variable: this is obsdate if recorded,
  #otherwise enterdate unless enterdate is before 1900
  obstaball[, eventdate :=
              as.IDate(ifelse(
                !is.na(obsdate),
                obsdate,
                ifelse(enterdate > "1900-01-01", enterdate, NA)
              ),
              format = "%d/%m/%Y")]
  
  setkey(obstaball, patid, eventdate)
  
  #saving this is before deleting any observations
  write_fst(obstaball, data_dir_CPRD("obstaball_raw_Dec13.fst"), 100)
}


#######################################################################

####  BASIC CLEANING OF OBSERVATION FILE DATA  ####

#######################################################################


# Adding in basic patient info
obstaball[patient, on = 'patid', `:=` (
  yob = i.yob,
  censordate = i.censordate,
  regenddate = i.regenddate,
  cprd_ddate = i.cprd_ddate
)]
obstaball[year(eventdate) < yob, .N] #3,101obs before year of birth
obstaball[eventdate > censordate, .N] #16072 obs after censordate
obstaball[eventdate - censordate > 30, .N] #8078
obstaball[eventdate - censordate > 365, .N] # 2603
obstaball[(cprd_ddate < eventdate & eventdate < regenddate) |
            (regenddate < eventdate & eventdate < cprd_ddate),
          .N] #3,666 obs fall in between regenddate & cprd_ddate

#Getting rid of obs before year of birth and after censor year
obstaball <- obstaball[!is.na(eventdate) &
                         yob <= year(eventdate) &
                         censordate >= eventdate]

#Removing unwanted cols
obstaball[, c("yob", "censordate", "regenddate", "cprd_ddate") := NULL]

# Looking at obstype
#1	Allergy; 2	Annotated Image; 3	Document; 4	Family history; 5	Immunisation
#6	Investigation; 7	Observation; 8	Referral; 9	Test Request; 10	Value;
#11	Care plan
obstaball[, obstypeid := as.numeric(obstypeid)]
obstaball[, .N, keyby = obstypeid]
#obstypeid  N
#1    14365
#2     1873
#3    80614
#4     2791
#5      877
#6  2633914
#7 17163358
#8    84825
#9    28456
#10 25848407

#not including "family history" records (N = 2791)
obstaball <- obstaball[obstypeid != 4,]
setkey(obstaball, patid, eventdate)


#In order to apply alogrithms, separating out into:
#single code ever recorded, multivisit criteria, test results,

#Single code ever recorded
obstabsingle <-
  obstaball[medcodeid %in% codelistsinglemedid$medcodeid,]
obstabsingle <- merge(
  obstabsingle,
  fullcodelistmedid[, .(medcodeid, disease_num)],
  by = "medcodeid",
  all.x = T,
  all.y = F,
  allow.cartesian = TRUE
)
obstabsingle[,.N] #N = 16,626,523


#Test codes
obstabtest <-
  obstaball[medcodeid %in% testcodes$medcodeid,] #N = 27,620,559
obstabtest[, .N]
obstabtest[, value := as.numeric(value)]
obstabtest[aurumdict, on = 'medcodeid', `:=`
           (descr = i.term)] #adding in medcodeid descr etc
obstabtest[numunit, on = 'numunitid',
           unittype := i.Description] #adding in unittype descr


#Getting rid of  records recorded before yob
obstabtest[patient, on = 'patid', yob := i.yob]
obstabtest[yob + 18 > year(eventdate),
           .N] #538,418 were recorded before 18
obstabtest[is.na(value) | value <= 0,
           .N] #864,876 obs have value unrecorded or <=0

obstabtest <- obstabtest[yob + 18 <= year(eventdate) &
                           #excluding tests before 18!is.na(value) &
                           value > 0, ] #excluding obs with a negative value
obstabtest[, yob := NULL]
write_fst(obstabtest, data_dir_CPRD("obstabtest_14Dec.fst"), 100)




#Multi visit criteria required
obstabmulti <-
  obstaball[medcodeid %in% codelistmultimedid$medcodeid,]
obstabmulti <- merge(
  obstabmulti,
  fullcodelistmedid[, .(medcodeid, disease_num)],
  by = "medcodeid",
  all.x = T,
  all.y = F,
  allow.cartesian = TRUE
)
obstabmulti[, .N] #N = 12,902,947
rm(obstaball)



#######################################################################

####  CLEANING OF TEST VALUE DATA - CKD  ####

#######################################################################


## Adding in CKD from test results

# GFR & Creatinine test results for CKD
#A patient is defined as having had CKD stage 3 or above at a specified date:
# IF egfr_ckdepi recorded on or before specified date, THEN
#IF egfr_ckdepi <60 ml/min on the most recent date (index date)
#before the specified date
#AND
#IF egfr_ckdepi <60 ml/min on any date greater than 90 days
#BEFORE the index date above
#THEN classify as having CKD3 or above
#ELSE the patient is not defined as having CKD stage 3 or above.

#Where egfr_ckdepi up to and including 31 Dec 2013 is defined as:
#  egfr_ckdepi = 141 * min(crea_gprd * 0.010746 / K, 1)^alpha
#* max(crea_gprd * 0.010746 / K, 1)^-1.209
#* 0.993^age * 1.018 [if female]  * 1.159 [if black]

#where:
#  alpha = -0.329 for females, -0.411 for males
#K = 0.7 for females, 0.9 for males

#Where egfr_ckdepi from and including 1 Jan 2014 is defined as:
#  egfr_ckdepi = 141 * min(crea_gprd * 0.010746 / K, 1)^alpha
#* max(crea_gprd * 0.0.011312/ K, 1)^-1.209
#* 0.993^age * 1.018 [if female]  * 1.159 [if black]

#where:
#  alpha = -0.329 for females, -0.411 for males
#K = 0.7 for females, 0.9 for males

#Where crea_gprd is defined as:
#  IF enttype = 165 [Serum creatinine]
#AND data1 [Operator] = 3 ["="] AND data2 [Value] > 0
#THEN crea_gprd = data2
#Assumes serum creatinine unit is micromol/L
#(https://www.caliberresearch.org/portal/show/crea_gprd)
#"Measurement of serum creatinine in primary care. Values of zero are discarded.
#The distribution of measurements with units recorded as mmol/L, mol/L or umol/L
#are all consistent with being recorded in micromol/L. Hence we assume the
#units are micromol/L regardless of the actual units stated."

## looking at serum creatinine. #The Kuan algorithm specifies serum creatinine

#Looking at code frequencies
#select & restrict to recorded values
creaobs <-
  obstabtest[medcodeid %in% c("380389013", "259054018", "457927010")
             & !is.na(value) & value > 0,]

creaobs[, .N] #4,403,816

creaobs[, .N, by = c("descr", "medcodeid")][order(-N)]

#DECIDED TO USE THESE THREE CODES:
#1: 380389013                                  Serum creatinine 4,400,928
#2: 259054018                              Serum creatinine NOS    2,108
#3: 457927010                  Corrected serum creatinine level     780




#Looking at unit frequencies
creaobs[, .N, by = unittype][order(-N)]

creaobs[, summary(value)]


#Normal ref range Adult F 45 - 84 M 59 - 104
#https://www.esht.nhs.uk/wp-content/uploads/2017/08/Clinical-Biochemistry-reference-ranges-handbook.pdf

#https://onlinelibrary-wiley-com.liverpool.idm.oclc.org/doi/epdf/10.1002/pds.2203
#this gives biological plausability limits of less than20 micromol/L (0.2 mg/dL)
#and  greater than 1800 micromol/L (20 mg/dL)
#Let's stick to these.
#1800 micromol/L in other valid units would still be below 1800:
creaobs[value < 1800 ,
        summary(value)]
#Regardless of units stated, the IQR is within the normal ref range
ggplot(creaobs[value < 1800], aes(x = value)) +
  geom_histogram(binwidth = 1)
#outliers are 1.5*IQR either side of lower & upper quartiles - 4.23%:
ggplot(creaobs[value < 1800], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )



micromol <- c("umol/l",
              "micmol/l",
              "micromol/l",
              "??mol/l",
              "umol//l",
              "mcmol/l")
creaobs[tolower(descr) %like% "serum" &
          tolower(unittype) %in% micromol, .N] #4,291,157 (97% of serum)

#looking at only those with the specified units
creaobs[value < 1800 &
          tolower(unittype) %in% micromol, summary(value)]
#1800 micromol/L in other valid units would still be below 1800
#Regardless of units stated, the IQR is within the normal ref range
ggplot(creaobs[value < 1800 &
                 tolower(unittype) %in% micromol, ], aes(x = value)) +
  geom_histogram(binwidth = 1)
#outliers are 1.5*IQR either side of lower & upper quartiles - 4.23%
ggplot(creaobs[value < 1800 &
                 tolower(unittype) %in% micromol, ], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
quantile(creaobs[tolower(unittype)  %in% micromol, value],
         c(0.005, .025, .25, .75, .975, .995))

#plotting only those within the 'plausible' limits
ggplot(creaobs[value < 1800 &
                 value > 20 & tolower(unittype) %in% micromol],
       aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
#outliers are 1.5*IQR either side of lower & upper quartiles
ggplot(creaobs[value < 1800 &
                 value > 20 & tolower(unittype) %in% micromol],
       aes(x = value)) +
  geom_histogram(binwidth = 1)
creaobs[value < 1800 &
          value > 20 & tolower(unittype) %in% micromol,
        .N] #4289848/4403816 - 97.4% all records with micromol units
summary(creaobs[value < 1800 &
                  value > 20 & tolower(unittype) %in% micromol,
                value])


#looking at those with different units
#For those with different units, most still in the range above.
#If they were actually recorded in mmol/l mol/l mg/dl they wouldn't be.
#So assuming that they are
creaobs[value < 1800 &
          !tolower(unittype) %in% micromol, summary(value)]
#1800 micromol/L in other valid units would still be below 1800
#Regardless of units stated, the IQR is within the normal ref range
ggplot(creaobs[value < 1800 &
                 !tolower(unittype) %in% micromol, ], aes(x = value)) +
  geom_histogram(binwidth = 1)
ggplot(creaobs[value < 1800 &
                 !tolower(unittype) %in% micromol, ], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
#outliers are 1.5*IQR either side of lower & upper quartiles - 4.23%
#Have quite similar distributions?
quantile(creaobs[!tolower(unittype)  %in% micromol, value],
         c(0.005, .025, .25, .75, .975, .995)) #more in the tails
#plotting only those within the 'plausible' limits for micromol units
ggplot(creaobs[value < 1800 &
                 value > 20 & !tolower(unittype) %in% micromol],
       aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
#outliers are 1.5*IQR either side of lower & upper quartiles
ggplot(creaobs[value < 1800 &
                 value > 20 & !tolower(unittype) %in% micromol],
       aes(x = value)) +
  geom_histogram(binwidth = 1)
creaobs[value < 1800 &
          value > 20 & !tolower(unittype) %in% micromol, .N]
#111273/112659 - 98.77% all records without micromol units
summary(creaobs[value < 1800 &
                  value > 20 & !tolower(unittype) %in% micromol,
                value])


#same range in mmol/l is 0.02 to 1.8
creaobs[tolower(unittype) %in% "mmol/l" &
          value < 1.8 & value > 0.02,
        .N] #there are 65 of these

#same range in mol/l is 0.00002 to 0.0018
creaobs[tolower(unittype) %in% "mmol/l" &
          value < 0.0018 & value > 0.00002,
        .N] #there are 0 of these

#same range in mg/dl is 0.00002 to 0.0018
creaobs[tolower(unittype) %in% "mg/dl"  &
          value < 20.34 & value > 0.226,
        .N] #there are 25 of these


#Just sticking to those within a plausible range, regardless of units
creaobs <- creaobs[value < 1800 & value > 20, ]


rm(micromol)



#adding in patient characteristics needed for formula
creaobs[patient, on = 'patid', `:=` (
  gender = i.gender,
  dob = i.dob,
  black = i.black,
  pracid = i.pracid
)]


#pmin/pmax takes the min/max of values in the vector
creaobs[, egfr := 141 * pmin(value * 0.010746 / 
                               ifelse(gender == "M", 0.9, 0.7), 1) ^
          ifelse(gender == "M",-0.411,-0.329) *
          pmax(
            value * ifelse(eventdate < "2014-01-01", 0.010746, 0.011312) /
              ifelse(gender == "M", 0.9, 0.7),
            1
          ) ^ -1.209 *
          0.993 ^ ((eventdate - dob) / 365.25) *
          ifelse(gender == "M", 1, 1.018) *
          ifelse(black == 1, 1.159, 1)]



# Using recorded gfr and gfr calc from creatinine

#combi  gfr
#select & restrict to eGFR readings (for CKD) <60
#using all the codes from aurum dict that incl "gfr"
eGFRresults <- obstabtest[medcodeid %in% testcodes[
  (tolower(term) %like% "gfr" | tolower(term) %like% "glomerular") &
    test == "ckd", medcodeid] & !is.na(value) & value > 0 & value < 60,
                          .(patid, obsid, value, eventdate, descr)]

#Using these 5:
#1:1942831000006114    eGFR using creatinine (CKD-EPI) per 1.73 square metres
#2:976481000006110                            GFR calculated abbreviated MDRD
#3:1540241000006111 GFRcalculated abbreviated MDRD adj for African Americ orign
#4:1866321000006117 Estimated GFR using CKD-Epi formula per 1.73 square metres
#5:1942821000006111     eGFR using cystatin C (CKD-EPI) per 1.73 square metres

eGFRresults[, .N, by = descr]
#                              GFR calculated abbreviated MDRD    606235
#       eGFR using creatinine (CKD-EPI) per 1.73 square metres    73429
#       eGFR using cystatin C (CKD-EPI) per 1.73 square metres    13
#   Estimated GFR using CKD-Epi formula per 1.73 square metres    334
# GFR calculated abbreviated MDRD adj for African Americ orign    3321

eGFRresults[, descr := NULL]

#Changing the creaobs DT to match the structure for binding
creaobs <-  creaobs[, .(patid, obsid, eventdate, value = egfr)]
eGFRresults <- rbind(eGFRresults, creaobs)
rm(creaobs)
setkey(eGFRresults, patid, eventdate)

eGFRresultsminmax <- eGFRresults[,
                                 .(max = max(value), min = min(value)),
                                 by = patid]
eGFRresultsminmax[, .N] #560,331 people have egfr results
eGFRresultsminmax[min < 60 &
                    max >= 60] #93,069 have at least one value >=60
#& at least one value <60
rm(eGFRresultsminmax)

#only interested in values less than 60
eGFR_CKD <- eGFRresults[value < 60,]

#looking at how many people had X numbers of abnormal records
eGFR_CKD[, numobs := .N, by = patid][, .N, keyby = numobs]
uniqueN(eGFR_CKD[,  patid]) #118,554 have at least one value <60

#need to have at least 2 values less than 60
eGFR_CKD <- eGFR_CKD[numobs > 1,]
uniqueN(eGFR_CKD[, patid]) #99,807 have at least 2 values <60
setkey(eGFR_CKD, patid, eventdate)

#calculating the number of days between the 1st and last value <60 recorded
time <-
  eGFR_CKD[, .(max = max(eventdate), min = min(eventdate)), by = patid]

#getting rid of patients who don't have at least 90 days
#between 1st and last recordings
eGFR_CKD <-
  eGFR_CKD[patid %in% time[as.IDate(max) - as.IDate(min) > 90,
                           patid],]
setkey(eGFR_CKD, patid, eventdate)
rm(time)

firstegrf <- eGFR_CKD[eGFR_CKD[, .I[1], keyby = c("patid")]$V1]
#adding in the date of the 1st <60 value recorded
eGFR_CKD <-
  eGFR_CKD[firstegrf, on = 'patid', `:=` (firstobsdate = i.eventdate)]
rm(firstegrf)

#want records 90+ days after the 1st <60 value recorded
eGFR_CKD <-
  eGFR_CKD[as.IDate(eventdate) - as.IDate(firstobsdate) > 90,]
setkey(eGFR_CKD, patid, eventdate)

#only need the earliest of recordings after the 90 time frame
eGFR_CKD <- eGFR_CKD[eGFR_CKD[, .I[1], keyby = c("patid")]$V1]

#Taking these observations from the original test result data so can combine
eGFR_CKD <-
  obstabtest[obsid %in% eGFR_CKD$obsid][, disease_num := 31]
#testdiagnoses <- rbind(testdiagnoses, eGFRresults)
uniqueN(eGFR_CKD[, patid])  #83,822 individuals
eGFR_CKD[, `:=` (unittype = NULL, descr = NULL)]

#looking at how test value diagnosis compares to medcodeid CKD recording
uniqueN(obstabsingle[disease_num == 31, patid]) #66,0100 have CKD from codes alone
uniqueN(obstabsingle[(disease_num == 31) &
                       patid %in% eGFR_CKD[, patid], patid])
#45,394 of these also have 'eligible' egfr values
eGFR_CKD[!patid %in% obstabsingle[(disease_num == 31), patid], .N]
#38,428 individuals have test values only and no diagnostic codes


#adding this in to obstabsingle
obstabsingle <- rbind(obstabsingle, eGFR_CKD)
rm(eGFRresults, eGFR_CKD)

#keeping only the 1st of CKD & ESRD
temp <- obstabsingle[disease_num == 31, ]
obstabsingle <- obstabsingle[!(disease_num == 31), ]
temp[, disease_num := 31]
setkey(temp, patid, eventdate)
temp <- temp[temp[, .I[1], keyby = c("patid")]$V1]
obstabsingle <- rbind(obstabsingle, temp)
rm(temp)

obstabsingle[disease_num == 31, .N] # 104,528 overall



#######################################################################

####  CLEANING OF CONDITIONS WITH SINGLE CODE EVER RECORDED  ####

#######################################################################



## Stable angina

#A patient is considered to have had stable angina IF
#they meet the criteria for any of the following
#1. Recorded diagnosis of stable angina in primary or secondary care
#2. Coronary revascularisation without unstable angina or
#myocardial infarction in the previous 30 days
#3. Primary care record of abnormal coronary angiogram or
#test showing evidence of myocardial ischaemia

stableangina <- obstabsingle[disease_num == 191]

stableangina[fullcodelistmedid[disease_num == 191, ],
             on = 'medcodeid', `:=` (category = i.category, descr =  i.descr)]

#stable angina/ischaemic chest pain
stableangina_def <- stableangina[category %in% c("Stable angina (4)",
                                                 "Chest pain, attributed to coronary causes (4)")]
stableangina_def <- rbind(stableangina_def,
                          stableangina[category %in% c(
                            "Results abnormal (3)",
                            "Results abnormal- indicative of ischaemia (T-wave changes or ST-depression) (2)"
                          )])

stableangina_proc <-
  stableangina[category %in% c("CABG Performed (2)",
                               "PCI performed (2)"), ]
stableangina_proc <-
  rbind(stableangina_proc[,-c("category", "descr")],
        obstabsingle[disease_num %in% c(100, 206) &
                       patid %in% stableangina_proc$patid])
setkey(stableangina_proc, patid, eventdate)
stableangina_proc[, `:=` (last_cond = (shift_bypid(disease_num, 1L, patid, 0L)),
                          last_date = as.IDate(
                            shift_bypid(eventdate, 1L, patid, 0L)))]
stableangina_proc <- stableangina_proc[!(disease_num == 191 &
                                           last_cond %in% c(100, 206) &
                                           eventdate - last_date < 30)]
stableangina <- rbind(stableangina_def[,-c("descr", "category")],
                      stableangina_proc[,-c("last_cond", "last_date")])
obstabsingle <- obstabsingle[disease_num != 191]
obstabsingle <- rbind(obstabsingle, stableangina)
setkey(obstabsingle, patid, eventdate)
rm(stableangina_proc, stableangina, stableangina_def)

###

## Unstable angina

#removing if have stable angina or MI recorded (this is our addition)
unstableangina <- obstabsingle[disease_num == 206]
unstableangina <- unstableangina[!patid %in%
                                   obstabsingle[disease_num %in% c(191, 100), ]]
obstabsingle <- obstabsingle[disease_num != 206]
obstabsingle <- rbind(obstabsingle, unstableangina)
setkey(obstabsingle, patid, eventdate)
rm(unstableangina)

###

## CHD NOS

#removing if have MI, stable angina or unstable angina

#CALIBER: No previous records meeting the criteria for
# stable angina OR unstable angina OR myocardial infarction
mi <- obstabsingle[disease_num == 39]
mi <-
  mi[!patid %in% obstabsingle[disease_num %in% c(191, 100, 206), patid]]
obstabsingle <- obstabsingle[disease_num != 39]
obstabsingle <- rbind(obstabsingle, mi)
setkey(obstabsingle, patid, eventdate)
rm(mi)

###

## Anaemia - others.

#CALIBER: don't want 'possible diagnosis' if also have any other anaemia,
#sickle cell or thalissemia recorded
anaemothcode <- fullcodelistmedid[disease_num == 8 &
                                    category ==
                                    "Possible Diagnosis of Other anaemias",
                                  medcodeid]
anaem <-
  obstabsingle[disease_num == 8 & medcodeid %in% anaemothcode]
anaem <- anaem[!patid %in% obstabsingle[disease_num %in%
                                          c(13, 61, 111, 184, 196, 211), patid]]
obstabsingle <-
  obstabsingle[!(disease_num == 8 & medcodeid %in% anaemothcode)]
obstabsingle <- rbind(obstabsingle, anaem)
setkey(obstabsingle, patid, eventdate)
rm(anaem, anaemothcode)

###

## Other haemolytic anaemias

#CALIBER: don't want 'possible diagnosis' if also have Thalassaemia or
#Sickle Cell Anaemia.
haemanaemothcode <- fullcodelistmedid[disease_num == 111 &
                                        category == "Possible Diagnosis of Other haemolytic anaemias",
                                      medcodeid]
haemanaem <-
  obstabsingle[disease_num == 111 & medcodeid %in% haemanaemothcode]
haemanaem <-
  haemanaem[!patid %in% obstabsingle[disease_num %in% c(184, 196),
                                     patid]]
obstabsingle <-
  obstabsingle[!(disease_num == 111 & medcodeid %in% haemanaem)]
obstabsingle <- rbind(obstabsingle, haemanaem)
setkey(obstabsingle, patid, eventdate)
rm(haemanaem, haemanaemothcode)

###

#Primary malignancy other

#1. Primary Malignancy Other organs diagnosis or history of diagnosis
#during a consultation
#OR
#2. Primary Malignancy Other organs possible diagnosis during a consultation IF
#NO record satisfying criteria for Primary Malignancy of any other organ in
#this document OR Haematological Malignancy (Hodgkin Lymphoma,
#Non-Hodgkin Lymphoma, Multiple Myeloma (Plasma Cell Malignancy), Leukaemia)
malignumbs <- diseasesum[system_num == 1 & disease_num != 141 &
                           (
                             Disease %like% "Primary" |
                               Disease %like% "Hodgkin" |
                               Disease %in% c("Leukaemia",
                                              "Plasma Cell Malignancy")
                           ) ,
                         disease_num]
maligotherposscodes <-
  fullcodelistmedid[disease_num == 141 &
                      category %in% "Possible Diagnosis of Primary Malignancy_Other Organs", 
                    medcodeid]
maligother <- obstabsingle[disease_num == 141 &
                             medcodeid %in% maligotherposscodes, ]
obstabsingle <- obstabsingle[!(disease_num == 141 &
                                 medcodeid %in% maligotherposscodes), ]
maligother <-
  maligother[!(patid %in% obstabsingle[disease_num %in% malignumbs,
                                       patid])]
obstabsingle <- rbind(obstabsingle, maligother)
rm(maligother, malignumbs, maligotherposscodes)

###

## Secondary malignancy other

#1. Secondary Malignancy Other organs diagnosis or
#history of diagnosis during a consultation
#OR
#2. Secondary Malignancy Other organs possible diagnosis during a consultation
#IF NO record satisfying criteria for Secondary Malignancy of any other organ
secmalignumbs <- diseasesum[system_num == 1 &
                              Disease %like% "Secondary" &
                              disease_num != 177, disease_num]
secmaligotherposscodes <-
  fullcodelistmedid[disease_num == 177 &
                      category %in% "Possible diagnosis of Secondary Malignancy_Other organs",
                    medcodeid]
secmaligother <- obstabsingle[disease_num == 177 &
                                medcodeid %in% secmaligotherposscodes, ]
obstabsingle <- obstabsingle[!(disease_num == 177 &
                                 medcodeid %in% secmaligotherposscodes), ]
secmaligother <-
  secmaligother[!(patid %in% obstabsingle[disease_num %in% 
                                            secmalignumbs, patid])]
obstabsingle <- rbind(obstabsingle, secmaligother)
rm(secmaligother, secmalignumbs, secmaligotherposscodes)


#### TAKING THE SINGLE CODE EVER RECORDED ####
setkey(obstabsingle, patid, eventdate)

#this gives first date of each disease, and gives NA for those with no diseases
firstobs <- obstabsingle[obstabsingle[, .I[1],
                                      by = c("patid", "disease_num")]$V1]



#### Other code amendments ####

## chronic sinusitits
#1. Chronic sinusitis diagnosis or history of diagnosis during a consultation
#OR
#2. There are at least 2 records satisfying the criteria for Possible diagnosis
#of Chronic sinusitis during a consultation more than 84 days apart.

#chronic codes for automatic inclusion
d32chronic <- fullcodelistmedid[disease_num == 32 &
                                  tolower(descr) %like% "chronic", medcodeid]
d32conditional <- fullcodelistmedid[disease_num == 32 &
                                      !medcodeid %in% d32chronic, medcodeid]
tmp <- obstabsingle[medcodeid %in% d32conditional,]
setkey(tmp, patid, eventdate)
tmp2 <- tmp[, .N, by = patid]
tmp2 <- tmp2[N == 1, patid]
tmp <- tmp[!patid %in% tmp2,]
tmp[, prevcode := as.IDate(shift_bypid(eventdate, 1L, patid, 0L))]
chronic <-
  tmp[eventdate - prevcode > 84 & prevcode != "1970-01-01", ]

#this gives the first time the second record is at least 84 days after
#the previous one. BUT what about if have 3 visits and the 1st and 3rd are
#84 days apart? I have ignored these
chronic <- chronic[chronic[, .I[1], by = c("patid")]$V1]
firstobs <-
  firstobs[disease_num != 32,] #removing all chronic sinusitus codes
#from firstobs
chronicrecords <- obstabsingle[obsid %in% chronic$obsid,]
firstobs <- rbind(firstobs, chronicrecords)
setkey(firstobs, patid, eventdate)
rm(tmp, tmp2, chronic, chronicrecords, d32chronic, d32conditional)

###

## Diabetes other

# Diabetes Other - if have both T1 and T2 records, change to DM other and
# take 1st date of the two diagnoses
anyDM <-
  firstobs[disease_num == 203  |
             disease_num == 204 | disease_num == 45,][, n := .N, by = patid] #adding a variable to see how many DM types have
anyDM[n == 1, .N] #those with only 1 DM type are fine
anyDM[n == 3, .N] #those with  3 DM types (T1+T2+other):  change to other
anyDM[n == 3, disease_num := 45]
#For those with 2 DM types: if both T1 & T2 -> other; otherwise change
#the 45 to either 203 or 204
anyDM[n == 2 , .N] #N=116380
anyDM[n == 2 &
        disease_num == 45, .N, by = patid] # all have one 45 code N= 58190
anyDM[n == 2 &
        disease_num == 203, .N, by = patid] #3987 have one 203 code
anyDM[n == 2 &
        disease_num == 204, .N, by = patid] #54203 have one 204 code
#for those with a 45 & a 203 code, changing the 45 code to 203:
anyDM[n == 2 & disease_num == 45 &
        patid %in% anyDM[disease_num == 203, patid], disease_num := 203]
#for those with a 45 & a 204 code, changing the 45 code to 204
anyDM[n == 2 & disease_num == 45 &
        patid %in% anyDM[disease_num == 204, patid], disease_num := 204]
#If there were some people with both T1&T2 cases, would change both to 45
anyDM <-
  anyDM[anyDM[, .I[1], by = c("patid", "disease_num")]$V1][, n := NULL]

#Removing the Diabetes codes from the main list and then reading these back in
firstobs <- firstobs[!disease_num %in% c(45, 203, 204)]
firstobs <- rbind(firstobs, anyDM)
rm(anyDM)
setkey(firstobs, patid, eventdate)

## Stroke NOS
#No record for subarachnoid haemorrhage, ischaemic stroke or intracerebral
#haemorrhage at any time on or before the specified date
#AND
#Primary care: Stroke NOS diagnosis or history of diagnosis during consultation
#taking specified stroke diagnoses:
stroke <-
  firstobs[disease_num == 87 |
             disease_num == 193 | disease_num == 85,]
stroke <-
  stroke[stroke[, .I[1], by = c("patid")]$V1] #only care about the 1st
strokenos <- firstobs[disease_num == 192,]  #taking NOS diagnoses
stroke <-
  stroke[patid %in% strokenos$patid,] #only interested in stroke
#diagnosis if have NOS diagnosis as well
strokenos <-
  strokenos[patid %in% stroke$patid,] #only interested in NOS
#diagnosis if have stroke diagnosis as well

stroke <- rbind(stroke, strokenos) #putting these together
setkey(stroke, patid, eventdate) #ordering them
#want to know which one came second:
stroke <- stroke[stroke[, .I[2], by = c("patid")]$V1]
strokedelete <-
  stroke[disease_num == 192, obsid] #if NOS diagnosis came second,
#this is what we want to delete
firstobs <-
  firstobs[!(obsid %in% strokedelete &
               disease_num == 192), ] #deleting
#these rows in the overall data table based on obsid
rm(stroke, strokenos, strokedelete)

###

## Keeping only the first record of these conditions

setkey(firstobs, patid, eventdate)
firstobs <-
  firstobs[firstobs[, .I[1], by = c("patid", "disease_num")]$V1]

write_fst(firstobs, data_dir_CPRD("obstabsinglefirst_14Dec.fst"), 100)




#######################################################################

####  CLEANING OF CONDITIONS WITH A MULTIPLE VISIT CRITERIA  ####

# We have defined these as three codes recorded

#######################################################################

setkey(obstabmulti, patid, eventdate) #ordering these by date for each patient

#Removing under 18 records for acne
acne <- obstabmulti[disease_num == 3,]
acne[patient, on = 'patid',  `:=` (dob = i.dob)]
acneDELETE <- acne[eventdate - dob < (365.24 * 18), obsid]
obstabmulti <-
  obstabmulti[!(obsid %in% acneDELETE & disease_num == 3), ]
rm(acne, acneDELETE)

#Removing under 18 records for Asthma
asthma <- obstabmulti[disease_num == 15,]
asthma[patient, on = 'patid',  `:=` (dob = i.dob)]
asthmaDELETE <- asthma[eventdate - dob < (365.24 * 18), obsid]
obstabmulti <-
  obstabmulti[!(obsid %in% asthmaDELETE & disease_num == 15), ]
rm(asthma, asthmaDELETE)

#Removing under 18 records for dermatitis
dermatitis <- obstabmulti[disease_num == 44,]
dermatitis[patient, on = 'patid',  `:=` (dob = i.dob)]
dermatitisDELETE <- dermatitis[eventdate - dob < (365.24 * 18), obsid]
obstabmulti <-
  obstabmulti[!(obsid %in% dermatitisDELETE & disease_num == 44), ]
rm(dermatitis, dermatitisDELETE)


#Removing under 18 records for Obstructive and reflux uropathy
uropathy <- obstabmulti[disease_num == 107,]
uropathy[patient, on = 'patid',  `:=` (dob = i.dob)]
uropathyDELETE <- uropathy[eventdate - dob < (365.24 * 18), obsid]
obstabmulti <-
  obstabmulti[!(obsid %in% uropathyDELETE & disease_num == 107), ]
rm(uropathy, uropathyDELETE)


#### TEST RESULTS####

## Raised total cholesterol (d157) - Kuan def

#1. IF the highest value EVER recorded for Total Cholesterol for a patient
#on or before the specified date is greater than:
#  a) serum: 5 mmol/L
#OR
#b) serum: 193.35 mg/dL
#OR
#c) plasma: 4.8544 mmol/L
#OR
#d) plasma: 187.7184 mg/dL

## looking at total chol

#Creating a data table for the test results will include as disease diagnoses
testdiagnoses <- data.table()


#select & restrict to recorded values
totcholobs <- obstabtest[medcodeid %in%
                           c(
                             "150921000006118",
                             "1484985012",
                             "2478443019",
                             "667191000006111",
                             "405621000000116",
                             "259252011"
                           )
                         & !is.na(value) & value > 0,]

totcholobs[, .N, by = c("descr", "medcodeid")][order(-N)]
#DECIDED TO INCL ONLY THESE CODES:
#                             descr       medcodeid       N
#1:               Serum cholesterol 150921000006118   2316817
#2:  Plasma total cholesterol level      1484985012   100693
#3:   Serum total cholesterol level      2478443019   216781
#4:       Fasting serum cholesterol 667191000006111   9093
#5: Serum fasting total cholesterol 405621000000116   1467
#6:           Serum cholesterol NOS       259252011   1773


#Looking at unit frequencies
totcholobs[, .N, by = unittype][order(-N)]
totcholobs[tolower(unittype) %like% "^mmol/l",
           .(medcodeid, value, descr, unittype)] # this gives mmol/l units
ggplot(totcholobs[tolower(unittype) %like% "^mmol/l"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(totcholobs[tolower(unittype) %like% "^mmol/l"], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(totcholobs[tolower(unittype) %like% "^mmol/l", value])

totcholobs[tolower(unittype) %like% "^mg/dl" |
             tolower(unittype) %like% "^mg/100",
           .N, by = unittype] # this gives mg/dl units
ggplot(totcholobs[tolower(unittype) %like% "^mg/dl" |
                    tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(totcholobs[tolower(unittype) %like% "^mg/dl" |
                    tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(totcholobs[tolower(unittype) %like% "^mg/dl" |
                     tolower(unittype) %like% "^mg/100", value])
totcholobs[tolower(unittype) %like% "^mg/dl" |
             tolower(unittype) %like% "^mg/100", .N, by = numrangehigh]
#These all seem to be in the mmol/L range


totcholobs[!tolower(unittype) %like% "^mmol/l" &
             !tolower(unittype) %like% "^mg/dl" &
             !tolower(unittype) %like% "^mg/100" &
             !is.na(unittype), .N] # 2,506 obs have an obscure unit
totcholobs[!tolower(unittype) %like% "^mmol/l" &
             !tolower(unittype) %like% "^mg/dl" &
             !tolower(unittype) %like% "^mg/100" &
             !is.na(unittype) & value > 50,]
summary(totcholobs[!tolower(unittype) %like% "^mmol/l" &
                     !tolower(unittype) %like% "^mg/dl" &
                     !tolower(unittype) %like% "^mg/100" &
                     !is.na(unittype), value])
ggplot(totcholobs[!tolower(unittype) %like% "^mmol/l" &
                    !tolower(unittype) %like% "^mg/dl" &
                    !tolower(unittype) %like% "^mg/100" &
                    !is.na(unittype)], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(totcholobs[!tolower(unittype) %like% "^mmol/l" &
                    !tolower(unittype) %like% "^mg/dl" &
                    !tolower(unittype) %like% "^mg/100" &
                    !is.na(unittype)], aes(x = value)) +
  geom_histogram()


totcholobs[is.na(unittype), .N] #59,769(/2646624; 2.3%) no unit recorded
ggplot(totcholobs[is.na(unittype)], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(totcholobs[is.na(unittype)], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(totcholobs[is.na(unittype), value])


#Reference range 0-5 mmol/L
#Limit to 0-50 mmol/L like the others?
totcholobs[value > 50,] #only 216 of these: surely most of these are in mmol/l,
#ignore the rest

#Keeping only those within the reference range
totcholobs <- totcholobs[value <= 50,]
ggplot(totcholobs, aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(totcholobs, aes(x = value)) +
  geom_histogram(binwidth = .1)
summary(totcholobs[, value])

#Applying Kuan rule
raisedcholobs <-
  totcholobs[(tolower(descr) %like% "serum" & value > 5) |
               (tolower(descr) %like% "plasma" &
                  value > 4.8544), ]

rm(totcholobs)

raisedcholobs[, .N, by = descr][order(-N)]
#               Serum cholesterol 1078103
#   Serum total cholesterol level   97377
#  Plasma total cholesterol level   51205
#       Fasting serum cholesterol    6210
#           Serum cholesterol NOS     1137
# Serum fasting total cholesterol     874

#looking at how many people had X numbers of raised records
raisedcholobs[, numobs := .N, by = patid][, .N, keyby = numobs]
raisedcholobs[numobs > 1, .N] # 1,139,151 - 92% of those with a raised
#totchol obs have at least 2 obs

setkey(raisedcholobs, patid, eventdate)

raisedcholobs[, `:=` (disease_num = 157, numobs = NULL)]
testdiagnoses <- rbind(testdiagnoses, raisedcholobs)
rm(raisedcholobs)

##

## Low HDL-C (d90) - Kuan def

#Primary care
#1. IF FEMALE the lowest value EVER recorded for HDL Cholesterol for a patient
#on or before the specified date is less than:
#  a) serum: 1.2 mmol/L
#OR
#b) serum: 46.404 mg/dL
#OR
#c) plasma: 1.1650 mmol/L
#OR
#d) plasma: 45.0524 mg/dL

#2. IF MALE the lowest value EVER recorded for HDL Cholesterol for a patient on
#or before the specified date is less than:
#  a) serum: 1 mmol/L
#OR
#b) serum: 38.67 mg/dL
#OR
#c) plasma: 0.9709 mmol/L
#OR
#d) plasma: 37.5437 mg/dL

# looking at low hdl chol
#select & restrict to recorded values
hdlcobs <-
  obstabtest[medcodeid %in% c(
    "259232010",
    "458313013",
    "259242012",
    "259243019",
    "855771000006118",
    "259567011",
    "259564016"
  ) &
    !is.na(value) & value > 0,]
hdlcobs[, .N, c("descr", "medcodeid")][order(-N)]
#                               descr       medcodeid       N
#           Serum HDL cholesterol level       259232010 2004611
#          Plasma HDL cholesterol level       458313013   92376
#   Serum fasting HDL cholesterol level       259242012    7982
#    Serum random HDL cholesterol level       259243019    4190
#                 HDL cholesterol level 855771000006118    3394
#  Plasma fasting HDL cholesterol level       259567011    1579
#   Plasma random HDL cholesterol level       259564016     202

#https://portal.caliberresearch.org/phenotypes/hdl-cholesterol
#We extracted all values recorded in mmol/L or mg/dL
#and filtered out values below 0 or above 50 mmol/L.
#Where multiple measurements were present in a single consultation
#(identified by the consid identifier)
#we calculated the average.


#Adding in patient info for the formula
hdlcobs[patient, on = 'patid', gender := i.gender]
hdlcobs[, summary(value)]
ggplot(hdlcobs, aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(hdlcobs, aes(x = value)) +
  geom_histogram(binwidth = 1)


#Looking at code frequencies
hdlcobs[, .N, by = c("descr", "medcodeid")][order(-N)]

#Looking at unit frequencies
hdlcobs[, .N, by = unittype][order(-N)]
hdlcobs[tolower(unittype) %like% "^mmol/l",
        .(medcodeid, value, descr, unittype)] # this gives mmol/l units
ggplot(hdlcobs[tolower(unittype) %like% "^mmol/l"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(hdlcobs[tolower(unittype) %like% "^mmol/l"], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(hdlcobs[tolower(unittype) %like% "^mmol/l", value])
hdlcobs[tolower(unittype) %like% "^mmol/l" &
          value > 50, .(value, numrangehigh)]

hdlcobs[tolower(unittype) %like% "^mg/dl" |
          tolower(unittype) %like% "^mg/100",
        .N, by = unittype] # this gives mg/dl units - 164
ggplot(hdlcobs[tolower(unittype) %like% "^mg/dl" |
                 tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(hdlcobs[tolower(unittype) %like% "^mg/dl" |
                 tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(hdlcobs[tolower(unittype) %like% "^mg/dl" |
                  tolower(unittype) %like% "^mg/100", value])
hdlcobs[(tolower(unittype) %like% "^mg/dl" |
           tolower(unittype) %like% "^mg/100") & value > 50,
        .(value, numrangehigh)] #these are basically all in the mmol/l range


hdlcobs[!tolower(unittype) %like% "^mmol/l" &
          !tolower(unittype) %like% "^mg/dl" &
          !tolower(unittype) %like% "^mg/100" & !is.na(unittype),
        .N] # 2,189 obs have an obscure unit
hdlcobs[is.na(unittype), .N] #31,342 don't have a unit recorded
ggplot(hdlcobs[!tolower(unittype) %like% "^mmol/l" &
                 !tolower(unittype) %like% "^mg/dl" &
                 !tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(hdlcobs[!tolower(unittype) %like% "^mmol/l" &
                 !tolower(unittype) %like% "^mg/dl" &
                 !tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(hdlcobs[!tolower(unittype) %like% "^mmol/l" &
                  !tolower(unittype) %like% "^mg/dl" &
                  !tolower(unittype) %like% "^mg/100", value])
hdlcobs[(
  !tolower(unittype) %like% "^mmol/l" &
    !tolower(unittype) %like% "^mg/dl" &
    !tolower(unittype) %like% "^mg/100"
) & value > 50,
.(value, numrangehigh, unittype)] #these are  all in the mmol/l range

hdlcobs[value > 50,] #only 105 of these - surely most of these are in mmol/l,
#ignore the rest

#Taking those within the reference range
hdlcobs <- hdlcobs[value <= 50,]
ggplot(hdlcobs, aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(hdlcobs, aes(x = value)) +
  geom_histogram(binwidth = 0.1)
summary(hdlcobs[, value])

#Applying the CALIBER algorithm
lowhdlcobs <- hdlcobs[(gender == "F"  &
                         tolower(descr) %like% "serum" &
                         value < 1.2) |
                        (gender == "M"  &
                           tolower(descr) %like% "serum" &
                           value < 1) |
                        (gender == "F" &
                           tolower(descr) %like% "plasma" &
                           value < 1.1650) |
                        (gender == "M" &
                           tolower(descr) %like% "plasma" &
                           value < 0.9709),]

rm(hdlcobs)

lowhdlcobs[, .N, by = descr][order(-N)]
# Serum HDL cholesterol level         363350
#Plasma HDL cholesterol level         18716
#Serum fasting HDL cholesterol level  1160
#Serum random HDL cholesterol level   775
#Plasma fasting HDL cholesterol level 301
#Plasma random HDL cholesterol level  34

#looking at how many people had X numbers of raised records
lowhdlcobs[, numobs := .N, by = patid][, .N, keyby = numobs]
setkey(lowhdlcobs, patid, eventdate)


lowhdlcobs[, `:=` (disease_num = 90,
                   gender = NULL,
                   numobs = NULL)]
testdiagnoses <- rbind(testdiagnoses, lowhdlcobs)
rm(lowhdlcobs)

###

## Raised LDL-C (d156) - Kuan def

# IF the highest value EVER recorded for LDL Cholesterol for a patient on
#or before the specified date is greater than:
#  a) serum: 3 mmol/L
#OR
#b) serum: 116.01 mg/dL
#OR
#c) plasma: 2.9126 mmol/L
#OR
#d) plasma: 112.6311 mg/dL

# looking at raised ldl chol
#select & restrict to recorded values
ldlcobs <-
  obstabtest[medcodeid %in% c(
    "259233017",
    "1488764011",
    "458314019",
    "259246010",
    "857141000006116",
    "259247018",
    "259570010",
    "259569014"
  )
  & !is.na(value) & value > 0, ]

#https://portal.caliberresearch.org/phenotypes/ldl-cholesterol
#Excludes below 0 or above 50 mmol/L
#"We extracted all values recorded in mmol/L or mg/dL
#and filtered out values below 0 or above 50 mmol/L.

#Looking at code frequencies
ldlcobs[, .N, by = c("descr", "medcodeid")][order(-N)]

#Looking at unit frequencies
ldlcobs[, .N, by = unittype][order(-N)]
ldlcobs[tolower(unittype) %like% "^mmol/l",
        .(medcodeid, value, descr, unittype)] # this gives mmol/l units
ggplot(ldlcobs[tolower(unittype) %like% "^mmol/l"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(ldlcobs[tolower(unittype) %like% "^mmol/l"], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(ldlcobs[tolower(unittype) %like% "^mmol/l" , value])
ldlcobs[(tolower(unittype) %like% "^mmol/l") & value > 50,
        .(value, numrangehigh, unittype)] #47 of these

ldlcobs[tolower(unittype) %like% "^mg/dl" |
          tolower(unittype) %like% "^mg/100",
        .N, by = unittype] # this gives mg/dl units - only 84
ggplot(ldlcobs[tolower(unittype) %like%  "^mg/dl" |
                 tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(ldlcobs[tolower(unittype) %like%  "^mg/dl" |
                 tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(ldlcobs[tolower(unittype)  %like% "^mg/dl" |
                  tolower(unittype) %like% "^mg/100", value])
ldlcobs[(tolower(unittype) %like% "^mg/dl" |
           tolower(unittype) %like% "^mg/100")
        & value > 50,
        .(value, numrangehigh, unittype)] #7 of these

ldlcobs[!tolower(unittype) %like% "^mmol/l" &
          !tolower(unittype) %like% "^mg/dl" &
          !tolower(unittype) %like% "^mg/100" & !is.na(unittype),
        .N] # 21,653 obs have an obscure unit
ldlcobs[is.na(unittype), .N] #22,262 don't have a unit recorded
ggplot(ldlcobs[!tolower(unittype) %like% "^mmol/l" &
                 !tolower(unittype) %like% "^mg/dl" &
                 !tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(ldlcobs[!tolower(unittype) %like% "^mmol/l" &
                 !tolower(unittype) %like% "^mg/dl" &
                 !tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(ldlcobs[!tolower(unittype) %like% "^mmol/l" &
                  !tolower(unittype) %like% "^mg/dl" &
                  !tolower(unittype) %like% "^mg/100", value])
ldlcobs[(
  !tolower(unittype) %like% "^mmol/l" &
    !tolower(unittype) %like% "^mg/dl" &
    !tolower(unittype) %like% "^mg/100"
) & value > 50,
.(value, numrangehigh, unittype)] #7 of these


ldlcobs[value > 50, .N] #only 61 greater than 50

#assuming all actually recorded in mmol/l. Capping at 50
ldlcobs <- ldlcobs[value <= 50, ]

ggplot(ldlcobs, aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(ldlcobs, aes(x = value)) +
  geom_histogram(binwidth = .1)
summary(ldlcobs[, value])


#Applying Kuan rule
raisedldlcobs <-
  ldlcobs[(tolower(descr) %like% "serum" & value > 3) |
            (tolower(descr) %like% "plasma" &
               value > 2.9126) |
            (!tolower(descr) %like% "serum" &
               !tolower(descr) %like% "plasma" &
               value > 3),]


rm(ldlcobs)

raisedldlcobs[, .N, by = descr][order(-N)]
#          Serum LDL cholesterol level 464974
#     Calculated LDL cholesterol level 109960
#         Plasma LDL cholesterol level  22999
#  Serum fasting LDL cholesterol level   4045
#                LDL cholesterol level   3533
#   Serum random LDL cholesterol level   1710
# Plasma fasting LDL cholesterol level    911
#  Plasma random LDL cholesterol level    258

#looking at how many people had X numbers of raised records
raisedldlcobs[, numobs := .N, by = patid][, .N, keyby = numobs]
setkey(raisedldlcobs, patid, eventdate)

raisedldlcobs[, `:=`(disease_num = 156, numobs = NULL)]
testdiagnoses <- rbind(testdiagnoses, raisedldlcobs)

rm(raisedldlcobs)


###

## Raised triglycerides (d158)

#Primary care
#1. IF the highest value EVER recorded for Triglycerides for a patient on or
#before the specified date is greater than:
#  a) serum: 1.7 mmol/L
#OR
#b) serum: 150.569 mg/dL
#OR
#c) plasma: 1.6521 mmol/L
#OR
#d) plasma: 146.3256 mg/dL

#looking at raised triglycerides
#https://portal.caliberresearch.org/phenotypes/triglycerides
#We extracted all values recorded in mmol/L (SUM lookup 96) or
#mg/dL (SUM lookup 82) and
#filtered out values below 0 or above 50 mmol/L.
#Where multiple measurements were present in a single consultation
#we calculated the average.



#select triglycerides codes
triglycobs <- obstabtest[medcodeid %in% c("145471000006115",
                                          "259263014",
                                          "854481000006110") &
                           !is.na(value) & value > 0, ]

#Looking at code frequencies
triglycobs[, .N, by = c("descr", "medcodeid")][order(-N)]
#     Serum triglycerides    145471000006115 1699110
# Fasting triglycerides     854481000006110   11643
# Serum triglycerides NOS         259263014    2918

#Looking at unit frequencies
triglycobs[, .N, by = unittype][order(-N)]

triglycobs[tolower(unittype) %like% "^mmol/l", .N]
ggplot(triglycobs[tolower(unittype) %like%  "^mmol/l"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(triglycobs[tolower(unittype) %like%  "^mmol/l"], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(triglycobs[tolower(unittype)  %like% "^mmol/l", value])
triglycobs[(tolower(unittype) %like% "^mmol/l") & value > 50,
           .(value, numrangehigh, unittype)] #120 of these

triglycobs[tolower(unittype) %like% "^mg/dl" |
             tolower(unittype) %like% "^mg/100", .N] #43
ggplot(triglycobs[tolower(unittype) %like%  "^mg/dl" |
                    tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(triglycobs[tolower(unittype) %like%  "^mg/dl" |
                    tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(triglycobs[tolower(unittype)  %like% "^mg/dl" |
                     tolower(unittype) %like% "^mg/100", value])
triglycobs[(tolower(unittype) %like% "^mg/dl" |
              tolower(unittype) %like% "^mg/100") & value > 50,
           .(value, numrangehigh, unittype)] #8 of these

triglycobs[tolower(unittype) %like% "^g/l" , .N] #10253
ggplot(triglycobs[tolower(unittype) %like%  "^g/l"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(triglycobs[tolower(unittype) %like%  "^g/l"], aes(x = value)) +
  geom_histogram(binwidth = .1)
summary(triglycobs[tolower(unittype)  %like% "^g/l", value])
triglycobs[(tolower(unittype) %like% "^g/l") & value > 50,
           .(value, numrangehigh, unittype)] #1 of these

triglycobs[!tolower(unittype) %like% "^mmol/l" &
             !tolower(unittype) %like% "^mg/dl" &
             !tolower(unittype) %like% "^g/l" &
             !tolower(unittype) %like% "^mg/100" &
             !tolower(unittype) %like% "^g/l" & !is.na(unittype),
           .N] # 1665 obs have an obscure unit
triglycobs[is.na(unittype), .N] #27,389 don't have a unit recorded...
ggplot(triglycobs[!tolower(unittype) %like% "^mmol/l" &
                    !tolower(unittype) %like% "^mg/dl" &
                    !tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(triglycobs[!tolower(unittype) %like% "^mmol/l" &
                    !tolower(unittype) %like% "^mg/dl" &
                    !tolower(unittype) %like% "^mg/100"], aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(triglycobs[!tolower(unittype) %like% "^mmol/l" &
                     !tolower(unittype) %like% "^mg/dl" &
                     !tolower(unittype) %like% "^mg/100", value])
triglycobs[(
  !tolower(unittype) %like% "^mmol/l" &
    !tolower(unittype) %like% "^mg/dl" &
    !tolower(unittype) %like% "^mg/100"
) & value > 50,
.(value, numrangehigh, unittype)] #8 of these

triglycobs[value > 50, .N] #only 136 greater than 50

#assuming all actually recorded in mmol/l or g/l if specified . Capping at 50
triglycobs <- triglycobs[value <= 50, ]

ggplot(triglycobs, aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(triglycobs, aes(x = value)) +
  geom_histogram(binwidth = .1)
summary(triglycobs[, value])


#Applying Kuan rule

raisedtriglycobs <- triglycobs[(!tolower(unittype) %like% "^g/l" &
                                  !tolower(descr) %like% "plasma" &
                                  value > 1.7) |
                                 (!tolower(unittype) %like% "^g/l" &
                                    tolower(descr) %like% "plasma" &
                                    value > 1.6521) |
                                 (tolower(unittype) %like% "^g/l" &
                                    !tolower(descr) %like% "plasma" &
                                    value > 1.50569) |
                                 (tolower(unittype) %like% "^g/l" &
                                    tolower(descr) %like% "plasma" &
                                    value > 1.463256) , ]

rm(triglycobs)
raisedtriglycobs[, .N, by = c("descr", "medcodeid")][order(-N)]

#     Serum triglycerides 145471000006115 564447
#   Fasting triglycerides 854481000006110   4858
# Serum triglycerides NOS       259263014   1197

#looking at how many people had X numbers of raised records
raisedtriglycobs[, numobs := .N, by = patid][, .N, keyby = numobs]
setkey(raisedtriglycobs, patid, eventdate)


raisedtriglycobs[, `:=` (disease_num = 158, numobs = NULL)]
testdiagnoses <- rbind(testdiagnoses, raisedtriglycobs)
rm(raisedtriglycobs)

###

## Obesity (d105) test results - BMI

#BMI > 30
#OR
#Weight (in kg)  / height (kg) ^2 > 30

# bmi test results
#selecting bmi results and restricting to sensible values
bmiobs <- obstabtest[medcodeid %in%
                       c("100716012", "923861000006112", "3484801000006114") &
                       !is.na(value) & value > 10 & value < 200,]
bmiobs[, .N]
ggplot(bmiobs, aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(bmiobs, aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(bmiobs[, value])
bmiobs[value > 80, #this is the kuan cutoff
       .(value, numrangehigh, unittype)] #1453 of these
ggplot(bmiobs[value > 80], aes(x = value)) +
  geom_histogram(binwidth = 1)

#https://portal.caliberresearch.org/phenotypes/body-mass-index
#We excluded measurements below 10 kg/m2 or above 80 kg/m2.

#Looking at code frequencies
bmiobs[, .N, by = descr][order(-N)]

#Looking at unit frequencies
bmiobs[, .N, by = unittype][order(-N)]

summary(bmiobs$value)

#Applying Kuan rule - >30
#Assumng the units are kg/m2 regardless
obesebmi <- bmiobs[value > 30,]

#looking at how many people had X numbers of raised records
obesebmi[, numobs := .N, by = patid][, .N, keyby = numobs]

#looking at which codes
obesebmi[, .N, by = c("descr", "medcodeid")][order(-N)]


#DECIDED TO USE ONLY THESE 2 CODES:
#                   descr        medcodeid       N
#1:       Body mass index        100716012 1193268
#2:       Body mass index  923861000006112      63
#3: BMI - Body mass index 3484801000006114       4


setkey(obesebmi, patid, eventdate)
obesebmi[, `:=` (numobs = NULL, disease_num = 105)]

testdiagnoses <- rbind(testdiagnoses, obesebmi)
rm(obesebmi)


# Obesity test results - weight & height
#Primary care (these are for GOLD database, Aurum doesnt have enttype etc)
#1. Obesity diagnosis or history of diagnosis or procedure during a consultation
#OR
#2. IF  enttype = 13 (Weight) available on or before specified date
#AND data3 not missing,
#BMI = data3. If BMI > 30, patient is defined as having had Obesity.
#OR
#3. If enttype = 13 (Weight) available on or before specified date
#AND data3 missing,
#BMI = data1 (enttype 13) /(data2 ^2) (enttype 14 = Height).
#If BMI > 30, patient is defined as having had Obesity.
#IF height not available on same eventdate as weight, use most recent
#height for age > 18 years.

# weight  test results
#selecting weight results
weightobs <- obstabtest[medcodeid %in% c("253677014",
                                         "253688015",
                                         "1910901000006116",
                                         "1910911000006118")
                        & !is.na(value) & value > 0, ]

#Looking at code frequencies
weightobs[, .N, by = descr][order(-N)]
#1:     O/E - weight 4412289
#2: O/E - weight NOS    1738
#3: Estimated weight     215
#4:  Reported weight      17

#Looking at unit frequencies
weightobs[, .N, by = unittype][order(-N)]
summary(weightobs$value)

#Applying Kuan rule
#For now lets just keep sensible units
weightunit <-
  c(
    "kg",
    "decimal stones",
    "kgs",
    "kilos",
    "kg.",
    "kilograms",
    "weight in kilos",
    "weight (stones/pounds)",
    "stones",
    "st",
    "stone",
    "llbs"
  )
weightobs[!(tolower(unittype) %in%  weightunit | is.na(unittype)),
          .N] #only 42 so remove
weightobs <-
  weightobs[(tolower(unittype) %in%  weightunit |
               is.na(unittype)),]


ggplot(weightobs[(tolower(unittype) %in% c("kg", "kgs", "kilos", "kg.",
                                           "kilograms", "weight in kilos"))],
       aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(weightobs[(tolower(unittype) %in% c("kg", "kgs", "kilos", "kg.",
                                           "kilograms", "weight in kilos"))],
       aes(x = value)) +  geom_histogram()
summary(weightobs[(tolower(unittype) %in% c("kg", "kgs", "kilos", "kg.",
                                            "kilograms", "weight in kilos")),
                  value])

ggplot(weightobs[(
  tolower(unittype) %in% c(
    "decimal stones",
    "weight (stones/pounds)",
    "stones",
    "st",
    "stone",
    "llbs"
  )
)],
aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(weightobs[(
  tolower(unittype) %in% c(
    "decimal stones",
    "weight (stones/pounds)",
    "stones",
    "st",
    "stone",
    "llbs"
  )
)],
aes(x = value)) + geom_histogram(binwidth = 1)
summary(weightobs[(
  tolower(unittype) %in% c(
    "decimal stones",
    "weight (stones/pounds)",
    "stones",
    "st",
    "stone",
    "llbs"
  )
), value])

weightobs[(
  tolower(unittype) %in% c(
    "decimal stones",
    "weight (stones/pounds)",
    "stones",
    "st",
    "stone",
    "llbs"
  )
) & value > 12.6, #this is the kuan cutoff in stone
.(value, numrangehigh, unittype)] #560 of these

ggplot(weightobs[(
  tolower(unittype) %in% c(
    "decimal stones",
    "weight (stones/pounds)",
    "stones",
    "st",
    "stone",
    "llbs"
  )
) & value > 12.6], aes(x = value)) +
  geom_histogram(binwidth = .1)


#Converting stones & lls to kgs
#1 stone = 6.35029318 kg; 1kg = 2.20462llbs
weightobs[tolower(unittype) %like% "stone" |
            tolower(unittype) == "st" ,
          `:=` (value = value * 6.35029318, unittype = "kg")]
weightobs[tolower(unittype) %in% "llbs",
          `:=` (value = value / 2.20462, unittype = "kg")]

# What about those without units??
weightobs[is.na(unittype), .(medcodeid, value, descr)]
weightobs[is.na(unittype), .N] #43319
ggplot(weightobs[is.na(unittype), ], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(weightobs[is.na(unittype), ], aes(x = value)) +
  geom_histogram()
summary(weightobs[is.na(unittype), value])



# Probably need a max weight... using 500kg
summary(weightobs$value)
summary(weightobs[value > 5 & value < 500, value])
summary(weightobs[tolower(unittype) %in%  weightunit  &
                    value > 5 & value < 500, value])
summary(weightobs[!tolower(unittype) %in%  weightunit  &
                    value > 5 & value < 500, value])
#these have very similar distributions
rm(weightunit)




# restricting values a bit - maybe need to restrict more?
weightobs <- weightobs[value > 5 & value < 500 &
                         !is.na(value) , ]
ggplot(weightobs, aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(weightobs, aes(x = value)) +
  geom_histogram(binwidth = 1)
summary(weightobs[, value])


weightobs[, .N, by = descr][order(-N)]
#              descr        medcodeid    N
#1:     O/E - weight        253677014 4410383
#2: O/E - weight NOS        253688015    1734
#3: Estimated weight 1910901000006116     215
#4:  Reported weight 1910911000006118      17



# Height results
# r height  test results
#selecting height results
heightobs <- obstabtest[medcodeid %in%
                          c("253669010",
                            "253676017",
                            "1910931000006112",
                            "1910921000006114")
                        & !is.na(value) & value > 0,]


#Looking at code frequencies
heightobs[, .N, by = c("descr", "medcodeid")][order(-N)]
#O/E - height         253669010           2493398
#O/E - height NOS     253676017           691
#Estimated height     1910921000006114    224
#Reported height      1910931000006112    22


#Looking at unit frequencies
heightobs[, .N, by = unittype][order(-N)]
#these are nearly all metres & cms
heightunit <- c("cm", "m", "metres", "cms", "meters") #senible unit


summary(heightobs$value)

#heights in cm
ggplot(heightobs[(tolower(unittype) %in% c("cm", "cms"))], aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(heightobs[(tolower(unittype) %in% c("cm", "cms"))], aes(x = value)) +
  geom_histogram()
summary(heightobs[(tolower(unittype) %in% c("cm", "cms")), value])
#looking at implausiable values
heightobs[(tolower(unittype) %in% c("cm", "cms")) &
            value > 250, .N] #1593 - will discard these as seem implausible
summary(heightobs[(tolower(unittype) %in% c("cm", "cms")) &
                    value > 250, value]) #
heightobs[(tolower(unittype) %in% c("cm", "cms")) &
            value < 50, .N] #7903 -
summary(heightobs[(tolower(unittype) %in% c("cm", "cms")) &
                    value < 50, value]) #most of these are likely to be in m





#heights in m
ggplot(heightobs[(tolower(unittype) %in% c("m", "metres", "meters"))], 
       aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(heightobs[(tolower(unittype) %in% c("m", "metres", "meters"))], 
       aes(x =  value)) +
  geom_histogram()
summary(heightobs[(tolower(unittype) %in% c("m", "metres", "meters")), value])
#looking at implausiable values
heightobs[(tolower(unittype) %in% c("m", "metres", "meters")) &
            value > 2.5, .N] #637 - a lot of these are likely to be in cm
summary(heightobs[(tolower(unittype) %in% c("m", "metres", "meters")) &
                    value > 2.5, value]) #
heightobs[(tolower(unittype) %in% c("m", "metres", "meters")) &
            value < .5, .N] #62 - will discard these as seem implausible
summary(heightobs[(tolower(unittype) %in% c("m", "metres", "meters")) &
                    value < .5, value]) #



heightobs[!tolower(unittype) %in% heightunit, ] #very few, but they look 
#plausible aas eith cm or m
rm(heightunit)

#For observations with weird or no unit , making assumption that those <3 are 
#in metres and those >2.5 are in cms. this is because no obs seem 
#to have ft/inch recorded

#assuming plausible height is 0.5-2.5m
#Going to assume that all values 0.5 < x <2.5 are in m and
#all 50 < x <250 are cm (approx world's tallest & min 0.5  world's shortest).
#excluding the rest

heightobs[, `:=`(value = ifelse(value < 2.5,
                                value, value / 100),
                 unittype = "m")]
heightobs <-
  heightobs[value < 2.5 &
              value > 0.5] # excluding heights above 2.5m & below 0.5m


ggplot(heightobs, aes(x = value)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(heightobs, aes(x = value)) +
  geom_histogram(binwidth = .02)
summary(heightobs[, value])



#These are the height codes:
heightobs[, .N, by = descr][order(-N)]
#     O/E - height  2490091
# O/E - height NOS     691
# Estimated height     221
#  Reported height      22


# Calculating BMI based on height & weight 
#only want obs of patients who have both height and weight
heightobs <- heightobs[patid %in% weightobs[, patid],]
weightobs <- weightobs[patid %in% heightobs[, patid],]


#add in height when recorded on same date as weight
weightobs[heightobs, on = c("patid", "eventdate"), height := i.value]

#combining weight and height so can add in previous height values if needed
weightobs[, type := "weight"]
heightobs[, `:=` (height = value, type = "height")]

#adding in height observations for weight obs without height on same day
weightobs <- rbind(weightobs,
                   heightobs[patid %in% weightobs[is.na(height), patid]])
rm(heightobs)
#carrying forward height obs
setkey(weightobs, patid, eventdate)
weightobs[, height := nafill(height, "locf"), by = patid]
#then carrying backwards if still empty:
weightobs[, height := nafill(height, "nocb"), by = patid] 

weightobs <-
  weightobs[type != "height",][, type := NULL] #getting rid of the height obs

weightobs[, bmi := value / height ^ 2]
ggplot(weightobs, aes(x = bmi)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(weightobs, aes(x = bmi)) +
  geom_histogram(binwidth = .02)
summary(weightobs[, bmi])


#excluding the extreme observations 
weightobs <- weightobs[bmi > 10 & bmi < 200]
weightobs[, bmi := value / height ^ 2]
ggplot(weightobs, aes(x = bmi)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  )
ggplot(weightobs, aes(x = bmi)) +
  geom_histogram(binwidth = .02)
summary(weightobs[, bmi])

obese <-
  weightobs[value / height ^ 2 > 30 &
              value / height ^ 2 < 200, ] [, disease_num := 105]
rm(weightobs)

setkey(obese, patid, eventdate)
obese[, height := NULL]
obese[, `:=` (value = bmi, descr = "calc_bmi")][, bmi := NULL]
obese[, .N, by = medcodeid][order(-N)]
testdiagnoses <- rbind(testdiagnoses, obese)
rm(obese)

testdiagnoses[, c("unittype", "descr") := NULL]

testdiagnoses <-
  rbind(testdiagnoses, obstabmulti[disease_num == 105, ])

setkey(testdiagnoses, patid, eventdate)



#two visit ever function
twovisits <- function(x, y) {
  #x: dataset to use; y = disease name
  allcodes <- #taking all the records with that disease
    x[disease_num == y, .(patid, obsid, eventdate, medcodeid)] 
  tmp <-#counting how many each patient has
    allcodes[, .N, by = patid] 
  tmp2 <-
    tmp[N >= 2, patid] #only interested in those patients with  2+ records
  multcodes <-
    allcodes[patid %in% tmp2,] #only interested in patients with 2+ records
  setkey(multcodes, patid, eventdate)
  multcodes[, prevcode := as.IDate(shift_bypid(eventdate, 1L, patid, 0L))] 
  #NB: 1st event per patient set at 1970-01-01
  
  multcodes <-
    multcodes[prevcode != eventdate,] #keeping only 1 obs record from each visit
  setkey(multcodes, patid, eventdate)
  #taking the first record which meets the multiple visit criteria:
  multcodes <-
    multcodes[multcodes[, .I[2], keyby = c("patid")]$V1] 
  diagnosisobs <-
    multcodes[, obsid] #pulling out the obsid for the records of interest
  diagnosis <-
    x[obsid %in% diagnosisobs,] #extract only the records meeting the criteria
}


testval <- unique(testdiagnoses[, disease_num])
library(foreach)
library(doParallel)
registerDoParallel(10L)


multtest <- data.table()
multtest <-
  foreach(i = testval, .combine = rbind) %do% {
    #this makes one big table combining all the disease together
    twovisits(testdiagnoses, i)
  }
rm(i, twovisits, testval)

multtest[, .N, keyby = disease_num]

write_fst(multtest, data_dir_CPRD("multtest_14Dec.fst"), 100)





#three visit function
threevisits <- function(x, y) {
  #x: dataset to use; y = disease name
  allcodes <- #taking all the records with that disease
    x[disease_num == y, .(patid, obsid, eventdate, medcodeid)] 
  tmp <-
    allcodes[, .N, by = patid] #counting how many each patient has
  tmp2 <-
    tmp[N >= 3, patid] #only interested in those patients with  3+ records
  multcodes <-
    allcodes[patid %in% tmp2,] #only interested in  patients with 3+ records
  setkey(multcodes, patid, eventdate)
  multcodes[, prevcode := as.IDate(shift_bypid(eventdate, 1L, patid, 0L))] 
  #NB: 1st event per patient set at 1970-01-01
  multcodes <-
    multcodes[prevcode != eventdate,] #keeping only 1 obs record from each visit
  multcodes[, prev2code := as.IDate(shift_bypid(eventdate, 2L, patid, 0L))] 
  #NB: 1st two events per patient set at 1970-01-01
  
  #extracting records where previous but 1 record happened within the last year 
  #& each record is from a different date:
  oneyeargap <-
    multcodes[eventdate - prev2code < 365,]   oneyeargap <-
    oneyeargap[oneyeargap[, .I[1], keyby = c("patid")]$V1] #taking the first 
  #record which meets the multiple visit criteria
  diagnosisobs <-
    oneyeargap[, obsid] #pulling out the obsid for the records of interest
  diagnosis <-
    x[obsid %in% diagnosisobs,] #extract only the records meeting the criteria
}


#including obesity in with the test values, so want to remove them here

obstabmulti <- obstabmulti[disease_num != 105, ]

setkey(codelistmultimedid, disease_num)
mv <- unique(codelistmultimedid[, disease_num])

library(foreach)
data.table::setDTthreads(1) #this is so that don't use all the processors

registerDoParallel(10L)
#Extracting the 3rd record within a year for diseases with multi visit criteria
multdiagnosis <- data.table()
multdiagnosis <-
  foreach(i = mv, .combine = rbind) %do% {
    #this makes one big table combining all the disease together
    threevisits(obstabmulti, i)
  }
rm(i, threevisits)

#automatically including certain diagnoses
autoinclcodes <-
  c(fullcodelistmedid[disease_num == 6 & #allergy
                        tolower(descr) %like% "chronic", medcodeid],
    fullcodelistmedid[disease_num == 2 & # abdominal hernia
                        tolower(descr) %like% "irreducible", medcodeid],
    fullcodelistmedid[disease_num == 48 & #Diaphragmatic hernia
                        tolower(descr) %like% "irreducible", medcodeid],
    fullcodelistmedid[disease_num == 114 & #Pericardial effusion
                        tolower(descr) %like% "chronic", medcodeid]) 
    
    
    autoincl <- obstabmulti[medcodeid %in% autoinclcodes, ]
    #rm(obstabmulti)
    
    #combining those with 3 visit criteria and autoincl
    multdiagnosis <-
      rbind(multdiagnosis, autoincl) 
    rm(autoincl, autoinclcodes)
    
    #Getting rid of any duplicates of diseases with autoincl criteria
    setkey(multdiagnosis, patid, eventdate)
    multdiagnosis <-
      multdiagnosis[multdiagnosis[, .I[1], 
                                  keyby = c("patid", "disease_num")]$V1]
    multdiagnosis[, .N, keyby = disease_num]
    
    write_fst(multdiagnosis, data_dir_CPRD("multdiagnosis_14Dec.fst"), 100)
    
    
    #######################################################################
    
    #### COMBINE ALL THE 'DIAGNOSES' TOGETHER ####
    
    #######################################################################
    
    #combining the diagnosis records for single read code and 
    #multi read code diseases
    obstab <- rbind(firstobs, multdiagnosis, multtest)
    setkey(obstab, patid, eventdate)
    
    obstab <-
      obstab[obstab[, .I[1], keyby = c("patid", "disease_num")]$V1]
    obstab[, .N] #2931338
    
    rm(firstobs, multdiagnosis, multtest)
    setkey(obstab, patid, eventdate)
    
    write_fst(obstab, data_dir_CPRD("obstab14Dec.fst"), 100)
    
    rm(list = ls())
    gc()