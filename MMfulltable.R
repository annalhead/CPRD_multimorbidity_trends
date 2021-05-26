
# Formatting data into a nice table

#This file is for creating a data.table with a row for each patient for each 
#year they are included in the study
#this is to make it easy to create prevalence & incidence graphs

#Our patient sample was derived from May 2020 Aurum database

#Inputs:
#1. patient info data.table (output of datacleaning file - patient.fst)
#2. observation records of conditions of interest  data.table (output of 
#datacleaning file - obstab.fst)
#3. summary disease list (diseasesummary.csv, available here:
#https://github.com/annalhead/CPRD_multimorbidity_codelists )

#######################################################################

#### SET UP #### 

#######################################################################
library(ggplot2)
library(Rcpp)
library(fst)
library(data.table)
library(CKutils) # This is one of CK packages. Can be installed from GitHub

sourceCpp("./aux_functions.cpp", cacheDir = ".cache")
data.table::setDTthreads(10) #this is so that don't use all the processors 


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

# Creating a function to add a col for 5yr age group: 
# x is the dataset to be used (must include yob), y is the year for calculating
# age

agegroup5yfn <- function(x) { 
  lev <-
    c(
      "15-19 years",
      "20-24 years",
      "25-29 years",
      "30-34 years",
      "35-39 years",
      "40-44 years",
      "45-49 years",
      "50-54 years",
      "55-59 years",
      "60-64 years",
      "65-69 years",
      "70-74 years",
      "75-79 years",
      "80-84 years",
      "85-89 years",
      "90plus years"
    )
  x[, agegrp5 := fcase(
    (age) >= 90L,           factor(lev[16], levels = lev),
    between(age, 85, 89), factor(lev[15], levels = lev),
    between(age, 80, 84), factor(lev[14], levels = lev),
    between(age, 75, 79), factor(lev[13], levels = lev),
    between(age, 70, 74), factor(lev[12], levels = lev),
    between(age, 65, 69), factor(lev[11], levels = lev),
    between(age, 60, 64), factor(lev[10], levels = lev),
    between(age, 55, 59), factor(lev[9], levels = lev),
    between(age, 50, 54), factor(lev[8], levels = lev),
    between(age, 45, 49), factor(lev[7], levels = lev),
    between(age, 40, 44), factor(lev[6], levels = lev),
    between(age, 35, 39), factor(lev[5], levels = lev),
    between(age, 30, 34), factor(lev[4], levels = lev),
    between(age, 25, 29), factor(lev[3], levels = lev),
    between(age, 20, 24), factor(lev[2], levels = lev),
    between(age, 15, 19), factor(lev[1], levels = lev)
  )
  ]
}

# Reading in the data  ----
patient <- read_fst(data_dir_CPRD("patient.fst"), 
                    columns = c("patid", "gender", "imd", 
                                "black", "yob", "dob", "reg1yr", "censordate", 
                                "censorreason", "region"), 
                    as.data.table = TRUE) 

#Observations of interest:
obstab <- read_fst(data_dir_CPRD("obstab14Dec.fst"), 
                   columns = c("patid", "eventdate", "disease_num"), 
                   as.data.table = TRUE) 

#Summary disease list
diseasesum <- fread(data_dir_DL("DiseaseSumm.csv"), 
                    header = TRUE, sep = ",", 
                    select = c("disease_num", "system_num"))

#make a list of relevant patients i.e. those registered in that year
tt <- data.table( # Auxiliary table
  y = 2004:2019,
  year = 2004:2019)

patient[, `:=`(
  year_reg1yr = year(reg1yr),
  year_censordate = year(censordate)
)]
pats <- patient[tt, 
                on = c("year_reg1yr <= y", "year_censordate >= y")
][year - yob >= 18, ]
#uniqueN(pats[, patid]) #N= 991250
#lost 8990 here because not within study period 

#patient[year(reg1yr) == 2020, .N] #8951 only registered in 2020

#don't want anybody with less than 1 day vetween reg1yr & censor
pats <- pats[censordate - reg1yr > 0] #removes 7 people

obstab[, year_eventdate := year(eventdate)]
obstab[, summary(year_eventdate)]
setkey(obstab, patid, eventdate)
incd_date <- obstab[obstab[,.I[2], by = c("patid")]$V1]
incd_date <- incd_date[!is.na(patid)]

cmm3plus_incdate <- obstab[obstab[,.I[3], by = c("patid")]$V1]
cmm3plus_incdate <- cmm3plus_incdate[!is.na(patid)]

tt <- pats[, .(patid = unique(patid)), keyby = year][, y := year] 
setkey(tt, patid, year)



# incidence ----
incd <- tt[obstab, on = c("patid", "y == year_eventdate")
][, y := NULL] 
incd <- incd[!is.na(year)]
setkey(incd, patid, year)
incd_w <- dcast(incd, patid + year ~ disease_num, fun.aggregate = length)
# View(head(incd, 100))
# View(head(incd_w, 100))
rm(incd)
# some disease incd appears multiple times in a year
incd_w[, lapply(.SD, max), .SDcols = as.character(1:212)]
for (j in as.character(1:212)) {
  # all incd denoted with 1L
  set(incd_w, NULL, j, fifelse(incd_w[[j]] > 0L, 1L, 0L))
}
#NB will get a warning for col 154(Ptosis) as removed
rm(j)

prvl <- tt[obstab, on = c("patid", "y >= year_eventdate")
][, y := NULL]
prvl <- prvl[!is.na(year)]
setkey(prvl, patid, year)
# nrow(prvl) # 17137608

# prvl does not contain completely healthy person-years. Let's add them
prvl <- prvl[tt, on = .(patid, year)]

# nrow(prvl) # 19882406
# Let's create disease_num = 0L to denote healthy
prvl[is.na(disease_num), disease_num := 0L] #disease_num = 0 -> healthy for 
#that year

prvl_w <- dcast(prvl, patid + year ~ disease_num, fun.aggregate = length)
# View(head(prvl, 100))
# View(head(prvl_w, 100))
rm(prvl)
prvl_w[, lapply(.SD, max), .SDcols = as.character(1:212)]
prvl_w[, `0` := 0L]
for (j in as.character(1:212)) { # all prvl denoted with 2L
  set(prvl_w, NULL, j, fifelse(prvl_w[[j]] > 0L, 2L, 0L))
}
#NB will get a warning for col 154(Ptosis) as removed
rm(j)
setnames(incd_w, as.character(1:212), paste0("i_", 1:212))
CKutils::absorb_dt(prvl_w, incd_w, on = c("patid", "year"))

for (j in as.character(1:212)) {
  i_j <- paste0("i_", j)
  absorb_incd(prvl_w[[j]], prvl_w[[i_j]])
}
prvl_w[, paste0("i_", 1:212) := NULL]
rm(incd_w, j, i_j)

# Housekeeping
setkey(prvl_w, patid, year) # IMPORTAND for C++ functions

prvl_w[, pid := rleid(patid)]
prvl_w[, pid := mk_new_simulant_markers(pid)] # aux col

# Ensure that for each patid an 1 is always followed by 2 and a 2 is 
#always followed by 2
invisible(prvl_w[, lapply(.SD, fill_gaps, pid), .SDcols = as.character(1:212)])

# Assuming that you have disease 3 & 5 that are cured after i.e. 10 years you 
#can use  invisible(prvl_w[, lapply(.SD, cured, 10L, pid), 
#.SDcols = as.character(c(4, 5))])
prvl_w[, pid := NULL]
prvl_w[, `154` := 0] #I removed #154 (ptosis), but putting all as 0s so as 
#not to have to redo all the numbering


# 1. incidence multimorbidity is when a disease row has
#   a. exactly one 2 and at least one 1
#   b. no 2 and at least two 1
# 2. prevalent multimorbidity is when a disease row has at least 2 twos
prvl_w[, count2s := Reduce("+", lapply(.SD, function(x)
  x == 2L)), .SDcols = as.character(1:212)]
prvl_w[, count1s := Reduce("+", lapply(.SD, function(x)
  x == 1L)), .SDcols = as.character(1:212)]
prvl_w[, mm := 0L] # mm = multimorbidity
prvl_w[(count2s == 1L & count1s >= 1L) | (count2s == 0L & count1s >= 2L), 
       mm := 1L]
prvl_w[count2s >= 2L, mm := 2L]

#Adding in the date of incident contions
prvl_w[incd_date, incdate := i.eventdate, on = "patid"]
prvl_w[cmm3plus_incdate, cmm3plus_incdate := i.eventdate, on = "patid"]
rm(incd_date, cmm3plus_incdate)


#Adding in the number of conditions 
prvl_w[, n_cond := count1s + count2s]


# Making the table for a sensitivity analysis with complex = 3+
prvl_w[, cmm3plus := 0L]
prvl_w[(count2s == 2L & count1s >= 1L) | 
         (count2s == 1L & count1s >= 2L) | 
         (count2s == 0L & count1s >= 3L), 
       cmm3plus := 1L]
prvl_w[count2s >= 3L, cmm3plus := 2L]
#check it makes sense: 
prvl_w[, table(cmm3plus, n_cond)]

prvl_w[patient, `:=`(
  age = year - year(i.dob),
  gender = i.gender, 
  region = i.region,
  black = i.black,
  imd = i.imd,
  reg1yr = i.reg1yr,
  censordate = i.censordate,
  censorreason = i.censorreason
), on = "patid"]


# Consistency check
temp <- prvl_w[year == 2019 & mm == 0L, patid]
obstab[patid %in% temp & year_eventdate <= 2019, .N, by = patid][N > 1]
temp <- prvl_w[year == 2004 & mm == 0L, patid]
obstab[patid %in% temp & year_eventdate <= 2004, .N, by = patid][N > 1]
temp <- prvl_w[year == 2004 & mm == 1L, patid]
obstab[patid %in% temp & year_eventdate <= 2004, .N, by = patid][, table(N)]
temp <- prvl_w[year == 2019 & mm == 1L, patid]
obstab[patid %in% temp & year_eventdate <= 2019, .N, by = patid][, table(N)]
temp <- prvl_w[year == 2004 & mm == 2L, patid]
obstab[patid %in% temp & year_eventdate <= 2004, .N, by = patid][, table(N)]
temp <- prvl_w[year == 2019 & mm == 2L, patid]
obstab[patid %in% temp & year_eventdate <= 2019, .N, by = patid][, table(N)]
rm(temp)


agegroup5yfn(prvl_w) 

#Store it 
write_fst(prvl_w, data_dir_CPRD("prvl_w_14Dec.fst"), 100) 

#now saved, getting rid of disease numbers
prvl_w <- prvl_w[, .(patid, year, gender, age, agegrp5, imd, region, reg1yr, 
                     censordate, censorreason, incdate, mm, n_cond, cmm3plus,
                     cmm3plus_incdate)]

# check prvl by year
ggplot(prvl_w[, sum(mm == 2L)/.N *100, keyby = c("year", "agegrp5")], 
       aes(year, V1, col = agegrp5)) + geom_line()




# Making the system table for CMM ----
obstab[diseasesum, on = 'disease_num', system_num := i.system_num]  
syst_obstab <- obstab[, -c("disease_num")]
rm(obstab, diseasesum)
setkey(syst_obstab, patid, eventdate)
syst_obstab <- syst_obstab[syst_obstab[, .I[1], #takes the nth line 
                                       keyby = c("patid", "system_num")]$V1]  
setkey(syst_obstab, patid, eventdate)
syst_incdate <- syst_obstab[syst_obstab[, .I[3], keyby = c("patid")]$V1]
syst_incdate <- syst_incdate[!is.na(patid)]

#Inc CMM ----
cmm_incd <- tt[syst_obstab, on = c("patid", "y == year_eventdate")
][, y := NULL] 
cmm_incd <- cmm_incd[!is.na(year)]
setkey(cmm_incd, patid, year)
cmm_incd_w <- dcast(cmm_incd, patid + year ~ system_num, 
                    fun.aggregate = length)
# View(head(cmm_incd, 100))
# View(head(cmm_incd_w, 100))
rm(cmm_incd)
# some disease incd appears multiple times in a year
cmm_incd_w[, lapply(.SD, max), .SDcols = as.character(1:15)]
for (j in as.character(1:15)) {
  # all incd denoted with 1L
  set(cmm_incd_w, NULL, j, fifelse(cmm_incd_w[[j]] > 0L, 1L, 0L))
}
rm(j)

# CMM prevalence ----
cmm_prvl <- tt[syst_obstab, on = c("patid", "y >= year_eventdate")
][, y := NULL]
cmm_prvl <- cmm_prvl[!is.na(year)]
setkey(cmm_prvl, patid, year)
# nrow(cmm_prvl) # 11603575
# cmm_prvl does not contain completely healthy person-years. Let's add them
cmm_prvl <- cmm_prvl[tt, on = .(patid, year)]

# nrow(cmm_prvl) # 14348373
# Let's create system_num = 0L to denote healthy
cmm_prvl[is.na(system_num), system_num := 0L] 

cmm_prvl_w <- dcast(cmm_prvl, patid + year ~ system_num, fun.aggregate = length)
# View(head(cmm_prvl, 100))
# View(head(cmm_prvl_w, 100))
rm(cmm_prvl)
cmm_prvl_w[, lapply(.SD, max), .SDcols = as.character(1:15)]
cmm_prvl_w[, `0` := 0L]
for (j in as.character(1:15)) { # all cmm_prvl denoted with 2L
  set(cmm_prvl_w, NULL, j, fifelse(cmm_prvl_w[[j]] > 0L, 2L, 0L))
}
rm(j)
setnames(cmm_incd_w, as.character(1:15), paste0("i_", 1:15))
CKutils::absorb_dt(cmm_prvl_w, cmm_incd_w, on = c("patid", "year"))
#rm(cmm_incd_w, tt, pats)


for (j in as.character(1:15)) {
  i_j <- paste0("i_", j)
  absorb_incd(cmm_prvl_w[[j]], cmm_prvl_w[[i_j]])
}
cmm_prvl_w[, paste0("i_", 1:15) := NULL]
rm(cmm_incd_w, pats, tt, j, i_j)


# Housekeeping
setkey(cmm_prvl_w, patid, year) # IMPORTAND for C++ functions

cmm_prvl_w[, pid := rleid(patid)]
cmm_prvl_w[, pid := mk_new_simulant_markers(pid)] # aux col

# Ensure that for each patid an 1 is always followed by 2 and a 2 is 
#always followed by 2
# Note fill_gaps alter x by reference so no assignment is necessary
#c(0L, rep(1L, 10))
#fill_gaps(c(0L, rep(1L, 10)), c(T, rep(F, 10)))
#c(0L, 2, rep(0L, 10))
#fill_gaps(c(0L, 2, rep(1L, 10)), c(T, rep(F, 11)))
invisible(cmm_prvl_w[, lapply(.SD, fill_gaps, pid), 
                     .SDcols = as.character(1:15)])

cmm_prvl_w[, pid := NULL]

# in principle:
# 1. incidence multimorbidity is when a disease row has
#   a. exactly one 2 and at least one 1
#   b. no 2 and at least two 1
# 2. prevalent multimorbidity is when a disease row has at least 2 twos
cmm_prvl_w[, count2s := Reduce("+", lapply(.SD, function(x) x == 2L)), 
           .SDcols = as.character(1:15)]
cmm_prvl_w[, count1s := Reduce("+", lapply(.SD, function(x) x == 1L)), 
           .SDcols = as.character(1:15)]
cmm_prvl_w[, cmm := 0L] # cmm = complex multimorbidity
cmm_prvl_w[(count2s == 2L & count1s >= 1L) | 
             (count2s == 1L & count1s >= 2L) | 
             (count2s == 0L & count1s >= 3L), cmm := 1L]
cmm_prvl_w[count2s >= 3L, cmm := 2L]

cmm_prvl_w[syst_incdate, cmm_incdate := i.eventdate, on = "patid"]
rm(syst_incdate)

# Consistency check
temp <- cmm_prvl_w[year == 2019 & cmm == 0L, patid]
syst_obstab[patid %in% temp & year(eventdate) <= 2019, .N, by = patid][N > 2]
temp <- cmm_prvl_w[year == 2004 & cmm == 0L, patid]
syst_obstab[patid %in% temp & year(eventdate) <= 2004, .N, by = patid][N > 2]
temp <- cmm_prvl_w[year == 2004 & cmm == 1L, patid]
syst_obstab[patid %in% temp & year(eventdate) <= 2004, .N, by = patid][
  , table(N)]
temp <- cmm_prvl_w[year == 2019 & cmm == 1L, patid]
syst_obstab[patid %in% temp & year(eventdate) <= 2019, .N, by = patid][
  , table(N)]
temp <- cmm_prvl_w[year == 2004 & cmm == 2L, patid]
syst_obstab[patid %in% temp & year(eventdate) <= 2004, .N, by = patid][
  , table(N)]
temp <- cmm_prvl_w[year == 2019 & cmm == 2L, patid]
syst_obstab[patid %in% temp & year(eventdate) <= 2019, .N, by = patid][
  , table(N)]
rm(syst_obstab, temp)


cmm_prvl_w[patient, `:=`(
  age = year - year(i.dob),
  gender = i.gender, 
  region = i.region,
  black = i.black,
  imd = i.imd,
  reg1yr = i.reg1yr,
  censordate = i.censordate, 
  censorreason = i.censorreason
), on = "patid"]
rm(patient)
agegroup5yfn(cmm_prvl_w) 

#temporarily writing it here
write_fst(cmm_prvl_w, data_dir_CPRD("cmm_prvl_w_14Dec.fst"), 100) 

#now saved, getting rid of disease numbers
cmm_prvl_w <- cmm_prvl_w[, 
                         .(patid, year, gender, age, agegrp5, imd, region, 
                           reg1yr, censordate,  censorreason, cmm_incdate, cmm)]



# combining bmm & cmm ----
combi_mm <- prvl_w[cmm_prvl_w, on = c("patid", "year"), 
                   `:=` (cmm = i.cmm, cmm_incdate = i.cmm_incdate)]
rm(prvl_w, cmm_prvl_w)

#Adding in a 10yr age group for graphs 
agegroup10yfn <- function(x) { 
  lev <- c("15-19 years", "20-29 years",
           "30-39 years", "40-49 years", 
           "50-59 years", "60-69 years", 
           "70-79 years", "80-89 years",
           "90plus years")
  x[, agegrp10 := fcase(
    (age) >= 90L,           factor(lev[9], levels = lev),
    between(age, 80, 89), factor(lev[8], levels = lev),
    between(age, 70, 79), factor(lev[7], levels = lev),
    between(age, 60, 69), factor(lev[6], levels = lev),
    between(age, 50, 59), factor(lev[5], levels = lev),
    between(age, 40, 49), factor(lev[4], levels = lev),
    between(age, 30, 39), factor(lev[3], levels = lev),
    between(age, 20, 29), factor(lev[2], levels = lev),
    between(age, 15, 19), factor(lev[1], levels = lev)
  )
  ]
}
agegroup10yfn(combi_mm)


#Adding in a 10yr age group for graphs - grouping under 30s & over 80s
combi_mm[, agegrp10_simple := agegrp10]
combi_mm[agegrp10 %in% c("15-19 years", "20-29 years"), 
         agegrp10_simple := "18-29 years"]
combi_mm[agegrp10 %in% c("80-89 years", "90plus years"), 
         agegrp10_simple := "80plus years"]
combi_mm[, agegrp10_simple := factor(
  agegrp10_simple,
  levels =  c(
    "18-29 years",
    "30-39 years",
    "40-49 years",
    "50-59 years",
    "60-69 years",
    "70-79 years",
    "80plus years"
  )
)]


#If incident MM is in 1st year in study, changing to prev if before reg1yr date
combi_mm[incdate < reg1yr & mm == 1, mm := 2]
combi_mm[cmm_incdate < reg1yr & cmm == 1, cmm := 2]
combi_mm[cmm3plus_incdate  < reg1yr & cmm3plus  == 1, cmm3plus  := 2]

write_fst(combi_mm, data_dir_CPRD("combi_mm_14Dec.fst"), 100)


#Patient file for descriptives needs to have only the patients 
#in the study period
patient <- read_fst(data_dir_CPRD("patient.fst"), 
                    columns = c("patid", "gender", "imd", 
                                "black", "yob", "dob", "reg1yr", "censordate", 
                                "censorreason", "region"), 
                    as.data.table = TRUE) 
study_pats <- patient[patid %in% combi_mm[, patid], ]
write_fst(study_pats, data_dir_CPRD("study_pats"), 100)

# time in study ----
study_pats[,  enterstudy := as.IDate(
  ifelse(year(reg1yr) < 2004, 
         as.IDate(paste0("20040101"), format = "%Y%m%d"), 
         reg1yr), format = "%Y%m%d")]
study_pats[enterstudy < censordate, (sum(censordate - enterstudy ))/365.24]

#Median IQR follow-uo
study_pats[enterstudy < censordate, summary((censordate - enterstudy)/365.24)]

#Total by gender
study_pats[, .N]
study_pats[, .(.N, .N/991243), by = gender ]        



# creating a table of study pop for weighting 
cprd_pop_tab <- combi_mm[gender != "I" & imd != "", 
                         .N, keyby = c("year", "gender", "imd", "agegrp5")]


cprd_pop_tab[, agegrp5 := factor(
  agegrp5,
  labels =  c(
    "15-19",
    "20-24",
    "25-29",
    "30-34",
    "35-39",
    "40-44",
    "45-49",
    "50-54",
    "55-59",
    "60-64",
    "65-69",
    "70-74",
    "75-79",
    "80-84",
    "85-89",
    "90+"
  )
)]
cprd_pop_tab[, gender := factor(gender, labels =  c("Men",
                                                    "Women"))]

combi_mm[, gender := factor(gender, 
                            labels =  c("Men", "Women", "Indeterminate"))]


#Downloading the ONS data on pop by age, sex, IMD quintile & year 
data_dir_lookup <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Dictionaries etc/", x)
ons_pop_tab <-  fread(data_dir_lookup("ons_pop_year_sex_age_imd.csv"), 
                      header = TRUE, sep = ",")

ons_pop_tab[, imd:= factor(imd)]
cprd_pop_tab[ons_pop_tab, on = c("year", "gender", "imd", "agegrp5"), 
             ons_pop := i.ons_pop]
cprd_pop_tab[, wt := ons_pop/N]
cprd_pop_tab[, agegrp5 := paste0(agegrp5, " years")]
cprd_pop_tab[agegrp5 == "90+ years", agegrp5 := "90plus years"  ]

combi_mm_wt <- combi_mm[gender != "Indeterminate" & imd != "",]
combi_mm_wt[cprd_pop_tab, 
            on = c("year", "gender", "imd", "agegrp5"), 
            wt := i.wt]
combi_mm_wt[, wt := wt*.N/sum(wt, na.rm = TRUE), by = year]

#checking that the weights add up to the number of people in each year. 
combi_mm_wt[, .(.N, sum(wt)), by = year]
write_fst(combi_mm_wt, data_dir_CPRD("combi_mm_weighted.fst"), 100)



### Making a table for all the graphs ----
cprd_reg <- fread(file = data_dir_lookup("Region.txt"), 
                  header = TRUE, sep = "\t") 

setnames(combi_mm, "region", "regionid")
cprd_reg[, regionid := as.factor(regionid)]
combi_mm[cprd_reg, on = 'regionid', region := i.Description]
combi_mm[, region := factor(
  region,
  levels = c(
    "London",
    "South West",
    "South Central",
    "South East Coast",
    "West Midlands",
    "East Midlands",
    "East of England",
    "North West",
    "Yorkshire And The Humber",
    "North East"
  ),
  labels = c(
    "London",
    "South West",
    "South Central",
    "South East Coast",
    "West Midlands",
    "East Midlands",
    "East of England",
    "North West",
    "Yorkshire And The Humber",
    "North East"
  )
)]
rm(cprd_reg)

#let's add the weights in for those who have them 
combi_mm[combi_mm_wt, on = c("patid", "year"), wt := i.wt]
rm(combi_mm_wt)

#Adding in stuff for BMM inc
combi_mm[, `:=` (
  t.start_incbmm = #start time at risk: 1st of year or reg1yr date - whichever 
    #later
    #those who turn 18 that year contribute only hlf a year
    as.IDate(ifelse(
      year(reg1yr) < year,
      as.IDate((paste0("0101", year)), format = "%d%m%Y"),
      reg1yr
    ), format = "Y%/m%/d%"),
  t.stop_incbmm = as.IDate(ifelse(
    year(censordate) > year,
    as.IDate((paste0("3112", (
      year
    ))), format = "%d%m%Y"),
    censordate
  ), format = "Y%/m%/d%")
)]
combi_mm[mm == 1, t.stop_incbmm := incdate]
combi_mm[mm == 2, c("t.start_incbmm", "t.stop_incbmm") := NA]
combi_mm[, bmm_daysar := t.stop_incbmm - t.start_incbmm + 1][
  bmm_daysar < 0 | is.na(bmm_daysar), bmm_daysar := 0 ]

#Adding in stuff for CMM3plus inc
combi_mm[, `:=` (
  t.start_inccmm3plus = #start time at risk:
    as.IDate(ifelse(
      year(reg1yr) < year,
      as.IDate((paste0("0101", year)), format = "%d%m%Y"),
      reg1yr
    ), format = "Y%/m%/d%"),
  t.stop_inccmm3plus = as.IDate(ifelse(
    year(censordate) > year,
    as.IDate((paste0("3112", (
      year
    ))), format = "%d%m%Y"),
    censordate
  ), format = "Y%/m%/d%")
)]
combi_mm[cmm3plus == 1, t.stop_inccmm3plus := cmm3plus_incdate]
combi_mm[cmm3plus == 2, c("t.start_inccmm3plus", "t.stop_inccmm3plus") := NA]
combi_mm[, cmm3plus_daysar := t.stop_inccmm3plus - t.start_inccmm3plus + 1][
  cmm3plus_daysar < 0 | is.na(cmm3plus_daysar), cmm3plus_daysar := 0]

#Adding in stuff for CMM inc
combi_mm[, `:=` (
  t.start_inccmm = #start time at risk:
    as.IDate(ifelse(
      year(reg1yr) < year,
      as.IDate((paste0("0101", year)), format = "%d%m%Y"),
      reg1yr
    ), format = "Y%/m%/d%"),
  t.stop_inccmm = as.IDate(ifelse(
    year(censordate) > year,
    as.IDate((paste0("3112", (
      year
    ))), format = "%d%m%Y"),
    censordate
  ), format = "Y%/m%/d%")
)]
combi_mm[cmm == 1, t.stop_inccmm := cmm_incdate]
combi_mm[cmm == 2, c("t.start_inccmm", "t.stop_inccmm") := NA]
combi_mm[, cmm_daysar := t.stop_inccmm - t.start_inccmm + 1][
  cmm_daysar < 0 | is.na(cmm_daysar), cmm_daysar := 0 ]


#Adding stuff for prevalence
combi_mm[, `:=` (
  t.start_prev = #start time at risk: 1st of year or reg1yr
    #date or mm incdate - whichever later
    #those who turn 18 that year contribute only hlf a year
    as.IDate(ifelse(
      year(reg1yr) < year,
      as.IDate((paste0("0101", year)),
               format = "%d%m%Y"),
      reg1yr
    ), format = "Y%/m%/d%"),
  t.stop_prev = as.IDate(ifelse(
    year(censordate) > year,
    as.IDate((paste0("3112", (
      year
    ))), format = "%d%m%Y"),
    censordate
  ), format = "Y%/m%/d%")
)]
combi_mm[, totaldays := t.stop_prev - t.start_prev + 1] #denom for each person
#creating bmm numerator for each person:
combi_mm[, bmm_days := ifelse(mm == 2,
                              t.stop_prev - t.start_prev + 1,
                              ifelse(mm == 1, t.stop_prev - incdate, 0))]

#creating3+ cmm numerator for each person
combi_mm[, cmm3plus_days := ifelse(
  cmm3plus == 2,
  t.stop_prev - t.start_prev + 1,
  ifelse(cmm3plus == 1, t.stop_prev - cmm3plus_incdate + 1, 0)
)]

#creating cmm numerator for each person
combi_mm[, cmm_days := ifelse(cmm == 2,
                              t.stop_prev - t.start_prev + 1,
                              ifelse(cmm == 1, t.stop_prev - cmm_incdate + 1, 
                                     0))]

setnames(combi_mm, "mm", "bmm")


# Adding in stuff for point prev on 1 July year Y
#Denominator
combi_mm[, pointprev_dnom := ifelse(
  t.start_prev <= as.IDate(paste0("0107", year), format = "%d%m%Y") & 
    t.stop_prev >= as.IDate(paste0("0107", year), format = "%d%m%Y"), 1,0)]
#BMM numerator 
combi_mm[, pointprev_bmm := ifelse(
  pointprev_dnom == 1 & 
    incdate <= as.IDate(paste0("0107", year), format = "%d%m%Y"), 1, 0)]
combi_mm[is.na(pointprev_bmm), pointprev_bmm := 0]
#CMM numerator 
combi_mm[, pointprev_cmm := ifelse(
  pointprev_dnom == 1 & 
    cmm_incdate <= as.IDate(paste0("0107", year), format = "%d%m%Y"), 1, 0)]
combi_mm[is.na(pointprev_cmm), pointprev_cmm := 0]

#keeping only 1 censored code per person
combi_mm[censorreason == 1 & year(censordate) != year, censorreason := "0"]

combi_mm <- combi_mm[, .(patid, year, gender, age, imd, regionid, region, 
                         censorreason, agegrp5, agegrp10_simple, wt,
                         bmm, cmm, cmm3plus, n_cond,
                         totaldays, bmm_days, cmm_days, cmm3plus_days,
                         bmm_daysar, incdate, 
                         cmm_daysar, cmm_incdate,
                         cmm3plus_daysar, cmm3plus_incdate, 
                         pointprev_dnom, pointprev_bmm, pointprev_cmm)]


combi_mm[, agegrp10_simple := factor(agegrp10_simple,
                                     labels =  c("18-29",
                                                 "30-39", "40-49",
                                                 "50-59", "60-69",
                                                 "70-79", "80+"))]
combi_mm[, gender := factor(gender, labels =  c("Men",
                                                "Women", "Indeterminate"))]

write_fst(combi_mm, data_dir_CPRD("combi_mm_detailed.fst"), 100)
