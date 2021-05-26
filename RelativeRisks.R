#Relative risks 

library(data.table)
library(fst)
library(glm2) # better algo to fit glm, same syntax
library(sandwich) # robust CI
library(MASS)
library(gamlss)

# In the data, the relevant variables are:
# -	bmm/cmm: 0 = no, 1 = incident, 2 = prevalent (in that year)
# -	totaldays: total days registered in the year
# -	bmm_days/cmm_days: days in the year with bmm/cmm
data_dir_CPRD <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Data May 2020/", x)

data_dir_lookup <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Dictionaries etc/", x)


mydt <- read_fst(data_dir_CPRD("combi_mm_detailed.fst"), as.data.table = TRUE)
names(mydt)

mydt <- mydt[imd != "" & gender != "Indeterminate" & region != ""]

# get robust CI using sandwich (just testing, do not use for publication)
# Note that while family = quasibinomial leads to a robust estimate of the 
#variance (αp^(1−p^)), it is not the same robust estimator as the sandwich 
#estimator (in which for each observation, the variance is estimated as 
#(yi−y^i)2 rather than by some assumed functional form)
robustCI <- function(model, type = "HC0", digits = 3) {
  cov_m1 <- sandwich::vcovHC(model, type = type)
  
  std_err <- sqrt(diag(cov_m1))
  
  q_val <- qnorm(0.975)
  hlp <- coef(model)/std_err
  r_est <- cbind(
    Estimate = coef(model)
    , "Robust SE" = std_err
    , z = (hlp)
    , "Pr(>|z|) " = 2 * pnorm(abs(hlp), lower.tail = FALSE)
    , LL = coef(model) - q_val  * std_err
    , UL = coef(model) + q_val  * std_err
    , RR = exp(coef(model))
    , LRR = exp(coef(model) - q_val  * std_err)
    , URR = exp(coef(model) + q_val  * std_err)
  )
  
  round(r_est, digits = digits)
}

# Incident bmm ----
tt <- mydt[bmm != "2" & imd != "" & totaldays > 0]
tt[, bmm := as.integer(as.character(bmm))]
tt[, hist(bmm/totaldays)]
tt[, mean(bmm/totaldays), keyby = .(year, agegrp10_simple)]
tt[, mean(bmm/totaldays), keyby = year][, plot(year, V1)]
tt[, mean(bmm/totaldays), keyby = age][, plot(age, V1)]

start_p <- sum(tt$bmm)/sum(tt$totaldays)
bmm_inc <-
  glm2(
    cbind(bmm, totaldays - bmm_days) ~ year + gender + region + 
      imd + agegrp10_simple,
    family = binomial(link = "log"),
    data = tt,
    start = c(log(start_p), rep(0, 22)), control = glm.control(trace = TRUE)
  )
summary(bmm_inc)
bmm_inc_rr_Wald <- format(round(
  cbind(exp(coef(bmm_inc)), exp(confint.default(bmm_inc))), 3), 
  scientific = FALSE) # Wald confidence interval:
#https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

bmm_inc_rr_robust <- robustCI(bmm_inc)
bmm_inc_rr_robust # use this one LRR & URR are lower and upper 95% CI for RR



stt <- na.omit(tt[, .(bmm = sum(bmm), 
                      totaldays = sum(totaldays - bmm_days + bmm),
                      bmm_days = sum(bmm_days)), 
                  by = .(year, gender, region, imd, agegrp10_simple)])
stt[, hist(bmm)]

bmm_inc2 <-
  glm2(
    cbind(bmm, totaldays - bmm) ~ 
      year + gender + region + imd + agegrp10_simple,
    family = binomial(link = "log"),
    data = stt,
    start = c(log(start_p), rep(0, 22)),
    control = glm.control(trace = TRUE)
  )
View(cbind(bmm_inc_rr_robust[, 7:9],
           robustCI(bmm_inc2)[, 7:9],
           round(cbind(exp(coef(bmm_inc2)),
                       exp(confint(bmm_inc2))), 3)))

# profile confidence intervals are very similar to robust ones
# Given the later method is much faster + as accurate will use it for all models
bmm_inc_rr2 <- round(cbind(exp(coef(bmm_inc2)), exp(confint(bmm_inc2))), 3)

bmm_inc3 <-
  glm.nb(
    bmm ~ offset(log(totaldays)) + 
      year + gender + region + imd + agegrp10_simple,
    data = stt
  )
bmm_inc_rr3 <- round(cbind(exp(coef(bmm_inc3)), exp(confint(bmm_inc3))), 3)
View(cbind(bmm_inc_rr2, bmm_inc_rr3, robustCI(bmm_inc3)[, 7:9]))

bmm_inc4 <-
  glm2(
    bmm ~ offset(log(totaldays)) + 
      year + gender + region + imd + agegrp10_simple,
    family = poisson(link = "log"),
    data = stt
  )
bmm_inc_rr4 <- round(cbind(exp(coef(bmm_inc4)), exp(confint(bmm_inc4))), 3)
View(cbind(bmm_inc_rr2, bmm_inc_rr4, robustCI(bmm_inc4)[, 7:9]))

bmm_inc5 <-
  glm2(
    bmm ~ offset(log(totaldays)) + 
      year + gender + region + imd + agegrp10_simple,
    family = quasipoisson(link = "log"),
    data = stt
  )
bmm_inc_rr5 <- round(cbind(exp(coef(bmm_inc5)), exp(confint(bmm_inc5))), 3)
View(cbind(bmm_inc_rr2, bmm_inc_rr5, robustCI(bmm_inc5)[, 7:9]))

#Comparing the nb & the qp
View(cbind(bmm_inc_rr3, bmm_inc_rr5, robustCI(bmm_inc5)[, 7:9]))
plot(bmm_inc_rr3)
plot(bmm_inc_rr5)
#Comparing the nb & the qp
AER::dispersiontest(bmm_inc4, trafo = 1) # qp a > 0, p < 0.05 hence 
#overdispersion
AER::dispersiontest(bmm_inc4, trafo = 2) # nb a > 0, p < 0.05 hence 
#overdispersion but (less)

# Hence, an important diagnostic is to plot (Yi - Î¼i)^2 against Î¼i.
# Function for variance-to-mean plot
# The graph plots the mean versus the variance and overlays the curves
# corresponding to the over-dispersed quasi-Poisson model, where the 
#variance is ÏÎ¼,
# and the negative binomial model, where the variance is Î¼(1+Î¼Ï2).

var2mean <- function(nb_model, qp_model, subtitle) {
  xb <- log(nb_model$fitted.values) # same as predict(nb_model). 
  #At linear predictor scale
  g <- cut(xb, breaks = quantile(xb, seq(0, 100, 5) / 100))
  m <- tapply(nb_model$y, g, mean)
  v <- tapply(nb_model$y, g, var)
  size <- tapply(nb_model$y, g, length)
  range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  
  mat <- matrix(1:2, 2)
  layout(mat, widths = c(3.5, 3.5), heights = c(8,6))
  
  plot(m,
       v,
       cex = range01(size),
       xlab = "Mean",
       ylab = "Variance",
       main = "Mean-Variance Relationship")
  
  mtext(subtitle, padj = -0.5)
  x <- seq(min(m), max(m), 0.02)
  phi <- summary(qp_model)$dispersion
  lines(x, x * phi, lty = "dashed")
  lines(x, x * (1 + x / nb_model$theta))
  legend(
    "topleft",
    lty = c("dashed", "solid"),
    legend = c("Q. Poisson", "Neg. Binom."),
    inset = 0.05
  )
  
  # weights plot
  w_nb <- tapply(nb_model$weights, g, mean)
  w_qp <- tapply(qp_model$weights, g, mean)
  
  plot(m, w_nb, type = "l",
       ylim = c(0, max(c(w_nb, w_qp))),
       xlab = "Mean",
       ylab = "weights",
       main = "Estimated weights as a function of the mean"
  )
  
  lines(m, w_qp, lty = "dashed")
}

var2mean(bmm_inc3, bmm_inc5, "BMM incidence")




# Prevalent bmm ----
tt <- mydt[imd != "" & totaldays > 0] 
tt[, bmm := as.integer(as.character(bmm))]
tt[bmm == 0L, bmm2 := 0L]
tt[bmm == 1L, bmm2 := bmm_days]
tt[bmm == 2L, bmm2 := totaldays]

tt[, hist(bmm2/totaldays)] 
tt[, mean(bmm2/totaldays), keyby = .(year, agegrp10_simple)]
tt[, mean(bmm2/totaldays), keyby = year][, plot(year, V1)]
tt[, mean(bmm2/totaldays), keyby = age][, plot(age, V1)]

stt <- tt[, .(bmm = sum(bmm2), totaldays = sum(totaldays),
              bmm_days = sum(bmm_days)), 
          by = .(year, gender, region, imd, agegrp10_simple)]

stt[, bmm_mean := mean(bmm, na.rm = T),
    keyby = .(year, gender, region, imd, agegrp10_simple)]
stt[, bmm_var := var(bmm, na.rm = T, use = "pairwise.complete.obs"), 
    keyby = .(year, gender, region, imd, agegrp10_simple)]
stt[, var(bmm), keyby = .(year, gender, imd, agegrp10_simple)]

tmptab <- stt[, .(bmm_mean = mean(bmm), bmm_var = var(bmm)), 
              keyby = .(year, gender, imd, agegrp10_simple)]

bmm_prv <-
  glm2(
    bmm ~ offset(log(totaldays)) + 
      year + gender + region + imd + agegrp10_simple,
    family = quasipoisson(link = "log"),
    data = stt
  )
bmm_prv_rr <- round(cbind(exp(coef(bmm_prv)), exp(confint(bmm_prv))), 3)
View(cbind(bmm_prv_rr, robustCI(bmm_prv)[, 7:9]))

dispersiontest(bmm_prvpois, trafo = 1)

bmm_prvpois <-
  glm2(
    bmm ~ offset(log(totaldays)) + 
      year + gender + region + imd + agegrp10_simple,
    family = poisson(link = "log"),
    data = stt
  )
#comparing with binomial
bmm_prv_nb <-
  glm.nb(
    bmm ~ offset(log(totaldays)) + 
      year + gender + region + imd + agegrp10_simple,
    data = stt
  )
bmm_prvnb_rr <- round(cbind(exp(coef(bmm_prv_nb)), 
                            exp(confint(bmm_prv_nb))), 3)
View(cbind(bmm_prv_rr, bmm_prvnb_rr, robustCI(bmm_prv_nb)[, 7:9]))
plot(bmm_prv)
plot(bmm_prv_nb)

var2mean(bmm_prv_nb, bmm_prv, "BMM prevalence")


# Case fatality bmm --------------------
# assuming censorreason = 1 means death
tt <- mydt[imd != "" & totaldays > 0 & bmm != "0"]
tt[, bmm := as.integer(as.character(bmm))]
tt[, censorreason := as.integer(as.character(censorreason))]

stt <- na.omit(tt[, .(censorreason = sum(censorreason), 
                      totaldays = sum(totaldays),
                      bmm_days = sum(bmm_days)), 
                  by = .(year, gender, region, imd, agegrp10_simple)])
stt[, hist(censorreason)]

bmm_cf <-
  glm2(
    censorreason ~ offset(log(bmm_days)) + 
      year + gender + region + imd + agegrp10_simple,
    family = quasipoisson(link = "log"),
    data = stt
  )
bmm_cf_rr <- round(cbind(exp(coef(bmm_cf)), exp(confint(bmm_cf))), 3)
View(cbind(bmm_cf_rr, robustCI(bmm_cf)[, 7:9]))

d <- fitDist(stt$censorreason, type = "counts",
             try.gamlss = TRUE, trace = TRUE)
d
bmm_cf2 <-
  gamlss(
    censorreason ~ offset(log(bmm_days)) + 
      year + gender + region + imd + agegrp10_simple,
    family = ZANBI(),
    data = stt
  )
bmm_cf_rr2 <- round(cbind(exp(coef(bmm_cf2)), exp(confint(bmm_cf2))), 3)
bmm_cf_rr2

#comparing with binomial
bmm_cf_nb <-
  glm.nb(
    censorreason ~ offset(log(bmm_days)) + 
      year + gender + region + imd + agegrp10_simple,
    data = stt
  )
bmm_cfnb_rr <- round(cbind(exp(coef(bmm_cf_nb)), exp(confint(bmm_cf_nb))), 3)
View(cbind(bmm_cf_rr, bmm_cfnb_rr, robustCI(bmm_cf_nb)[, 7:9]))
plot(bmm_cf)
plot(bmm_cf_nb)

AIC(bmm_cf_nb, bmm_cf2) 
var2mean(bmm_cf_nb, bmm_cf, "BMM case-fatality")


# Incident cmm ----------------------------
tt <- mydt[cmm != "2" & imd != "" & totaldays > 0] 
tt[, cmm := as.integer(as.character(cmm))]
original <- copy(tt)
tt[, hist(cmm/totaldays)] 
tt[, mean(cmm/totaldays), keyby = .(year, agegrp10_simple)][, summary(V1)]
tt[, mean(cmm/totaldays), keyby = year][, plot(year, V1)]
tt[, mean(cmm/totaldays), keyby = age][, plot(age, V1)]

stt <- tt[, .(cmm = sum(cmm), totaldays = sum(totaldays - cmm_days + cmm),
              cmm_days = sum(cmm_days)), 
          by = .(year, gender, region, imd, agegrp10_simple)]

cmm_inc <-
  glm2(
    cmm ~ offset(log(totaldays)) + 
      year + gender + region + imd + agegrp10_simple,
    family = quasipoisson(link = "log"),
    data = stt
  )
cmm_inc_rr <- round(cbind(exp(coef(cmm_inc)), exp(confint(cmm_inc))), 3)
View(cbind(cmm_inc_rr, robustCI(cmm_inc)[, 7:9]))


#comparing with nb
cmm_inc_nb <-
  glm.nb(
    cmm ~ offset(log(totaldays)) + 
      year + gender + region + imd + agegrp10_simple,
    data = stt
  )
cmm_inc_nb_rr <- round(cbind(exp(coef(cmm_inc_nb)), 
                             exp(confint(cmm_inc_nb))), 3)
View(cbind(cmm_inc_rr, cmm_inc_nb_rr, robustCI(cmm_inc_nb)[, 7:9]))
plot(cmm_inc)
plot(cmm_inc_nb)
var2mean(cmm_inc_nb, cmm_inc, "CMM incidence")


# Prevalent cmm ----------------------------
tt <- mydt[imd != "" & totaldays > 0] 
tt[, cmm := as.integer(as.character(cmm))]
tt[cmm == 0L, cmm2 := 0L]
tt[cmm == 1L, cmm2 := cmm_days]
tt[cmm == 2L, cmm2 := totaldays]

tt[, hist(cmm2/totaldays)] 
tt[, mean(cmm2/totaldays), keyby = .(year, agegrp10_simple)]
tt[, mean(cmm2/totaldays), keyby = year][, plot(year, V1)]
tt[, mean(cmm2/totaldays), keyby = age][, plot(age, V1)]

stt <- tt[, .(cmm = sum(cmm2), totaldays = sum(totaldays),
              cmm_days = sum(cmm_days)), 
          by = .(year, gender, region, imd, agegrp10_simple)]

cmm_prv <-
  glm2(
    cmm ~ offset(log(totaldays)) + 
      year + gender + region + imd + agegrp10_simple,
    family = quasipoisson(link = "log"),
    data = stt
  )
cmm_prv_rr <- round(cbind(exp(coef(cmm_prv)), exp(confint(cmm_prv))), 3)
View(cbind(cmm_prv_rr, robustCI(cmm_prv)[, 7:9]))

#comparing with nb
cmm_prv_nb <-
  glm.nb(
    cmm ~ offset(log(totaldays)) + 
      year + gender + region + imd + agegrp10_simple,
    data = stt
  )
cmm_prv_nb_rr <- round(cbind(exp(coef(cmm_prv_nb)), 
                             exp(confint(cmm_prv_nb))), 3)
View(cbind(cmm_prv_rr, cmm_prv_nb_rr, robustCI(cmm_prv_nb)[, 7:9]))
plot(cmm_prv)
plot(cmm_prv_nb)

var2mean(cmm_prv_nb, cmm_prv, "CMM prevalence")

# Case fatality cmm --------------------
# assuming censorreason = 1 means death
tt <- mydt[imd != "" & totaldays > 0 & cmm != "0"] 
tt[, cmm := as.integer(as.character(cmm))]
tt[, censorreason := as.integer(as.character(censorreason))]

stt <- na.omit(tt[, .(censorreason = sum(censorreason), 
                      totaldays = sum(totaldays),
                      cmm_days = sum(cmm_days)), 
                  by = .(year, gender, region, imd, agegrp10_simple)])
stt[, hist(censorreason)] 

cmm_cf <-
  glm2(
    censorreason ~ offset(log(cmm_days)) + 
      year + gender + region + imd + agegrp10_simple,
    family = quasipoisson(link = "log"),
    data = stt
  )
cmm_cf_rr <- round(cbind(exp(coef(cmm_cf)), exp(confint(cmm_cf))), 3)
View(cbind(cmm_cf_rr, robustCI(cmm_cf)[, 7:9]))

#comparing with binomial
cmm_cf_nb <-
  glm.nb(
    censorreason ~ offset(log(cmm_days)) + 
      year + gender + region + imd + agegrp10_simple,
    data = stt
  )
cmm_cfnb_rr <- round(cbind(exp(coef(cmm_cf_nb)), exp(confint(cmm_cf_nb))), 3)
View(cbind(cmm_cf_rr, cmm_cfnb_rr, robustCI(cmm_cf_nb)[, 7:9]))
plot(cmm_cf)
plot(cmm_cf_nb)

var2mean(cmm_cf_nb, cmm_cf, "CMM case-fatality")

