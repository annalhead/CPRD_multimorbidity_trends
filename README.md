# CPRD_multimorbidity_trends

This repository shares the code for the analysis and figures for our study of trends and inequalities in multimorbidity between 2004 - 2019 in England, using CPRD Aurum data. The results of our analysis are published in The Lancet Healthy Longevity as _Inequalities in incident and prevalent multimorbidity in England, 2004–19: a population-based, descriptive study_ [https://doi.org/10.1016/S2666-7568(21)00146-X](https://doi.org/10.1016/S2666-7568(21)00146-X)
 
* Please note - this code is for the work within my UoLiverpool/SPHR PhD studentship (as is all other code on this github). This is not the phenotyping or analysis code used in other projects I am involved in (e.g. work with The Health Foundation or the IMPACTncd microsimulation model), and does not include implementation of the Cambridge Multimorbidity Score. I am working on making this code available, but in the meantime please email me (ahead@liverpool.ac.uk) with any questions about code for other projects.

We have used 2 definitions of multimorbidity:
Basic Multimorbidity - 2 or more conditions (from a list of 211, details here: https://github.com/annalhead/CPRD_multimorbidity_codelists)
Complex Multimorbidity - 3 or more conditions (as above) affecting at least 3 body systems (as defined by Harrison et al. https://bmjopen.bmj.com/lookup/doi/10.1136/bmjopen-2013-004694)
We have also included a sensitivity analysis where complex multimorbidity is defined as 3 or more of any condition (regardless of body system) 

The files in this repository:
1. Datacleaning_neat.R - this reads in the raw CPRD data and extracts the observations of interest based on the codelists stored in https://github.com/annalhead/CPRD_multimorbidity_codelists 
2. MMfulltable.R - this takes the output of datacleaning_neat.R and makes a table for the analysis and figures: each individual has 1 row for every year they are included in the study population. Files 3- all use the output from this file
3. DescriptiveGraphs.Rmd - study descriptives
4. Prevalencegraphs.Rmd - crude & standardised prevalence + measures of inequality (SII & RII)
5. Incidencegraphs.Rmd - crude & standardised incidence + measures of inequality (SII & RII)
6. CFgraphs.Rmd - crude & standardised case fatality + measures of inequality (SII & RII)
7. AgeOnsetGraph.Rmd - median age of onset 
8. RelativeRisks.R - relative risks (quassi-poisson regression) for incidence, prevalence and case fatality.

Files for sensitivity analyses are still to be added. 
