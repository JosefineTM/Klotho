##########################
# CHS first look
#
##########################

rm(list = ls())

library(tidyverse)
library(haven)
library(gt)
library(gtsummary)

# Import data
source('paths.R')

data_soma <- read.csv(file.path(protein_data_path, 'data_soma_klotho.csv'))
str(data_soma)


# Import other CHS data using Krisitnes code: 
# '/Volumes/auditing-groupdirs/SUN-IFSV-CHS-data/Kristine/CHS Sphingo/Kristine/Cognition and dementia/Code'


# all_data_clean <- read_sas(file.path(data_path_CHS, "all_data_clean.sas7bdat"))

all_data_dem <- read_sas (file.path(data_path_CHS, "all_data_dem.sas7bdat"))
# 
# all_data_cog <- read_sas(file.path(data_path_CHS, "all_data_cogedit2.sas7bdat"))
# 
# all_sl_cog <- read_sas(file.path(data_path_CHS, "all_sl_cog.sas7bdat"))


all_data_dem_out <- all_data_dem %>%
  select(c(idno, year, baseage, gend01, race01, white, clinic, grade01, starts_with('base'), 
           demclass_updated, demtype, incad, incmvad, incmvad, ttoincdem)) %>%
  mutate(
    year = as.factor(year),
    gend01 = case_when(
      gend01 == 0 ~ 'Female', 
      gend01 == 1 ~ 'Male', 
    ), 
    white = case_when(
      white == 0 ~ 'White', 
      white == 1 ~ 'Non-white'
      ), 
    race01 = case_when(
      race01 == 1 ~ 'White', 
      race01 == 2 ~ 'Black', 
      race01 == 3 ~ 'Native american', 
      race01 == 4 ~ 'Asian', 
      race01 == 5 ~ 'Other'
    ), 
    education_group = case_when(
      grade01 <= 11                  ~ "Less than high school",
      grade01 %in% c(12, 13)         ~ "High school",
      grade01 %in% c(14, 15, 16)     ~ "High school",        # vocational
      grade01 >= 17 & grade01 <= 98  ~ "College and above",
      grade01 == 99                  ~ NA_character_,
      TRUE                           ~ NA_character_
    ), 
    basesmoke = case_when(
      basesmoke == 1 ~ 'Never', 
      basesmoke == 2 ~ 'Former', 
      basesmoke == 3 ~ 'Current'
    ), 
    demclass_updated = case_when(
      demclass_updated == 0 ~ 'Normal',
      demclass_updated == 1 ~ 'Incident',
      demclass_updated == 2 ~ 'MCI', 
      demclass_updated == 3 ~ 'Prevalent',
      TRUE ~ 'Unknown'
    ), 
    demtype = case_when(
      demtype == 0 ~ 'Normal', 
      demtype == 1 ~ 'AD',
      demtype == 2 ~ 'Mixed',
      demtype == 3 ~ 'VaD',
      demtype == 9 ~ 'Unknown',
    ), 
    basehlth = case_when(
      basehlth == 1 ~ 'Excellent', 
      basehlth == 2 ~ 'Very good', 
      basehlth == 3 ~ 'Good', 
      basehlth == 4 ~ 'Fair', 
      basehlth == 5 ~ 'Poor', 
    ), 
    # baseapoe = case_when(
    #   baseapoe == 1 ~ "e2/e2",
    #   baseapoe == 2 ~ "e2/e3",
    #   baseapoe == 3 ~ "e2/e4",
    #   baseapoe == 4 ~ "e3/e3",
    #   baseapoe == 5 ~ "e3/e4",
    #   baseapoe == 6 ~ "e4/e4",
    #   TRUE      ~ NA_character_
    # )
    ) %>%
  as.data.frame() %>%
  rename(ID = idno) 

str(all_data_dem_out %>% select(c(baseapoe)))

write_csv(all_data_dem_out, file.path(save_path, 'data_bl_klotho.csv'))


# Make summary table

variables <- all_data_dem_out %>% 
  select(-c(ID, gend01, white, grade01, baseincp, baseinc, basebmio,baseheight, 
            baseweight, basekcal, basealcohc, baseestroc, baseestrocm, basedepscr, 
            basediabada, basechd, basestk, basetrig, baseglu, baseins, basealbp, 
            basealb, baseapoem, basehtnmed, basehtn, basesttn)
         ) %>% 
  colnames()
  
table1 <- all_data_dem_out |>
  tbl_summary(
    by = gend01,
    include = all_of(variables),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(
      year = "Year",
      baseage = 'Age',
      race01 = 'Race',
      `CHS CLINIC` = 'CHS clinic number',
      basemm3 = '3MSE score (baseline)',
      basedsst = 'DST score (baseline)',
      basebmi = 'BMI',
      basewaist = 'Waist circumference (cm)',
      basealcoh = 'Alchol consumption (#/week)',
      basesmoke = 'Smoking status',
      baseestroc = 'Takes estrogen or other female hormones', 
      basehlth = 'Health status', 
      basesysbp = 'SBP',
      basecrp = 'CRP',
      basechol = 'Cholesterol', 
      baseldl = 'LDL', 
      basehdl = 'HDL',
      baseapoe = 'APOE genotype', 
      basemci = 'MCI at baseline', 
      demclass_updated = 'Dementia class', 
      demtype = 'Dementia type',
      incad = 'Incident AD', 
      incmvad = 'Incident VaD', 
      education_group = 'Education'
    ),
    missing_text = "(Missing)"
  ) 

table1

gtsave(as_gt(table1), file.path(save_path, "table1.docx"))
gtsave(as_gt(table1), file.path(save_path, "table1.png"))




# basealcohc = AMT OF ALCOH CONSUMMED PER WK


#Dementia status

#Prevalent all-cause dementia
table(all_data_dem_out$demclass_updated, useNA = "always")
all_data_dem_out$demclass_updated <- as.factor(all_data_dem$demclass_updated)

#Incidence all-cause dementia
table(all_data_dem$incdem, useNA = "always")
typeof(all_data_dem$incdem)
all_data_dem$incdem <- as.factor(all_data_dem$incdem)
plot(all_data_dem$incdem)

#Incidence dementia subtypes
table(all_data_dem$demtype, useNA = "always")
# 0 = normal
# 1 = AD
# 2 = mixed
# 3 = VaD
# 9 = unknown
typeof(all_data_dem$demtype)
all_data_dem$demtype <- as.factor(all_data_dem$demtype)
plot(all_data_dem$demtype)

#Time to incidence dementia from date of year5/year7 visit or first MRI (days)
with(all_data_dem, summary(ttoincdem))
hist(all_data_dem$ttoincdem)











# -----------------------------------------------------

#BMI (year 5 or 7) - (con)
with(all_data_dem, summary(basebmi))
hist(all_data_dem$basebmi)
qqnorm(scale(all_data_dem$basebmi), xlim=c(-4,4), ylim=c(-4,4), main='Baseline BMI')
abline(0,1)
typeof(all_data_dem$basebmi)

#Physical activity (kcal) - (con) !
with(all_data_dem, summary(basekcal))
hist(all_data_dem$basekcal)
qqnorm(scale(all_data_dem$basekcal), xlim=c(-4,4), ylim=c(-4,4), main='Baseline physical activity (kcal/week)')
abline(0,1)
typeof(all_data_dem$basekcal)

#Alcohol consumption (current) (cat)
table(all_data_dem$basealcohc, useNA = "always")
plot(all_data_dem$basealcohc)
typeof(all_data_dem$basealcohc)
all_data_dem$basealcohc <- as.factor(all_data_dem$basealcohc)

#Depression score (con) - !
summary(all_data_dem$basedepscr, useNA = "always")
hist(all_data_dem$basedepscr)
typeof(all_data_dem$basedepscr)

#Prevalent diabetes (cat)
table(all_data_dem$basediabada, useNA = "always")
plot(all_data_dem$basediabada)
typeof(all_data_dem$basediabada)
all_data_dem$basediabada <- as.factor(all_data_dem$basediabada)

#Prevalent Coronary heart disease (cat)
table(all_data_dem$basechd, useNA = "always")
plot(all_data_dem$basechd)
typeof(all_data_dem$basechd)
all_data_dem$basechd <- as.factor(all_data_dem$basechd)

#Prevalent stroke (cat)
table(all_data_dem$basestk, useNA = "always")
plot(all_data_dem$basestk)
typeof(all_data_dem$basestk)
all_data_dem$basestk <- as.factor(all_data_dem$basestk)

#Hypertension med-user (cat)
table(all_data_dem$basehtnmed, useNA = "always")
plot(all_data_dem$basehtnmed)
typeof(all_data_dem$basehtn)
all_data_dem$basehtnmed <- as.factor(all_data_dem$basehtnmed)

#Lipid-lowering med-user (cat)
table(all_data_dem$basesttn, useNA = "always")
plot(all_data_dem$basesttn)
typeof(all_data_dem$basesttn)
all_data_dem$basesttn <- as.factor(all_data_dem$basesttn)

#Baseline estrogen therapy (cat)
table(all_data_dem$baseestroc, useNA = "always")
plot(all_data_dem$baseestroc)
typeof(all_data_dem$baseestroc)
all_data_dem$baseestroc <- as.factor(all_data_dem$baseestroc)

#APOE carrier (cat)
table(all_data_dem$baseapoe, useNA = "always")
plot(all_data_dem$baseapoe)
typeof(all_data_dem$baseapoe)
all_data_dem$baseapoe <- as.factor(all_data_dem$baseapoe)

#HDL - 10 NA!
with(all_data_dem, summary(basehdl))
hist(all_data_dem$basehdl)
qqnorm(scale(all_data_dem$basehdl), xlim=c(-4,4), ylim=c(-4,4), main='Baseline HDL levels')
abline(0,1)
typeof(all_data_dem$basehdl)

#LDL - 25 NA!
with(all_data_dem, summary(baseldl))
hist(all_data_dem$baseldl)
qqnorm(scale(all_data_dem$baseldl), xlim=c(-4,4), ylim=c(-4,4), main='Baseline LDL levels')
abline(0,1)
typeof(all_data_dem$baseldl)

#Baseline MCI status (MM < 88 = 1) - 3 NA
table(all_data_dem$basemci, useNA = "always")
plot(all_data_dem$basemci)
typeof(all_data_dem$basemci)
all_data_dem$basemci <- as.factor(all_data_dem$basemci)

#Baseline income
table(all_data_dem$baseinc, useNA = "always")
plot(all_data_dem$baseinc)
typeof(all_data_dem$baseinc)
all_data_dem$baseinc <- as.factor(all_data_dem$baseinc)

