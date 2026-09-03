##########################
# CHS first look
#
##########################

rm(list = ls())

library(tidyverse)
library(haven)
library(gt)
library(gtsummary)
library(glue)

# Import data
source('paths.R')

data_soma <- read.csv(file.path(protein_data_path, 'data_soma_klotho.csv'))
str(data_soma)

data_clin <- read.csv(file.path(protein_data_path, 'data_bl_klotho.csv'))
str(data_clin)

save_path = file.path(proj_path, 'analysis', '00_first_look')


# --------------------#
#                     #
#         SEX         #
#                     #
# --------------------#

n_missing = data_soma %>%left_join(data_clin %>% select(ID, gend01), by = "ID") %>% filter(is.na(gend01)) %>% nrow

# Table 1 with T-statistic
table1 <- data_soma %>%
  left_join(data_clin %>% select(ID, gend01), by = "ID") %>%
  tbl_summary(
    by = gend01,
    include = data_soma %>%
      select(where(is.numeric), -c(ID)) %>%
      colnames(),
    digits = all_continuous() ~ 2
  ) %>%
  add_p(
    test = all_continuous() ~ "t.test",
    pvalue_fun = ~style_pvalue(.x)
  ) %>%
  as_gt() %>%
  tab_footnote(
    footnote = glue("Missing values (N = {n_missing}) were excluded from the analysis"),
    locations = cells_title(groups = "title")  # 👈 attach somewhere
  )
  

table1

gtsave(table1 , file.path(save_path, 'klotho_sex_strat.png'))
gtsave(table1, file.path(save_path, 'klotho_sex_strat.docx'))


data_soma %>%
  pivot_longer(cols = -c(ID),
               names_to = 'Protein', 
               values_to = 'Intensity') %>%
  left_join(data_clin %>% select(ID, gend01), by = "ID") %>%
  ggplot(aes(x = Protein, y = Intensity)) + 
  geom_boxplot(aes(colour = gend01))


data_soma %>%
  pivot_longer(cols = -c(ID),
               names_to = 'Protein', 
               values_to = 'Intensity') %>%
  left_join(data_clin %>% select(ID, baseage), by = "ID") %>%
  mutate(age_groups = case_when(
    baseage <= 69 ~ '65-69', 
    baseage >= 70 & baseage < 75 ~ '70-74', 
    baseage >= 75 & baseage < 80 ~ '75-79', 
    baseage >= 80 & baseage < 85 ~ '80-84', 
    baseage >= 85 & baseage < 90 ~ '85-89', 
    baseage >= 90 & baseage < 90 ~ '90-94', 
    baseage >= 95 ~ '>94'
    
  )) %>%
  ggplot(aes(x = Protein, y = Intensity)) + 
  geom_boxplot(aes(colour = age_groups))



library(ggplot2)
library(dplyr)

data_soma %>%
  left_join(data_clin %>% select(ID, baseage), by = "ID") %>%
  ggplot(aes(x = baseage, y = KL.15384.15)) +
  geom_point(na.rm = TRUE) +                     # 👈 add points
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE)


library(ggpmisc)

data_soma %>%
  left_join(data_clin %>% select(ID, baseage), by = "ID") %>%
  ggplot(aes(x = baseage, y = KL.15384.15)) +
  geom_point(na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE
  )
