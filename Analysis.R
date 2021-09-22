library(readxl)
library(dplyr)
library(rstatix)
library(naniar)
library(epitools)
library(lubridate)
library(survival)
library(survminer)
library(brglm)
library(pROC)


df <- read_excel('Data.xlsx',
                 na = c('not avaliable','Not assessed','not assessed','UNKNOWN','Not avaliable'))

df <- df %>%
  mutate_if(is.character,as.factor)

levels(df$`IGHV gene status`) <- c('Mutated','Mutated','Unmutated','Unmutated','Unmutated')

levels(df$`TP53 mutation status`) <- c('MUTATED', 'MUTATED', 'UNMUTATED')

levels(df$`del17p                         (last assessment)`) <- c('NEGATIVE','NEGATIVE','POSITIVE','POSITIVE')

levels(df$Obesity) <- c('NO','NO','YES','YES')

levels(df$Smoking) <- c('Current smoker','Ex-smoker','NEVER','NEVER','NEVER')

levels(df$`Was the patient on treatment for CLL at the time of COVID-19?`) <- c('NO','NO', 'NO', 'YES', 'YES')

levels(df$`Was the patient on treatment with corticosteroids for CLL or other disease`) <- c('NO','NO','NO','YES')

levels(df$`CLL treatment status`) <- c('TREATED','TREATED','UNTREATED')

levels(df$Hospitalisation.) <- c('Home treatment','Home treatment','Hospitalization','Hospitalization')

levels(df$`Oxygen supplementation`) <- c('High flow','IV','NIV','NO','NO','NO', NA, 'Regular')

levels(df$ICU) <- c('NO','NO','NO','YES','YES')

levels(df$Antiviral) <- c('NO','NO', NA, 'YES', 'YES')

levels(df$`Hydroxycloroquine or similar`) <- c('NO','NO',NA,'YES','YES')

levels(df$Azithromycin) <-  c('NO','NO',NA,'YES','YES')

levels(df$Steroids) <- c('NO','NO','NO',NA,'YES','YES','YES')

levels(df$`Anti-IL6/IL6R`) <- c('NO','NO',NA,'YES')

levels(df$`Convalescent/ Hyperimmune plasma`) <- c('NO', NA, 'YES')

levels(df$`Complications of COVID-19 infection`) <- c('DIC','None','None','None','Other','Other','VTE')

levels(df$`Infection outcome`) <- c('DEATH','DEATH','Resolution','Still under medical care/confinement')

df$Hospital <- ifelse(df$`Measures taken for the management of COVID-19` == 'Confinement at home only', 'NO','YES')

df$cirs <- ifelse(df$`CIRS Score` %in% c('0','1','2','3','4','5','6'), '<=6','>6')

df$cirs[which(is.na(df$`CIRS Score`))] <- NA

df$tp53_del17p_comb <- ifelse((df$`TP53 mutation status` == 'MUTATED' | df$`del17p                         (last assessment)` == 'POSITIVE'), 'YES','NO')

df$`Treated in the last 12 months`[which(df$`CLL treatment status` == 'UNTREATED')] <- 'NO'

levels(df$`For treated patients, total lines of treatment`) <- c('>4','1','1','2','3','4',NA)


df$disease <- ifelse(df$`Measures taken for the management of COVID-19` %in% c("Confinement at home only", "Hospitalization without need of oxygen"), 'NonSevere','Severe')

df$disease[which(is.na(df$`Measures taken for the management of COVID-19`))] <- NA


df$age65 <- ifelse(df$`Age at COVID-19` >= 65, '>=65','<65')
df$age75 <- ifelse(df$`Age at COVID-19` >= 75, '>=75', '<75')

df$time <- difftime(df$End_time, df$Start_time, units = 'days')
df$time[which(df$time < 0)] <- NA

df$outcome <- ifelse(df$`Infection outcome` == 'DEATH', 1, 0)

df$survival <- Surv(time = df$time/30, event= df$outcome)


df$Wave <- ifelse(df$Start_time < '2020-06-30', 'First','Second')


df$country <- NA

df$country[which(df$Site %in% c('ASTI',
                                'AL',
                                'Bergamo',
                                '1512',
                                'OSR',
                                'Rome',
                                'PV',
                                '956',
                                'new956',
                                'Modena',
                                'Trento',
                                '782',
                                '618',
                                '736',
                                'UDINE',
                                '951',
                                '874',
                                '1239',
                                '496',
                                '883',
                                'Brescia',
                                'Emat. Torino',
                                'emat. Torino',
                                'Cuneo AO S Croce e Carle',
                                '1326',
                                'Cardarelli-Napoli',
                                'CAGLIARI-1043',
                                'LECCE UOC Ematologia',
                                'PESARO',
                                'Policlinico Tor Vergata',
                                'BergaMO'))] <- 'Italy'


df$country[which(df$Site %in% c('1271',
                                'Hospital General Universitario (Valencia, Spain)',
                                'Cordoba',
                                '1303',
                                '918',
                                '1221',
                                'Princesa',
                                '1224',
                                '564',
                                '694',
                                '565',
                                'Hospital Universitario de Burgos, Burgos, Spain',
                                '314',
                                'HUMV'))] <- 'Spain'

df$country[which(is.na(df$country))] <- 'Other'


df_treat <- df %>%
  filter(`CLL treatment status` == 'TREATED')

df_outcome <- df %>%
  filter(`Infection outcome` != "Still under medical care/confinement") %>%
  droplevels()

df_outcome_treat <- df_outcome %>%
  filter(`CLL treatment status` == 'TREATED')


df_severe <- df %>%
  filter(disease == 'Severe') %>%
  droplevels()


df_severe_outcome <- df_severe %>%
  filter(disease == 'Severe', `Infection outcome` != 'Still under medical care/confinement') %>%
  droplevels()

################################################################################



# Table 1
# Baseline patient characteristics
# All patients

df %>%
  get_summary_stats(`Age at COVID-19`, type = 'common')

n_miss(df$`Age at COVID-19`)

df %>%
  freq_table(Gender)

n_miss(df$Gender)

df %>%
  freq_table(Diagnosis)

n_miss(df$Diagnosis)


df %>%
  freq_table(Obesity)

n_miss(df$Obesity)
pct_miss(df$Obesity)

df %>%
  freq_table(Smoking)

n_miss(df$Smoking)
pct_miss(df$Smoking)

df %>%
  freq_table(Hypogammaglobulinemia)

n_miss(df$Hypogammaglobulinemia)
pct_miss(df$Hypogammaglobulinemia)

df %>%
  get_summary_stats(`CIRS Score`, type = 'common')


df %>%
  get_summary_stats(`Number of comorbidities at the time of suspected/documented COVID-19`, type = 'common')

n_miss(df$`Number of comorbidities at the time of suspected/documented COVID-19`)



df %>%
  freq_table(`Other Respiratory`)

n_miss(df$`Other Respiratory`)
pct_miss(df$`Other Respiratory`)

df %>%
  freq_table(Asthma)

n_miss(df$Asthma)
pct_miss(df$Asthma)

df %>%
  freq_table(COPD)

n_miss(df$COPD)
pct_miss(df$COPD)

df %>%
  freq_table(`Cardiac failure`)

n_miss(df$`Cardiac failure`)
pct_miss(df$`Cardiac failure`)

df %>%
  freq_table(Arrythmias)

n_miss(df$Arrythmias)
pct_miss(df$Arrythmias)

df %>%
  freq_table(`Coronary artery disease`)

n_miss(df$`Coronary artery disease`)
pct_miss(df$`Coronary artery disease`)

df %>%
  freq_table(`Other Cardiovascular`)

n_miss(df$`Other Cardiovascular`)
pct_miss(df$`Other Cardiovascular`)


df %>%
  freq_table(Hypertension)

n_miss(df$Hypertension)
pct_miss(df$Hypertension)



df %>%
  freq_table(Diabetes)

n_miss(df$Diabetes)
pct_miss(df$Diabetes)

df %>%
  freq_table(`  Chronic renal disease`)

n_miss(df$`  Chronic renal disease`)
pct_miss(df$`  Chronic renal disease`)

df %>%
  freq_table(`Other hematological malignancies`)

n_miss(df$`Other hematological malignancies`)
pct_miss(df$`Other hematological malignancies`)

df %>%
  freq_table(`  Other non hem malignancies (excluding non melanoma skin cancer)`)

n_miss(df$`  Other non hem malignancies (excluding non melanoma skin cancer)`)
pct_miss(df$`  Other non hem malignancies (excluding non melanoma skin cancer)`)


# Table 2.
# Information about CLL-directed therapy 

df %>%
  freq_table(`CLL treatment status`)

n_miss(df$`CLL treatment status`)

df %>%
  freq_table(`Treated in the last 12 months`)


n_miss(df$`Treated in the last 12 months`)
pct_miss(df$`Treated in the last 12 months`)


df %>%
  freq_table(`Was the patient on treatment for CLL at the time of COVID-19?`)

n_miss(df_treat$`Was the patient on treatment for CLL at the time of COVID-19?`)
pct_miss(df_treat$`Was the patient on treatment for CLL at the time of COVID-19?`)


n_miss(df_treat$`How did you manage the CLL treatment`)
pct_miss(df_treat$`How did you manage the CLL treatment`)

df %>%
  filter(`Was the patient on treatment for CLL at the time of COVID-19?` == 'YES') %>%
  freq_table(`How did you manage the CLL treatment`, na.rm = T)

df_treat %>%
  freq_table(`For treated patients, total lines of treatment`)

n_miss(df_treat$`For treated patients, total lines of treatment`)
pct_miss(df_treat$`For treated patients, total lines of treatment`)

df_treat_at_covid <- df %>%
  filter(`Was the patient on treatment for CLL at the time of COVID-19?` == 'YES')

df_treat_at_covid %>%
  freq_table(`Type of treatment at the time of COVID-19`)

n_miss(df_treat_at_covid$`Type of treatment at the time of COVID-19`)
pct_miss(df_treat_at_covid$`Type of treatment at the time of COVID-19`)

# Supplemental Table 2.

freq_table(df$Fever)

n_miss(df$Fever)
pct_miss(df$Fever)


freq_table(df$Dyspnea)

n_miss(df$Dyspnea)
pct_miss(df$Dyspnea)


freq_table(df$Cough)

n_miss(df$Cough)
pct_miss(df$Cough)


freq_table(df$Fatigue)

n_miss(df$Fatigue)
pct_miss(df$Fatigue)

freq_table(df$Headache)

n_miss(df$Headache)
pct_miss(df$Headache)


freq_table(df$`GI symptoms`)

n_miss(df$`GI symptoms`)
pct_miss(df$`GI symptoms`)


freq_table(df$`Anosmia/Ageusia`)

n_miss(df$`Anosmia/Ageusia`)
pct_miss(df$`Anosmia/Ageusia`)


freq_table(df$`Myalgias/Arthalgias`)
n_miss(df$`Myalgias/Arthalgias`)
pct_miss(df$`Myalgias/Arthalgias`)

freq_table(df$Other)
n_miss(df$Other)
pct_miss(df$Other)

df %>%
  get_summary_stats(`CRP x times the uln`, type = 'common')

n_miss(df$`CRP x times the uln`)
pct_miss(df$`CRP x times the uln`)

df %>%
  get_summary_stats(`D-DIMERS x times the uln`, type = 'common')

n_miss(df$`D-DIMERS x times the uln`)
pct_miss(df$`D-DIMERS x times the uln`)

df$ALC <- as.numeric(df$`ALC       (x10^9/l) (peak value during COVID-19 infection)`)


df %>%
  get_summary_stats(ALC, type = 'common')


n_miss(df$`ALC       (x10^9/l) (peak value during COVID-19 infection)`)
pct_miss(df$`ALC       (x10^9/l) (peak value during COVID-19 infection)`)

# Table 3.

df %>%
  freq_table(`Measures taken for the management of COVID-19`)

n_miss(df$`Measures taken for the management of COVID-19`)
pct_miss(df$`Measures taken for the management of COVID-19`)

df %>%
  freq_table(disease)

n_miss(df$disease)
pct_miss(df$disease)


df %>%
  freq_table(Hospital)

n_miss(df$Hospital)
pct_miss(df$Hospital)


df %>%
  freq_table(Antiviral)

n_miss(df$Antiviral)
pct_miss(df$Antiviral)

df %>%
  freq_table(`Hydroxycloroquine or similar`)

n_miss(df$`Hydroxycloroquine or similar`)
pct_miss(df$`Hydroxycloroquine or similar`)

df %>%
  freq_table(Azithromycin)

n_miss(df$Azithromycin)
pct_miss(df$Azithromycin)

df %>%
  freq_table(Steroids)

n_miss(df$Steroids)
pct_miss(df$Steroids)


df %>%
  freq_table(`Anti-IL6/IL6R`)

n_miss(df$`Anti-IL6/IL6R`)
pct_miss(df$`Anti-IL6/IL6R`)


df %>%
  freq_table(`Convalescent/ Hyperimmune plasma`)

n_miss(df$`Convalescent/ Hyperimmune plasma`)
pct_miss(df$`Convalescent/ Hyperimmune plasma`)

df %>%
  freq_table(`COVID-19 PNEUMONIA`)

n_miss(df$`COVID-19 PNEUMONIA`)
pct_miss(df$`COVID-19 PNEUMONIA`)


df %>%
  freq_table(`Complications of COVID-19 infection`)

n_miss(df$`Complications of COVID-19 infection`)
pct_miss(df$`Complications of COVID-19 infection`)


df_vte <- df %>%
  filter(`Complications of COVID-19 infection` == 'VTE')

df_vte %>%
  freq_table(PE)

n_miss(df_vte$PE)
pct_miss(df_vte$PE)


df %>%
  freq_table(`Infection outcome`)

n_miss(df$`Infection outcome`)
pct_miss(df$`Infection outcome`)


# Supplemental Table 3.
# Disease severity in all patients.


df %>%
  freq_table(age65, disease)

chisq.test(table(df$age65, df$disease))


df %>%
  freq_table(age75, disease)

chisq.test(table(df$age75, df$disease))


df %>%
  freq_table(Gender, disease)

chisq.test(table(df$Gender, df$disease))


df %>%
  freq_table(`IGHV gene status`,disease)

chisq.test(table(df$disease, df$`IGHV gene status`))


df %>%
  freq_table(`del13q                         (last assessment)`,disease)

chisq.test(table(df$disease, df$`del13q                         (last assessment)`))


df %>%
  freq_table(`del11q                        (last assessment)`, disease)

chisq.test(table(df$disease, df$`del11q                        (last assessment)`))


df %>%
  freq_table(`trisomy 12                 (last assessment)`, disease)

chisq.test(table(df$disease, df$`trisomy 12                 (last assessment)`))


df %>%
  freq_table(`del17p                         (last assessment)`, disease)

chisq.test(table(df$disease, df$`del17p                         (last assessment)`))


df %>%
  freq_table(`TP53 mutation status`, disease)

chisq.test(table(df$disease, df$`TP53 mutation status`))


df %>%
  freq_table(tp53_del17p_comb, disease)

chisq.test(table(df$tp53_del17p_comb, df$disease))

df %>%
  freq_table(cirs, disease)

chisq.test(table(df$disease, df$cirs))



df %>%
  freq_table(`Other Respiratory`, disease)

chisq.test(table(df$disease, df$`Other Respiratory`))

df %>%
  freq_table(Asthma, disease)

chisq.test(table(df$disease, df$Asthma))

df %>%
  freq_table(COPD, disease)

chisq.test(table(df$disease, df$COPD))

df %>%
  freq_table(`Cardiac failure`, disease)

chisq.test(table(df$disease, df$`Cardiac failure`))


df %>%
  freq_table(Arrythmias, disease)

chisq.test(table(df$disease, df$Arrythmias))

df %>%
  freq_table(`Coronary artery disease`, disease)

chisq.test(table(df$disease, df$`Coronary artery disease`))

df %>%
  freq_table(`Other Cardiovascular`, disease)

chisq.test(table(df$disease, df$`Other Cardiovascular`))

df %>%
  freq_table(Hypertension, disease)

chisq.test(table(df$disease, df$Hypertension))

df %>%
  freq_table(Diabetes, disease)

chisq.test(table(df$disease, df$Diabetes))

df %>%
  freq_table(`  Chronic renal disease`, disease)

chisq.test(table(df$disease, df$`  Chronic renal disease`))


df %>%
  freq_table(`Other hematological malignancies`, disease)

fisher.test(table(df$disease, df$`Other hematological malignancies`))


df %>%
  freq_table(`  Other non hem malignancies (excluding non melanoma skin cancer)`, disease)

chisq.test(table(df$disease, df$`  Other non hem malignancies (excluding non melanoma skin cancer)`))


df %>%
  freq_table(Obesity, disease)

chisq.test(table(df$disease, df$Obesity))


df %>%
  freq_table(Smoking, disease)

chisq.test(table(df$disease, df$Smoking))


df %>%
  freq_table(Hypogammaglobulinemia, disease)

chisq.test(table(df$disease, df$Hypogammaglobulinemia))

df %>%
  freq_table(`CLL treatment status`, disease)

chisq.test(table(df$disease, df$`CLL treatment status`))


df %>%
  freq_table(`Was the patient on treatment for CLL at the time of COVID-19?`, disease)

chisq.test(table(df$disease, df$`Was the patient on treatment for CLL at the time of COVID-19?`))


df %>%
  freq_table(`Treated in the last 12 months`, disease)

chisq.test(table(df$disease, df$`Treated in the last 12 months`))

df_treat %>%
  freq_table(`Treated in the last 12 months`, disease)

chisq.test(table(df_treat$`Treated in the last 12 months`, df_treat$disease))

# multivariate in disease severity
newdf <- df %>%
  dplyr::select(disease, `Age at COVID-19`, cirs, COPD, `Coronary artery disease`,
                Diabetes, `  Chronic renal disease`, Hypogammaglobulinemia) %>%
  mutate_if(is.character, as.factor) %>%
  na.omit() %>%
  droplevels()


model <- glm(disease ~ `Age at COVID-19` + cirs + COPD + `Coronary artery disease` +
               Diabetes + `  Chronic renal disease` + Hypogammaglobulinemia, 
             family = 'binomial', data = newdf)

summary(model)


model1 <- glm(disease ~ `Age at COVID-19` + cirs + `Coronary artery disease` +
                Diabetes + `  Chronic renal disease` + Hypogammaglobulinemia, 
              family = 'binomial', data = newdf)

summary(model1)

model2 <- glm(disease ~ `Age at COVID-19` + cirs + `Coronary artery disease`+
                `  Chronic renal disease` + Hypogammaglobulinemia, 
              family = 'binomial', data = newdf)

summary(model2)

model3 <- glm(disease ~ `Age at COVID-19` + cirs + `Coronary artery disease`+
                Hypogammaglobulinemia, 
              family = 'binomial', data = newdf)

summary(model3)

model4 <- glm(disease ~ `Age at COVID-19` + `Coronary artery disease`+
                Hypogammaglobulinemia, 
              family = 'binomial', data = newdf)

summary(model4)

exp(coef(model4))
exp(confint(model4))

pred <- predict(model4, type = 'response')

class <- as.factor(ifelse(pred > 0.5, 'Severe','NonSevere'))

table(class, newdf$disease)

roc(newdf$disease ~ pred, plot = T, print.auc = T)

# Supplemental table 4.
# Severity only for BTKi patients.

btki <- df %>%
  filter(`Type of treatment at the time of COVID-19` %in% c("BTKi", "BTKi+Venetoclax","BTKi+Venetoclax+ anti-CD20"))


btki %>%
  freq_table(age65, disease)

chisq.test(table(btki$age65, btki$disease))

btki %>%
  freq_table(age75, disease)

chisq.test(table(btki$age75, btki$disease))


btki %>%
  freq_table(Gender, disease)

chisq.test(table(btki$Gender, btki$disease))


btki %>%
  freq_table(`IGHV gene status`,disease)

chisq.test(table(btki$disease, btki$`IGHV gene status`))

btki %>%
  freq_table(`del13q                         (last assessment)`,disease)

chisq.test(table(btki$disease, btki$`del13q                         (last assessment)`))

btki %>%
  freq_table(`del11q                        (last assessment)`, disease)

chisq.test(table(btki$disease, btki$`del11q                        (last assessment)`))

btki %>%
  freq_table(`trisomy 12                 (last assessment)`, disease)

chisq.test(table(btki$disease, btki$`trisomy 12                 (last assessment)`))


btki %>%
  freq_table(`del17p                         (last assessment)`, disease)

chisq.test(table(btki$disease, btki$`del17p                         (last assessment)`))


btki %>%
  freq_table(`TP53 mutation status`, disease)

chisq.test(table(btki$disease, btki$`TP53 mutation status`))


btki %>%
  freq_table(tp53_del17p_comb, disease)

chisq.test(table(btki$tp53_del17p_comb, btki$disease))

btki %>%
  freq_table(cirs, disease)

chisq.test(table(btki$disease, btki$cirs))

btki %>%
  freq_table(`Other Respiratory`, disease)

chisq.test(table(btki$disease, btki$`Other Respiratory`))

btki %>%
  freq_table(Asthma, disease)

fisher.test(table(btki$disease, btki$Asthma))

btki %>%
  freq_table(COPD, disease)

chisq.test(table(btki$disease, btki$COPD))

btki %>%
  freq_table(`Cardiac failure`, disease)

fisher.test(table(btki$disease, btki$`Cardiac failure`))


btki %>%
  freq_table(Arrythmias, disease)

chisq.test(table(btki$disease, btki$Arrythmias))

btki %>%
  freq_table(`Coronary artery disease`, disease)

chisq.test(table(btki$disease, btki$`Coronary artery disease`))

btki %>%
  freq_table(`Other Cardiovascular`, disease)

chisq.test(table(btki$disease, btki$`Other Cardiovascular`))

btki %>%
  freq_table(Hypertension, disease)

chisq.test(table(btki$disease, btki$Hypertension))

btki %>%
  freq_table(Diabetes, disease)

chisq.test(table(btki$disease, btki$Diabetes))

btki %>%
  freq_table(`  Chronic renal disease`, disease)

chisq.test(table(btki$disease, btki$`  Chronic renal disease`))


btki %>%
  freq_table(`Other hematological malignancies`, disease)

fisher.test(table(btki$disease, btki$`Other hematological malignancies`))


btki %>%
  freq_table(`  Other non hem malignancies (excluding non melanoma skin cancer)`, disease)

chisq.test(table(btki$disease, btki$`  Other non hem malignancies (excluding non melanoma skin cancer)`))


btki %>%
  freq_table(Obesity, disease)

chisq.test(table(btki$disease, btki$Obesity))


btki %>%
  freq_table(Smoking, disease)

fisher.test(table(btki$disease, btki$Smoking))


btki %>%
  freq_table(Hypogammaglobulinemia, disease)

chisq.test(table(btki$disease, btki$Hypogammaglobulinemia))

btki %>%
  freq_table(Antiviral, disease)

chisq.test(table(btki$disease, btki$Antiviral))


btki %>%
  freq_table(`Hydroxycloroquine or similar`, disease)

chisq.test(table(btki$disease, btki$`Hydroxycloroquine or similar`))


btki %>%
  freq_table(Azithromycin, disease)

chisq.test(table(btki$disease, btki$Azithromycin))

btki %>%
  freq_table(Steroids, disease)

chisq.test(table(btki$disease, btki$Steroids))

btki %>%
  freq_table(`Anti-IL6/IL6R`, disease)

chisq.test(table(btki$disease, btki$`Anti-IL6/IL6R`))


btki %>%
  freq_table(`CLL treatment status`, disease)

chisq.test(table(btki$disease, btki$`CLL treatment status`))

btki %>%
  freq_table(`Treated in the last 12 months`, disease)

chisq.test(table(btki$disease, btki$`Treated in the last 12 months`))

# Multivariate analysis
model <- brglm(as.factor(disease) ~ `Age at COVID-19` + Hypogammaglobulinemia +
               `Coronary artery disease`,
             data = btki, family = 'binomial')

summary(model)

exp(coef(model))

exp(confint(model))

model1 <- brglm(as.factor(disease) ~ `Age at COVID-19` + Hypogammaglobulinemia,
               data = btki, family = 'binomial')

summary(model1)

exp(coef(model1))

exp(confint(model1))

# Correlation of Disease Severity with Infection outcome
# All patients

df %>%
  freq_table(disease, `Infection outcome`)

df %>%
  freq_table(`Infection outcome`, disease)

# Supplemental Table 5
# Baseline characteristics between the two waves.

df_severe_outcome %>%
  group_by(Wave) %>%
  get_summary_stats(`Age at COVID-19`)

df_severe_outcome %>%
  group_by(Wave) %>%
  get_summary_stats(`Number of comorbidities at the time of suspected/documented COVID-19`)


df_severe_outcome %>%
  freq_table(Wave, age65)

chisq.test(table(df_severe_outcome$age65, df_severe_outcome$Wave))


df_severe_outcome %>%
  freq_table(Wave,age75)

chisq.test(table(df_severe_outcome$age75, df_severe_outcome$Wave))


df_severe_outcome %>%
  freq_table(Wave,Gender)

chisq.test(table(df_severe_outcome$Gender, df_severe_outcome$Wave))


df_severe_outcome %>%
  freq_table(Wave,`IGHV gene status`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`IGHV gene status`))


df_severe_outcome %>%
  freq_table(Wave,`del13q                         (last assessment)`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`del13q                         (last assessment)`))


df_severe_outcome %>%
  freq_table(Wave,`del11q                        (last assessment)`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`del11q                        (last assessment)`))


df_severe_outcome %>%
  freq_table(Wave,`trisomy 12                 (last assessment)`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`trisomy 12                 (last assessment)`))


df_severe_outcome %>%
  freq_table(Wave,`del17p                         (last assessment)`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`del17p                         (last assessment)`))


df_severe_outcome %>%
  freq_table(Wave,`TP53 mutation status`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`TP53 mutation status`))


df_severe_outcome %>%
  freq_table(Wave,tp53_del17p_comb)

chisq.test(table(df_severe_outcome$tp53_del17p_comb, df_severe_outcome$Wave))


df_severe_outcome %>%
  freq_table(Wave,cirs)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$cirs))


df_severe_outcome %>%
  freq_table(Wave,`Other Respiratory`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`Other Respiratory`))

df_severe_outcome %>%
  freq_table(Wave,Asthma)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$Asthma))

df_severe_outcome %>%
  freq_table(Wave,COPD)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$COPD))

df_severe_outcome %>%
  freq_table(Wave,`Cardiac failure`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`Cardiac failure`))


df_severe_outcome %>%
  freq_table(Wave,Arrythmias)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$Arrythmias))

df_severe_outcome %>%
  freq_table(Wave,`Coronary artery disease`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`Coronary artery disease`))

df_severe_outcome %>%
  freq_table(Wave,`Other Cardiovascular`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`Other Cardiovascular`))

df_severe_outcome %>%
  freq_table(Wave,Hypertension)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$Hypertension))

df_severe_outcome %>%
  freq_table(Wave,Diabetes)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$Diabetes))

df_severe_outcome %>%
  freq_table(Wave,`  Chronic renal disease`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`  Chronic renal Wave`))


df_severe_outcome %>%
  freq_table(Wave,`Other hematological malignancies`)

fisher.test(table(df_severe_outcome$Wave, df_severe_outcome$`Other hematological malignancies`))


df_severe_outcome %>%
  freq_table(Wave,`  Other non hem malignancies (excluding non melanoma skin cancer)`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`  Other non hem malignancies (excluding non melanoma skin cancer)`))


df_severe_outcome %>%
  freq_table(Wave,Obesity)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$Obesity))


df_severe_outcome %>%
  freq_table(Wave,Smoking)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$Smoking))


df_severe_outcome %>%
  freq_table(Wave,Hypogammaglobulinemia)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$Hypogammaglobulinemia))



df_severe_outcome %>%
  freq_table(Wave,`CLL treatment status`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`CLL treatment status`))

df_severe_outcome %>%
  freq_table(Wave,`Was the patient on treatment for CLL at the time of COVID-19?`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`Was the patient on treatment for CLL at the time of COVID-19?`))


df_severe_outcome %>%
  freq_table(Wave,`Treated in the last 12 months`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`Treated in the last 12 months`))

df_severe_outcome %>%
  filter(`CLL treatment status` == 'TREATED') %>%
  freq_table(Wave,`Type of last treatment  in the last 12 months`)


df_severe_outcome %>%
  freq_table(Wave,Antiviral)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$Antiviral))



df_severe_outcome %>%
  freq_table(Wave,`Hydroxycloroquine or similar`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`Hydroxycloroquine or similar`))


df_severe_outcome %>%
  freq_table(Wave,Azithromycin)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$Azithromycin))



df_severe_outcome %>%
  freq_table(Wave,Steroids)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$Steroids))


df_severe_outcome %>%
  freq_table(Wave,`Anti-IL6/IL6R`)

chisq.test(table(df_severe_outcome$Wave, df_severe_outcome$`Anti-IL6/IL6R`))


df_severe_outcome %>%
  freq_table(Wave, `Infection outcome`)

chisq.test(df_severe_outcome$Wave, df_severe_outcome$`Infection outcome`)

# Supplemental Table 6.
# Comparison between the two waves, after controlling for age

## For all patients

model <- glm(`Infection outcome` ~ Wave + `Age at COVID-19`,
             data = df_outcome, family = 'binomial')


summary(model)

exp(coef(model))


## For severe patients

model <- glm(`Infection outcome` ~ Wave + `Age at COVID-19`,
             data = df_severe_outcome, family = 'binomial')


summary(model)

exp(coef(model))


# Table 4.
# Risk factors of infection outcome for patients with severe COVID-19.

df_severe_outcome %>%
  freq_table(age65, `Infection outcome`)

chisq.test(table(df_severe_outcome$age65, df_severe_outcome$`Infection outcome`))

df_severe_outcome %>%
  freq_table(age75, `Infection outcome`)

chisq.test(table(df_severe_outcome$age75, df_severe_outcome$`Infection outcome`))


df_severe_outcome %>%
  freq_table(Gender, `Infection outcome`)

chisq.test(table(df_severe_outcome$Gender, df_severe_outcome$`Infection outcome`))


df_severe_outcome %>%
  freq_table(`IGHV gene status`,`Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`IGHV gene status`))

df_severe_outcome %>%
  freq_table(`del13q                         (last assessment)`,`Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`del13q                         (last assessment)`))

df_severe_outcome %>%
  freq_table(`del11q                        (last assessment)`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`del11q                        (last assessment)`))

df_severe_outcome %>%
  freq_table(`trisomy 12                 (last assessment)`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`trisomy 12                 (last assessment)`))


df_severe_outcome %>%
  freq_table(`del17p                         (last assessment)`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`del17p                         (last assessment)`))


df_severe_outcome %>%
  freq_table(`TP53 mutation status`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`TP53 mutation status`))


df_severe_outcome %>%
  freq_table(tp53_del17p_comb, `Infection outcome`)

chisq.test(table(df_severe_outcome$tp53_del17p_comb, df_severe_outcome$`Infection outcome`))


df_severe_outcome %>%
  freq_table(cirs, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$cirs))


df_severe_outcome %>%
  freq_table(`Other Respiratory`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`Other Respiratory`))

df_severe_outcome %>%
  freq_table(Asthma, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$Asthma))

df_severe_outcome %>%
  freq_table(COPD, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$COPD))

df_severe_outcome %>%
  freq_table(`Cardiac failure`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`Cardiac failure`))


df_severe_outcome %>%
  freq_table(Arrythmias, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$Arrythmias))

df_severe_outcome %>%
  freq_table(`Coronary artery disease`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`Coronary artery disease`))

df_severe_outcome %>%
  freq_table(`Other Cardiovascular`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`Other Cardiovascular`))

df_severe_outcome %>%
  freq_table(Hypertension, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$Hypertension))

df_severe_outcome %>%
  freq_table(Diabetes, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$Diabetes))

df_severe_outcome %>%
  freq_table(`  Chronic renal disease`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`  Chronic renal disease`))


df_severe_outcome %>%
  freq_table(`Other hematological malignancies`, `Infection outcome`)

fisher.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`Other hematological malignancies`))


df_severe_outcome %>%
  freq_table(`  Other non hem malignancies (excluding non melanoma skin cancer)`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`  Other non hem malignancies (excluding non melanoma skin cancer)`))


df_severe_outcome %>%
  freq_table(Obesity, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$Obesity))


df_severe_outcome %>%
  freq_table(Smoking, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$Smoking))


df_severe_outcome %>%
  freq_table(Hypogammaglobulinemia, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$Hypogammaglobulinemia))


df_severe_outcome %>%
  freq_table(`CLL treatment status`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`CLL treatment status`))



df_severe_outcome %>%
  freq_table(`Was the patient on treatment for CLL at the time of COVID-19?`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`Was the patient on treatment for CLL at the time of COVID-19?`))


df_severe_outcome %>%
  freq_table(`Treated in the last 12 months`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`Treated in the last 12 months`))



df_severe_outcome %>%
  freq_table(Antiviral, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$Antiviral))


df_severe_outcome %>%
  freq_table(`Hydroxycloroquine or similar`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`Hydroxycloroquine or similar`))


df_severe_outcome %>%
  freq_table(Azithromycin, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$Azithromycin))

df_severe_outcome %>%
  freq_table(Steroids, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$Steroids))

df_severe_outcome %>%
  freq_table(`Anti-IL6/IL6R`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Infection outcome`, df_severe_outcome$`Anti-IL6/IL6R`))

## Multivariate analysis

df_severe_outcome$`Infection outcome` <- relevel(df_severe_outcome$`Infection outcome`,
                                                 ref = 'Resolution')

model <- glm(`Infection outcome` ~ `Age at COVID-19` + `IGHV gene status` +
               cirs + `Coronary artery disease` + `  Chronic renal disease` +
               `Cardiac failure` + `Treated in the last 12 months`,
             family = 'binomial',
             data = df_severe_outcome)

model <- glm(`Infection outcome` ~ `Age at COVID-19` +
               `Cardiac failure` + `Treated in the last 12 months`,
             family = 'binomial',
             data = df_severe_outcome)

summary(model)

exp(coef(model))

exp(confint(model))


# OUT start
df_severe_outcome %>%
  freq_table(`Treated  in the last 3 months`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Treated  in the last 3 months`, df_severe_outcome$`Infection outcome`))


df_severe_outcome %>%
  freq_table(`Treated  in the last 6 months`, `Infection outcome`)

chisq.test(table(df_severe_outcome$`Treated  in the last 6 months`, df_severe_outcome$`Infection outcome`))



df_severe_outcome_treat <- df_severe_outcome %>%
  filter(`CLL treatment status` == 'TREATED')

df_severe_outcome_treat %>%
  freq_table(`Treated in the last 12 months`, `Infection outcome`)

chisq.test(table(df_severe_outcome_treat$`Infection outcome`, df_severe_outcome_treat$`Treated in the last 12 months`))
# OUT end




# Supplemental table 7.
## Repeat the analysis only for the Sites which participated in the first study.

desired <- c('OSR', 'Bergamo',
             'Brescia',
             '1512',
             '694',
             '1271',
             '1221',
             '1303',
             '1224','1623','956','145',
             '564','ASTI',
             'AL','618','PV','314','Modena',
             '883','560',
             '1677','1239','Rome',
             '1640','Trento',
             '1326','1671','638',
             '951','648','Cologne')


desired_df <- df_severe_outcome %>%
  filter(Site %in% desired)


desired_df %>%
  group_by(Wave) %>%
  get_summary_stats(`Age at COVID-19`)

desired_df %>%
  group_by(Wave) %>%
  get_summary_stats(`Number of comorbidities at the time of suspected/documented COVID-19`)


desired_df %>%
  freq_table(Wave, age65)

chisq.test(table(desired_df$age65, desired_df$Wave))


desired_df %>%
  freq_table(Wave,age75)

chisq.test(table(desired_df$age75, desired_df$Wave))


desired_df %>%
  freq_table(Wave,Gender)

chisq.test(table(desired_df$Gender, desired_df$Wave))


desired_df %>%
  freq_table(Wave,`IGHV gene status`)

chisq.test(table(desired_df$Wave, desired_df$`IGHV gene status`))


desired_df %>%
  freq_table(Wave,`del13q                         (last assessment)`)

chisq.test(table(desired_df$Wave, desired_df$`del13q                         (last assessment)`))


desired_df %>%
  freq_table(Wave,`del11q                        (last assessment)`)

chisq.test(table(desired_df$Wave, desired_df$`del11q                        (last assessment)`))


desired_df %>%
  freq_table(Wave,`trisomy 12                 (last assessment)`)

chisq.test(table(desired_df$Wave, desired_df$`trisomy 12                 (last assessment)`))


desired_df %>%
  freq_table(Wave,`del17p                         (last assessment)`)

chisq.test(table(desired_df$Wave, desired_df$`del17p                         (last assessment)`))


desired_df %>%
  freq_table(Wave,`TP53 mutation status`)

chisq.test(table(desired_df$Wave, desired_df$`TP53 mutation status`))


desired_df %>%
  freq_table(Wave,tp53_del17p_comb)

chisq.test(table(desired_df$tp53_del17p_comb, desired_df$Wave))


desired_df %>%
  freq_table(Wave,cirs)

chisq.test(table(desired_df$Wave, desired_df$cirs))


desired_df %>%
  freq_table(Wave,`Other Respiratory`)

chisq.test(table(desired_df$Wave, desired_df$`Other Respiratory`))

desired_df %>%
  freq_table(Wave,Asthma)

chisq.test(table(desired_df$Wave, desired_df$Asthma))

desired_df %>%
  freq_table(Wave,COPD)

chisq.test(table(desired_df$Wave, desired_df$COPD))

desired_df %>%
  freq_table(Wave,`Cardiac failure`)

chisq.test(table(desired_df$Wave, desired_df$`Cardiac failure`))


desired_df %>%
  freq_table(Wave,Arrythmias)

chisq.test(table(desired_df$Wave, desired_df$Arrythmias))

desired_df %>%
  freq_table(Wave,`Coronary artery disease`)

chisq.test(table(desired_df$Wave, desired_df$`Coronary artery disease`))

desired_df %>%
  freq_table(Wave,`Other Cardiovascular`)

chisq.test(table(desired_df$Wave, desired_df$`Other Cardiovascular`))

desired_df %>%
  freq_table(Wave,Hypertension)

chisq.test(table(desired_df$Wave, desired_df$Hypertension))

desired_df %>%
  freq_table(Wave,Diabetes)

chisq.test(table(desired_df$Wave, desired_df$Diabetes))

desired_df %>%
  freq_table(Wave,`  Chronic renal disease`)

chisq.test(table(desired_df$Wave, desired_df$`  Chronic renal Wave`))


desired_df %>%
  freq_table(Wave,`Other hematological malignancies`)

fisher.test(table(desired_df$Wave, desired_df$`Other hematological malignancies`))


desired_df %>%
  freq_table(Wave,`  Other non hem malignancies (excluding non melanoma skin cancer)`)

chisq.test(table(desired_df$Wave, desired_df$`  Other non hem malignancies (excluding non melanoma skin cancer)`))


desired_df %>%
  freq_table(Wave,Obesity)

chisq.test(table(desired_df$Wave, desired_df$Obesity))


desired_df %>%
  freq_table(Wave,Smoking)

chisq.test(table(desired_df$Wave, desired_df$Smoking))


desired_df %>%
  freq_table(Wave,Hypogammaglobulinemia)

chisq.test(table(desired_df$Wave, desired_df$Hypogammaglobulinemia))



desired_df %>%
  freq_table(Wave,`CLL treatment status`)

chisq.test(table(desired_df$Wave, desired_df$`CLL treatment status`))

desired_df %>%
  freq_table(Wave,`Was the patient on treatment for CLL at the time of COVID-19?`)

chisq.test(table(desired_df$Wave, desired_df$`Was the patient on treatment for CLL at the time of COVID-19?`))


desired_df %>%
  freq_table(Wave,`Treated in the last 12 months`)

chisq.test(table(desired_df$Wave, desired_df$`Treated in the last 12 months`))

desired_df %>%
  filter(`CLL treatment status` == 'TREATED') %>%
  freq_table(Wave,`Type of last treatment  in the last 12 months`)


desired_df %>%
  freq_table(Wave,Antiviral)

chisq.test(table(desired_df$Wave, desired_df$Antiviral))



desired_df %>%
  freq_table(Wave,`Hydroxycloroquine or similar`)

chisq.test(table(desired_df$Wave, desired_df$`Hydroxycloroquine or similar`))


desired_df %>%
  freq_table(Wave,Azithromycin)

chisq.test(table(desired_df$Wave, desired_df$Azithromycin))



desired_df %>%
  freq_table(Wave,Steroids)

chisq.test(table(desired_df$Wave, desired_df$Steroids))


desired_df %>%
  freq_table(Wave,`Anti-IL6/IL6R`)

chisq.test(table(desired_df$Wave, desired_df$`Anti-IL6/IL6R`))

desired_df_outcome <- desired_df %>%
  filter(`Infection outcome` != 'Still under medical care/confinement')

desired_df_outcome %>%
  freq_table(Wave, `Infection outcome`)

chisq.test(desired_df_outcome$Wave, desired_df_outcome$`Infection outcome`)

# Supplemental Table 8.
# Risk factors of infection outcome for patients with severe COVID-19,
# only for centers that participated in both studies .

desired_df_outcome %>%
  freq_table(age65, `Infection outcome`)

chisq.test(table(desired_df_outcome$age65, desired_df_outcome$`Infection outcome`))

desired_df_outcome %>%
  freq_table(age75, `Infection outcome`)

chisq.test(table(desired_df_outcome$age75, desired_df_outcome$`Infection outcome`))


desired_df_outcome %>%
  freq_table(Gender, `Infection outcome`)

chisq.test(table(desired_df_outcome$Gender, desired_df_outcome$`Infection outcome`))


desired_df_outcome %>%
  freq_table(`IGHV gene status`,`Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`IGHV gene status`))

desired_df_outcome %>%
  freq_table(`del13q                         (last assessment)`,`Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`del13q                         (last assessment)`))

desired_df_outcome %>%
  freq_table(`del11q                        (last assessment)`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`del11q                        (last assessment)`))

desired_df_outcome %>%
  freq_table(`trisomy 12                 (last assessment)`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`trisomy 12                 (last assessment)`))


desired_df_outcome %>%
  freq_table(`del17p                         (last assessment)`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`del17p                         (last assessment)`))


desired_df_outcome %>%
  freq_table(`TP53 mutation status`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`TP53 mutation status`))


desired_df_outcome %>%
  freq_table(tp53_del17p_comb, `Infection outcome`)

chisq.test(table(desired_df_outcome$tp53_del17p_comb, desired_df_outcome$`Infection outcome`))


desired_df_outcome %>%
  freq_table(cirs, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$cirs))


desired_df_outcome %>%
  freq_table(`Other Respiratory`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`Other Respiratory`))

desired_df_outcome %>%
  freq_table(Asthma, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$Asthma))

desired_df_outcome %>%
  freq_table(COPD, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$COPD))

desired_df_outcome %>%
  freq_table(`Cardiac failure`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`Cardiac failure`))


desired_df_outcome %>%
  freq_table(Arrythmias, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$Arrythmias))

desired_df_outcome %>%
  freq_table(`Coronary artery disease`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`Coronary artery disease`))

desired_df_outcome %>%
  freq_table(`Other Cardiovascular`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`Other Cardiovascular`))

desired_df_outcome %>%
  freq_table(Hypertension, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$Hypertension))

desired_df_outcome %>%
  freq_table(Diabetes, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$Diabetes))

desired_df_outcome %>%
  freq_table(`  Chronic renal disease`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`  Chronic renal disease`))


desired_df_outcome %>%
  freq_table(`Other hematological malignancies`, `Infection outcome`)

fisher.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`Other hematological malignancies`))


desired_df_outcome %>%
  freq_table(`  Other non hem malignancies (excluding non melanoma skin cancer)`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`  Other non hem malignancies (excluding non melanoma skin cancer)`))


desired_df_outcome %>%
  freq_table(Obesity, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$Obesity))


desired_df_outcome %>%
  freq_table(Smoking, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$Smoking))


desired_df_outcome %>%
  freq_table(Hypogammaglobulinemia, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$Hypogammaglobulinemia))


desired_df_outcome %>%
  freq_table(`CLL treatment status`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`CLL treatment status`))



desired_df_outcome %>%
  freq_table(`Was the patient on treatment for CLL at the time of COVID-19?`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`Was the patient on treatment for CLL at the time of COVID-19?`))


desired_df_outcome %>%
  freq_table(`Treated in the last 12 months`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`Treated in the last 12 months`))



desired_df_outcome %>%
  freq_table(Antiviral, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$Antiviral))


desired_df_outcome %>%
  freq_table(`Hydroxycloroquine or similar`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`Hydroxycloroquine or similar`))


desired_df_outcome %>%
  freq_table(Azithromycin, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$Azithromycin))

desired_df_outcome %>%
  freq_table(Steroids, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$Steroids))

desired_df_outcome %>%
  freq_table(`Anti-IL6/IL6R`, `Infection outcome`)

chisq.test(table(desired_df_outcome$`Infection outcome`, desired_df_outcome$`Anti-IL6/IL6R`))

## Multivariate analysis

desired_df_outcome$`Infection outcome` <- relevel(desired_df_outcome$`Infection outcome`,
                                                 ref = 'Resolution')

desired_df_outcome$`CLL treatment status` <- relevel(desired_df_outcome$`CLL treatment status`,
                                                     ref = 'UNTREATED')

model <- glm(`Infection outcome` ~ `Age at COVID-19` + cirs +
               `  Chronic renal disease` +
               `Cardiac failure` + `CLL treatment status`,
             family = 'binomial',
             data = desired_df_outcome)

summary(model)

model1 <- glm(`Infection outcome` ~ `Age at COVID-19` + cirs +
               `  Chronic renal disease` +
               `CLL treatment status`,
             family = 'binomial',
             data = desired_df_outcome)

summary(model1)

model2 <- glm(`Infection outcome` ~ `Age at COVID-19` + cirs +
                `CLL treatment status`,
              family = 'binomial',
              data = desired_df_outcome)

summary(model2)

model3 <- glm(`Infection outcome` ~ `Age at COVID-19` +
                `CLL treatment status`,
              family = 'binomial',
              data = desired_df_outcome)

summary(model3)

exp(coef(model3))

exp(confint(model3))


# Supplemental Figure 1.
## Number of cases per month

cases <- as.data.frame(table(df$`Month of Diagnosis`))
colnames(cases) <- c('Month of Diagnosis', 'Frequency')
cases$Month <- factor(cases$`Month of Diagnosis`, levels = c('1/20.', '2/20.','3/20.','4/20.','5/20.','6/20.','7/20.',
                                                             '8/20.','9/20.','10/20.','11/20.','12/20.','1/21.','2/21.','3/21.'))

ggplot(cases, aes(x = Month, y = Frequency)) +
  geom_line(group = 1, colour = 'navy') +
  theme_classic() 

# Table 5.
# Risk factors of OS.

newdf <- df_severe %>%
  dplyr::select(survival, Gender, `IGHV gene status`,
                `del13q                         (last assessment)`,
                `del11q                        (last assessment)`,
                `del17p                         (last assessment)`,
                `trisomy 12                 (last assessment)`,
                `TP53 mutation status`,
                cirs,
                `Other Respiratory`,
                Asthma, COPD, `Cardiac failure`,
                Arrythmias, `Coronary artery disease`,
                `Other Cardiovascular`,
                Hypertension,
                Diabetes,
                `  Chronic renal disease`,
                `Other hematological malignancies`,
                `  Other non hem malignancies (excluding non melanoma skin cancer)`,
                Obesity, Smoking, Hypogammaglobulinemia,
                `CLL treatment status`,
                `Was the patient on treatment for CLL at the time of COVID-19?`,
                `Treated in the last 12 months`)

newdf$cirs <- as.factor(newdf$cirs)

for (i in 2:ncol(newdf)){
  filename <- paste(i,'.jpeg')
  pred <- newdf[,i]
  km <- survfit(survival ~ pred, data = newdf)
  p <- ggsurvplot(fit=km,
                  xlab='Months',
                  ylab = 'Overall survival probability',
                  surv.median.line = "hv",
                  conf.int = FALSE,
                  pval=TRUE,
                  risk.table = FALSE,
                  palette = c('blue','red','green'),
                  legend.labs = paste(colnames(newdf)[i], levels(newdf[,i]))) 
  jpeg(file = filename)
  print(p)
  dev.off()
}

# Multivariate for OS 
# Only for severe patients


coxmodel <- coxph(survival ~ df_severe_outcome$`Age at COVID-19` +
                    df_severe_outcome$`CLL treatment status` +
                    df_severe_outcome$`Cardiac failure`,
                  data = df_severe_outcome)

summary(coxmodel)$coefficients

exp(confint(coxmodel))

# Figure 4 and Supplemental Figure 4.

df_outcome$`Treatment categories (ANTI-CD20) or Untreated` <- relevel(df_outcome$`Treatment categories (ANTI-CD20) or Untreated`, ref = 'Untreated')

df_outcome$`Infection outcome` <- relevel(df_outcome$`Infection outcome`, ref = 'Resolution')


### Correlation of Infection outcome between patients who received
### BTKi at the time of COVID-19
### Venetoclax at the time of COVID-19
### Chemoimmunotherapy in the last 12 months

my_data <- df_outcome %>%
  dplyr::select(`Infection outcome`, `Treatment categories or Untreated`) %>%
  filter(`Treatment categories or Untreated` != 'Untreated') %>%
  droplevels() %>%
  na.omit()

table(my_data$`Treatment categories or Untreated`)

model1 <- glm(`Infection outcome` ~  `Treatment categories or Untreated`,
              family = 'binomial',
              data = my_data)

model2 <- glm(`Infection outcome` ~ 1,
              family = 'binomial',
              data = my_data)



anova(model1, model2, test = 'LRT')


### Correlation of Infection outcome between patients who received
### BTKi at the time of COVID-19
### Venetoclax at the time of COVID-19
### Anti-CD20 in the last 12 months

my_data <- df_outcome %>%
  dplyr::select(`Infection outcome`, `Treatment categories (ANTI-CD20) or Untreated`) %>%
  filter(`Treatment categories (ANTI-CD20) or Untreated` != 'Untreated') %>%
  droplevels() %>%
  na.omit()

table(my_data$`Treatment categories (ANTI-CD20) or Untreated`)


model1 <- glm(`Infection outcome` ~ `Treatment categories (ANTI-CD20) or Untreated`,
              family = 'binomial',
              data = my_data)

model2 <- glm(`Infection outcome` ~ 1,
              family = 'binomial',
              data = my_data)



anova(model1, model2, test = 'LRT')


my_data <- df_outcome %>%
  filter(`Treatment categories (ANTI-CD20) or Untreated` != 'Untreated') %>%
  droplevels() %>%
  dplyr::select(`Infection outcome`, `Treatment categories (ANTI-CD20) or Untreated`) %>%
  na.omit()

model1 <- glm(`Infection outcome` ~ `Treatment categories (ANTI-CD20) or Untreated`,
              family = 'binomial',
              data = my_data)

model2 <- glm(`Infection outcome` ~ 1,
              family = 'binomial',
              data = my_data)


anova(model1, model2, test = 'LRT')

exp(coef(model1))

exp(confint(model1))



#HERE

#########################################################################################################


# Correlations of Infection outcome with clinico-biological characteristics,
# COVID-19 disease severity, CLL treatment status, type of treatment,
# Management of treatment
## For all patients

df_outcome %>%
  freq_table(age65, `Infection outcome`)

chisq.test(table(df$age65, df$`Infection outcome`))


df_outcome %>%
  freq_table(age75, `Infection outcome`)

chisq.test(table(df$age75, df$`Infection outcome`))

df_outcome %>%
  freq_table(Gender, `Infection outcome`)

chisq.test(table(df_outcome$Gender, df_outcome$`Infection outcome`))


df_outcome %>%
  freq_table(`IGHV gene status`,`Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`IGHV gene status`))


df_outcome %>%
  freq_table(`del13q                         (last assessment)`,`Infection outcome`)


chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`del13q                         (last assessment)`))

df_outcome %>%
  freq_table(`del11q                        (last assessment)`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`del11q                        (last assessment)`))


df_outcome$`trisomy 12                 (last assessment)`

df_outcome %>%
  freq_table(`trisomy 12                 (last assessment)`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`trisomy 12                 (last assessment)`))


df_outcome %>%
  freq_table(`del17p                         (last assessment)`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`del17p                         (last assessment)`))


df_outcome %>%
  freq_table(`TP53 mutation status`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`TP53 mutation status`))


df_outcome %>%
  freq_table(tp53_del17p_comb, `Infection outcome`)

chisq.test(table(df_outcome$tp53_del17p_comb, df_outcome$`Infection outcome`))

df_outcome %>%
  freq_table(`Other Respiratory`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`Other Respiratory`))

df_outcome %>%
  freq_table(Asthma, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$Asthma))

df_outcome %>%
  freq_table(COPD, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$COPD))

df_outcome %>%
  freq_table(`Cardiac failure`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`Cardiac failure`))


df_outcome %>%
  freq_table(Arrythmias, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$Arrythmias))

df_outcome %>%
  freq_table(`Coronary artery disease`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`Coronary artery disease`))

df_outcome %>%
  freq_table(`Other Cardiovascular`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`Other Cardiovascular`))

df_outcome %>%
  freq_table(Hypertension, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$Hypertension))

df_outcome %>%
  freq_table(Diabetes, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$Diabetes))

df_outcome %>%
  freq_table(`  Chronic renal disease`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`  Chronic renal disease`))


df_outcome %>%
  freq_table(`Other hematological malignancies`, `Infection outcome`)

fisher.test(table(df_outcome$`Infection outcome`, df_outcome$`Other hematological malignancies`))


df_outcome %>%
  freq_table(`  Other non hem malignancies (excluding non melanoma skin cancer)`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`  Other non hem malignancies (excluding non melanoma skin cancer)`))


df_outcome %>%
  freq_table(Obesity, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$Obesity))


df_outcome %>%
  freq_table(Smoking, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$Smoking))


df_outcome %>%
  freq_table(Hypogammaglobulinemia, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$Hypogammaglobulinemia))



df_outcome %>%
  freq_table(`CLL treatment status`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`CLL treatment status`))


df_outcome %>%
  freq_table(`Treated in the last 12 months`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`Treated in the last 12 months`))



df_outcome %>%
  freq_table(Antiviral, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$Antiviral))



df_outcome %>%
  freq_table(`Hydroxycloroquine or similar`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`Hydroxycloroquine or similar`))


df_outcome %>%
  freq_table(Azithromycin, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$Azithromycin))



df_outcome %>%
  freq_table(Steroids, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$Steroids))


df_outcome %>%
  freq_table(`Anti-IL6/IL6R`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`Anti-IL6/IL6R`))


df_outcome %>%
  freq_table(cirs, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$cirs))


df_outcome %>%
  freq_table(`Was the patient on treatment for CLL at the time of COVID-19?`, `Infection outcome`)

chisq.test(table(df_outcome$`Infection outcome`, df_outcome$`Was the patient on treatment for CLL at the time of COVID-19?`))


df_outcome %>%
  group_by(`Infection outcome`) %>%
  get_summary_stats(`Age at COVID-19`)


age <- glm(`Infection outcome` ~ `Age at COVID-19`, family = 'binomial', data = df_outcome)

summary(age)
exp(coef(age))


df_outcome %>%
  group_by(`Infection outcome`) %>%
  get_summary_stats(`Number of comorbidities at the time of suspected/documented COVID-19`)

comorbidities <- glm(`Infection outcome` ~ `Number of comorbidities at the time of suspected/documented COVID-19`, family = 'binomial', data = df_outcome)

summary(comorbidities)
exp(coef(comorbidities))


df_outcome_treat %>%
  freq_table(`Treated in the last 12 months`, `Infection outcome`)

chisq.test(table(df_outcome_treat$`Treated in the last 12 months`, df_outcome_treat$`Infection outcome`))

########################################################################################################


# Supplemental Figure 2.

km <- survfit(survival ~ 1, data = df)

ggsurvplot(data = df, fit=km, xlab='Months',
           ylab = 'Overall survival probability',
           surv.median.line = "hv",
           conf.int = FALSE, pval = F,
           censor = F,
           palette = 'blue',
           legend = 'none',
           legend.title = '',
           legend.labs = '',
           ggtheme = theme_classic(), risk.table = T,
           break.x.by = 1,
           xlim = c(0,8))


# Supplemental Figure 3.

km <- survfit(survival ~ 1, data = df_severe_outcome)

ggsurvplot(data = df_severe_outcome, fit=km, xlab='Months',
           ylab = 'Overall survival probability',
           surv.median.line = "hv",
           conf.int = FALSE, pval = F,
           censor = F,
           palette = 'blue',
           legend = 'none',
           legend.title = '',
           legend.labs = '',
           ggtheme = theme_classic(), risk.table = T,
           break.x.by = 1,
           xlim = c(0,8))

# Supplemental Figure 4.

km <- survfit(survival ~ df_severe_outcome$`Treatment categories (ANTI-CD20) or Untreated`,
              data = df_severe_outcome)

ggsurvplot(fit=km,
           xlab='Months',
           ylab = 'Overall survival probability',
           surv.median.line = "hv",
           conf.int = FALSE,
           pval=TRUE,
           censor = F,
           risk.table = T,
           legend.labs = c('Anti-CD20','BTKi', 'Untreated', 'Venetoclax')) 


# Figure 1.
# Overall survival in patients with severe COVID-19:
# comparison between treated and untreated patients
km <- survfit(survival ~ df_severe_outcome$`CLL treatment status`, data = df_severe_outcome)

ggsurvplot(fit=km,
           xlab='Months',
                ylab = 'Overall survival probability',
           surv.median.line = "hv",
           conf.int = FALSE,
           pval=TRUE,
           censor = F,
           risk.table = T,
           palette = c('blue','red'),
           legend.labs = c('Treated', 'Untreated'),
           legend.title = 'CLL treatment status',
           break.x.by = 1,
           xlim = c(0,8)) 



# Figure 2.
# Overall survival in patients with severe COVID-19 according to treatment:
# comparisons between
# BTKi (At time of COVID-19),
# Venetoclax (At time of COVID-19),
# Chemoimmunotherapy in last 12 months,
# Untreated


km <- survfit(survival ~ df_severe_outcome$Comparison4,
              data = df_severe_outcome)

ggsurvplot(fit=km,
           xlab='Months',
           ylab = 'Overall survival probability',
           surv.median.line = "hv",
           conf.int = FALSE,
           pval=TRUE,
           risk.table = T,
           censor = F,
           legend.labs = c('BTKi', 'Chemoimmunotherapy', 'Untreated', 'Venetoclax'),
           break.x.by = 1,
           xlim = c(0,8),
           legend.title = element_blank()) 


# Table 6.

freq_table(df_severe_outcome$`Type of treatment at the time of COVID-19`)

btki_at_covid <- df_severe_outcome %>%
  filter(`Type of treatment at the time of COVID-19` == 'BTKi')

btki_at_covid %>%
  freq_table(`How did you manage the CLL treatment`)


df_severe_outcome %>%
  freq_table(`Treatment strategy only for those with BTKi or Untreated`,
             `Infection outcome`)

chisq.test(df_severe_outcome$`Treatment strategy only for those with BTKi or Untreated`,
           df_severe_outcome$`Infection outcome`)
