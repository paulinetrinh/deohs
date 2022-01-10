library(tidyverse)
library(breakaway)
library(magrittr)
library(phyloseq)

readRDS("/home/gradstudent/ptrinh88/Iquitos_fecal.rds") -> iquitos
# update metadata 
readRDS("/home/gradstudent/ptrinh88/Iquitos_metadata.rds") -> metadata
readRDS("/home/gradstudent/ptrinh88/Iquitos_metadata_diet.RDS") -> metadata2
readRDS("/home/gradstudent/ptrinh88/Iquitos_metadata_diet_updateApr28.RDS") -> metadata2

metadata2 %>% 
  subset(Time == "18-Feb") %>% 
  mutate(tertile_CHO_feb =`CHO.Energy%.overall.feb` %>% as.numeric() %>% ntile(.,3)) %>%
  mutate(tertile_PRO_feb = `PRO.Energy%.overall.feb` %>% as.numeric() %>% ntile(.,3)) %>% 
  mutate(tertile_FAT_feb = `FAT.Energy%.overall.feb` %>% as.numeric() %>% ntile(.,3)) %>% 
  mutate(tertile_FIBER_feb = `total_KCAL_feb_dietary_fiber_Energy%` %>% as.numeric() %>% ntile(.,3)) %>% 
  mutate(tertile_CHO_feb_factor = ifelse(tertile_CHO_feb == 1, "low",
                                         ifelse(tertile_CHO_feb == 2, "middle", 
                                                ifelse(tertile_CHO_feb == 3, "high",NA)))) %>% 
  mutate(tertile_PRO_feb_factor = ifelse(tertile_PRO_feb == 1, "low",
                                         ifelse(tertile_PRO_feb == 2, "middle", 
                                                ifelse(tertile_PRO_feb == 3, "high",NA)))) %>% 
  mutate(tertile_FAT_feb_factor = ifelse(tertile_FAT_feb == 1, "low",
                                         ifelse(tertile_FAT_feb == 2, "middle", 
                                                ifelse(tertile_FAT_feb == 3, "high",NA)))) %>% 
  mutate(tertile_FIBER_feb_factor = ifelse(tertile_FIBER_feb == 1, "low",
                                        ifelse(tertile_FIBER_feb == 2, "middle", 
                                               ifelse(tertile_FIBER_feb == 3, "high",NA))))-> feb_metadata 

metadata2 %>% 
  subset(Time == "18-Jul") %>% 
  mutate(tertile_CHO_jul =`CHO.Energy%.overall.july` %>% as.numeric() %>% ntile(.,3)) %>%
  mutate(tertile_PRO_jul = `PRO.Energy%.overall.july` %>% as.numeric() %>% ntile(.,3)) %>% 
  mutate(tertile_FAT_jul = `FAT.Energy%.overall.july` %>% as.numeric() %>% ntile(.,3)) %>% 
  mutate(tertile_FIBER_jul = `total_KCAL_july_dietary_fiber_Energy%` %>% as.numeric() %>% ntile(.,3)) %>% 
  mutate(tertile_CHO_jul_factor = ifelse(tertile_CHO_jul == 1, "low",
                                         ifelse(tertile_CHO_jul == 2, "middle", 
                                                ifelse(tertile_CHO_jul == 3, "high",NA)))) %>% 
  mutate(tertile_PRO_jul_factor = ifelse(tertile_PRO_jul == 1, "low",
                                         ifelse(tertile_PRO_jul == 2, "middle", 
                                                ifelse(tertile_PRO_jul == 3, "high",NA)))) %>% 
  mutate(tertile_FAT_jul_factor = ifelse(tertile_FAT_jul == 1, "low",
                                         ifelse(tertile_FAT_jul == 2, "middle", 
                                                ifelse(tertile_FAT_jul == 3, "high",NA)))) %>% 
  mutate(tertile_FIBER_jul_factor = ifelse(tertile_FIBER_jul == 1, "low",
                                           ifelse(tertile_FIBER_jul == 2, "middle", 
                                                  ifelse(tertile_FIBER_jul == 3, "high",NA))))-> jul_metadata 

sample_data(metadata) -> test
sample_names(test) <- metadata$id
sample_data(iquitos) <- test
sample_data(iquitos)

sample_data(metadata2) -> test
sample_names(test) <- metadata2$id.x
sample_data(iquitos) <- test
sample_data(iquitos) 

############################ 
# Analysis for the February baseline + diet information 
#############################
Feb_phylo <- subset_samples(iquitos, Time == "18-Feb")
#SP <- sample_data(feb_metadata)
#sample_names(SP) <- feb_metadata$id.x
#sample_data(Feb_phylo) <- SP

richness_iquitos_Feb <- Feb_phylo %>% breakaway
plot(richness_iquitos_Feb, physeq = Feb_phylo, color = "tertile_CHO_feb_factor") -> CHO.Energy.breakaway.plot

sample_data(Feb_phylo)
metadata <- Feb_phylo %>% sample_data %>% as_tibble %>% mutate("sample_names" = Feb_phylo %>% sample_names)
Feb_breakaway <- metadata %>%
  left_join(summary(richness_iquitos_Feb),
            by = "sample_names")

png(filename = "CHO_Intake_breakaway.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Feb_breakaway, aes(x = tertile_CHO_feb_factor, y = estimate, color = tertile_CHO_feb_factor, ymin=lower, ymax=upper)) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness in ALL, Tertiles of CHO Intake")
dev.off()

png(filename = "PRO_Intake_breakaway.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Feb_breakaway, aes(x = tertile_PRO_feb_factor, y = estimate, color = tertile_PRO_feb_factor, ymin=lower, ymax=upper)) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness in ALL, Tertiles of PROTEIN Intake")
dev.off()

png(filename = "FAT_Intake_breakaway.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Feb_breakaway, aes(x = tertile_FAT_feb_factor, y = estimate, color = tertile_FAT_feb_factor, ymin=lower, ymax=upper)) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness in ALL, Tertiles of FAT Intake")
dev.off()

png(filename = "Fiber_Intake_breakaway_FEB.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Feb_breakaway, aes(x = FIBER.gr._feb, y = estimate, ymin=lower, ymax=upper)) + 
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness in ALL, Tertiles of FIBER Intake, FEB")
dev.off()


### Richness plots 
png(filename = "Fiber_Intake_breakaway_FEB_Apr28.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Feb_breakaway, aes(x = tertile_FIBER_feb_factor, y = estimate, color = tertile_FIBER_feb_factor, ymin=lower, ymax=upper)) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness in ALL, Tertiles of FIBER Intake, FEB")
dev.off()

# devtools::install_github("adw96/breakaway")
### Look at statistical significance of these differences 
library(breakaway)
Feb_breakaway %>% drop_na(tertile_CHO_feb_factor) -> CHO
CHO_breakaway <- betta(chats = CHO$estimate,
                                         ses = CHO$error,
                                         X = model.matrix(~tertile_CHO_feb_factor + Sex.x + Age.x, data = CHO))
#Estimates Standard Errors p-values
#(Intercept)                  253.25302        10.29388    0.000
#tertile_CHO_feb_factorlow    -57.45614        17.32718    0.001
#tertile_CHO_feb_factormiddle -35.21135        18.09564    0.052
# lower richness in the low CHO group compared to high CHO (where fiber comes in?)

Feb_breakaway %>% drop_na(tertile_PRO_feb_factor) -> PRO
PRO_breakaway <- betta(chats = PRO$estimate,
                       ses = PRO$error,
                       X = model.matrix(~tertile_PRO_feb_factor + Sex.x + Age.x, data = PRO))
#Estimates Standard Errors p-values
#(Intercept)                  191.30105        10.29388    0.000
#tertile_PRO_feb_factorlow     49.00474        17.32323    0.005
#tertile_PRO_feb_factormiddle  40.13741        18.10153    0.027
# Lower protein intake higher diversity 

Feb_breakaway %>% drop_na(tertile_FAT_feb_factor) -> FAT
FAT_breakaway <- betta(chats = FAT$estimate,
                       ses = FAT$error,
                       X = model.matrix(~tertile_FAT_feb_factor + Sex.x + Age.x, data = FAT))
#Estimates Standard Errors p-values
#(Intercept)                  192.02006973      10.2938762    0.000
#tertile_FAT_feb_factorlow     69.24151030      17.3285133    0.000
#tertile_FAT_feb_factormiddle   8.23598235      18.0971034    0.649
#Sex.x                         12.76289433      18.9784341    0.501
#Age.x                         -0.04400159       0.2815033    0.876
# lower fat intake, higher diversity nb h

Feb_breakaway %>% drop_na(tertile_FIBER_feb_factor) -> FIBER
FIBER_breakaway <- betta(chats = FIBER$estimate,
                       ses = FIBER$error,
                       X = model.matrix(~tertile_FIBER_feb_factor + Sex.x + Age.x, data = FIBER))
#Estimates Standard Errors p-values
#(Intercept)                    254.3976060      10.2938763    0.000
#tertile_FIBER_feb_factorlow    -14.9988090      17.3267593    0.387
#tertile_FIBER_feb_factormiddle -56.1394201      18.0949805    0.002
#Sex.x                           15.6044734      18.9784342    0.411
#Age.x                           -0.4727859       0.2815033    0.093


Feb_breakaway %>% drop_na(FIBER.gr._feb) -> FIBER
FIBER_breakaway <- betta(chats = FIBER$estimate,
                         ses = FIBER$error,
                         X = model.matrix(~FIBER.gr._feb + Sex.x + Age.x, data = FIBER))
# Estimates Standard Errors p-values
# (Intercept)   221.5680832      10.2954110    0.000
# FIBER.gr._feb   0.4680409       0.6317803    0.459
# Sex.x          11.0157646      18.9776779    0.562
# Age.x          -0.3482831       0.2815590    0.216

########################################################
#######JULY LOOK AT FIBER # #########
##############################################################
Jul_phylo <- subset_samples(iquitos, Time == "18-Jul")
SP <- sample_data(jul_metadata)
sample_names(SP) <- jul_metadata$id.x
sample_data(Jul_phylo) <- SP

richness_iquitos_Jul <- Jul_phylo %>% breakaway
metadata <- Jul_phylo %>% sample_data %>% as_tibble %>% mutate("sample_names" = Jul_phylo %>% sample_names)
Jul_breakaway <- metadata %>%
  left_join(summary(richness_iquitos_Jul),
            by = "sample_names")
png(filename = "Fiber_Intake_breakaway_JULY.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Jul_breakaway, aes(x = tertile_FIBER_jul_factor, y = estimate, color = tertile_FIBER_jul_factor, ymin=lower, ymax=upper)) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness in ALL, Tertiles of FIBER Intake, JULY")
dev.off()
Jul_breakaway %>% drop_na(tertile_FIBER_jul_factor) -> FIBER_jul
FIBER_breakaway_jul <- betta(chats = FIBER_jul$estimate,
                         ses = FIBER_jul$error,
                         X = model.matrix(~tertile_FIBER_jul_factor + Sex.x + Age.x, data = FIBER_jul))
#Estimates Standard Errors p-values
#(Intercept)                    237.6847724       8.9776453    0.000
#tertile_FIBER_jul_factorlow    -25.7439705      15.1105010    0.088
#tertile_FIBER_jul_factormiddle -37.8316032      15.7834044    0.017
#Sex.x                           15.7852135      16.5535404    0.340
#Age.x                           -0.3148284       0.2455006    0.200

Jul_breakaway %>% drop_na(FIBER.gr._july) -> FIBER_jul
FIBER_breakaway_jul <- betta(chats = FIBER_jul$estimate,
                             ses = FIBER_jul$error,
                             X = model.matrix(~FIBER.gr._july + Sex.x + Age.x, data = FIBER_jul))
# Estimates Standard Errors p-values
# (Intercept)    215.4503799       8.9749413    0.000
# FIBER.gr._july   0.2090292       0.2838060    0.461
# Sex.x           14.3304515      16.5491240    0.387
# Age.x           -0.4555650       0.2454348    0.063

# shows a trend towards higher fiber --> higher species richness but not significant

##################################################################
##### Diversity over Time comparing Increase in Diet vs. No Increase/Decrease 
##############################################################
# Need Time*Covariate of Interest + adjustment for confounders 
# need to create a new variable that subtracts 
# "CHO.Energy%.overall.july" "PRO.Energy%.overall.july" "FAT.Energy%.overall.july"
# "CHO.Energy%.overall.feb" "PRO.Energy%.overall.feb" "FAT.Energy%.overall.feb" 
metadata2 %>% 
  mutate(CHO_diff = `CHO.Energy%.overall.july` - `CHO.Energy%.overall.feb`, 
         PRO_diff = `PRO.Energy%.overall.july` -`PRO.Energy%.overall.feb`, 
         FAT_diff = `FAT.Energy%.overall.july` - `FAT.Energy%.overall.feb`, 
         FIBER_diff = `total_KCAL_july_dietary_fiber_Energy%` - `total_KCAL_feb_dietary_fiber_Energy%`) %>% 
  mutate(CHO_increase = ifelse(CHO_diff > 0, 1, 0), 
         PRO_increase = ifelse(PRO_diff > 0 , 1, 0), 
         FAT_increase = ifelse(FAT_diff > 0, 1, 0), 
         FIBER_increase = ifelse(FIBER_diff > 0, 1,0)) -> metadata3

metadata2 %>% 
  mutate(FIBER_diff = `FIBER(gr)_july` - `FIBER(gr)_feb`, 
         FIBER_increase = ifelse(FIBER_diff > 0, 1,0)) -> metadata3

sample_data(metadata3) -> test
sample_names(test) <- metadata3$id.x
sample_data(iquitos) <- test
sample_data(iquitos)
richness_iquitos_all <- iquitos %>% breakaway

sample_data(iquitos)
metadata <- iquitos %>% sample_data %>% as_tibble %>% mutate("sample_names" = iquitos %>% sample_names)
Iquitos_breakaway <- metadata %>%
  left_join(summary(richness_iquitos_all),
            by = "sample_names")
png(filename = "CHO_Increase_breakaway.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Iquitos_breakaway, aes(x = Time, y = estimate, color = as.factor(CHO_increase), group = person, ymin=lower, ymax=upper)) + 
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .25) + 
  geom_line(aes(linetype=as.factor(CHO_increase)))+
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness between Timepoints, CHO-Increase")
dev.off()

Iquitos_breakaway %>% filter(!is.na(CHO_increase)) -> CHO_binary
time_fixed_person_random_CHO_Increase <- betta_random(chats = CHO_binary$estimate,
                                         ses = CHO_binary$error,
                                         X = model.matrix(~Time*CHO_increase + Age.y, data = CHO_binary),
                                         groups=CHO_binary$person)
CHO_binary %>%
  group_by(Time,CHO_increase) %>%
  summarise(mean(estimate), sd(estimate))
# So it looks like on average everyone is decreasing diversity, but the group that is decreasing less quickly is the CHO people 
# Time   CHO_increase `mean(estimate)` `sd(estimate)`
# <fct>         <dbl>            <dbl>          <dbl>
#   1 18-Feb            0             221.           61.2
# 2 18-Feb            1             222.           60.7
# 3 18-Jul            0             204.           39.5
# 4 18-Jul            1             219.           61.9

#### PRO_increase
png(filename = "PRO_Increase_breakaway.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Iquitos_breakaway, aes(x = Time, y = estimate, color = as.factor(PRO_increase), group = person, ymin=lower, ymax=upper)) + 
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .25) + 
  geom_line(aes(linetype=as.factor(PRO_increase)))+
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness between Timepoints, PRO-Increase")
dev.off()

Iquitos_breakaway %>% filter(!is.na(PRO_increase)) -> PRO_binary
time_fixed_person_random_PRO_Increase <- betta_random(chats = PRO_binary$estimate,
                                                      ses = PRO_binary$error,
                                                      X = model.matrix(~Time*PRO_increase + Age.y, data = PRO_binary),
                                                      groups=PRO_binary$person)
PRO_binary %>%
  group_by(Time,PRO_increase) %>%
  summarise(mean(estimate), sd(estimate)) # those with increase in protein decrease diversity more than those without 
# Time   PRO_increase `mean(estimate)` `sd(estimate)`
# <fct>         <dbl>            <dbl>          <dbl>
#   1 18-Feb            0             205.           59.0
# 2 18-Feb            1             239.           57.6
# 3 18-Jul            0             195.           48.4
# 4 18-Jul            1             228.           52.3

### FAT_increase
png(filename = "FAT_Increase_breakaway.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Iquitos_breakaway, aes(x = Time, y = estimate, color = as.factor(FAT_increase), group = person, ymin=lower, ymax=upper)) + 
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .25) + 
  geom_line(aes(linetype=as.factor(FAT_increase)))+
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness between Timepoints, FAT-Increase")
dev.off()
Iquitos_breakaway %>% filter(!is.na(FAT_increase)) -> FAT_binary
time_fixed_person_random_FAT_Increase <- betta_random(chats = FAT_binary$estimate,
                                                      ses = FAT_binary$error,
                                                      X = model.matrix(~Time*FAT_increase + Age.y, data = FAT_binary),
                                                      groups=FAT_binary$person)
FAT_binary %>%
  group_by(Time,FAT_increase) %>%
  summarise(mean(estimate), sd(estimate)) # those with increase in fat decrease diversity on average more than those without 
# Time   FAT_increase `mean(estimate)` `sd(estimate)`
# <fct>         <dbl>            <dbl>          <dbl>
#   1 18-Feb            0             222.           59.6
# 2 18-Feb            1             221.           62.9
# 3 18-Jul            0             218.           60.7
# 4 18-Jul            1             203.           37.6

#### FIBER INCREASE 
png(filename = "FIBER_Increase_breakaway.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Iquitos_breakaway, aes(x = Time, y = estimate, color = as.factor(FIBER_increase), group = person, ymin=lower, ymax=upper)) + 
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .25) + 
  geom_line(aes(linetype=as.factor(FIBER_increase)))+
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness between Timepoints, FIBER-Increase")
dev.off()
Iquitos_breakaway %>% filter(!is.na(FIBER_increase)) -> FIBER_binary
time_fixed_person_random_FIBER_Increase <- betta_random(chats = FIBER_binary$estimate,
                                                      ses = FIBER_binary$error,
                                                      X = model.matrix(~Time*FIBER_increase + Age.y + Sex.y, data = FIBER_binary),
                                                      groups=FIBER_binary$person)
FIBER_binary %>%
  group_by(Time,FIBER_increase) %>%
  summarise(mean(estimate), sd(estimate))

####################################################
###### FIBER Analysis with Shannon Diversity #######
####################################################
library(DivNet)
iquitos_genus <- iquitos %>%
  tax_glom(taxrank="Genus")
iquitos_dv <- DivNet::divnet(iquitos_genus, X = NULL)
#combined_shannon <- metadata %>%
#  left_join(iquitos_dv$shannon %>% summary,
#            by = "sample_names")
#### Feb 7, 2021 update 
combined_shannon <- metadata3 %>%
  left_join(iquitos_dv$shannon %>% summary,
            by = "sample_names")
combined_shannon

metadata3 <- iquitos %>% sample_data %>% as_tibble %>% mutate("sample_names" = id.x)

combined_shannon <- metadata3 %>%
  left_join(iquitos_dv$shannon %>% summary,
            by = "sample_names")
combined_shannon

### FIBER 
shannon_plot_fiber <- shannon_plot %>% 
                      mutate(fiber = ifelse(Time == "18-Feb",`total_KCAL_feb_dietary_fiber_Energy%`,`total_KCAL_july_dietary_fiber_Energy%`))
png(filename = "Fiber_PrePost.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(shannon_plot_fiber, aes(x = Time, y = fiber, color = Time)) + 
  geom_boxplot(show.legend = FALSE) + 
  geom_point(position=position_dodge(width=0.5),show.legend = FALSE) +
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  theme_bw() + 
  ylab("% Fiber Intake")
dev.off()
# Subset weight data before treatment
before <- subset(shannon_plot_fiber, Time == "18-Feb", fiber,
                 drop = TRUE)
# subset weight data after treatment
after <- subset(shannon_plot_fiber, Time == "18-Jul", fiber,
                 drop = TRUE)
res <- t.test(before, after, paired = TRUE)
res
res <- t.test(fiber ~ Time, data = shannon_plot_fiber, paired = TRUE)
res
res$p.value

combined_shannon %>% 
    filter(Time %in% c("18-Feb","18-Jul")) -> shannon_plot
png(filename = "Shannon_All_PrePost.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(shannon_plot, aes(x = Time, y = estimate, ymin = lower, ymax = upper, color = Time)) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  theme_bw() + 
  ylab("Shannon Diversity Estimate") + 
  ggtitle("Shannon Diversity Estimates Pre-Post")
dev.off()

png(filename = "Shannon_All_PrePost_boxplot.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(shannon_plot, aes(x = Time, y = estimate, ymin = lower, ymax = upper, color = Time)) + 
  geom_boxplot() + 
  #geom_point(position=position_dodge(width=0.5)) +
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  theme_bw() + 
  ylab("Shannon Diversity Estimate") + 
  ggtitle("Shannon Diversity Estimates Pre-Post")
dev.off()

png(filename = "Shannon_All_PrePost_boxplot_points.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(shannon_plot, aes(x = Time, y = estimate, ymin = lower, ymax = upper, color = Time)) + 
  geom_boxplot() + 
  geom_point(position=position_dodge(width=1)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  theme_bw() + 
  ylab("Shannon Diversity Estimate") + 
  ggtitle("Shannon Diversity Estimates Pre-Post")
dev.off()


##### Microbiome comparison Diarrhea vs. no Diarrhea 
shannon_plot %>% 
  filter(Time == "18-Jul") %>% 
  filter(Diarrhea.x %in% c("TRUE","FALSE")) -> july_shannon 
png(filename = "Shannon_Jul_Diarrhea.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(july_shannon, aes(x = Diarrhea.y, y = estimate, ymin = lower, ymax = upper)) + 
  geom_boxplot() + 
  geom_point(position=position_dodge(width=1)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  theme_bw() + 
  ylab("Shannon Diversity Estimate") + 
  xlab("Diarrhea Status") + 
  scale_x_discrete(labels = c("No Diarrhea","Diarrhea")) + 
  ggtitle("Shannon Diversity Estimates Diarrhea vs. No Diarrhea: July")
dev.off() 
Shannon_Diarrhea_July <- betta(chats = july_shannon$estimate,
                               ses = july_shannon$error,
                               X = model.matrix(~Diarrhea.y + Sex.y + Age.y, data = july_shannon))

shannon_plot %>% 
  filter(Time == "18-Feb") %>% 
  filter(Diarrhea.x %in% c("TRUE","FALSE")) -> february_shannon 
png(filename = "Shannon_Feb_Diarrhea.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(february_shannon, aes(x = Diarrhea.x, y = estimate, ymin = lower, ymax = upper)) + 
  geom_boxplot() + 
  geom_point(position=position_dodge(width=1)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  theme_bw() + 
  ylab("Shannon Diversity Estimate") + 
  xlab("Diarrhea Status") + 
  scale_x_discrete(labels = c("No Diarrhea","Diarrhea")) + 
  ggtitle("Shannon Diversity Estimates Diarrhea vs. No Diarrhea: February")
dev.off() 
Shannon_Diarrhea_Feb <- betta(chats = february_shannon$estimate,
                           ses = february_shannon$error,
                           X = model.matrix(~Diarrhea.x + Sex.x + Age.x, data = february_shannon))
february_shannon$fiber <- february_shannon$`total_KCAL_feb_dietary_fiber_Energy%`
shannon_plot %>%
  filter(Time == "18-Feb") %>% 
    filter(!is.na(FIBER.gr._feb)) -> Fiber_Feb # There are 7 people missing fiber information 
Shannon_FiberContinuous_Feb <- betta(chats = Fiber_Feb$estimate,
                              ses = Fiber_Feb$error,
                              X = model.matrix(~FIBER.gr._feb + Sex.x + Age.x, data = Fiber_Feb))
Shannon_FiberContinuous_Feb$table
# Estimates Standard Errors p-values
# (Intercept)    3.028310158     0.059354740    0.000
# FIBER.gr._feb -0.010371190     0.003642567    0.004
# Sex.x         -0.217970940     0.109448615    0.046
# Age.x          0.003390195     0.001623028    0.037

png(filename = "Shannon_Continuous_Fiber_Feb_Apr28.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Fiber_Feb, aes(x = FIBER.gr._feb, y = estimate, ymin = lower, ymax = upper)) + 
  geom_point(position=position_dodge(width=1)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  theme_bw() + 
  ylab("Shannon Diversity Estimate") + 
  xlab("% Fiber Intake") + 
  ggtitle("Shannon Diversity Estimates Fiber: February")
dev.off()

png(filename = "Shannon_Tertiles_Fiber_Feb.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Feb_shannon, aes(x = as.factor(tertile_FIBER_feb), y = estimate, ymin = lower, ymax = upper, color = as.factor(tertile_FIBER_feb))) + 
  geom_point(position=position_dodge(width=1), show.legend = FALSE) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.06, position=position_dodge(width=0.5), show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  theme_bw() + 
  ylab("Shannon Diversity Estimate") + 
  xlab("Fiber Intake") + 
  scale_x_discrete(labels = c("Low","Middle","High")) + 
  ggtitle("Shannon Diversity Estimates Fiber: February")
dev.off()

combined_shannon %>% filter(!is.na(FIBER_increase)) -> FIBER_binary_shannon

# This is looking at the increase in fiber group vs. no increase over time in their shannon diversity 
time_fixed_person_random_FIBER_Increase_Shannon <- betta_random(chats = FIBER_binary$estimate,
                                                        ses = FIBER_binary$error,
                                                        X = model.matrix(~Time*FIBER_increase + Age.y + Sex.y, data = FIBER_binary_shannon),
                                                        groups=FIBER_binary$person)

# look at the cross-sectional differences in shannon between the tertiles of groups 
FIBER_binary_shannon %>% 
  subset(Time == "18-Feb") %>% 
  mutate(tertile_FIBER_feb = `total_KCAL_feb_dietary_fiber_Energy%` %>% as.numeric() %>% ntile(.,3)) %>% 
  mutate(tertile_FIBER_feb_factor = ifelse(tertile_FIBER_feb == 1, "low",
                                           ifelse(tertile_FIBER_feb == 2, "middle", 
                                                  ifelse(tertile_FIBER_feb == 3, "high",NA))))-> Feb_shannon 
png(filename = "Fiber_Intake_Shannon_FEB.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(Feb_shannon, aes(x = tertile_FIBER_feb_factor, y = estimate, color = tertile_FIBER_feb_factor, ymin=lower, ymax=upper)) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Shannon Diversity in ALL, Tertiles of FIBER Intake, FEB")
dev.off()

Feb_shannon %>% drop_na(tertile_FIBER_feb_factor) -> FIBER_feb_shannon
FIBER_shannon_feb <- betta(chats = FIBER_feb_shannon$estimate,
                           ses = FIBER_feb_shannon$error,
                           X = model.matrix(~tertile_FIBER_feb_factor + Sex.x + Age.x, data = FIBER_feb_shannon))
#Estimates Standard Errors p-values
#(Intercept)                     2.928360864     0.061013011    0.000
#tertile_FIBER_feb_factorlow    -0.048611247     0.102702057    0.636
#tertile_FIBER_feb_factormiddle -0.086099754     0.107272331    0.422
#Sex.x                          -0.150002704     0.112512623    0.182
#Age.x                           0.002350914     0.001668469    0.159


FIBER_binary_shannon %>% 
  subset(Time == "18-Jul") %>% 
  mutate(tertile_FIBER_jul = `total_KCAL_july_dietary_fiber_Energy.` %>% as.numeric() %>% ntile(.,3)) %>% 
  mutate(tertile_FIBER_jul_factor = ifelse(tertile_FIBER_jul == 1, "low",
                                           ifelse(tertile_FIBER_jul == 2, "middle", 
                                                  ifelse(tertile_FIBER_jul == 3, "high",NA))))-> July_shannon 

png(filename = "Fiber_Intake_Shannon_JULY.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(July_shannon, aes(x = tertile_FIBER_jul_factor, y = estimate, color = tertile_FIBER_jul_factor, ymin=lower, ymax=upper)) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position=position_dodge(width=0.5)) + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Shannon Diversity in ALL, Tertiles of FIBER Intake, JULY")
dev.off()
July_shannon %>% drop_na(tertile_FIBER_jul_factor) -> FIBER_jul_shannon
FIBER_shannon_jul <- betta(chats = FIBER_jul_shannon$estimate,
                             ses = FIBER_jul_shannon$error,
                             X = model.matrix(~tertile_FIBER_jul_factor + Sex.x + Age.x, data = FIBER_jul_shannon))
#Estimates Standard Errors p-values
#(Intercept)                     2.91195721     0.055584063    0.000
#tertile_FIBER_jul_factorlow    -0.21764875     0.093568950    0.020
#tertile_FIBER_jul_factormiddle -0.13880191     0.097733697    0.156
#Sex.x                          -0.02803924     0.102465896    0.784
#Age.x                           0.00178367     0.001520276    0.241

#########################################################################
# Bacteroides vs Firmicutues change over time with CLR + paired t-tests 
#########################################################################
iquitos_phylum <- iquitos %>%
  tax_glom(taxrank="Phylum")
# ASV 8 is Firmicutes 
# ASV 10 is Bacteroidetes 
# ASV 22 is Actinobacteria 
# ASV3 is Proteobacteria 
phyloseq::transform_sample_counts(iquitos_phylum, function(x) x / sum(x) ) %>% 
  phyloseq::otu_table() %>% 
  as.data.frame() %>% 
  mutate(id = rownames(.)) -> phylum_RA

iquitos_phylum %>% 
  sample_data() %>% 
  as.data.frame() %>% 
  mutate(id = rownames(.)) -> meta_iquitos

all_wide <- phylum_RA %>% right_join(meta_iquitos, by = "id")
all_wide %>% 
  filter(Time == "18-Feb") -> wide_feb
colMeans(wide_feb[,1:17], na.rm = TRUE) -> means
apply(wide_feb[,1:17],2,sd) -> sd

wide_feb %>% 
  filter(Age.x > 17) -> adults_wide_feb 
colMeans(adults_wide_feb[,1:17], na.rm = TRUE) -> means_adults
apply(adults_wide_feb[,1:17],2,sd) -> sd_adults

wide_feb %>% 
  filter(Age.x < 18) -> children_wide_feb 
colMeans(children_wide_feb[,1:17], na.rm = TRUE) -> means_children
apply(children_wide_feb[,1:17],2,sd) -> sd_children
february_shannon %>% 
  mutate(adults = ifelse(Age.x > 17, 1, 0)) -> february_shannon
Shannon_adults_children_Feb <- breakaway::betta(chats = february_shannon$estimate,
                              ses = february_shannon$error,
                              X = model.matrix(~adults, data = february_shannon))

all_wide %>% 
  filter(Time == "18-Jul") -> wide_jul
colMeans(wide_jul[,1:17], na.rm = TRUE) -> means2

all_long <- phylum_RA %>%
  gather(taxa, abundance,ASV3:ASV3593) %>%
  filter(abundance > 0) %>%
  right_join(meta_iquitos, by= "id")

clr <- function(x) {
  (x %>% log) - (x %>% log %>% mean)
}
all_long %>%
  mutate("clr" = clr(abundance)) %>%
  group_by(taxa, Time) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) -> mean_clr

all_long %>%
  mutate("clr" = clr(abundance)) %>% 
  filter(taxa == "ASV8") -> Firmicutes
t_test_F <- t.test(clr ~ Time, data = Firmicutes, paired = TRUE )

all_long %>%
  mutate("clr" = clr(abundance)) %>% 
  filter(taxa == "ASV10") -> Bacteroides
t_test_B <- t.test(clr ~ Time, data = Bacteroides, paired = TRUE )

all_long %>%
  mutate("clr" = clr(abundance)) %>% 
  filter(taxa == "ASV22") -> Actinobacteria
t_test_A <- t.test(clr ~ Time, data = Actinobacteria, paired = TRUE )


july_clr_Firm <- all_long %>%
  mutate("clr_july" = clr(abundance)) %>%
  filter(taxa == "ASV8" & Time == "18-Jul")

feb_clr_Firm <-  all_long %>%
  mutate("clr_feb" = clr(abundance)) %>%
  filter(taxa == "ASV8" & Time == "18-Feb")
  
july_clr_Firm %>% right_join(feb_clr_Firm, by = "person") -> Firm_analysis 

Firm_paired_t <- t.test(Firm_analysis$clr_july, Firm_analysis$clr_feb, paired = TRUE) 



july_clr_Bac <- all_long %>%
  mutate("clr_july" = clr(abundance)) %>%
  filter(taxa == "ASV10" & Time == "18-Jul")

feb_clr_Bac <-  all_long %>%
  mutate("clr_feb" = clr(abundance)) %>%
  filter(taxa == "ASV10" & Time == "18-Feb")

july_clr_Bac %>% right_join(feb_clr_Bac, by = "person") -> Bac_analysis 

Bac_paired_t <- t.test(Bac_analysis$clr_july, Bac_analysis$clr_feb, paired = TRUE) 


july_clr_Act <- all_long %>%
  mutate("clr_july" = clr(abundance)) %>%
  filter(taxa == "ASV22" & Time == "18-Jul")

feb_clr_Act <-  all_long %>%
  mutate("clr_feb" = clr(abundance)) %>%
  filter(taxa == "ASV22" & Time == "18-Feb")

july_clr_Act %>% right_join(feb_clr_Act, by = "person") -> Act_analysis 

Act_paired_t <- t.test(Act_analysis$clr_july, Act_analysis$clr_feb, paired = TRUE) 
#####################################################

# Let's look at adults first 
adults = subset_samples(iquitos, Age.x > 17)
# It looks like there are 46 samples == 23 adults 
# 20 unique households, 3 have two adults sampled   
children = subset_samples(iquitos, Age.x < 18)
# 36 samples == 18 children 
# 16 unique households, 2 have two children sampled 
Feb_phylo <- subset_samples(iquitos, Time == "18-Feb")
February_all$sample_names <- February_all$id
SP <- sample_data(February_all)
sample_names(SP) <- February_all$id
sample_data(Feb_phylo) <- SP
# yay figured it out. 
### Updating to have entire metadata with change in vegetables 
All_metadata$changeHill <- All_metadata$FGhill.y - All_metadata$FGhill.x
All_metadata$changeHill_bin <- ifelse(All_metadata$changeHill > 0, 1, 0)
All_metadata$changeHome <- All_metadata$FGhome.y - All_metadata$FGhome.x
All_metadata$changeHome_bin <- ifelse(All_metadata$changeHome > 0, 1, 0)
All_metadata$changeMedHill <- All_metadata$HerbMedHill.y - All_metadata$HerbMedHill.x
All_metadata$changeMedHill_bin <- ifelse(All_metadata$changeMedHill > 0, 1, 0)
All_metadata$changeMedHome <- All_metadata$HerbMedHome.y - All_metadata$HerbMedHome.x
All_metadata$changeMedHome_bin <- ifelse(All_metadata$changeMedHome > 0, 1, 0)
All_metadata <- All_metadata[-c(83,84,85),]
sampledata_all <- sample_data(All_metadata)
sample_names(sampledata_all) <- All_metadata$id
sample_data(iquitos) <- sampledata_all

sample_data(iquitos_genus) <- sampledata_all

Feb_genus <- subset_samples(iquitos_genus, Time == "18-Feb")
SP <- sample_data(February_all)
sample_names(SP) <- February_all$id
sample_data(Feb_genus) <- SP # agglomerated at the Genus level 

# plot the breakaway estimates for adults and children 
richness_iquitos_adults <- adults %>% breakaway
plot(richness_iquitos_adults, physeq = adults, color = "Time") -> adults_richness_plot
summary(richness_iquitos_adults) %>% as_tibble #Can make this look nicer in ggplot with the tibble 
adults_metadata <- adults %>% sample_data %>% as_tibble %>% mutate("sample_names" = adults %>% sample_names)
adults_richness <- adults_metadata %>%
  left_join(summary(richness_iquitos_adults),
            by = "sample_names")

# plot for children 
richness_iquitos_children <- children %>% breakaway
plot(richness_iquitos_children, physeq = children, color = "Time") -> children_richness_plot
summary(richness_iquitos_children) %>% as_tibble #Can make this look nicer in ggplot with the tibble 
children_metadata <- children %>% sample_data %>% as_tibble %>% mutate("sample_names" = children %>% sample_names)
children_richness <- children_metadata %>%
  left_join(summary(richness_iquitos_children),
            by = "sample_names")


# These are breakaway estimates for the entire population 
richness_iquitos <- iquitos %>% breakaway
plot(richness_iquitos, physeq = iquitos, color = "Time") -> iquitos_breakaway_richness_plot
sample_data(iquitos)
metadata <- iquitos %>% sample_data %>% as_tibble %>% mutate("sample_names" = iquitos %>% sample_names)
combined_richness <- metadata %>%
  left_join(summary(richness_iquitos),
            by = "sample_names")

combined_richness <- All_metadata %>%
  left_join(summary(richness_iquitos),
            by = "sample_names")

summary(richness_iquitos) %>% as_tibble #Can make this look nicer in ggplot with the tibble 

## For plotting use combined_richness 
library(ggplot2)
# for adults 
png(filename = "adults_richness.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(adults_richness, aes(x = person, y = estimate, color = Time, ymin=lower, ymax=upper)) + 
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position = "dodge") + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness in Adults Before & After")
dev.off()
# for children 
png(filename = "children_richness.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(children_richness, aes(x = person, y = estimate, color = Time, ymin=lower, ymax=upper)) + 
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position = "dodge") + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Species Richness in Children Before & After")
dev.off()

# for everyone 
png(filename = "iquitos_richness.png", width = 20, height = 8, units = 'in', res = 300)
ggplot(combined_richness, aes(x = person, y = estimate, color = Time, ymin=lower, ymax=upper)) + 
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position = "dodge") + 
  theme(axis.text.x = element_text(angle = 65))
dev.off()
# Let's do some inferential modeling. 
# First pull out sample data into tibble 
metadata <- iquitos %>% sample_data %>% as_tibble %>% mutate("sample_names" = iquitos %>% sample_names)
readRDS("/home/gradstudent/ptrinh88/February_all.rds") -> February_all
read.table("/home/gradstudent/ptrinh88/February_all.csv", sep = ",", header = TRUE) -> February_all
read.table("/home/gradstudent/ptrinh88/All2.csv", sep = ",", header = TRUE) -> All_metadata

February_all$sample_names <- February_all$id

combined_richness <- metadata %>%
  left_join(summary(richness_iquitos),
            by = "sample_names")

time_fixed_person_random <- betta_random(chats = combined_richness$estimate,
                                       ses = combined_richness$error,
                                       X = model.matrix(~Time * Sex.x, data = combined_richness),
                                       groups=combined_richness$person)
time_fixed_person_random$table

combined_richness$adult <- ifelse(combined_richness$Age.x > 17, "Adult", "Child")
filter(combined_richness, Time == "18-Feb") -> Feb 
adult_fixed_household_random <- betta_random(chats = Feb$estimate,
                                         ses = Feb$error,
                                         X = model.matrix(~adult, data = Feb),
                                         groups=Feb$Household)
adult_fixed_household_random$table

adult_fixed <- betta(chats = Feb$estimate,
                       ses = Feb$error,
                         X = model.matrix(~adult, data = Feb)
                                             )
adult_fixed$table
# So here I'm fitting betta with person as a random effect and time as my fixed effect. 
# basically what's happening here is that there is a decrease in richness of about 5 taxa from July to February p = 0.092 
# am I fitting a LMM here? 
iquitos_phylum <- iquitos %>%
  tax_glom(taxrank="Phylum")
# What if we want to look at Shannon or another metric. Let's use DivNet 
iquitos_genus <- iquitos %>%
  tax_glom(taxrank="Genus")
iquitos_divnet <- DivNet::divnet(iquitos_genus, X = NULL)
# iquitos_divnet is object with diversity estimates 
combined_shannon <- metadata %>%
  left_join(iquitos_divnet$shannon %>% summary,
            by = "sample_names")
All_metadata$sample_names <- All_metadata$id
combined_shannon <- All_metadata %>%
  left_join(iquitos_divnet$shannon %>% summary,
            by = "sample_names")

# update metadata with garden and animal information 
# Let's make new variables for FGhill, FGhome, HerbMedHill, HerbMedHome 
February_all$FGhill_bin <- ifelse(February_all$FGhill %in% c("0","."), 0, 1)
February_all$FGhome_bin <- ifelse(February_all$FGhome == 0, 0, 1)
February_all$HerbMedHill_bin <- ifelse(February_all$HerbMedHill == 0, 0, 1)
February_all$HerbMedHome_bin <- ifelse(February_all$HerbMedHome == 0, 0, 1)
February_all[is.na(February_all)]<-0
February_all$Ani_tot <- February_all$Dogs + February_all$Cats + February_all$Chickens + February_all$Pigs + February_all$Parrots
February_all$AnyAni <- ifelse(February_all$Ani_tot == 0, 0, 1)
# Fixed the Animal Variables Files! 

combined_shannon <- February_all %>%
  left_join(iquitos_divnet$shannon %>% summary,
            by = "sample_names")
# Look at differences in shannon between FGhill, FGhome, HerbMedHill, HerbMedHome 
### FGhill-adults
adults_FGhill_shannon <- betta(chats = adults$estimate,
                                                        ses = adults$error,
                                                        X = model.matrix(~FGhill_bin, data = adults))
adults_FGhill_shannon$table
## FG home-adults
adults_FGhome_shannon <- betta(chats = adults$estimate,
                               ses = adults$error,
                               X = model.matrix(~FGhome_bin, data = adults))
adults_FGhome_shannon$table
##HerbMedHill-adults
adults_HerbMedHill_shannon <- betta(chats = adults$estimate,
                               ses = adults$error,
                               X = model.matrix(~HerbMedHill_bin, data = adults))
adults_HerbMedHill_shannon$table
##HerbMedHome-adults
adults_HerbMedHome_shannon <- betta(chats = adults$estimate,
                                    ses = adults$error,
                                    X = model.matrix(~HerbMedHome_bin, data = adults))
adults_HerbMedHome_shannon$table
#######################################################
### FGhill-children
children_FGhill_shannon <- betta(chats = children$estimate,
                               ses = children$error,
                               X = model.matrix(~FGhill_bin, data = children))
children_FGhill_shannon$table
## FG home-children
children_FGhome_shannon <- betta(chats = children$estimate,
                               ses = children$error,
                               X = model.matrix(~FGhome_bin, data = children))
children_FGhome_shannon$table
##HerbMedHill-children
children_HerbMedHill_shannon <- betta(chats = children$estimate,
                                    ses = children$error,
                                    X = model.matrix(~HerbMedHill_bin, data = children))
children_HerbMedHill_shannon$table
##HerbMedHome-children
children_HerbMedHome_shannon <- betta(chats = children$estimate,
                                    ses = children$error,
                                    X = model.matrix(~HerbMedHome_bin, data = children))
children_HerbMedHome_shannon$table
  
###############################################################
#Let's look at any animals vs no animals exposure.
# In Adults
adults_animal_shannon <- betta(chats = adults$estimate,
                               ses = adults$error,
                               X = model.matrix(~AnyAni, data = adults))
adults_animal_shannon$table
# In Children 
children_animal_shannon <- betta(chats = children$estimate,
                               ses = children$error,
                               X = model.matrix(~AnyAni, data = children))
children_animal_shannon$table

combined_shannon
combined_shannon$adult <- ifelse(combined_shannon$Age.x > 17, 1, 0)
# Is there a difference between two time points for all? 
time_fixed_person_random_shannon <- betta_random(chats = combined_shannon$estimate,
                                       ses = combined_shannon$error,
                                       X = model.matrix(~Time + adult, data = combined_shannon),
                                       groups=combined_shannon$person)
time_fixed_person_random_shannon$table
#################################################
##### Change in microbiome models #############
time_fixed_FGHill_shannon <- betta_random(chats = combined_shannon$estimate,
                                                 ses = combined_shannon$error,
                                                 X = model.matrix(~ Time * changeHill_bin + Age.x*Time, data = combined_shannon),
                                                 groups=combined_shannon$person)
time_fixed_FGHill_shannon$table

time_fixed_FGHome_shannon <- betta_random(chats = combined_shannon$estimate,
                                          ses = combined_shannon$error,
                                          X = model.matrix(~ Time * changeHome_bin + Age.x*Time, data = combined_shannon),
                                          groups=combined_shannon$person)
time_fixed_FGHome_shannon$table

time_fixed_MedHill_shannon <- betta_random(chats = combined_shannon$estimate,
                                          ses = combined_shannon$error,
                                          X = model.matrix(~ Time * changeMedHill_bin + Age.x*Time, data = combined_shannon),
                                          groups=combined_shannon$person)
time_fixed_MedHill_shannon$table

time_fixed_MedHome_shannon <- betta_random(chats = combined_shannon$estimate,
                                           ses = combined_shannon$error,
                                           X = model.matrix(~ Time * changeMedHome_bin + Age.x*Time, data = combined_shannon),
                                           groups=combined_shannon$person)
time_fixed_MedHome_shannon$table
#############################################################

# appears to be a difference let's break it out by group 
# -	Is there a difference between two time points within adults? 
filter(combined_shannon, adult == 1) -> adults
adults_time_fixed_person_random_shannon <- betta_random(chats = adults$estimate,
                                                 ses = adults$error,
                                                 X = model.matrix(~Time, data = adults),
                                                 groups=adults$person)
time_fixed_person_random_shannon$table
# -	Is there a difference between two time points within children? 
filter(combined_shannon, adult == 0) -> children 
children_time_fixed_person_random_shannon <- betta_random(chats = children$estimate,
                                                        ses = children$error,
                                                        X = model.matrix(~Time, data = children),
                                                        groups=children$person)

# Well here based on Shannon we're saying that there is a significantly lower shannon index 
# Is there a difference in Shannon diversity between adults vs. children accounting for household clustering?
filter(combined_shannon, Time == "18-Feb") -> baseline_all 
adultvchildren_shannon <- betta(chats = baseline_all$estimate,
                                                 ses = baseline_all$error,
                                                 X = model.matrix(~adult, data = baseline_all))
adultvchildren_shannon$table
###############################################################
# Difference by diarrhea in the overall diversity 
baseline_all$diarrhea <- ifelse(baseline_all$Diarrhea.x==TRUE, 1, 0)
baseline_all %>% filter(!person == "59112B") -> baseline

diarrhea_shannon <- betta(chats = baseline$estimate,
                                ses = baseline$error,
                                X = model.matrix(~Diarrhea.x, data = baseline))
diarrhea_shannon$table

filter(baseline, adult == 1) -> baseline_adults 
adults_diarrhea_shannon <- betta(chats = baseline_adults$estimate,
                          ses = baseline_adults$error,
                          X = model.matrix(~Diarrhea.x, data = baseline_adults))
adults_diarrhea_shannon$table

filter(baseline, adult == 0) -> baseline_children
children_diarrhea_shannon <- betta(chats = baseline_children$estimate,
                                 ses = baseline_children$error,
                                 X = model.matrix(~Diarrhea.x, data = baseline_children))
children_diarrhea_shannon$table
##############################################################
adults_genus = subset_samples(iquitos_genus, Age.x > 17)
children_genus = subset_samples(iquitos_genus, Age.x < 18)

adults_divnet <- DivNet::divnet(adults_genus, X = NULL)
children_divnet <- DivNet::divnet(children_genus, X = NULL)

time_fixed_person_random_shannon <- betta_random(chats = combined_shannon$estimate,
                                                 ses = combined_shannon$error,
                                                 X = model.matrix(~Time, data = combined_shannon),
                                                 groups=combined_shannon$person)
adult_shannon <- filter(combined_shannon, Age.x>17)
adult_shannon <- adults_metadata %>%
  left_join(adults_divnet$shannon %>% summary,
            by = "sample_names")
children_shannon <- filter(combined_shannon, Age.x<18)
children_shannon <- children_metadata %>%
  left_join(children_divnet$shannon %>% summary,
            by = "sample_names")
# Shannon diversity plots for children and adults from DivNet 
# for adults 
png(filename = "adults_shannon.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(adult_shannon, aes(x = person, y = estimate, color = Time, ymin=lower, ymax=upper)) + 
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position = "dodge") + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Shannon Diversity in Adults Before & After")
dev.off()
# for children 
png(filename = "children_shannon.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(children_shannon, aes(x = person, y = estimate, color = Time, ymin=lower, ymax=upper)) + 
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5, position = "dodge") + 
  theme(axis.text.x = element_text(angle = 65)) + 
  ggtitle("Shannon Diversity in Children Before & After")
dev.off()


###########################
# Let's look at Simpson # 
###########################
combined_simpson <- metadata %>%
  left_join(iquitos_divnet$simpson %>% summary,
            by = "sample_names")
time_fixed_person_random_simpson <- betta_random(chats = combined_simpson$estimate,
                                                 ses = combined_simpson$error,
                                                 X = model.matrix(~Time, data = combined_simpson),
                                                 groups=combined_simpson$person)
time_fixed_person_random_simpson$table
# Based on simpson there's no difference between the simpson diversity between July and February p - 0.648

###############################################################################
#### Using Corncob for Genus-level differential abundance and variability #####
###############################################################################

adults_genus = subset_samples(Feb_genus, Age.x >= 18)
children_genus = subset_samples(Feb_genus, Age.x < 18)

adults_genus = subset_samples(iquitos_genus, Age.x >= 18)
children_genus = subset_samples(iquitos_genus, Age.x < 18)

library(corncob)
### Change over time only 
# Jointly testing for differential abundance and variability across Time, controlling for Age on abundance and dispersion 
set.seed(1)
ChangeOverTime <- differentialTest(formula = ~ Time + Age.x, 
                                   phi.formula = ~ Time + Age.x, 
                                   formula_null = ~ Age.x, 
                                   phi.formula_null = ~ Time + Age.x, 
                                   test = "Wald", boot = FALSE, 
                                   data = iquitos_genus, 
                                   fdr_cutoff = 0.05)

ChangeOverTime_personfixed <- differentialTest(formula = ~ Time + person, 
                                   phi.formula = ~ Time + person, 
                                   formula_null = ~ person, 
                                   phi.formula_null = ~ Time + person, 
                                   test = "Wald", boot = FALSE, 
                                   data = iquitos_genus, 
                                   fdr_cutoff = 0.05)

ChangeOverTime$significant_taxa
otu_to_taxonomy(OTU = ChangeOverTime$significant_taxa, data = iquitos_genus)
ChangeOverTime$significant_models

# Testing for differential abundance across Time controlling for Age on abundance 
set.seed(1)
ChangeOverTime2 <- differentialTest(formula = ~ Time + Age.x, 
                                   phi.formula = ~ 1, 
                                   formula_null = ~ Age.x, 
                                   phi.formula_null = ~ 1, 
                                   test = "Wald", boot = FALSE, 
                                   data = iquitos_genus, 
                                   fdr_cutoff = 0.05)

ChangeOverTime2$significant_taxa
otu_to_taxonomy(OTU = ChangeOverTime2$significant_taxa, data = iquitos_genus)
ChangeOverTime2$significant_models

### Looking at changes in taxa over time by these factors, stratified by age 
set.seed(1)
FGhill_corncob <- differentialTest(formula = ~ Time + changeHill_bin + Time*changeHill_bin + Time*Age.x, 
                                         phi.formula = ~ Time + changeHill_bin + Time*changeHill_bin + Time*Age.x, 
                                         formula_null = ~ Time + changeHill_bin + Time*Age.x, 
                                         phi.formula_null = ~ Time + changeHill_bin + Time*changeHill_bin + Time*Age.x, 
                                         test = "Wald", boot = FALSE, 
                                         data = iquitos_genus, 
                                         fdr_cutoff = 0.05)
FGhill_corncob$significant_taxa
otu_to_taxonomy(OTU = FGhill_corncob$significant_taxa, data = iquitos_genus)
FGhill_corncob$significant_models

set.seed(1)
FGhome_corncob <- differentialTest(formula = ~ Time + changeHome_bin + Time*changeHome_bin + Time*Age.x, 
                                   phi.formula = ~ Time + changeHome_bin + Time*changeHome_bin + Time*Age.x, 
                                   formula_null = ~ Time + changeHome_bin + Time*Age.x, 
                                   phi.formula_null = ~ Time + changeHome_bin + Time*changeHome_bin + Time*Age.x, 
                                   test = "Wald", boot = FALSE, 
                                   data = iquitos_genus, 
                                   fdr_cutoff = 0.05)
FGhome_corncob$significant_taxa
otu_to_taxonomy(OTU = FGhome_corncob$significant_taxa, data = iquitos_genus)
FGhome_corncob$significant_models

set.seed(1)
Medhill_corncob <- differentialTest(formula = ~ Time + changeMedHill_bin + Time*changeMedHill_bin + Time*Age.x, 
                                   phi.formula = ~ Time + changeMedHill_bin + Time*changeMedHill_bin + Time*Age.x, 
                                   formula_null = ~ Time + changeMedHill_bin + Time*Age.x, 
                                   phi.formula_null = ~ Time + changeMedHill_bin + Time*changeMedHill_bin + Time*Age.x, 
                                   test = "Wald", boot = FALSE, 
                                   data = iquitos_genus, 
                                   fdr_cutoff = 0.05)
Medhill_corncob$significant_taxa
otu_to_taxonomy(OTU = FGhill_corncob$significant_taxa, data = iquitos_genus)
Medhill_corncob$significant_models

set.seed(1)
Medhome_corncob <- differentialTest(formula = ~ Time + changeMedHome_bin + Time*changeMedHome_bin + Time*Age.x, 
                                    phi.formula = ~ Time + changeMedHome_bin + Time*changeMedHome_bin + Time*Age.x, 
                                    formula_null = ~ Time + changeMedHome_bin + Time*Age.x, 
                                    phi.formula_null = ~ Time + changeMedHome_bin + Time*changeMedHome_bin + Time*Age.x, 
                                    test = "Wald", boot = FALSE, 
                                    data = iquitos_genus, 
                                    fdr_cutoff = 0.05)
Medhome_corncob$significant_taxa
otu_to_taxonomy(OTU = Medhome_corncob$significant_taxa, data = iquitos_genus)
Medhome_corncob$significant_models

####### Looking at Diarrhea vs. Non Diarrhea 
set.seed(1)
all_diarrhea_corncob <- differentialTest(formula = ~ Diarrhea.x, 
                                            phi.formula = ~ Diarrhea.x, 
                                            formula_null = ~ 1, 
                                            phi.formula_null = ~ Diarrhea.x, 
                                            test = "Wald", boot = FALSE, 
                                            data = Feb_genus, 
                                            fdr_cutoff = 0.05)
all_diarrhea_corncob$significant_taxa
otu_to_taxonomy(OTU = all_diarrhea_corncob$significant_taxa, data = Feb_genus)
all_diarrhea_corncob$significant_models

### FOR ADULTS 
set.seed(1)
adults_diarrhea_corncob <- differentialTest(formula = ~ Diarrhea.x, 
                                phi.formula = ~ Diarrhea.x, 
                                formula_null = ~ 1, 
                                phi.formula_null = ~ Diarrhea.x, 
                                test = "Wald", boot = FALSE, 
                                data = adults_genus, 
                                fdr_cutoff = 0.05)
adults_diarrhea_corncob$significant_taxa
adults_diarrhea_corncob$significant_models
otu_to_taxonomy(OTU = adults_diarrhea_corncob$significant_taxa, data = adults_genus)

set.seed(1)
adults_dv_corncob <- differentialTest(formula = ~ Diarrhea.x,
                                phi.formula = ~ Diarrhea.x,
                                formula_null = ~ Diarrhea.x,
                                phi.formula_null = ~ 1,
                                data = adults_genus,
                                test = "LRT", boot = FALSE,
                                fdr_cutoff = 0.05)
adults_dv_corncob$significant_taxa
otu_to_taxonomy(OTU = adults_dv_corncob$significant_taxa, data = adults_genus)

### FOR CHILDREN 
set.seed(1)
children_diarrhea_corncob <- differentialTest(formula = ~ Diarrhea.x, 
                                            phi.formula = ~ Diarrhea.x, 
                                            formula_null = ~ 1, 
                                            phi.formula_null = ~ Diarrhea.x, 
                                            test = "Wald", boot = FALSE, 
                                            data = children_genus, 
                                            fdr_cutoff = 0.05)
children_diarrhea_corncob$significant_taxa
otu_to_taxonomy(OTU = children_diarrhea_corncob$significant_taxa, data = Feb_genus)
children_diarrhea_corncob$significant_models

### Amy's analysis

richness_iquitos %>% summary %>% print(n=Inf)
richness_iquitos %>% build_frequency_count_tables

summary(richness_iquitos) %>% as_tibble 
combined_richness %>%
  ggplot(aes(col = Time, x = sample_names, y = estimate)) +
  geom_point()

combined_richness %>%
  group_by(Time) %>%
  summarise(mean(estimate), sd(estimate))

pivot_longer
library(tidyverse)
iquitos_divnet$shannon %>% summary
combined_shannon %>%
  group_by(Time) %>%
  summarise(mean(estimate), sd(estimate))
# log ratios 
# does eating garden food 

library(corncob)
library(phyloseq)
fecal_phylum <- iquitos %>%
  tax_glom("Phylum")

refseq(fecal_phylum) -> sequences

sequences[3]
set.seed(1)
dv_analysis_phylum_test <- differentialTest(formula = ~ Time + person,
                                       phi.formula = ~ Time + person,
                                       formula_null = ~ person,
                                       phi.formula_null = ~ Time + person,
                                       data = iquitos_phylum,
                                       test = "LRT", boot = FALSE,
                                       fdr_cutoff = 0.05)

set.seed(1)
dv_analysis_phylum <- differentialTest(formula = ~ Time + Age.x,
                                phi.formula = ~ Time + Age.x,
                                formula_null = ~ Age.x,
                                phi.formula_null = ~ Time + Age.x,
                                data = iquitos_phylum,
                                test = "LRT", boot = FALSE,
                                fdr_cutoff = 0.05)
dv_analysis_phylum$significant_taxa

otu_to_taxonomy(OTU = dv_analysis_phylum$significant_taxa, data = iquitos_phylum)
dv_analysis_phylum$significant_models

set.seed(1)
da_analysis <- differentialTest(formula = ~ Time, 
                                phi.formula = ~ Time, 
                                formula_null = ~ 1, 
                                phi.formula_null = ~ Time, 
                                test = "Wald", boot = FALSE, 
                                data = iquitos_genus, 
                                fdr_cutoff = 0.05)
da_analysis$significant_taxa

set.seed(1)
dv_analysis <- differentialTest(formula = ~ Time,
                                phi.formula = ~ Time,
                                formula_null = ~ Time,
                                phi.formula_null = ~ 1,
                                data = iquitos_genus,
                                test = "LRT", boot = FALSE,
                                fdr_cutoff = 0.05)
dv_analysis$significant_taxa

otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = iquitos_genus)


fecal_species <- iquitos %>%
  tax_glom("Species")
set.seed(1)
da_analysis_species <- differentialTest(formula = ~ Time, 
                                    phi.formula = ~ Time, 
                                    formula_null = ~ 1, 
                                    phi.formula_null = ~ Time, 
                                    test = "Wald", boot = FALSE, 
                                    data = fecal_species, 
                                    fdr_cutoff = 0.05)
da_analysis_species$significant_taxa

set.seed(1)
da_analysis_ASV <- differentialTest(formula = ~ Time, 
                                phi.formula = ~ Time, 
                                formula_null = ~ 1, 
                                phi.formula_null = ~ Time, 
                                test = "Wald", boot = FALSE, 
                                data = iquitos, 
                                fdr_cutoff = 0.05)
da_analysis_ASV$significant_taxa



### Let's look at Figure 2: Overall makeup of microbiome phyla 

# PCOA Of diarrhea vs. microbiomes in Adults and Children 
# Adults 
unifrac.adults <- phyloseq::distance(adults, method = "wunifrac")
library(vegan)
sampledf <- data.frame(sample_data(fecal_FEB1))
adonis(unifrac.iquitos ~ Diarrhea, data = sampledf) # PERMANOVA
# according to PERMANOVA there is no difference in the centroids between those 
# who have diarrhea and those who do not. P = 0.879
adonis(distance(fecal_FEB, method = "wunifrac") ~ Sex, data = sampledf)

# plot PCOA
ordu = ordinate(fecal_FEB1, "PCoA", "unifrac", weighted=TRUE)
p = plot_ordination(fecal_FEB1, ordu, color="Diarrhea")
p = p + geom_point(size=3) + ggtitle("MDS/PCoA on weighted-UniFrac distance by Diarrhea, Iquitos")
# it doesn't look like there's any different between having diarrhea and not having diarrhea

unfrac.all <- phyloseq::distance(iquitos, method = "wunifrac")

## Figure 2: Stacked barplots by phyla 
fecal_phylum <- iquitos %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  dplyr::filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  dplyr::arrange(Phylum)                                      # Sort data frame alphabetically by phylum

fecal_phylum <- iquitos %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  dplyr::mutate(plot_phyla = ifelse(Abundance > 0.02, Phylum, NA)) %>%                         # Filter out low abundance taxa
  dplyr::arrange(Phylum)  

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)


adults_phylum = filter(fecal_phylum, Age.x >= 18)
# It looks like there are 46 samples == 23 adults 
# 20 unique households, 3 have two adults sampled   
children_phylum = filter(fecal_phylum, Age.x < 18) 

# Plot February timepoint only 

png("FEB_phylum_barplot.png", units="in", width=10, height=4, res=300)
ggplot(FEB_Phylum, aes(x = person, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text.x = element_text(angle = 45)) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Iquitos Bacterial Communities in February") 
dev.off()

png("diarrhea_phylum.png", units="in", width=10, height=4, res=300)
ggplot(fecal_phylum, aes(x = person, y = Abundance, fill = Phylum)) + 
  facet_grid(Diarrhea~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text.x = element_text(angle = 45)) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Iquitos Bacterial Communities Diarrhea vs. Non Diarrhea") 
dev.off()

###### Comparison of microbbiomes between Diarrhea vs. Non-Diarrhea 
fecal_phylum %>% 
  filter(Time == "18-Feb" & Diarrhea.x != "NA") -> Feb_phylum 

# new facet labels for Diarrhea 
diarrhea.labs <- c("Diarrhea", "No Diarrhea")
names(diarrhea.labs) <- c("TRUE", "FALSE")
png("diarrhea_phylum.png", units="in", width=10, height=4, res=300)
ggplot(Feb_phylum, aes(x = person, y = Abundance, fill = Phylum)) + 
  facet_grid(~ Diarrhea.x, labeller = labeller(Diarrhea.x = diarrhea.labs)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Iquitos Bacterial Communities Diarrhea vs. Non Diarrhea") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()

Feb_phylum %>% filter(Diarrhea.x == FALSE) -> NoDiarrhea_plot
NoDiarrhea <- ggplot(NoDiarrhea_plot, aes(x = person, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text.x = element_text(angle = 45)) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("No Diarrhea") + 
  theme_bw() + 
  theme(legend.position = "top")


Feb_phylum %>% filter(Diarrhea.x == TRUE) -> Diarrhea_plot
Diarrhea <- ggplot(Diarrhea_plot, aes(x = person, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text.x = element_text(angle = 45)) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Diarrhea") + 
  theme_bw()

grid.arrange(Diarrhea,NoDiarrhea,ncol = 2)

install.packages("tableone")
vars = "Diarrhea.x"
factorVars <- c("Diarrhea.x")
vars <- c("Diarrhea.x")
metadata$diarrhea_y <- ifelse(metadata$Diarrhea.y == "TRUE", 1, 0)

metadata %>% 
  filter(Time == "18-Jul") -> July 
metadata %>% 
  filter(Time == "18-Feb") -> Feb 
table(July$Diarrhea.y)
table(Feb$Diarrhea.x)
vars <- c("total_KCAL_feb_dietary_fiber_Energy%")
table <- tableone::CreateTableOne(vars = vars, strata = "Diarrhea.x", data = feb_metadata)

feb_metadata %>% filter(Age.x >= 18) -> adults 

table_a <- tableone::CreateTableOne(vars = "total_KCAL_feb_dietary_fiber_Energy%", strata = "Diarrhea.x", data = adults)

feb_metadata %>% filter(Age.x < 18) -> children 
table_c <- tableone::CreateTableOne(vars = "total_KCAL_feb_dietary_fiber_Energy%", strata = "Diarrhea.x", data = children)
table_a2 <- tableone::CreateTableOne(vars = "total_KCAL_feb_dietary_fiber_Energy%", data = adults)
table_c2 <- tableone::CreateTableOne(vars = "total_KCAL_feb_dietary_fiber_Energy%", data = children)
table2 <- tableone::CreateTableOne(vars = "total_KCAL_feb_dietary_fiber_Energy%", data = feb_metadata)
####### LMM analysis with Diarrhea as the outcome 
View(metadata)
if (metadata$Time == "18-Feb"){
  metadata$Diarrhea <- metadata$Diarrhea.x
}else{
  metadata$Diarrhea <- metadata$Diarrhea.y
}

metadata$Diarrhea_bin <- ifelse(metadata$Diarrhea=="TRUE", 1, 0)

library(lme4)
#Model 1: 
m_glmm_1 <- glmer(Diarrhea_bin ~ Time*changeHill_bin + Time*Age.x + (1|person), 
                  data = metadata, family = binomial, nAGQ = 20)

summary(m_glmm_1)

m_glmm_2 <- glmer(Diarrhea_bin ~ Time*changeHome_bin + Time*Age.x + (1|person), 
                  data = metadata, family = binomial, nAGQ = 20)

summary(m_glmm_2)

m_glmm_3 <- glmer(Diarrhea_bin ~ Time*changeMedHill_bin + Time*Age.x + (1|person), 
                  data = metadata, family = binomial, nAGQ = 20)

summary(m_glmm_3)


m_glmm_4 <- glmer(Diarrhea_bin ~ Time*changeMedHome_bin + Time*Age.x + (1|person), 
                  data = metadata, family = binomial, nAGQ = 20)

summary(m_glmm_4)


### LMM with Diarrheas as outcome to diet 
metadata2 %>% 
  mutate(Diarrhea = ifelse(Time == "18-Feb", Diarrhea.x, Diarrhea.y)) -> metadata2

metadata2$Diarrhea_bin <- ifelse(metadata2$Diarrhea=="TRUE", 1, 0)

library(lme4)
#Model 1: 
diarrhea_time <- glmer(Diarrhea_bin ~ Time*Age.x + (1|person), 
                  data = metadata2, family = binomial (link = "logit"), nAGQ = 20)
summary(diarrhea_time)
table_d <- tableone::CreateTableOne(vars = "Diarrhea", strata = "Time", factorVars = "Diarrhea", data = metadata2)

#Model 1: 
metadata3 %>% 
  mutate(Diarrhea = ifelse(Time == "18-Feb", Diarrhea.x, Diarrhea.y)) %>% 
  filter(!is.na(Household)) -> metadata3

metadata3$Diarrhea_bin <- ifelse(metadata3$Diarrhea=="TRUE", 1, 0)

fiber_time <- glmer(Diarrhea_bin ~ FIBER_increase*Time + Age.x + (1|person), 
                       data = metadata3, family = binomial (link = "logit"), nAGQ = 20)
summary(fiber_time)

#######################################################
# Adding in a table about the July characteristics ####
#######################################################

metadata %>% dplyr::filter(Time == "18-Jul") -> July
July %>% dplyr::distinct(Household, .keep_all=TRUE) -> distinct_July

factorVars <- c("FGhill.y", "FGhome.y", "HerbMedHill.y", "HerbMedHome.y")
vars <- c("FGhill.y", "FGhome.y", "HerbMedHill.y", "HerbMedHome.y")
table <- tableone::CreateTableOne(vars = vars, data = July, factorVars=factorVars)


library(DivNet)
library(tidyverse)
data(Lee)
DivNet::divnet(Lee) -> test

