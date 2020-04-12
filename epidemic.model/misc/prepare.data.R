  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(data.table)
  
  # device
  if(1){
    indir = "~/Box\ Sync/2020/social_mixing_ncov"
  }
  
  
  ### UK ###
  age = 0:100
  thresholds = seq(0, 100, 10); thresholds[length(thresholds)] = 101
  age_cat =  c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-100")
  df = data.table(age = age, AGE_BIN = findInterval(age, thresholds, all.inside = T)) %>%
    mutate(part.age = age_cat[AGE_BIN])
  
  # contact
  polymod_participants = read.table(file.path(indir, "data", "participants_final.txt"), header = T)
  polymod_contacts = read.table(file.path(indir, "data", "contacts_final.txt"), header = T)
  contact = polymod_participants %>%
    select(local_id, participant_age, country) %>%
    subset(!is.na(participant_age) & country == "GB") %>%
    merge(df, by.x = "participant_age", by.y = "age") %>%
    merge(polymod_contacts, by = "local_id") %>%
    mutate(agej = ifelse(is.na(cnt_age_mean), mean(c(cnt_age_l, cnt_age_r)), cnt_age_mean)) %>%
    merge(df, by.x = "agej", by.y = "age") %>%
    rename(part.age = part.age.x, cont.age = part.age.y) %>%
    select(part.age, cont.age) %>%
    group_by(part.age, cont.age) 
  contact$part.age = factor(contact$part.age, levels = age_cat)
  contact$cont.age = factor(contact$cont.age, levels = age_cat)
  contact.tab <- with(contact,
                           as.data.frame(table(part.age = part.age, cont.age = cont.age),
                                         responseName = "y"))
  contact_N = polymod_participants %>%
    select(local_id, participant_age, country) %>%
    subset(!is.na(participant_age) & country == "GB") %>%
    merge(df, by.x = "participant_age", by.y = "age") %>%
    group_by(part.age) 
  contact_N$part.age = factor(contact_N$part.age, levels = age_cat)
  contact.tab.N <- with(contact_N,
                      as.data.frame(table(part.age = part.age),
                                    responseName = "N"))
  contact.tab = merge(contact.tab, contact.tab.N, by = "part.age")
  
  #demography
  pop.data <- read.csv(file = file.path(indir, "data", "demo_UK_2008.csv")) %>%
    group_by(age) %>%
    summarise(pop = sum(pop))
  pop.data = rbind(pop.data, data.table(age = 91:100, pop = round(630818/10)))
  demographic = pop.data %>%
    merge(df, by = "age") %>%
    group_by(part.age) %>%
    rename(cont.age = part.age) %>%
    summarise(T = sum(pop))
  demographic$cont.age.midpoint = seq(4.5,94.5,10)
  demographic$cont.age = factor(demographic$cont.age, levels = age_cat)
  contact.tab = merge(contact.tab, demographic, by = "cont.age")
  contact.tab = contact.tab %>%
    mutate(U = ifelse(N != 0,N*T/1e8, 1)) # Divide by 1e8 for better numerical properties (inflates contact rate c)
  save(contact.tab, file = file.path(indir, "analyses", "contact.tab_UK10.rda"))



age = 0:100
thresholds = seq(0, 100, 5); thresholds[length(thresholds)] = 101
age_cat =  c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", 
             "80-84", "85-89", "90-94", "95-100")
df = data.table(age = age, AGE_BIN = findInterval(age, thresholds, all.inside = T)) %>%
  mutate(part.age = age_cat[AGE_BIN])

# contact
polymod_participants = read.table(file.path(indir, "data", "participants_final.txt"), header = T)
polymod_contacts = read.table(file.path(indir, "data", "contacts_final.txt"), header = T)
contact = polymod_participants %>%
  select(local_id, participant_age, country) %>%
  subset(!is.na(participant_age) & country == "GB") %>%
  merge(df, by.x = "participant_age", by.y = "age") %>%
  merge(polymod_contacts, by = "local_id") %>%
  mutate(agej = ifelse(is.na(cnt_age_mean), mean(c(cnt_age_l, cnt_age_r)), cnt_age_mean)) %>%
  merge(df, by.x = "agej", by.y = "age") %>%
  rename(part.age = part.age.x, cont.age = part.age.y) %>%
  select(part.age, cont.age) %>%
  group_by(part.age, cont.age) 
contact$part.age = factor(contact$part.age, levels = age_cat)
contact$cont.age = factor(contact$cont.age, levels = age_cat)
contact.tab <- with(contact,
                    as.data.frame(table(part.age = part.age, cont.age = cont.age),
                                  responseName = "y"))
contact_N = polymod_participants %>%
  select(local_id, participant_age, country) %>%
  subset(!is.na(participant_age) & country == "GB") %>%
  merge(df, by.x = "participant_age", by.y = "age") %>%
  group_by(part.age) 
contact_N$part.age = factor(contact_N$part.age, levels = age_cat)
contact.tab.N <- with(contact_N,
                      as.data.frame(table(part.age = part.age),
                                    responseName = "N"))
contact.tab = merge(contact.tab, contact.tab.N, by = "part.age")

#demography
pop.data <- read.csv(file = file.path(indir, "data", "demo_UK_2008.csv")) %>%
  group_by(age) %>%
  summarise(pop = sum(pop))
pop.data = rbind(pop.data, data.table(age = 91:100, pop = round(630818/10)))
demographic = pop.data %>%
  merge(df, by = "age") %>%
  group_by(part.age) %>%
  rename(cont.age = part.age) %>%
  summarise(T = sum(pop))
demographic$cont.age.midpoint = seq(2,97.5,5)
demographic$cont.age = factor(demographic$cont.age, levels = age_cat)
contact.tab = merge(contact.tab, demographic, by = "cont.age")
contact.tab = contact.tab %>%
  mutate(U = ifelse(N != 0,N*T, 1)) # Divide by 1e6 for better numerical properties (inflates contact rate c)
save(contact.tab, file = file.path(indir, "analyses", "contact.tab_UK5.rda"))







#outbreak
tmp = read.table(file.path(indir, "data", "andre_estimates_21_02.txt"))
outbreak = tmp %>%
  mutate(Z.under5 = V1 + V2, 
            Z.5to24 = V3 + V4, 
            Z.25to64 = V5 + V6, 
            Z.over65 = V7) %>%
  select(Z.under5, Z.5to24, Z.25to64, Z.over65)
save(outbreak, file = file.path(indir, "analyses", "outbreak_UK_H1N1.rda"))




### ncov
## UK 
### 5 YEARS BANDWIDTH
age = 0:100
thresholds = seq(0, 100, 5); thresholds[length(thresholds)] = 101
age_cat =  c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64",
                    "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-100")
age_cat_ref = data.table(AGE_C = age, AGE_BIN = findInterval(age, thresholds, all.inside = T)) %>%
  mutate(AGE_i = age_cat[AGE_BIN])
### 4 AGE GROUPS
part.data <-read.table(file.path(indir, "data", "participants_final.txt"), header = T) %>%
  subset(country == "GB")
cont.data <- read.table(file.path(indir, "data", "contacts_final.txt"), header = T) %>%
  subset(country == "GB")
# Select relevant columns
part.data  <- subset(part.data, select = c("local_id", "participant_age", "dayofweek"))
cont.data  <- subset(cont.data, select = c("local_id", "cnt_age_l", "cnt_age_r", "cnt_touch"))
# Clean part.data
names(part.data) <- gsub(names(part.data), pattern = "_", replacement = "." )
names(part.data) <- gsub(names(part.data), pattern = "participant", replacement = "part")
part.data <- within(part.data, {
  # Give dayofweek informative labels and call it part.day
  part.day <- factor(dayofweek, levels = 0:6, labels = c("sunday", "monday", "tuesday", "wednesday", "thursday", "friday", "saturday"))
  # Remove old variabels
  rm(dayofweek)
})
# Replace _ by . and cnt by cont
names(cont.data) <- gsub(names(cont.data), pattern = "_", replacement = ".")
names(cont.data) <- gsub(names(cont.data), pattern = "cnt", replacement = "cont")
cont.data <- within(cont.data, {
  # Make cont.age.l numeric and make NA if onbekend
  cont.age.l <- as.numeric(ifelse(cont.age.l == "onbekend", yes = NA, no = as.character(cont.age.l)))
})
set.seed(1)
# 1) Take care of the age range
cont.data <- within(cont.data, {
  # Start with empty cont.age
  cont.age <- NA
  # All ok? Then cont.age is cont.age.l
  ix <- (!is.na(cont.age.l) & is.na(cont.age.r)) | (!is.na(cont.age.l) & !is.na(cont.age.r) & cont.age.r<=cont.age.l)
  cont.age[ix] <- cont.age.l[ix]
  # Missing cont.age.l and cont.age.r reported? Then cont.age is cont.age.r
  ix <- is.na(cont.age.l) & !is.na(cont.age.r)
  cont.age[ix] <- cont.age.r[ix]
  # Age given in a range? Then sample uniformly from range
  ix <- !is.na(cont.age.l) & !is.na(cont.age.r) & cont.age.r > cont.age.l
  cont.age[ix] <- mapply(FUN = function(l, r) sample(l:r, size = 1), l = cont.age.l[ix], r = cont.age.r[ix])
  # Remove old variabels
  rm(cont.age.l, cont.age.r, ix)
})
polymod.data <- merge(part.data, cont.data, all = TRUE)
polymod.data <- droplevels(subset(polymod.data,
                                  !(is.na(part.age) | part.age > max(age) | is.na(part.day) |
                                      is.na(cont.age) | cont.age > max(age))))

# Make age a categorical variable (needed for cross tabulation in next section)
polymod.data <- within(polymod.data, {
  part.age <- factor(part.age, levels = age)
  cont.age <- factor(cont.age, levels = age)
})
polymod.data = merge(polymod.data, age_cat_ref, by.x = "part.age", by.y = "AGE_C")
contact$AGE_Cj = cont.data
contact = contact %>%
  merge(df, by.x = "AGE_Cj", by.y = "AGE_C") %>%
  rename(AGE_i = AGE_i.x, AGE_j = AGE_i.y) %>%
  group_by(AGE_i, AGE_j) %>%
  summarise(n_contact = n())
contact_N = polymod_contacts %>%
  select(part_id, part_age, country, month) %>%
  subset(!is.na(part_age) & country == "United Kingdom" & !is.na(month)) %>%
  merge(df, by.x = "part_age", by.y = "AGE_C") %>%
  group_by(AGE_i) %>%
  summarise(N = n()) 
contact$AGE_i = factor(contact$AGE_i, levels = thresholds_age)
contact_N$AGE_i = factor(contact_N$AGE_i, levels = thresholds_age)
save(contact, contact_N, file = file.path(indir, "analyses", "contact_UK_5.rda"))

## Prepare demographic dataset 
pop = read.table(file.path(indir, "data", "age_sizes.txt"))[,1]
demographic = data.table(pop = pop, AGE_C = 1:90) %>%
  merge(df, by = "AGE_C") %>%
  group_by(AGE_i) %>%
  summarise(T = sum(pop)) %>%
  merge(tmp, by = "AGE_i")
tmp = data.table(AGE_i = thresholds_age, 
                 AGE_MIDPOINT = c(mean(c(0,4)), 
                             mean(c(5,9)), 
                             mean(c(10,14)), 
                             mean(c(15,19)), 
                             mean(c(20,24)), 
                             mean(c(25,29)),
                             mean(c(30,34)), 
                             mean(c(35,39)),
                             mean(c(40,44)), 
                             mean(c(45,49)),
                             mean(c(50,54)), 
                             mean(c(55,59)),
                             mean(c(60,64)), 
                             mean(c(65,69)),
                             mean(c(70,74)), 
                             mean(c(75,79)),
                             mean(80,90)))
demographic$AGE_i = factor(demographic$AGE_i, levels = thresholds_age)
    save(demographic, file = file.path(indir, "analyses", "demographic_UK_5.rda"))

## Prepare outbreak dataset
outbreak = read.csv(file.path(indir, "data", "COVID19_2020_open_line_list - outside_Hubei.csv"))
#correct an error in the data
outbreak[which(outbreak$date_death_or_discharge == "02.02.2021"),"date_death_or_discharge"] = "02.02.2020"
df2 = outbreak %>%
  mutate(tt = lubridate::dmy(date_confirmation), 
         AGE_C = as.numeric(as.character(age)), 
         t = lubridate::dmy(date_death_or_discharge)) %>%
  subset(!is.na(AGE_C) & !is.na(tt) & country == "United Kingdom") %>%
  merge(df, by = "AGE_C") %>%
  group_by(t, AGE_BIN) %>%
  summarise(R_new = n()) 
tmp = outbreak %>%
  mutate(t = lubridate::dmy(date_confirmation), 
         AGE_C = as.numeric(as.character(age))) %>%
  subset(!is.na(AGE_C) & !is.na(t) & country == "United Kingdom") #%>%
  select(t, AGE_C) %>%
  merge(df, by = "AGE_C") %>%
  group_by(t, AGE_BIN, AGE_C) %>%
  summarise(I_new = n()) %>%
  merge(df2, by = c("t", "AGE_BIN"), all.x = T) %>%
  mutate_at("R_new", ~ replace(., is.na(.),0)) %>%
  group_by(AGE_BIN, AGE_C) %>%
  mutate(Z = cumsum(I_new - R_new)) %>%
  select(t, AGE_BIN, Z, AGE)
save(outbreak, file = file.path(indir, "analyses", "outbreak_10.rda"))



### CHINA ###
### 5 YEARS BANDWIDTH
thresholds = seq(0, 89, 5); thresholds[length(thresholds)] = 99
thresholds_age =  c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64",
                    "65-69", "70-74", "75-79", "80+")
df = data.table(AGE_C = 0:99, AGE_BIN = findInterval(0:99, thresholds, all.inside = T)) %>%
  mutate(AGE = thresholds_age[AGE_BIN])

## Prepare contact data
contact = read.csv(file.path(indir, "data", "rspb20140268supp2.csv"))
contact = contact %>%
  select("age", "c.age.under5", "c.age.6to19", "c.age.20to64", "c.age.over65") %>%
  mutate(AGE_i = ifelse(age == "0-4", "under5", 
                        ifelse(age == "5-9" | age == "10-14" | age == "15-19", "6to19", 
                               ifelse(age == "65-69" | age == "70-74" | age == "75-79" | age == "80+", "over65", 
                                      "20to64")))) %>%
  rename(AGE = age)
contact$AGE = factor(contact$AGE, levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64",
                                       "65-69", "70-74", "75-79", "80+"))
contact$AGE_i = factor(contact$AGE_i, levels = c("under5", "6to19", "20to64", "over65"))
save(contact, file = file.path(indir, "analyses", "contact_5.rda"))

## Prepare demographic dataset 
demographic = read.csv2(file.path(indir, "data", "UNdata_Export_20200302_164046625.txt"))
demographic = demographic %>%
  subset(Sex == "Both Sexes" & Source.Year == "2012" & Area == "Total" & !is.na(as.numeric(as.character(Age)))) %>%
  mutate(AGE_C = as.numeric(as.character(Age))) %>%
  select(AGE_C, Value) %>%
  merge(df, by = "AGE_C") %>%
  group_by(AGE_BIN, AGE) %>%
  summarise(T = sum(Value), 
            age_midpoint = mean(AGE_C)) 
save(demographic, file = file.path(indir, "analyses", "demographic_5.rda"))


## Prepare outbreak dataset
outbreak = read.csv(file.path(indir, "data", "COVID19_2020_open_line_list - outside_Hubei.csv"))
#correct a typo
outbreak[which(outbreak$date_death_or_discharge == "02.02.2021"),"date_death_or_discharge"] = "02.02.2020"
df2 = outbreak %>%
  mutate(tt = lubridate::dmy(date_confirmation), 
         AGE_C = as.numeric(as.character(age)), 
         t = lubridate::dmy(date_death_or_discharge)) %>%
  subset(!is.na(AGE_C) & !is.na(tt) & country == "China") %>%
  merge(df, by = "AGE_C") %>%
  group_by(t, AGE_BIN) %>%
  summarise(R_new = n()) 
outbreak = outbreak %>%
  mutate(t = lubridate::dmy(date_confirmation), 
         AGE_C = as.numeric(as.character(age))) %>%
  subset(!is.na(AGE_C) & !is.na(t) & country == "China") %>%
  select(t, AGE_C)%>%
  merge(df, by = "AGE_C") %>%
  group_by(t, AGE_BIN, AGE) %>%
  summarise(I_new = n()) %>%
  merge(df2, by = c("t", "AGE_BIN"), all.x = T) %>%
  mutate_at("R_new", ~replace(., is.na(.),0)) %>%
  group_by(AGE_BIN, AGE) %>%
  mutate(Z = cumsum(I_new - R_new)) %>%
  select(t, AGE_BIN, Z, AGE)
save(outbreak, file = file.path(indir, "analyses", "outbreak_5.rda"))


### 10 YEARS BANDWIDTH
thresholds = c(0, 5,seq(10, 59, 10), 60, 65, 70,80,99);# thresholds[length(thresholds)] = 99
thresholds_age =  c("0-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-64", "65-69", "70-79", "80+")
df = data.table(AGE_C = 0:99, AGE_BIN = findInterval(0:99, thresholds, all.inside = T)) %>%
  mutate(AGE = thresholds_age[AGE_BIN])

## Prepare contact data
contact = read.csv(file.path(indir, "data", "rspb20140268supp2.csv"))
contact = contact %>%
  select("age", "c.age.under5", "c.age.6to19", "c.age.20to64", "c.age.over65") %>%
  mutate(AGE_i = ifelse(age == "0-4", "under5", 
                        ifelse(age == "5-9" | age == "10-14" | age == "15-19", "6to19", 
                               ifelse(age == "65-69" | age == "70-74" | age == "75-79" | age == "80+", "over65", 
                                      "20to64"))),
         AGE = ifelse(age == "0-4", "0-4", 
                      ifelse(age == "5-9", "5-9", 
                          ifelse(age == "10-14" | age == "15-19", "10-19", 
                                 ifelse(age == "20-24" | age == "25-29", "20-29", 
                                        ifelse(age == "30-34" | age == "35-39", "30-39", 
                                               ifelse(age == "40-44" | age == "45-49", "40-49",
                                                      ifelse(age == "50-54" | age == "55-59", "50-59", 
                                                             ifelse(age == "60-64", "60-64",
                                                                    ifelse(age == "65-69", "65-69", 
                                                                         ifelse(age == "70-74" | age == "75-79", "70-79", "80+"))))))))))) %>%
  select(AGE, c.age.under5, c.age.6to19, c.age.20to64, c.age.over65, AGE_i)
contact$AGE = factor(contact$AGE, levels = c("0-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-64", "65-69", "70-79", "80+"))
contact$AGE_i = factor(contact$AGE_i, levels = c("under5", "6to19", "20to64", "over65"))
save(contact, file = file.path(indir, "analyses", "contact_10.rda"))

## Prepare demographic dataset 
demographic = read.csv2(file.path(indir, "data", "UNdata_Export_20200302_164046625.txt"))
demographic = demographic %>%
  subset(Sex == "Both Sexes" & Source.Year == "2012" & Area == "Total" & !is.na(as.numeric(as.character(Age)))) %>%
  mutate(AGE_C = as.numeric(as.character(Age))) %>%
  select(AGE_C, Value) %>%
  merge(df, by = "AGE_C") %>%
  group_by(AGE_BIN, AGE) %>%
  summarise(T = sum(Value), 
            AGE_MIDPOINT = mean(AGE_C)) 
save(demographic, file = file.path(indir, "analyses", "demographic_10.rda"))

## Prepare outbreak dataset
outbreak = read.csv(file.path(indir, "data", "COVID19_2020_open_line_list - outside_Hubei.csv"))
#correct an error in the data
outbreak[which(outbreak$date_death_or_discharge == "02.02.2021"),"date_death_or_discharge"] = "02.02.2020"
df2 = outbreak %>%
  mutate(tt = lubridate::dmy(date_confirmation), 
         AGE_C = as.numeric(as.character(age)), 
         t = lubridate::dmy(date_death_or_discharge)) %>%
  subset(!is.na(AGE_C) & !is.na(tt) & country == "China") %>%
  merge(df, by = "AGE_C") %>%
  group_by(t, AGE_BIN) %>%
  summarise(R_new = n()) 
outbreak = outbreak %>%
  mutate(t = lubridate::dmy(date_confirmation), 
         AGE_C = as.numeric(as.character(age))) %>%
  subset(!is.na(AGE_C) & !is.na(t) & country == "China") %>%
  select(t, AGE_C) %>%
  merge(df, by = "AGE_C") %>%
  group_by(t, AGE_BIN, AGE) %>%
  summarise(I_new = n()) %>%
  merge(df2, by = c("t", "AGE_BIN"), all.x = T) %>%
  mutate_at("R_new", ~ replace(., is.na(.),0)) %>%
  group_by(AGE_BIN, AGE) %>%
  mutate(Z = cumsum(I_new - R_new)) %>%
  select(t, AGE_BIN, Z, AGE)
save(outbreak, file = file.path(indir, "analyses", "outbreak_10.rda"))


### need to finish to code this up
### 4 AGE GROUPS
thresholds = c(0, 5, 20, 65, 99);# thresholds[length(thresholds)] = 99
thresholds_age =  c("0-4", "5-19", "20-64", "65+")
df = data.table(AGE_C = 0:99, AGE_BIN = findInterval(0:99, thresholds, all.inside = T)) %>%
  mutate(AGE = thresholds_age[AGE_BIN])

## Prepare contact data
contact = read.csv(file.path(indir, "data", "rspb20140268supp2.csv"))
contact = contact %>%
  select("age", "c.age.under5", "c.age.6to19", "c.age.20to64", "c.age.over65") %>%
  mutate(AGE_i = ifelse(age == "0-4", "under5", 
                        ifelse(age == "5-9" | age == "10-14" | age == "15-19", "6to19", 
                               ifelse(age == "65-69" | age == "70-74" | age == "75-79" | age == "80+", "over65", 
                                      "20to64"))),
         AGE = ifelse(age == "0-4", "0-4", 
                      ifelse(age == "5-9" | age == "10-14" | age == "15-19", "5-19", 
                              ifelse(age == "65-69" | age == "70-74" | age == "75-79" | age == "80+", "65+", "20-64")))) %>%
  select(AGE, c.age.under5, c.age.6to19, c.age.20to64, c.age.over65, AGE_i)
contact$AGE = factor(contact$AGE, levels = c("0-4", "5-19", "20-64", "65+"))
contact$AGE_i = factor(contact$AGE_i, levels = c("under5", "6to19", "20to64", "over65"))
save(contact, file = file.path(indir, "analyses", "contact_10.rda"))

## Prepare demographic dataset 
demographic = read.csv2(file.path(indir, "data", "UNdata_Export_20200302_164046625.txt"))
demographic = demographic %>%
  subset(Sex == "Both Sexes" & Source.Year == "2012" & Area == "Total" & !is.na(as.numeric(as.character(Age)))) %>%
  mutate(AGE_C = as.numeric(as.character(Age))) %>%
  select(AGE_C, Value) %>%
  merge(df, by = "AGE_C") %>%
  group_by(AGE_BIN, AGE) %>%
  summarise(T = sum(Value), 
            AGE_MIDPOINT = mean(AGE_C)) 
save(demographic, file = file.path(indir, "analyses", "demographic_10.rda"))

## Prepare outbreak dataset
outbreak = read.csv(file.path(indir, "data", "COVID19_2020_open_line_list - outside_Hubei.csv"))
#correct an error in the data
outbreak[which(outbreak$date_death_or_discharge == "02.02.2021"),"date_death_or_discharge"] = "02.02.2020"
df2 = outbreak %>%
  mutate(tt = lubridate::dmy(date_confirmation), 
         AGE_C = as.numeric(as.character(age)), 
         t = lubridate::dmy(date_death_or_discharge)) %>%
  subset(!is.na(AGE_C) & !is.na(tt) & country == "China") %>%
  merge(df, by = "AGE_C") %>%
  group_by(t, AGE_BIN) %>%
  summarise(R_new = n()) 
outbreak = outbreak %>%
  mutate(t = lubridate::dmy(date_confirmation), 
         AGE_C = as.numeric(as.character(age))) %>%
  subset(!is.na(AGE_C) & !is.na(t) & country == "China") %>%
  select(t, AGE_C) %>%
  merge(df, by = "AGE_C") %>%
  group_by(t, AGE_BIN, AGE) %>%
  summarise(I_new = n()) %>%
  merge(df2, by = c("t", "AGE_BIN"), all.x = T) %>%
  mutate_at("R_new", ~ replace(., is.na(.),0)) %>%
  group_by(AGE_BIN, AGE) %>%
  mutate(Z = cumsum(I_new - R_new)) %>%
  select(t, AGE_BIN, Z, AGE)
save(outbreak, file = file.path(indir, "analyses", "outbreak_10.rda"))
