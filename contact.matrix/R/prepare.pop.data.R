## PREPARE DEMOGRAPHIC DATA

# COUNTRIES:
## Austria (AT), Belgium (BE), Switzerland (CH), Germany (DE), Denemark (DK), Spain (ES), Finland (FI), France (FR), Italy (IT), Luxembourg (LU)
## Netherlands (NL), Norway (NO), Poland (PL), Sweden (SE), United Kingdom (UK)

# SOURCE
  # `eurostat 2011 Census database`
  # https://ec.europa.eu/eurostat/web/population-and-housing-census/census-data/2011-census

indir = "~/git/social.mixing.ncov/contact.matrix"

demography = read.csv(file.path(indir, "data", "csvoutput_HC55_2020_04_12_18_44.csv"))
demography = demography %>%
  mutate(age = as.numeric(ifelse(AGE == "Y_GE100", 100,
                      ifelse(AGE == "Y_LT1", 0,
                             sub("Y","",AGE))))) %>%
  rename(country = GEO, pop = VALUE) %>%
  select(country, pop, age)
demography <- demography[order(demography$country, demography$age),]

write.csv(demography, file.path(indir, "data", "pop.data.csv"), row.names = F)



### CONTACT SURVEY ###
# Mossong et al. (2008): Belgium, Germany, Finland, Great Britain, Italy, Luxembourg, the Netherlands, Poland
# Stromgren et al. (2012): Sweden
# Beraud et al. (2015): France
