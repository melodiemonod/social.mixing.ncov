## PREPARE DEMOGRAPHIC DATA

# SOURCE
  # United Nations: Department of Economic and Social Affairs
  # © August 2019 by United Nation

outdir = "~/git/social.mixing.ncov/contact.matrix"
link = "https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/EXCEL_FILES/5_Interpolated/WPP2019_INT_F03_1_POPULATION_BY_AGE_ANNUAL_BOTH_SEXES.xlsx"

library(readxl)
library(httr)
library(data.table)
library(dplyr)

countries.of.interest = c("Austria", "Belgium", "Switzerland", "Germany", "Denmark", "Spain", "Finland", "France", "Italy", "Luxembourg", "Netherlands", "Norway",
                          "Poland", "Sweden", "United Kingdom")

GET(link, write_disk(tf <- tempfile(fileext = ".xlsx")))
df <- read_excel(tf)
demography = df[-(1:11),]
names(demography) <- as.matrix(demography[1, ])
demography <- demography[-1, -c(1:2, 4:7)]

demography = demography %>%
  rename(country = `Region, subregion, country or area *`, year = `Reference date (as of 1 July)`) %>%
  subset(country %in% countries.of.interest & year == 2019) %>%
  reshape2::melt(id.vars = c("country", "year")) %>%
  mutate(age = as.numeric(as.character(variable)),
         pop = as.numeric(as.character(value))*1000) %>%
  select(country, age, pop)
  
demography <- demography[order(demography$country, demography$age),]

write.csv(demography, file.path(outdir, "data", "pop.data.csv"), row.names = F)


