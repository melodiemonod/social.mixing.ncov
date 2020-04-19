library(data.table)

demog_WPP <- read.csv("~/git/usa_ifrs/country_inputs.csv")
countries.of.interest = c("Norway", "Austria", "Switzerland", "Denmark", "Spain", "France", "Sweden", "Belgium", "Germany", "Italy", "United Kingdom", "Greece",
                            "Netherlands", "Portugal")

mapping.country = data.table(Matrix = 1:15, country.close = c("United Kingdom", "Belgium", "Germany", "Finland", "Poland", "Luxembourg", "Netherlands", 
                                            "Italy", "France", "Hong Kong", "Russia", "China", "Peru", "India", "Zimbabwe"))
mapping = select(demog_WPP, Country_or_region, Matrix) %>%
  rename(country = Country_or_region) %>%
  merge(mapping.country, by = "Matrix") %>%
  subset(country %in% countries.of.interest) %>%
  select(country, country.close) %>%
  mutate(country.close.abb = ifelse(country.close == "United Kingdom", "UK", 
                                    ifelse(country.close == "Belgium", "BE",
                                           ifelse(country.close == "Germany", "DE",
                                                  ifelse(country.close == "Italy", "IT", 
                                                         ifelse(country.close == "Netherlands", "NL", "FR"))))))
mapping$country = plyr::mapvalues(mapping$country, from = "United Kingdom", to = "United_Kingdom")
mapping$country.close = plyr::mapvalues(mapping$country.close, from = "United Kingdom", to = "United_Kingdom")
write.csv(mapping, file = "~/git/social.mixing.ncov/contact.matrix/analyses/mapping.csv", row.names = F)
  