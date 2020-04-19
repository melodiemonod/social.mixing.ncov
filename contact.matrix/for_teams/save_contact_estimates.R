save_contact_estimates = function(countries_CAR, countries_Perm, indir)
{
  
  countries <- c(
    "Denmark",
    "Italy",
    "Germany",
    "Spain",
    "United_Kingdom",
    "France",
    "Norway",
    "Belgium",
    "Austria", 
    "Sweden",
    "Switzerland",
    "Greece",
    "Portugal",
    "Netherlands"
  )
  
  indir = "~/git/social.mixing.ncov/contact.matrix"
  
  contact_tab = list(); contact_tab_agg = list(); contact_tab_long = list()
  
  mapping = read.csv("~/git/social.mixing.ncov/contact.matrix/analyses/mapping.csv")
  
  for(country in countries){
    
    country.name = country
    country.used <- subset(mapping, country == country.name)$country.close
    
    country.used.abb = ifelse(country.used =="Belgium", "BE", 
                          ifelse(country.used == "Germany", "DE",
                                 ifelse(country.used == "Italy", "IT",
                                        ifelse(country.used == "United_Kingdom", "GB", 
                                               ifelse(country.used == "France", "FR", 
                                                      ifelse(country.used == "Netherlands", "NL", NA))))))
    
    load(file.path(indir, "results", paste0("polymod.tab.bin_", country.used.abb, ".rda")))
    polymod.tab_long = polymod.tab %>%
      select(part.age, cont.age, m)
    contact_tab_long[[country.name]] = polymod.tab_long
    polymod.tab = polymod.tab %>%
      reshape2::dcast(part.age ~ cont.age, value.var = "m")
    polymod.tab = polymod.tab[,-1]
    contact_tab[[country.name]] = polymod.tab

    load(file.path(indir, "results", paste0("contact.tab.agg_", country.used.abb, ".rda")))
    contact.tab.agg = contact.tab.agg %>%
      reshape2::dcast(part.age.cat ~ cont.age.cat, value.var = "m")
    contact.tab.agg = contact.tab.agg[,-1]
    contact_tab_agg[[country.name]] = contact.tab.agg
  }
  
  save(contact_tab, file = file.path(indir, "for_teams", "contact.intensities.estimates_160420.rda"))
  save(contact_tab_agg, file = file.path(indir, "for_teams", "contact.intensities.estimates.agg_160420.rda"))
  save(contact_tab_long, file = file.path(indir, "for_teams", "contact.intensities.estimates.long_160420.rda"))
}
