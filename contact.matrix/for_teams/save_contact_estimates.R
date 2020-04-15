save_contact_estimates = function(countries_CAR, countries_Perm, indir)
{
  countries_CAR = c("BE", "DE", "IT", "GB")
  countries_Perm = c("Austria", "Switzerland", "Denmark", "Spain", "France", "Sweden")
  
  df = NULL
  for(country in countries_CAR){
    load(file.path(indir, "results", paste0("polymod.tab.bin_", country, ".rda")))
    country.name = ifelse(country == "BE", "Belgium", 
                          ifelse(country == "DE", "Germany",
                                 ifelse(country == "IT", "Italy",
                                        ifelse(country == "GB", "United Kingdom", NA))))
    polymod.tab = polymod.tab %>%
      rename(pop = T) %>%
      subset(part.age < 80 & cont.age < 80) %>%
      mutate(country = country.name, 
             method = "CAR") %>%
      select(part.age, cont.age, pop, m, c, country, method)
    df = rbind(df, polymod.tab)
  }
  for(country in countries_Perm){
    load(file.path(indir, "results_perm", paste0("contact.matrix_",  str_replace_all(string=country, pattern=" ", repl=""), ".rda")))
    contact.tab$method = "Perm"
    df = rbind(df,contact.tab)
  }
  write_csv(df, path = file.path(indir, "for_teams", "contact.intensities.estimates_150420.csv"))
}
