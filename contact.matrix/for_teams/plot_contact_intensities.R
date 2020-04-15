## For Teams

plot_contact_intensities = function(country, path_to_estimates){
  library(dplyr)
  library(data.table)
  library(RColorBrewer)
  library(stringr)
  
  `%notin%` <- Negate(`%in%`)
  
  possible.countries = c("Austria", "Belgium", "Switzerland", "Germany", "Denmark", "Spain", "France", "Italy", "Sweden", "United Kingdom")
  
  if(country %notin% possible.countries) stop(paste("no estimate for this country - country should take one of those values:", paste(possible.countries, collapse=", ")))
  
  country.of.interest = country
  
  df = read.csv(file.path(path_to_estimates, "contact.intensities.estimates_150420.csv")) %>%
    subset(country == country.of.interest)
  
  age = sort(unique(df$part.age)); n.age = length(age)
  age.cat <- cut(age, breaks = seq(min(age),(max(age+1)),5), include.lowest = TRUE, right = FALSE); n.age.cat = length(levels(age.cat))
  
  ## Settings
  cols <- colorRampPalette(brewer.pal(name = "YlOrRd", n = 9))(n = max(age)) # Colors
  euro.levs <- as.vector(outer(c(1, 2, 5), 10^(-3:3)))  
  
  ticks <- seq(from = 0, to = 80, by = 10)                              # Tick locations
  z <- log(df$c)
  z.range <- c(-20, -10)
  
  pdf(file.path(path_to_estimates, paste0("c_", str_replace_all(string=country, pattern=" ", repl=""), ".pdf")),width=7,height=7,paper='special')
  plot.new()
  plot.window(
    xlim = c(min(age) - 0.5, max(age) + 0.5),
    ylim = c(min(age) - 0.5, max(age) + 0.5))
  image  (age, age, matrix(        z[0*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
  axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  box(lwd = 0.5)
  mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
  mtext("Contacts",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
  dev.off()
}
