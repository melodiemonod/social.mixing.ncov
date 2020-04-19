clean_contact_intensities = function(country, path_to_estimates, indir)
{
  tic = Sys.time()
  
  age = 0:79; n.age = length(age)
  age.cat <- cut(age, breaks = seq(0,80,5), include.lowest = TRUE, right = FALSE); n.age.cat = length(levels(age.cat))
  
  # Load estimate from perm et al.
  country.of.interest = ifelse(country == "United Kingdom", "United Kingdom of Great Britain", country)
  index = ifelse(substring(country.of.interest, 1, 3) < "Mor", 1, 2)
  col_names = ifelse(index == 1, TRUE, FALSE)
  df <- suppressWarnings(suppressMessages(read_excel(file.path(path_to_estimates, paste0("MUestimates_all_locations_", index, ".xlsx")), 
                                                     sheet = country.of.interest, col_names = col_names)))
  
  # Load population counts
  country.of.interest = country
  pop.data.wide <- read.csv(file = file.path(indir, "data", paste0("pop.data.csv"))) %>%
    subset(country == country.of.interest & age < 80)
  
  # Aggregate population numbers  
  pop.data.wide <- cbind(pop.data.wide, age.cat = age.cat)
  pop.data.agg <- aggregate(cbind(pop) ~ age.cat, FUN = sum, data = pop.data.wide)
  
  # estimate by 5 y age bands
  contact.tab.agg = data.table(part.age.cat = rep(levels(age.cat), n.age.cat), cont.age.cat = rep(levels(age.cat), each = n.age.cat), 
                   m = as.vector(as.matrix(df))) %>%
    merge(pop.data.agg, by.x = "cont.age.cat", by.y = "age.cat") %>%
    mutate(c = m/pop) %>%
    select(part.age.cat, cont.age.cat, m, c, pop)
  contact.tab.agg$country = country
  save(contact.tab.agg, file = file.path(indir, "results_perm", paste0("contact.matrix.agg_",  str_replace_all(string=country, pattern=" ", repl=""), ".rda")))
  
  # obtain estimate by 1 y age band
  m.age = matrix(nrow = n.age, ncol = n.age, 0); 
  for(x in (age+1)){
    age.x = x; age_group.x = (floor(age/5)+1)[x]
    m.age[,age.x] = rep((df[,age_group.x]*pop.data.wide$pop[age.x]/pop.data.agg$pop[age_group.x])[,1], each = 5)
  }
  
  contact.tab = data.table(part.age = rep(age, n.age), cont.age = rep(age, each = n.age), 
                   m = as.vector(m.age)) %>%
    merge(pop.data.wide, by.x = "cont.age", by.y = "age") %>%
    mutate(c = m/pop) %>%
    select(part.age, cont.age, m, c, country, pop)
  save(contact.tab, file = file.path(indir, "results_perm", paste0("contact.matrix_",  str_replace_all(string=country, pattern=" ", repl=""), ".rda")))
  
  # print time difference
  toc = Sys.time()
  print(paste("clean contact intensities", "---", "for country", country, "---", round(as.numeric(toc-tic), digits = 4), "seconds"))
}

make_figures_perm = function(country, indir)
{
  tic = Sys.time()
  
  library(RColorBrewer)
  
  age = 0:79; n.age = length(age)
  age.cat <- cut(age, breaks = seq(0,80,5), include.lowest = TRUE, right = FALSE); n.age.cat = length(levels(age.cat))
  
  ## Settings
  cols <- colorRampPalette(brewer.pal(name = "YlOrRd", n = 9))(n = max(age)) # Colors
  euro.levs <- as.vector(outer(c(1, 2, 5), 10^(-3:3)))  
  
  # one y age band
  load(file.path(indir, "results_perm", paste0("contact.matrix_",  str_replace_all(string=country, pattern=" ", repl=""), ".rda")))
  
  ticks <- seq(from = 0, to = 80, by = 10)                              # Tick locations
  z <- log(contact.tab$c)
  z.range <- c(-20, -10)
  
  pdf(file.path(indir, "figures_perm", paste0("c_", str_replace_all(string=country, pattern=" ", repl=""), ".pdf")),width=7,height=7,paper='special')
  plot.new()
  plot.window(
    xlim = c(min(age) - 0.5, max(age) + 0.5),
    ylim = c(min(age) - 0.5, max(age) + 0.5))
  image  (age, age, matrix(        z[0*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
  #contour(age, age, matrix(1e6*exp(z[0*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
  axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  box(lwd = 0.5)
  mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
  mtext("Contacts",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
  dev.off()
  
  z <- contact.tab$m
  z.range.m <- c(0, 8)
  
  pdf(file.path(indir, "figures_perm", paste0("m_", str_replace_all(string=country, pattern=" ", repl=""), ".pdf")),width=7,height=7,paper='special')
  plot.new()
  plot.window(
    xlim = c(min(age) - 0.5, max(age) + 0.5),
    ylim = c(min(age) - 0.5, max(age) + 0.5))
  image  (age, age, matrix(        z[0*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range.m, col = cols, useRaster = TRUE, add = TRUE)
  #contour(age, age, matrix(1e6*exp(z[0*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
  axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  box(lwd = 0.5)
  mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
  mtext("Contacts",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
  dev.off()
  
  #aggregated by 5 y age bands
  load(file.path(indir, "results_perm", paste0("contact.matrix.agg_",  str_replace_all(string=country, pattern=" ", repl=""), ".rda")))
  ##settings
  contact.tab.agg$part.age.cat = factor(contact.tab.agg$part.age.cat, 
                                        levels = unique(age.cat))
  ticks <- sort(unique(as.numeric(contact.tab.agg$part.age.cat)))
  age = ticks; n.age = length(ticks)
  
  z <- log(contact.tab.agg$c)
     
  pdf(file.path(indir, "figures_perm", paste0("c.agg_", country, ".pdf")),width=7,height=7,paper='special')
  plot.new()
  plot.window(
    xlim = c(min(age) - 0.5, max(age) + 0.5),
    ylim = c(min(age) - 0.5, max(age) + 0.5))
  image  (age, age, matrix(        z[0*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
  contour(age, age, matrix(1e6*exp(z[0*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
  axis(side = 1, at = ticks, labels = sort(unique(contact.tab.agg$part.age.cat)), lwd = 0.5)
  axis(side = 2, at = ticks, labels = sort(unique(contact.tab.agg$part.age.cat)), lwd = 0.5)
  box(lwd = 0.5)
  mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
  mtext("Contacts",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
  dev.off()
  
  z <- contact.tab.agg$m
  
  pdf(file.path(indir, "figures_perm", paste0("m.agg_", country, ".pdf")),width=7,height=7,paper='special')
  plot.new()
  plot.window(
    xlim = c(min(age) - 0.5, max(age) + 0.5),
    ylim = c(min(age) - 0.5, max(age) + 0.5))
  image  (age, age, matrix(        z[0*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range.m, col = cols, useRaster = TRUE, add = TRUE)
  contour(age, age, matrix(1e6*exp(z[0*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
  axis(side = 1, at = ticks, labels = sort(unique(contact.tab.agg$part.age.cat)), lwd = 0.5)
  axis(side = 2, at = ticks, labels = sort(unique(contact.tab.agg$part.age.cat)), lwd = 0.5)
  box(lwd = 0.5)
  mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
  mtext("Contacts",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
  dev.off()
  # print time difference
  toc = Sys.time()
  print(paste("make figures", "---", "for country", country, "---", round(as.numeric(toc-tic), digits = 4), "seconds"))
}

obtain_contact_estimates_perm = function(country, path_to_estimates, indir = "~/")
  {
  `%notin%` <- Negate(`%in%`)
  
  possible.countries = c("Austria", "Belgium", "Switzerland", "Germany", "Denmark", "Spain", "Finland", "France", "Italy", "Luxembourg", "Netherlands",
                         "Poland", "Sweden", "United Kingdom")
  if(country %notin% possible.countries) stop(paste("no estimate for this country - country should take one of those values:", paste(possible.countries, collapse=", ")))
  
  clean_contact_intensities(country, path_to_estimates, indir)
  make_figures_perm(country, indir)
}