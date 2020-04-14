## PARTS OF THIS CODE IS TAKEN FROM "Contact-patterns" repository BELONGING TO kassteele
## https://github.com/kassteele/Contact-patterns

prepare_contactsurvey_and_demographic_data = function(country, indir)
  {
  tic = Sys.time()

  library(tidyverse)

  country.of.interest = country

  # POLYMOD SURVEY (MOSSONG ET AL. 2008)
  load(file.path(indir, "data", "participants_final.rda"))
  load(file.path(indir, "data", "contacts_final.rda"))
  part.data <- part.data %>%
    subset(country == country.of.interest)
  cont.data <- cont.data %>%
    subset(country == country.of.interest)
  # DEMOGRAPHIC DATA
  country.of.interest = ifelse(country == "BE", "Belgium", 
         ifelse(country == "DE", "Germany", 
                ifelse(country == "FI", "Finland", 
                       ifelse(country == "GB", "United Kingdom", 
                              ifelse(country == "IT", "Italy",
                                     ifelse(country == "LU", "Luxembourg", 
                                            ifelse(country == "NL", "Netherlands",
                                                   ifelse(country == "PL", "Poland", NA))))))))
  pop.data <- read.csv(file = file.path(indir, "data", paste0("pop.data.csv"))) %>%
    subset(country == country.of.interest)

   # Select relevant columns
  part.data  <- subset(part.data, select = c("local_id", "participant_age", "dayofweek"))
  cont.data  <- subset(cont.data, select = c("local_id", "cnt_age_l", "cnt_age_r", "cnt_touch"))

  # save
  save(part.data, cont.data, pop.data, file = file.path(indir, "analyses", paste0("contact.demography_", country, ".rda")))

  # print time difference
  toc = Sys.time()
  print(paste("prepare contact survey and demographic data", "---", round(as.numeric(toc-tic), digits = 4), "seconds"))
}

read_and_clean = function(country, age_range = c(0,100), indir)
  {
  tic = Sys.time()

  # Read data
  load(file.path(indir, "analyses", paste0("contact.demography_", country, ".rda")))

  #
  # Data pre-processing
  #

  # Clean pop.data
  pop.data <- within(pop.data, {
    # Rename pop to w
    T <- pop
    # Rename age to cont.age (needed for merge with polymod data later)
    cont.age <- age
    # Remove old variabels
    rm(age, pop)
  })

  # Clean part.data
  names(part.data) <- gsub(names(part.data), pattern = "_", replacement = "." )
  names(part.data) <- gsub(names(part.data), pattern = "participant", replacement = "part")
  part.data <- within(part.data, {
    # Give dayofweek informative labels and call it part.day
    part.day <- factor(dayofweek, levels = 0:6, labels = c("sunday", "monday", "tuesday", "wednesday", "thursday", "friday", "saturday"))
    # Remove old variabels
    rm(dayofweek)
  })

  # Modify cont.data
  # Replace _ by . and cnt by cont
  names(cont.data) <- gsub(names(cont.data), pattern = "_", replacement = ".")
  names(cont.data) <- gsub(names(cont.data), pattern = "cnt", replacement = "cont")
  cont.data <- within(cont.data, {
    # Make cont.age.l numeric and make NA if onbekend
    cont.age.l <- as.numeric(ifelse(cont.age.l == "onbekend", yes = NA, no = as.character(cont.age.l)))
  })

  #
  # Construct polymod.data for analysis
  #

  # Construct contact age in cont.data
  # Two problems: 1) Contact age is sometimes given as a range
  #               2) Digit preferencing in contact age

  # Set seed for reproducibility

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

  # # 2) Correction for digit preferencing in cont.age
  # # Extract age from cont.data
  age <- cont.data$cont.age
  # age.tab is a table with number of contacts for the non-pref ages
  age.tab <- table(age[age>=20 & !is.element(age, seq(20, 80, 5))])
  # Fit smooth spline through number of contacts for the non-pref ages
  age.smt <- smooth.spline(x = as.numeric(names(age.tab)), y = age.tab, cv = TRUE)
  # For each pref age, calculate the difference between the peak and the spline
  # Then redistribute the difference over the neighbouring ages -2:2
  for (age.i in seq(from = 20, to = 80, by = 5)) {
    ix1 <- age == age.i & !is.na(age)
    n.obs <- sum(ix1)                             # Observed number at age i
    n.exp <- round(predict(age.smt, x = age.i)$y) # Expected number at age i
    ix2 <- sample(1:n.obs, size = n.obs - n.exp)
    correction <- try(sample(-2:2, size = n.obs - n.exp, replace = TRUE))
    age[ix1][ix2] <- age[ix1][ix2] + correction
  }
  # Put it pack in cont.data
  cont.data$cont.age <- age

  # Remove unused variables
  rm(age, age.smt, age.tab, age.i, correction, ix1, ix2, n.exp, n.obs)

  #
  # Construct polymod.data from part.data and cont.data
  #

  # Set age range
  age <- seq(age_range[1], age_range[2], 1)
  n <- length(age)

  # Merge part.data and cont.data into polymod.data
  polymod.data <- merge(part.data, cont.data, all = TRUE)

  # Omit missing age and age > max(age)
  polymod.data <- droplevels(subset(polymod.data,
                                    !(is.na(part.age) | part.age > max(age) | is.na(part.day) |
                                        is.na(cont.age) | cont.age > max(age))))

  # Make age a categorical variable (needed for cross tabulation in next section)
  polymod.data <- within(polymod.data, {
    part.age <- factor(part.age, levels = age)
    cont.age <- factor(cont.age, levels = age)
  })

  # Also in pop.data, keep age <= max(age)
  pop.data <- subset(pop.data, cont.age <= max(age))

  # Reorder pop.data
  pop.data <- with(pop.data, pop.data[order(cont.age), ])

  # Save
  save(polymod.data, file = file.path(indir, "analyses", paste0("polymod.data.bin_", country, ".rda")))
  save(pop.data, file = file.path(indir, "analyses", paste0("pop.data.bin_", country, ".rda")))

  # print time difference
  toc = Sys.time()
  print(paste("clean data", "---", round(as.numeric(toc-tic), digits = 4), "seconds"))
}

cross_tabulate = function(country, indir)
  {
  tic = Sys.time()

  # Read
  load(file.path(indir, "analyses", paste0("polymod.data.bin_", country, ".rda")))
  load(file.path(indir, "analyses", paste0("pop.data.bin_", country, ".rda")))

  #
  # Cross tabulation
  #

  # Cross tabulate number of participants for each combination of age and sex
  # t = number of participants
  polymod.part.tab <- with(subset(polymod.data, subset = !duplicated(local.id)),
                           as.data.frame(table(part.age = part.age),
                                         responseName = "N"))

  # Cross tabulate number of contacts for each combination of age and sex
  # y = number of participants x number of contacts
  polymod.cont.tab <- with(polymod.data,
                           as.data.frame(table(part.age = part.age, cont.age = cont.age),
                                         responseName = "y"))

  # Make age an integer again (was factor in previous section)
  polymod.part.tab <- within(polymod.part.tab, {
    part.age <- as.numeric(as.character(part.age))
  })
  polymod.cont.tab <- within(polymod.cont.tab, {
    part.age <- as.numeric(as.character(part.age))
    cont.age <- as.numeric(as.character(cont.age))
  })

  # Merge polymod.cont.tab, polymod.part.tab and pop.data
  # polymod.tab is the tabulated version of polymod.data
  polymod.tab <- merge(polymod.part.tab, merge(polymod.cont.tab, pop.data, all = TRUE), all = TRUE)

  # Reorder polymod.tab
  polymod.tab <- polymod.tab[c("part.age", "cont.age", "y", "N", "T")]
  polymod.tab <- with(polymod.tab, polymod.tab[order(cont.age, part.age), ])

  # Calculate denominator U = number of participants N x population T
  # Divide by 1e6 for better numerical properties (inflates contact rate c)
  # If N or T = 0, then set U   = 1 and y = NA (record will not contribute to likelihood)
  polymod.tab <- within(polymod.tab, {
    U <- N*T/1e6
    y <- ifelse(U == 0, yes = NA, no = y)
    U <- ifelse(U == 0, yes = 1, no = U)
  })

  #
  # Save result
  #

  # Save
  save(polymod.tab, file = file.path(indir, "results", paste0("polymod.tab.bin_", country, ".rda")))

  # print time difference
  toc = Sys.time()
  print(paste("cross tabulate data", "---", round(as.numeric(toc-tic), digits = 4), "seconds"))
}

estimate_contact_intensities = function(country, indir)
  {
  tic = Sys.time()

  # Load packages
  library(INLA)

  # Source functions
  source(file = file.path(indir, "functions", "construct.recordID.R"))
  source(file = file.path(indir, "functions", "construct.nodeID.R"))
  source(file = file.path(indir, "functions", "construct.Rmat.R"))
  source(file = file.path(indir, "functions", "construct.Dmat.R"))

  #
  # Read data
  #

  # Read
  load(file = file.path(indir, "results", paste0("polymod.tab.bin_", country, ".rda")))

  #
  # Model contact patterns with INLA
  #

  # Get ages and number of age classes
  age <- unique(polymod.tab$part.age)
  n.age <- length(age)

  # Add node and record IDs to polymod.tab.list
  polymod.tab <- within(polymod.tab, {
    record.id <- with(construct.recordID(n = n.age, sex = F), c( rec))
    node.id   <- with(construct.nodeID  (n = n.age, sex = F), c(node))
  })

  # Construct structure matrix R
  ord <- 2
  R <- construct.Rmat(n = n.age, order = ord, sex = F)$R

  # Run model
  polymod.mod <- inla(
    y ~ 1 + f(node.id, model = "generic0", Cmatrix = R,
              # log-Gamma(1, 0.0001) prior on log-precision
              hyper = list(prec = list(prior = "loggamma", param = c(1, 0.0001))),
              rankdef = 3*ord^2, constr = TRUE, diagonal = 0.001),
    E = U,
    family = "nbinomial",
    data = polymod.tab,
    # Normal(0, 0.001) prior on intercept
    control.fixed = list(mean.intercept = 0, prec.intercept = 0.001),
    # Normal(0, 0.001) prior on log-dispersion parameter
    control.family = list(
      hyper = list(theta = list(prior = "gaussian", param = c(0, 0.001)))),
    # Integration strategy eb for faster computation
    control.inla = list(int.strategy = "eb"),
    # Compute linear predictor
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(
      waic = TRUE,    # Compute WAIC
      cpo = TRUE,     # Compute cross-validated predictive measures
      config = TRUE)) # Enable posterior sampling

  # Show summary
  summary(polymod.mod)

  #
  # Post processing
  #

  # Add expected c and m to polymod.tab
  polymod.tab <- within(polymod.tab, {
    # Compute expected linear predictor
    linpred <- polymod.mod$summary.linear.predictor[, "mean"]
    linpred.sd <- polymod.mod$summary.linear.predictor[, "sd"]
    # c = exp(linear predictor) = contact rate
    # Divide c by 1e6 to go back to original scale (deflates contact rate c)
    c <- exp(linpred)/1e6
    # m = w*c = contact intensity
    m <- T*c
    # Remove linpred
    rm(linpred,linpred.sd)
  })

  # Generate samples from approximated joint posterior distribution to obtain the empirical standard deviation
  n.samples <- 1000
  inla.sam <- inla.posterior.sample(n = n.samples, result = polymod.mod)
  # Extract samples for linear predictor and put them in a matrix
  # Size: nrow(polymod.tab) x n.samples
  linpred.sam <- sapply(X = inla.sam, FUN = function(x) x$latent[grepl(pattern = "Predictor", x = rownames(x$latent))])
  polymod.tab$c.sd = apply(exp(linpred.sam)/1e6, 1, sd)


  #
  # Save results
  #
  #suppressWarnings(suppressMessages(save(polymod.mod, file = file.path(indir, "results", paste0("polymod.mod.bin_", country, ".rda")))))
  save(polymod.tab, file = file.path(indir, "results", paste0("polymod.tab.bin_", country, ".rda")))

  # print time difference
  toc = Sys.time()
  print(paste("estimate contact intensities", "---", round(as.numeric(toc-tic), digits = 4), "seconds"))
}

aggregate_contact_intensities = function(country, age_range, age_bands, indir)
  {
  tic = Sys.time()

  load(file = file.path(indir, "results", paste0("polymod.tab.bin_", country, ".rda")))
  load(file = file.path(indir, "analyses", paste0("pop.data.bin_", country, ".rda")))

  # Get ages and number of age classes
  age <- unique(polymod.tab$part.age)
  n.age <- length(age)

  # Define age cateogries
  age.cat <- cut(age, breaks = seq(age_range[1],age_range[2],age_bands), include.lowest = TRUE, right = FALSE)

  # Reshape pop.data in wide format.
  pop.data.wide <- cbind(
    data.frame(age = age),
    matrix(pop.data$T, nrow = n.age, ncol = 1, dimnames = list(NULL, c("T"))))
  pop.data.wide <- within(pop.data.wide, T <- T)

  # Create a dataframe contact.data with n.age^2 rows.
  contact.data <- cbind(
    expand.grid(part.age = age, cont.age = age),
    matrix(polymod.tab$c, nrow = n.age^2, ncol = 1, dimnames = list(NULL, c("c"))),
    matrix(polymod.tab$m, nrow = n.age^2, ncol = 1, dimnames = list(NULL, c("m"))))

  #
  # Aggregate over ages (equation 6.2, suppl mat AOAS paper)
  #

  # Add age categories to contact.data and pop.data
  contact.data <- cbind(contact.data, expand.grid(part.age.cat = age.cat, cont.age.cat = age.cat))
  pop.data.wide <- cbind(pop.data.wide, age.cat = age.cat)

  # Aggregate population numbers
  pop.data.agg <- aggregate(cbind(T) ~ age.cat, FUN = sum, data = pop.data.wide)

  # Aggegrate contact intensities over ages
  record.id.part <- match(x = contact.data$part.age, table = pop.data.wide$age)
  record.id.part.agg <- match(x = contact.data$part.age.cat, table = pop.data.agg$age.cat)
  contact.tab.agg <- within(contact.data, {
    m   <- pop.data.wide[record.id.part, "T" ]*m  /pop.data.agg[record.id.part.agg, "T" ]
  })
  contact.tab.agg <- aggregate(cbind(m) ~ part.age.cat + cont.age.cat, FUN = sum, data = contact.tab.agg)

  # Calculate contact rates
  record.id.cont.agg <- match(x = contact.tab.agg$cont.age.cat, table = pop.data.agg$age.cat)
  contact.tab.agg <- within(contact.tab.agg, {
    T  <- pop.data.agg[record.id.cont.agg, "T" ]
    c   <- m  / T
  })

  # Reorder columns
  contact.tab.agg <- contact.tab.agg[, c("part.age.cat", "cont.age.cat", "T", "m", "c")]

  # save
  save(contact.tab.agg, file = file.path(indir, "results", paste0("contact.tab.agg_", country, ".rda")))

  # print time difference
  toc = Sys.time()
  print(paste("aggregate contact intensities", "---", round(as.numeric(toc-tic), digits = 4), "seconds"))
}

make_figures = function(country, age_range, age_bands, indir)
  {
  tic = Sys.time()

  library(RColorBrewer)

  # load results
  load(file.path(indir, "results", paste0("polymod.tab.bin_", country, ".rda")))

  ## Settings
  cols <- colorRampPalette(brewer.pal(name = "YlOrRd", n = 9))(n = age_range[2]) # Colors
  euro.levs <- as.vector(outer(c(1, 2, 5), 10^(-3:3)))                  # "euro" levels for contour lines
  ticks <- seq(from = age_range[1], to = age_range[2], by = 10)                              # Tick locations

  ## crude estimate
  # Set variable and range
  age <- unique(polymod.tab$part.age)
  n.age <- length(age)

  z <- with(polymod.tab, log(1 + y/(U)))
  z.range <- range(z, na.rm = TRUE)

  pdf(file.path(indir, "figures", paste0("c.crude_", country, ".pdf")),width=7,height=7,paper='special')
  plot.new()
  plot.window(
    xlim = c(min(age) - 0.5, max(age) + 0.5),
    ylim = c(min(age) - 0.5, max(age) + 0.5))
  image(age, age, matrix(z[0*n.age^2 + 1:n.age^2], n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
  axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  box(lwd = 0.5)
  mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
  mtext("Contacts",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
  dev.off()

  ## Smooth estimate
  z <- log(polymod.tab$c)
  z.range <- range(z)

  pdf(file.path(indir, "figures", paste0("c.smooth_", country, ".pdf")),width=7,height=7,paper='special')
  plot.new()
  plot.window(
    xlim = c(min(age) - 0.5, max(age) + 0.5),
    ylim = c(min(age) - 0.5, max(age) + 0.5))
  image  (age, age, matrix(        z[0*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
  contour(age, age, matrix(1e6*exp(z[0*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
  axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  box(lwd = 0.5)
  mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
  mtext("Contacts",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
  dev.off()

  ## Smooth estimate aggregated
  if(!is.null(age_bands)){
    load(file.path(indir, "results", paste0("contact.tab.agg_", country, ".rda")))

    ticks <- sort(unique(as.numeric(contact.tab.agg$part.age.cat)))
    age = ticks; n.age = length(ticks)

    z <- log(contact.tab.agg$c)

    pdf(file.path(indir, "figures", paste0("c.smooth.agg_", country, ".pdf")),width=7,height=7,paper='special')
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
  }

  ## Standard deviation of smooth estimate
  age <- unique(polymod.tab$part.age)
  n.age <- length(age)

  z <- log(polymod.tab$c.sd)
  z.range <- range(z)

  pdf(file.path(indir, "figures", paste0("sd.c.smooth_", country, ".pdf")),width=7,height=7,paper='special')
  plot.new()
  plot.window(
    xlim = c(min(age) - 0.5, max(age) + 0.5),
    ylim = c(min(age) - 0.5, max(age) + 0.5))
  image  (age, age, matrix(        z[0*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
  contour(age, age, matrix(1e6*exp(z[0*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
  axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  box(lwd = 0.5)
  mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
  mtext("Contacts",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
  dev.off()

  # print time difference
  toc = Sys.time()
  print(paste("make figures", "---", round(as.numeric(toc-tic), digits = 4), "seconds"))
}

obtain_contact_estimates = function(country, age_range = c(0,100), age_bands = NULL, indir= "~/")
 {
  `%notin%` <- Negate(`%in%`)

  if(max(age_range) > 100) stop("Max of age range should be lower or equal to 100")
  if(min(age_range) < 0) stop("Min of age range should be higher or equal to 0")
  if(!is.null(age_bands)){
    if(age_bands > 100) stop("Age bands should be lower or equal to 100")
    if(age_bands < 1) stop("Age bands should be higher or equal to 1")
    if(max(age_range) < age_bands) stop("age bands larger than the maximum of age range")
    if(round(age_bands) != age_bands) stop("age_range should be an integer")
  }
  if(country %notin% c("BE", "DE", "FI", "GB", "IT", "LU", "NL", "PL")) stop("country should take one of those values: BE, DE, FI, GB, IT, LU, NL or PL")

  prepare_contactsurvey_and_demographic_data(country, indir)
  read_and_clean(country, age_range, indir)
  cross_tabulate(country, indir)
  estimate_contact_intensities(country, indir)
  if(!is.null(age_bands)) aggregate_contact_intensities(country, age_range, age_bands, indir)
  make_figures(country, age_range, age_bands, indir)
 }
