indir = "~/git/social.mixing.ncov/contact.matrix"

### Call "obtain.contact.estimates" to obtain contact estimates.

## Inputs:
# - country: BE (Belgium), DE (Germany), FI (Finland), GB (Great Britain), IT (Italy), LU (Luxembourg), NL (the Netherlands) or PL (Poland)
# - age_range: age range for contact matrix, by default 0-100. 
# - age_bands (OPTIONAL): obtain estimates aggregated by age bands > 1

## Outputs:
# in ./results/
# polymod.tab.bin_country.rda:
  #   part.age: a
  #   cont.age: b
  #   y_ab: total number of contact between participants in a with people in b in contact survey
  #   N_a: Number of participants in a
  #   T_b: Population in b
  #   m_ab: average number of contact between one individual in a and people in b
  #   c_ab: probability that two randomly selected individuals in a and b have contact
  #   c.sd_ab = empirical sd of c_ab
# contact.tab.agg_country.rda (OPTIONAL): Same as above aggregated by age bands
# polymod.mod.bin_GB.rda: INLA model
# in ./figures/
# c.crude_country.pdf: Crude estimate of contact matrix (i.e., Y/U)
# c.smooth_country.pdf: Smooth estimates estimated with INLA by one year age bands
# c.smooth.agg_country.pdf (OPTIONAL): Same as above aggregated by age bands
# sd.c.smooth_country.pdf: Standard deviaton of smooth estimates estimated with Monte Carlo by one year age bands

## Load functions ##
source(file.path(indir, "R", "obtain.contact.estimates.R"))

## Call functions ##
# example
obtain_contact_estimates(country = "NL", age_range = c(0,100), age_bands = 5, indir = indir)
obtain_contact_estimates(country = "FR", age_range = c(0,100), age_bands = 5, indir = indir)

