### OBTAIN CONTACT MATRIX ESTIMATE BY 1 Y AGE BANDS ###
library(readxl)
library(httr)
library(data.table)
library(dplyr)

indir = "~/git/social.mixing.ncov/contact.matrix"

countries.of.interest = c("Austria", "Belgium", "Switzerland", "Germany", "Denmark", "Spain", "Finland", "France", "Italy", "Luxembourg", "Netherlands",
                          "Poland", "Sweden", "United Kingdom")


### Call "obtain.contact.estimates.perm" to obtain contact estimates.

## Inputs:
# - country: "Austria", "Belgium", "Switzerland", "Germany", "Denmark", "Spain", "Finland", "France", "Italy", "Luxembourg", "Netherlands", 
# "Poland", "Sweden", "United Kingdom"

## Outputs:
# in ./results_perm/
# contact.matrix_country.rda: with a and b \in {0, ..., 79}
#   part.age: a
#   cont.age: b
#   T_b: Population in b
#   m_ab: average number of contact between one individual in a and people in b
#   c_ab: probability that two randomly selected individuals in a and b have contact
# contact.tab.agg_country.rda: Same as above aggregated by 5 years age bands 
# in ./figures_perm/
# c_country.pdf: Perm estimates estimated by one year age bands
# c.agg_country.pdf: Same as above aggregated by 5 years age bands

## Load functions ##
source(file.path(indir, "R", "obtain.contact.estimates.perm.R"))

## Load estimates from Perm et al. ##
link = "https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=info:doi/10.1371/journal.pcbi.1005697.s002"
GET(link, write_disk(tf <- tempfile(fileext = ".zip")))
unzip(zipfile = tf, exdir = dirname(tf))
path_to_estimates = file.path(dirname(tf), "contact_matrices_152_countries")

for(country in countries.of.interest){
  obtain_contact_estimates_perm(country = country, path_to_estimates = path_to_estimates, indir = indir)
}

