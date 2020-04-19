
cont.data <- read.csv("~/git/social.mixing.ncov/contact.matrix/data/2015_Beraud_France_contact_common.csv") %>%
  mutate(cnt_age_l = ifelse(!is.na(cnt_age_exact), cnt_age_exact, cnt_age_est_min)) %>%
  rename(local_id = part_id, cnt_age_r = cnt_age_est_max, cnt_touch = phys_contact)
save(cont.data, file = "~/git/social.mixing.ncov/contact.matrix/data/contacts_final_france.rda")

part.data <- read.csv("~/git/social.mixing.ncov/contact.matrix/data/2015_Beraud_France_participant_common.csv") %>%
  rename(local_id = part_id, participant_age = part_age)
save(part.data, file = "~/git/social.mixing.ncov/contact.matrix/data/participants_final_france.rda")
