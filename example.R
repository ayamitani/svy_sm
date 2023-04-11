### APPLICATION ###

data(arthritis)
head(arthritis)

source(svy_sm.R)

# take the baseline observation for arthritis
baseoa <- arthritis %>%
  group_by(id) %>%
  dplyr::slice(1) %>%
  dplyr::select(-time) %>%
  ungroup()

# select 3 patients from each outcome category (K=5)
odsoa <- baseoa %>%
  na.omit(y) %>%
  dplyr::group_by(y) %>%
  dplyr::mutate(num_rows=n()) %>%
  dplyr::sample_n(10) %>%
  dplyr::ungroup() %>%    
  dplyr::mutate(ods = 1) %>%
  dplyr::select(id, ods)
ods_sample <- left_join(baseoa, odsoa, by = "id") %>%
  dplyr::mutate(ods = replace_na(ods, 0)) %>%
  group_by(y) %>%
  dplyr::mutate(sampwt = mean(ods),
                invsampwt = 1/sampwt) %>%
  ungroup() %>%
  dplyr::filter(ods == 1) %>%
  as.data.frame()

# survey design
odsdes <- svydesign(id = ~1,
                    strata = ~ y, 
                    weights = ~ invsampwt, 
                    data = ods_sample)

# use above function
wsmout <- svy_sm(-as.numeric(y) ~ sex + age + trt + baseline, design = odsdes)
wsmout