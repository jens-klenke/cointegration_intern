library(haven)

null_dist <- read_dta(here::here('00_data/NullDistr.dta'))

save(null_dist,file = here::here('00_data/null_dist.rda'))

load(here::here('00_data/null_dist.rda'))
