# script to generate seasonality profiles for STP and Mali
# using Imperial's umbrella package

library(umbrella)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
rm(list=ls());gc()

spatial_bf <- st_read("./data/bf_spatial/bfregion.shp")
spatial_bf$ISO = "BFA"
spatial_bf$Country <- "Burkina Faso"
spatial_bf<- spatial_bf[, c("geometry", "Country", "ISO")]

spatial_kenya <- st_read("./data/kenya_spatial/ken_admbnda_adm0_iebc_20191031.shp")
spatial_kenya$ISO = "KEN"
spatial_kenya$Country <- "Kenya"
spatial_kenya<- spatial_kenya[, c("geometry", "Country", "ISO")]

start_date <- "2017-01-01"
end_date <- "2019-12-31"

daily_rain_raw <- pull_daily_rainfall(sf = spatial_bf, start_date = start_date, end_date = end_date)
# Process the raw data
daily_rain  <- daily_rain_raw %>%
  # Convert to long and format
  pivot_longer(-c(ISO, Country),
               names_to = "date",
               values_to = "rainfall",
               names_prefix = "X",
               names_pattern = "(.*)_precipitation") %>%
  mutate(date = as.Date(as.character(readr::parse_number(.data$date)), format = "%Y%m%d"),
         year = lubridate::year(.data$date),
         day_of_year = lubridate::yday(.data$date),
         t = lubridate::yday(.data$date) / 365,
         rainfall = as.numeric(rainfall)) %>%
  # Replace missing data with 0
  replace_na(replace = list(rainfall = 0)) %>%
  # Remove any leap year addtional days
  dplyr::filter(.data$day_of_year < 366)

rain <- daily_rain %>%
  group_by(ISO, Country) %>%
  summarise(
    data = list(data.frame(date, day_of_year, t, rainfall)),
    model = list(fit_fourier(data[[1]]$rainfall, data[[1]]$t)),
    profile = list(fourier_predict(model[[1]]$coefficients, t = 1:365 / 365, floor = model[[1]]$floor))
  )

# Visualise output
rd <- rain %>%
  select(ISO, Country, data) %>%
  tidyr::unnest(cols = c("data"))
pf <- rain %>%
  select(ISO, Country, profile) %>%
  tidyr::unnest(cols = c("profile"))

p1 <- ggplot() +
  geom_point(data = rd, aes(x = t, y = rainfall, col = Country), alpha = 1, size = 0.1) +
  geom_line(data = pf, aes(x = t, y = profile)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  facet_grid( ~ Country)

# write profiles to CSV
prof <- rain$profile[[1]]$profile
write.csv(prof, "./data/profile_bf.csv", col.names=F, row.names=F)


### KENYA ###
daily_rain_raw <- pull_daily_rainfall(sf = spatial_kenya, start_date = start_date, end_date = end_date)
# Process the raw data
daily_rain  <- daily_rain_raw %>%
  # Convert to long and format
  pivot_longer(-c(ISO, Country),
               names_to = "date",
               values_to = "rainfall",
               names_prefix = "X",
               names_pattern = "(.*)_precipitation") %>%
  mutate(date = as.Date(as.character(readr::parse_number(.data$date)), format = "%Y%m%d"),
         year = lubridate::year(.data$date),
         day_of_year = lubridate::yday(.data$date),
         t = lubridate::yday(.data$date) / 365,
         rainfall = as.numeric(rainfall)) %>%
  # Replace missing data with 0
  replace_na(replace = list(rainfall = 0)) %>%
  # Remove any leap year addtional days
  dplyr::filter(.data$day_of_year < 366)

rain <- daily_rain %>%
  group_by(ISO, Country) %>%
  summarise(
    data = list(data.frame(date, day_of_year, t, rainfall)),
    model = list(fit_fourier(data[[1]]$rainfall, data[[1]]$t)),
    profile = list(fourier_predict(model[[1]]$coefficients, t = 1:365 / 365, floor = model[[1]]$floor))
  )

# Visualise output
rd <- rain %>%
  select(ISO, Country, data) %>%
  tidyr::unnest(cols = c("data"))
pf <- rain %>%
  select(ISO, Country, profile) %>%
  tidyr::unnest(cols = c("profile"))

p2 <- ggplot() +
  geom_point(data = rd, aes(x = t, y = rainfall, col = Country), alpha = 1, size = 0.1) +
  geom_line(data = pf, aes(x = t, y = profile)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  facet_grid( ~ Country)

# write profiles to CSV
prof <- rain$profile[[1]]$profile
write.csv(prof, "./data/profile_kenya.csv", col.names=F, row.names=F)

library(gridExtra)
library(grid)

grid.arrange(p1,p2, nrow=2)

