library("Rcpp")
library("here")
library("dplyr")

# ----

sourceCpp(here("Study 1", "data_proc", "frameGrid.cpp"))

generateFrames(74, 7074)

# ----

sourceCpp(here("Study 1", "data_proc", "upsample_events.cpp"))

df <- data.frame(onset_time = c(1.5, 3.2, 5.7),
                 offset_time = c(2.9, 4.1, 6.8))

convertToFrames(df$onset_time, df$offset_time)

# ----

sourceCpp(here("Study 1", "data_proc", "temporal_profile.cpp"))

df <- data.frame(
  onset = c(1.5, 3.2, 5.7, 7.8),
  offset = c(2.9, 4.1, 6.2, 9.3),
  category = c("A", "B", "A", "B")
)

generateTemporalProfile(df, "offset", "A", -120, 50, "B")
