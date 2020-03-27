johansen <- function(x, type = "eigen", lags = 1, trend = "const"){
  df_vec <- tsDyn::VECM(x, lag = lags, r = 2,
                        include = trend, # Bezeichnung Trend anpassen
                        estim = 'ML')  # lag -1 zu vecrank
  test.stat <- summary(tsDyn::rank.test(df_vec, type = 'eigen'))$eigen[1]
  print(test.stat)
}

df <- read_csv(here::here('00_data/lutkepohl.csv'))
df <- df %>%
  dplyr::select(linvestment,
                lincome,
                lconsumption)
johansen(df)
?tsDyn::VECM
