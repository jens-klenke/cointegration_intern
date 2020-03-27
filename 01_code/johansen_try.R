load(here::here('00_data/null_dist.rda'))

trendtype <- 2
nvar <- 3

basecase <- 44 * (trendtype - 1) + 4 * (nvar -2) 


case <- basecase + 2


border <- 10000*0.975

null_dist[border, case]

#### tsDyn try

install.packages('tsDyn')

df <- read_csv(here::here('00_data/lutkepohl.csv'))
df <- df%>%
    dplyr::select(linvestment,
                  lincome, 
                  lconsumption)

df_vec <- tsDyn::VECM(df, lag = 1, r = 2, 
                      include = "both",  estim = 'ML')  # lag -1 zu vecrank 


summary(tsDyn::rank.test(df_vec, type = 'eigen'))


