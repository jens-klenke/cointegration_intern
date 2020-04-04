plot.bh.test <- function(obj) {

  # Obtain Parameters from the Bayerhanck object

  K <- obj$K
  if (identical(obj$trend, "none"))
    trendtype <- 1
  if (identical(obj$trend, "const"))
    trendtype <- 2
  if (identical(obj$trend, "trend"))
    trendtype <- 3

  bh.stat <- obj$bh.test
  stat.type <- obj$test.typ
  crit <- obj$crit



  basecase <- 44 * (trendtype - 1) + 4 * (nvar - 2)

  load("null_dist.rda")
  i <- basecase

  df_gg <- null_dist%>%
    dplyr::select(as.name(paste0('var', i)))%>%
    dplyr::mutate(y = rep(1/nrow(null_dist), nrow(null_dist)))
  colnames(df_gg) <- c('x', 'y')

  ggplot(data = df_gg, aes(x = x))+
    stat_ecdf(geom = "step")+
    geom_vline(xintercept = bh.stat, linetype = "dashed",
               color = "red", size = 1)+
    annotate("text", x = (bh.stat*1.05), y = 0.5, label = paste('B-H-S \n', round( bh.stat,2)) ,
             colour = 'red')+
    labs( x = "\n \n Bayer-Hanck-Statistic", y = "F(Bayer-Hanck-Statistic) \n",
          title = 'Empirical Cumulative Distribution Function')+
    theme_minimal()+
    theme(plot.margin = unit(c(1,1,1,1),"cm"),
          plot.title = element_text(hjust = 0.5))

#  ggplot(data = df_gg) +
#    geom_density(aes(x = x))+
 #   geom_vline(xintercept = b_h_stat_2, linetype = "dashed",
#               color = "red", size = 1)+
#    annotate("text", x = (b_h_stat_2*1.05), y = (max(density(df_gg$x)$y)/2),
#             label = paste('B-H-S \n', b_h_stat_2) , colour = 'red')+
 #   labs( x = "\n \n Bayer-Hanck-Statistic", y = "Density \n",
  #        title = 'Kernal Density')+
#    theme_minimal()+
#    theme(plot.margin = unit(c(1,1,1,1),"cm"),
#          plot.title = element_text(hjust = 0.5))
}

