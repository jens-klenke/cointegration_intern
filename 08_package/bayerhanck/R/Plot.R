plot.bh.test <- function(object, theme = "light") {

  load("null_dist.rda")
  i <- object$basecase

  df_gg <- null_dist %>%
    dplyr::select(paste0('var', i))%>%
    dplyr::mutate(y = rep(1/nrow(null_dist), nrow(null_dist)))
  colnames(df_gg) <- c('x', 'y')

  gg.bh <- ggplot2::ggplot(data = df_gg, ggplot2::aes(x = x)) +
    ggplot2::geom_vline(xintercept = object$bh.test, linetype = "dotted",
               color = "red", size = 1) +
    ggplot2::annotate("text", x = (object$bh.test*1.05), y = 0.5,
                      label = paste('B-H-S \n', round(object$bh.test, 2)),
                      colour = 'red') +
    ggplot2::labs(x = "\n \n Bayer-Hanck-Statistic", y = "F(Bayer-Hanck-Statistic) \n") +
    ggplot2::theme(plot.margin = ggplot2::unit(c(1, 1, 1, 1),"cm"),
          plot.title = ggplot2::element_text(hjust = 0.5))

  if (identical(theme, "light")) {
    gg.bh +
      ggplot2::stat_ecdf(geom = "step") +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.major.y = element_line(size = 0.3),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.x = element_line(size = 0.3),
                     panel.grid.minor.x = element_blank(),
                     axis.line = element_line(size = 0.3))
  } else if (identical(theme, "dark")) {
    gg.bh +
      ggplot2::stat_ecdf(geom = "step", color = "white") +
      ggplot2::theme(plot.background = element_rect(fill = "#1B2B37", colour = "#1B2B37"),
                     panel.background = element_rect(fill = "#1B2B37"),
                     panel.grid.major.y = element_line(size = 0.3, colour = "#546069"),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.x = element_line(size = 0.3, colour = "#546069"),
                     panel.grid.minor.x = element_blank(),
                     axis.text = element_text(colour = "#FFFFFF"),
                     axis.line = element_line(size = 0.3, colour = "#546069"),
                     axis.title = element_text(colour = "#FFFFFF"))
  }

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

