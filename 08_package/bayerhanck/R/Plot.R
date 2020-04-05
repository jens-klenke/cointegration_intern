plot.bh.test <- function(object, theme = "dark") {

  load("bayerhanck_cv.rda")

  #getting the distribution
  if (identical(object$trend, "const"))
    cv <- Null_Distr_E_J
  if (identical(object$test.type, "all"))
    cv <- Null_Distr_B_ECR_J_E

  if (identical(object$trend, "none")){
    trendtype = 1
  } else if (identical(object$trend, "const")){
    trendtype = 2
  } else if (identical(object$trend, "trend")){
    trendtype = 3
  }




  x <- cv[object$K, trendtype,]
  y <- rep((1/length(x)), length(x))

  df_gg <- data.frame(x,y)


  gg.bh <- ggplot2::ggplot(data = df_gg, ggplot2::aes(x = x)) +
    ggplot2::geom_vline(xintercept = object$bh.test, linetype = "dotted",
               color = "red", size = 1) +
    ggplot2::labs(x = "\n Bayer-Hanck-Statistic", y = "F(Bayer-Hanck-Statistic) \n") +
    ggplot2::theme(legend.position="right") #+
  #  ggplot2::annotate("text", x = (object$bh.test*(-1.1)), y = 0.5,
  #                    label = paste('B-H-S \n', round(object$bh.test, 2)),
  #                    colour = 'red')

  if (identical(theme, "light")) {
    gg.bh +
      ggthemes::theme_economist() +
      ggplot2::stat_ecdf(geom = "step") +
<<<<<<< HEAD
      ggplot2::theme(panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "#546069"),
                     plot.background = ggplot2::element_rect(fill = "#FFFFFF", colour = "#FFFFFF"),
                     panel.background = ggplot2::element_rect(fill = "#FFFFFF"))
  } else if (identical(theme, "dark")) {
    gg.bh +
      ggplot2::stat_ecdf(geom = "step", color = "white") +
      ggthemes::theme_economist() +
=======
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_line(size = 0.3),
                     panel.grid.minor.y = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_line(size = 0.3),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(size = 0.3))
  } else if (identical(theme, "dark")) {
    gg.bh +
      ggplot2::stat_ecdf(geom = "step", color = "white") +
>>>>>>> 5f0d74def01172f99a67cd856df3eb83eb377670
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "#1B2B37", colour = "#1B2B37"),
                     panel.background = ggplot2::element_rect(fill = "#1B2B37"),
                     panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "#546069"),
                     panel.grid.minor.y = ggplot2::element_blank(),
<<<<<<< HEAD
                     axis.text = ggplot2::element_text(colour = "#FFFFFF"),
                     axis.line = ggplot2::element_line(size = 0.3, colour = "#546069"),
                     axis.title = ggplot2::element_text(colour = "#FFFFFF"),
                     axis.ticks.x = ggplot2::element_line(color = "#FFFFFF"))
=======
                     panel.grid.major.x = ggplot2::element_line(size = 0.3, colour = "#546069"),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     axis.text = ggplot2::element_text(colour = "#FFFFFF"),
                     axis.line = ggplot2::element_line(size = 0.3, colour = "#546069"),
                     axis.title = ggplot2::element_text(colour = "#FFFFFF"))
>>>>>>> 5f0d74def01172f99a67cd856df3eb83eb377670
  }

}

