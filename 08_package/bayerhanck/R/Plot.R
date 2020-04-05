plot.bh.test <- function(object, theme = "dark") {

  load("bayerhanck_cv.rda")

  #getting the distribution
  if (identical(object$trend, "const"))
    cv <- Null_Distr_E_J
  if (identical(object$test.type, "all"))
    cv <- Null_Distr_B_ECR_J_E

  if (identical(object$trend, "none")) {
    trendtype = 1
  } else if (identical(object$trend, "const")) {
    trendtype = 2
  } else if (identical(object$trend, "trend")) {
    trendtype = 3
  }

  x <- cv[object$K, trendtype, ]
  y <- rep((1/length(x)), length(x))

  df_gg <- data.frame(x, y)

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
      ggthemes::theme_economist() +
      ggplot2::stat_ecdf(geom = "step") +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "#FFFFFF", colour = "#FFFFFF"),
                     panel.background = ggplot2::element_rect(fill = "#FFFFFF"),
                     panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "#546069"),
                     axis.line = ggplot2::element_line(size = 0.3, colour = "#546069"))
  } else if (identical(theme, "dark")) {
    gg.bh +
      ggthemes::theme_economist() +
      ggplot2::stat_ecdf(geom = "step", color = "white") +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "#1B2B37", colour = "#1B2B37"),
                     panel.background = ggplot2::element_rect(fill = "#1B2B37"),
                     panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "#546069"),
                     axis.text = ggplot2::element_text(colour = "#FFFFFF"),
                     axis.line = ggplot2::element_line(size = 0.3, colour = "#546069"),
                     axis.title = ggplot2::element_text(colour = "#FFFFFF"),
                     axis.ticks.x = ggplot2::element_line(colour = "#FFFFFF"))
  }
}

