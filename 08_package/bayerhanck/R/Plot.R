#' @export
plot.bh.test <- function(x, theme = "dark", ...) {
x
  #getting the distribution
  if (identical(x$test.type, "eg-j"))
    cv <- Null_Distr_E_J
  if (identical(x$test.type, "all"))
    cv <- Null_Distr_B_ECR_J_E

  if (identical(x$trend, "none")) {
    trendtype = 1
  } else if (identical(x$trend, "const")) {
    trendtype = 2
  } else if (identical(x$trend, "trend")) {
    trendtype = 3
  }

  z <- cv[x$K, trendtype, ]
  y <- rep((1/length(z)), length(z))

  df_gg <- data.frame(z, y)

  gg.bh <- ggplot2::ggplot(data = df_gg, ggplot2::aes(x = z)) +
    ggplot2::geom_vline(xintercept = x$bh.test, linetype = "dotted",
               color = "red", size = 1) +
    ggplot2::labs(x = "\n \n Bayer-Hanck-Statistic", y = "F(Bayer-Hanck-Statistic) \n") +
    ggplot2::theme(plot.margin = ggplot2::unit(c(1, 1, 1, 1),"cm"),
          plot.title = ggplot2::element_text(hjust = 0.5))

  if (identical(theme, "light")) {
    suppressWarnings(   gg.bh +
      ggthemes::theme_economist() +
      ggplot2::stat_ecdf(geom = "step") +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "#FFFFFF", colour = "#FFFFFF"),
                     panel.background = ggplot2::element_rect(fill = "#FFFFFF"),
                     panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "#546069"),
                     axis.line = ggplot2::element_line(size = 0.3, colour = "#546069")))
  } else if (identical(theme, "dark")) {
    suppressWarnings(   gg.bh +
      ggthemes::theme_economist() +
      ggplot2::stat_ecdf(geom = "step", color = "white") +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "#1B2B37", colour = "#1B2B37"),
                     panel.background = ggplot2::element_rect(fill = "#1B2B37"),
                     panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "#546069"),
                     axis.text = ggplot2::element_text(colour = "#FFFFFF"),
                     axis.line = ggplot2::element_line(size = 0.3, colour = "#546069"),
                     axis.title = ggplot2::element_text(colour = "#FFFFFF"),
                     axis.ticks.x = ggplot2::element_line(colour = "#FFFFFF")))
  }
}

