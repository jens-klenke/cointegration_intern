source(here::here('01_code/packages/packages.R'))

parallel::makeCluster(5)

data_case_1 <- Data %>%
  dplyr::filter(case == 1)

test_all_data_case_1_0.2 <- Data_100k %>%
  dplyr::filter(case == 1, 
                p_value_Fisher_all <= 0.2)

test_all_data_case_1 <- Data_100k %>%
  dplyr::filter(case == 1)

rm(Data, Data_100k)


#------------------
# Function
#------------------
metric_fun <- function(object, new_data, new_data_0.2, dep_VAR){
  RMSE <- object$results$RMSE
  
  PRED <- predict(object, new_data)
  RMSE_out <- sqrt((sum((new_data[, dep_VAR] - PRED)^2)/length(PRED)))
  PRED_0.2 <- predict(object, new_data_0.2)
  RMSE_0.2 <- sqrt((sum((new_data_0.2[, dep_VAR] - PRED_0.2)^2)/length(PRED_0.2)))

  call <- object$call$form[3]
  
  metric <- c(as.character(call), as.numeric(RMSE), as.numeric(RMSE_out), as.numeric(RMSE_0.2))
  
  return(metric)
}


#------------------
# Models
#------------------
CV_control <- caret::trainControl(method = "cv", number = 10) 

case_1_all_lm <- caret::train(p_value_Fisher_all ~ stat_Fisher_all + k, 
                              data = data_case_1,
                              method = "lm", 
                              trControl = CV_control)

case_1_all_lm_log <- caret::train(p_value_Fisher_all ~ log(stat_Fisher_all) + k, 
                              data = data_case_1,
                              method = "lm", 
                              trControl = CV_control)

case_1_all_poly_3 <- caret::train(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + k, 
                                  data = data_case_1,
                                  method = "lm", 
                                  trControl = CV_control)

case_1_all_poly_6 <- caret::train(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + k, 
                                  data = data_case_1,
                                  method = "lm", 
                                  trControl = CV_control)

case_1_all_lm_I <- caret::train(p_value_Fisher_all ~ stat_Fisher_all + I(1/k), 
                              data = data_case_1,
                              method = "lm", 
                              trControl = CV_control)

case_1_all_lm_log_I <- caret::train(p_value_Fisher_all ~ log(stat_Fisher_all) + I(1/k), 
                                  data = data_case_1,
                                  method = "lm", 
                                  trControl = CV_control)

case_1_all_lm_log_no_k <- caret::train(p_value_Fisher_all ~ log(stat_Fisher_all), 
                                  data = data_case_1,
                                  method = "lm", 
                                  trControl = CV_control)

case_1_all_lm <- caret::train(p_value_Fisher_all ~ stat_Fisher_all + k, 
                              data = data_case_1,
                              method = "gam", 
                              trControl = CV_control)

# model_metrics <- NULL
model_metrics <- rbind(model_metrics,
                      # metric_fun(case_1_all_lm, test_all_data_case_1, test_all_data_case_1_0.2, "p_value_Fisher_all"), 
                      # metric_fun(case_1_all_lm_log, test_all_data_case_1, test_all_data_case_1_0.2, "p_value_Fisher_all"),
                      # metric_fun(case_1_all_poly_3, test_all_data_case_1, test_all_data_case_1_0.2, "p_value_Fisher_all"),
                      # metric_fun(case_1_all_poly_6, test_all_data_case_1, test_all_data_case_1_0.2, "p_value_Fisher_all"),
                      # metric_fun(case_1_all_lm_I, test_all_data_case_1, test_all_data_case_1_0.2, "p_value_Fisher_all"), 
                      # metric_fun(case_1_all_lm_log_I, test_all_data_case_1, test_all_data_case_1_0.2, "p_value_Fisher_all"), 
                      # metric_fun(case_1_all_lm_log_no_k, test_all_data_case_1, test_all_data_case_1_0.2, "p_value_Fisher_all")
                      )

# colnames(model_metrics) <- c("Formula", "CV RMSE", "Out of sample RMSE", "Out of sample RMSE p <= 0.2")
save(model_metrics, file = "model_metrics.RData")

