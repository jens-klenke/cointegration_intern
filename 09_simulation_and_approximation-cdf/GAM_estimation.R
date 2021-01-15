#---- GAM ----

# Generalized Additive Models
models_gam <- list(
    
    # functional form: poly(t, p) + (1/k)
    gam_case.1_ns_3 <-  mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 3) + ns(I(1/k), 3),
                                  data = data_case_1),
    
    gam_case.1_ns_5 <-  mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 5) + ns(I(1/k), 5),
                                  data = data_case_1),
    
    gam_case.1_ns_7 <-  mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 7) + ns(I(1/k), 7),
                                  data = data_case_1),
    
    gam_case.1_ns_9 <-  mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 9) + ns(I(1/k), 9),
                                  data = data_case_1),
    
    # functional form: poly(t, p) + (1/k) + poly(t, p)*(1/k)
    gam_case.1_ns_3_3 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 3) + I(1/k) + ns(k, 3)*(1/k),
                                   data = data_case_1),
    
    gam_case.1_ns_5_5 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 5) + I(1/k) + ns(k, 5)*(1/k),
                                   data = data_case_1),
    
    gam_case.1_ns_7_7 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 7) + I(1/k) + ns(k, 7)*(1/k),
                                   data = data_case_1),
    
    gam_case.1_ns_9_9 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 9) + I(1/k) + ns(k, 9)*(1/k),
                                   data = data_case_1),
    
    # functional form: poly(t, p) + log(k) + poly(t, p)*log(k)
    gam_case.1_ns_3_log <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 3) + log(k) + ns(k, 3)*log(k),
                                     data = data_case_1),
    
    gam_case.1_ns_5_log <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 5) + log(k) + ns(k, 5)*log(k),
                                     data = data_case_1),
    
    gam_case.1_ns_7_log <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 7) + log(k) + ns(k, 7)*log(k),
                                     data = data_case_1),
    
    gam_case.1_ns_9_log <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 9) + log(k) + ns(k, 9)*log(k),
                                     data = data_case_1),
    
    # functional form: poly(t, p) + k + (1/k)
    gam_case.1_p_3_k_1 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 3) + k + I(1/k),
                                    data = data_case_1),
    
    gam_case.1_p_5_k_1 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 5) + k + I(1/k),
                                    data = data_case_1),
    
    gam_case.1_p_7_k_1 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 7) + k + I(1/k),
                                    data = data_case_1),
    
    gam_case.1_p_9_k_1 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 9) + k + I(1/k),
                                    data = data_case_1),
    
    # functional form: ns(log(t, p)) + log(k) 
    gam_case.1_ns_3_log_ad <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 3) + log(k),
                                        data = data_case_1),
    
    gam_case.1_ns_5_log_ad <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 5) + log(k),
                                        data = data_case_1),
    
    gam_case.1_ns_7_log_ad <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 7) + log(k),
                                        data = data_case_1),
    
    gam_case.1_ns_9_log_ad <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 9) + log(k),
                                        data = data_case_1),
    
    # functional form: ns(log(t, p)) * log(k) + (1/k)
    gam_case.1_ns_3_log_mu <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 3) * log(k),
                                        data = data_case_1),
    
    gam_case.1_ns_5_log_mu <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 5) * log(k),
                                        data = data_case_1),
    
    gam_case.1_ns_7_log_mu <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 7) * log(k),
                                        data = data_case_1),
    
    gam_case.1_ns_9_log_mu <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 9) * log(k),
                                        data = data_case_1),
    
    # bc 
    gam_case.1_bc_poly_log_m_10_poly_m_sqrt <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J_bc, 9) * log(k) + 
                                                             ns(stat_Fisher_E_J_bc, 9)*sqrt(k), data = data_case_1),
    
    # by ~ bc 
    gam_case.1_bc_poly_log_m_10_poly_m_sqrt <- mgcv::gam(p_value_Fisher_E_J_bc ~ ns(stat_Fisher_E_J_bc, 3) * log(k) + 
                                                             ns(stat_Fisher_E_J_bc, 3)*sqrt(k), data = data_case_1)
    
    
)


# save metrics
model_metrics_E_J <- bind_model_metrics(models_gam %>% purrr::map_dfr(metric_fun))
