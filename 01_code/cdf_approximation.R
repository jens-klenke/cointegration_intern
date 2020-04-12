load('C:/Users/HP/Dropbox/stata/Simulation_04.04.RData')


k_1_t_1 <- Null_Distr_B_ECR_J_E[1,1,]

k_1_t_1[c(length(k_1_t_1) - 1, length(k_1_t_1))] <- 73.27065

p_value <- (seq(1:length(k_1_t_1))/length(k_1_t_1))

dist <- data.frame(chi_sq = k_1_t_1,
                   p_value = p_value)

mod <- mgcv::gam(p_value ~ chi_sq ,data = dist)

summary(mod)

