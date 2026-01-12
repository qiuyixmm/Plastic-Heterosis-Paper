## define a function
calc_t_test_one_sided <- function(P_RC_samples, P_M, P_F) {
  N <- length(P_RC_samples) 
  P_RC_mean <- mean(P_RC_samples)
  H_percent <- ((P_RC_mean - (P_M + P_F) / 2) / ((P_M + P_F) / 2)) * 100 
  
  variance_P_RC <- sum((P_RC_samples - P_RC_mean)^2) / (N - 1) 
  denominator <- sqrt((2 * variance_P_RC / (N - 1)) / ((P_M + P_F) * sqrt(N)))  
  
  t_value <- H_percent / denominator
  
  df <- N - 1 
  p_value <- 1 - pt(t_value, df) 
  
  list(
    H_percent = round(H_percent, 2),
    t_value = round(t_value, 3),
    p_value = format.pval(p_value, digits = 3),
    significant = ifelse(p_value < 0.05, "significant", "insignificant")
  )
}


################# fed ad libitum
Peking_avg <- mean(c(58, 60, 63, 64, 67, 44, 47, 51, 57, 57)) 
Muscovy_avg <- mean(c(64, 72, 72, 76, 70, 93, 78, 91, 76))

#### For mule ducks
mule <- c(58, 56, 61, 58, 56, 65, 80, 65, 63, 52)
mule.result <- calc_t_test_one_sided(mule, Peking_avg, Muscovy_avg)
print(mule.result)  # H_percent: -8.14; t_value: -45.678; p_value: 1

#### For hinny ducks
hinny <- c(56, 66, 64, 66, 56, 70, 61, 54, 72, 66)
hinny.result <- calc_t_test_one_sided(hinny, Peking_avg, Muscovy_avg)
print(hinny.result)  # H_percent: -5.6; t_value: -39.699; p_value: 1


################# overfeeding
Peking_avg <- mean(c(122, 455, 519, 189, 338, 343, 362, 366, 385, 442))
Muscovy_avg <- mean(c(472, 580, 617, 422, 499, 551, 513, 599, 540, 555))

#### For mule ducks
mule <- c(654, 630, 639, 594, 713, 733, 745, 643, 539, 572)
mule.result <- calc_t_test_one_sided(mule, Peking_avg, Muscovy_avg)
print(mule.result)  # H_percent: 45.72; t_value: 75.31; p_value: 3.24e-14

#### For hinny ducks
hinny <- c(485, 448, 605, 524, 512, 661, 674, 632, 452, 685)
hinny.result <- calc_t_test_one_sided(hinny, Peking_avg, Muscovy_avg)
print(hinny.result)  # H_percent: 28.04; t_value: 33.655; p_value: 4.45e-11
