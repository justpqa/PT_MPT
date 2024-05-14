library(Rcpp)
library(parallel)
library(dplyr)
library(tidyr)
sourceCpp('~/Stats 322 S24/Submit/PT_MPT.cpp')
# THIS CODE IS THE VERSION OF THE RMD FILE BUT USING THE IMPLEMENTATION IN RCPP INSTEAD
mu_test <- 0
sigma_test <- 1
theta_list_test1 <- matrix(c(0.4, 0.6, 0, 0), nrow = 1)
theta_list_test2 <- matrix(c(0.4, 0.3, 0.6, 0.7, 0, 0.6, 0, 0.4), nrow = 2)
pdf_pt_prior_cpp(0, mu_test, sigma_test, theta_list_test1) # should be close to 0.4
cdf_pt_prior_cpp(0, mu_test, sigma_test, theta_list_test1) # should be close to 0.5
cdf_pt_prior_cpp(qnorm(1/4), mu_test, sigma_test, theta_list_test2) # should be close to 0.12
cdf_pt_prior_cpp(qnorm(7/8), mu_test, sigma_test, theta_list_test2)
pt_update_cpp(c(1/4, 1/2, -1/4, -1/4), 0, 1, 2)
data <- (mtcars$mpg)
# randomly chosen
# mu_0 <- 15
# sigma_0 <- 5
# how can it improved a informed parans
mu_0 <- mean(data)
sigma_0 <- sd(data)
normal <- data.frame(x = rnorm(1000, mean = mu_0, sd = sigma_0)) %>%
  mutate(y = dnorm(x, mean = mu_0, sd = sigma_0))
p <- ggplot() +
  geom_histogram(data = data.frame(data), aes(x = data), bins = 30) +
  ylab("Count")
p + geom_line(data = normal, aes(x = x, y = y *75, color = "Normal")) +
  scale_y_continuous(sec.axis = sec_axis(~ ./75, name = "Density")) +
  scale_color_manual(name = "Legend", 
                     values = c("Normal" = "red"))
pt_1_post <- pt_update_cpp(data, mu_0, sigma_0, J = 1)
pt_1 <- data.frame(x = rnorm(1000, mean = mu_0, sd = sigma_0))
# Now we draw 10 different distribution from PT prior
for (i in 1:10) {
  # get the column name to add into the data
  curr_col_name <- paste("y_", as.character(i), sep = "") 
  # now draw the random theta
  curr_theta_list <- draw_theta_cpp(TRUE, pt_1_post, 1, 1)
  # Now we fill the column
  curr_col <- unlist(lapply(pt_1$x, function(x) {return(pdf_pt_prior_cpp(x, mu_0, sigma_0, curr_theta_list))}))
  pt_1[,curr_col_name] <- curr_col
}
pt_1 <- pt_1 %>%
  pivot_longer(cols = starts_with("y"), names_to = "y", values_to = "Density")
p <- ggplot() +
  geom_histogram(data = data.frame(data), aes(x = data), bins = 30) +
  ylab("Count")
p + geom_line(data = pt_1, aes_string(x = "x", y = paste("Density", "*75",sep = ""), color = "y"), linetype = "dashed") + 
  scale_y_continuous(sec.axis = sec_axis(~ ./75, name = "Density"))
pt_1 <- data.frame(x = rnorm(1000, mean = mu_0, sd = sigma_0), 
                   y = numeric(1000))
for (i in 1:1000) {
  if (i %% 10 == 0) {
    # print(i)
  }
  pt_1[i, 2] <- pdf_pt_prior_mixture_c_cpp(pt_1[i, 1], mu_0, sigma_0, 1, 1, 10000, TRUE, pt_1_post)
}
p <- ggplot() +
  geom_histogram(data = data.frame(data), aes(x = data), bins = 30) +
  ylab("Count")
p + geom_line(data = normal, aes(x = x, y = y *75, color = "Normal")) +
  geom_line(data = pt_1, aes(x = x, y = y *75, color = "PT 1 level"))  + 
  scale_y_continuous(sec.axis = sec_axis(~ ./75, name = "Density")) +
  scale_color_manual(name = "Legend", 
                     values = c("Normal" = "red", 
                                "PT 1 level" = "blue")) + 
  geom_vline(xintercept = mosaic::fav_stats(data)$Q1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = mosaic::fav_stats(data)$median, linetype = "dashed", color = "black") +
  geom_vline(xintercept = mosaic::fav_stats(data)$Q3, linetype = "dashed", color = "black")
pt_2_post <- pt_update_cpp(data, mu_0, sigma_0, J = 2)
pt_2 <- data.frame(x = rnorm(1000, mean = mu_0, sd = sigma_0))
# Now we draw 10 different distribution from PT prior
for (i in 1:10) {
  # get the column name to add into the data
  curr_col_name <- paste("y_", as.character(i), sep = "") 
  # now draw the random theta
  curr_theta_list <- draw_theta_cpp(TRUE, pt_2_post, 1, 2)
  # Now we fill the column
  curr_col <- unlist(lapply(pt_2$x, function(x) {return(pdf_pt_prior_cpp(x, mu_0, sigma_0, curr_theta_list))}))
  pt_2[,curr_col_name] <- curr_col
}
pt_2 <- pt_2 %>%
  pivot_longer(cols = starts_with("y"), names_to = "y", values_to = "Density")
p <- ggplot() +
  geom_histogram(data = data.frame(data), aes(x = data), bins = 30) +
  ylab("Count")
p + geom_line(data = pt_2, aes_string(x = "x", y = paste("Density", "*75",sep = ""), color = "y"), linetype = "dashed") + 
  scale_y_continuous(sec.axis = sec_axis(~ ./75, name = "Density"))
pt_2 <- data.frame(x = rnorm(1000, mean = mu_0, sd = sigma_0), 
                   y = numeric(1000))
for (i in 1:1000) {
  if (i %% 10 == 0) {
    # print(i)
  }
  pt_2[i, 2] <- pdf_pt_prior_mixture_c_cpp(pt_2[i, 1], mu_0, sigma_0, 1, 2, 10000, TRUE, pt_2_post)
}
p <- ggplot() +
  geom_histogram(data = data.frame(data), aes(x = data), bins = 30) +
  ylab("Count")
p + geom_line(data = normal, aes(x = x, y = y *75, color = "Normal")) +
  geom_line(data = pt_2, aes(x = x, y = y *75, color = "PT 2 level"))  + 
  scale_y_continuous(sec.axis = sec_axis(~ ./75, name = "Density")) +
  scale_color_manual(name = "Legend", 
                     values = c("Normal" = "red", 
                                "PT 2 level" = "blue")) + 
  geom_vline(xintercept = mosaic::fav_stats(data)$Q1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = mosaic::fav_stats(data)$median, linetype = "dashed", color = "black") +
  geom_vline(xintercept = mosaic::fav_stats(data)$Q3, linetype = "dashed", color = "black")
# We started with 1 level first
mpt_1 <- data.frame(x = rnorm(1000, mean = mu_0, sd = sigma_0))
# Our prior for mu and sigma
d_mu_test <- function(x) {
  return(dnorm(x, mu_0, sigma_0 / 10))
}
r_mu_test <- function() {
  return(rnorm(1, mu_0, sigma_0 / 10))
}
d_sigma_test <- function(x) {
  return(dnorm(x, sigma_0, sigma_0 / 10))
}
r_sigma_test <- function() {
  return(rnorm(1, sigma_0, sigma_0 / 10))
}
# Draw 10 random theta
for (i in 1:10) {
  print(i)
  # get the column name to add into the data
  curr_col_name <- paste("y_", as.character(i), sep = "") 
  # now draw the random theta
  curr_theta_list <- draw_theta_cpp(TRUE, pt_1_post, 1, 1)
  # Now we fill the column
  curr_col <- unlist(mclapply(mpt_1$x, 
                            function(x) {return(pdf_mpt_prior_cpp(x, 
                                                                  mu_0, d_mu_test, r_mu_test, 
                                                                  sigma_0, d_sigma_test, r_sigma_test,
                                                                  curr_theta_list,
                                                                  10000,
                                                                  TRUE, data,
                                                                  1, 1))},
                            mc.cores = 4))
  mpt_1[,curr_col_name] <- curr_col
}
mpt_1 <- mpt_1 %>%
  pivot_longer(cols = starts_with("y"), names_to = "y", values_to = "Density")
p <- ggplot() +
  geom_histogram(data = data.frame(data), aes(x = data), bins = 30) +
  ylab("Count")
p + geom_line(data = mpt_1, aes_string(x = "x", y = paste("Density", "*75",sep = ""), color = "y"), linetype = "dashed") + 
  scale_y_continuous(sec.axis = sec_axis(~ ./75, name = "Density"))
mpt_2 <- data.frame(x = rnorm(1000, mean = mu_0, sd = sigma_0))
# Draw 10 random theta
for (i in 1:10) {
  print(i)
  # get the column name to add into the data
  curr_col_name <- paste("y_", as.character(i), sep = "") 
  # now draw the random theta
  curr_theta_list <- draw_theta_cpp(TRUE, pt_2_post, 1, 2)
  # Now we fill the column
  curr_col <- unlist(mclapply(mpt_2$x, 
                              function(x) {return(pdf_mpt_prior_cpp(x, 
                                                                    mu_0, d_mu_test, r_mu_test, 
                                                                    sigma_0, d_sigma_test, r_sigma_test,
                                                                    curr_theta_list,
                                                                    10000,
                                                                    TRUE, data,
                                                                    1, 1))},
                              mc.cores = 4))
  mpt_2[,curr_col_name] <- curr_col
}
mpt_2 <- mpt_2 %>%
  pivot_longer(cols = starts_with("y"), names_to = "y", values_to = "Density")
p <- ggplot() +
  geom_histogram(data = data.frame(data), aes(x = data), bins = 30) +
  ylab("Count")
p + geom_line(data = mpt_2, aes_string(x = "x", y = paste("Density", "*75",sep = ""), color = "y"), linetype = "dashed") + 
  scale_y_continuous(sec.axis = sec_axis(~ ./75, name = "Density"))
