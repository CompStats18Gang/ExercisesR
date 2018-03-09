# Computational Statistics, SS18
# Exercise 1: Linear Regression
# Daniel Regenass, IACETH


# parameter
intercept_true <- 1.0
slope_true <- 2.0
sigma_true <- 5.0
number_iterations <- 100

# x axis
x <- seq(1,40,1)

# fix random seed for reproducability
set.seed(42)

# empty vector (not list) for slope coefficients
beta_2 <- vector(length=number_iterations)



####################################################################################
# Exercises a) b) and d) are packed into functions such
# that they may be recycled for exercise 4.

# a) iteration to create 100 different realizations for y
# error distribution may be changed to solve exercises
random_iterator <- function(x, error_distribution) {
  for (i in 1:number_iterations) {
    if (error_distribution == "normal") {
      stochastic_errors <- sigma_true * rnorm(length(x))
    } else if (error_distribution == "chi_squared") {
      stochastic_errors <- sigma_true * (1 - rchisq(length(x), df = 1)) / sqrt(2)
    } else if (error_distribution == "quadratic") {
      stochastic_errors <- rnorm(length(x), mean = x^2 / 40 - 1, sd = 1)
    } else if (error_distribution == "toeplitz") {
      Sigma <- toeplitz(c(seq(from = 1, to = 0, by = -0.1), rep(0, 29)))
      stochastic_errors <- MASS::mvrnorm(n=1, mu=rep(0, length(x)), Sigma=Sigma)
    } else if (error_distribution == "reduced_sd") {
      stochastic_errors <- rnorm(length(x), mean = 0, sd = x / 20)
    }
    y <- intercept_true + slope_true * x + stochastic_errors
    regression <- lm(y~x)
    beta_2[[i]] <- regression$coefficients[2]
  }
  
  #a) plot one (i.e. the 100th) realization
  plot(x,y)
  abline(regression)
  
  # tukey-anscombe plot
  plot(regression, which=1)
  
  return(beta_2)
}


# Exercise b)
print_mean_sd <- function(beta_2) {
 #b) mean of estimated slopes:
 slopes_mean <- mean(beta_2)
 print(slopes_mean)

 #b) sd of estimated slopes:
 slopes_sd <- sd(beta_2)
 print(slopes_sd)
}

# Exercise d) Histogram plot:
plot_histogram <- function(beta_2) {
 hist(beta_2, freq=FALSE)
 lines(seq(1.8, 2.3, by=0.01), dnorm(seq(1.8, 2.3, by=0.01), mean=2, sd=0.05))
}

####################################################################################


# 3a) use function with normally distributed residuals
beta_2_normal <- random_iterator(x, "normal")

# 3b) and print mean and sd
print_mean_sd(beta_2_normal)


#3c) Theoretical variance:
var_beta_th <- solve(t(x)%*%x) 
var_beta_th <- sigma_true^2.* var_beta_th
sd_beta_th <- sqrt(var_beta_th)
print(t(x)%*%x)
print(sd_beta_th)
print(var_beta_th)
print(var(beta_2))

#3d) Histogram plot
plot_histogram(beta_2_normal)


# Ex 4) Test different error distibrutions

#4a)errors not normally, but chi-squared distributed
beta_2_chi2 <- random_iterator(x, "chi_squared")
# sample mean and sd
print_mean_sd(beta_2_chi2)
# histogram plot
plot_histogram(beta_2_chi2)

#4b) quadratic dependency of residuals on x
beta_2_quadratic <- random_iterator(x, "quadratic")
# sample mean and sd
print_mean_sd(beta_2_quadratic)
# histogram plot
plot_histogram(beta_2_quadratic)

#4c) who the f*** is toeplitz
beta_2_toeplitz <- random_iterator(x, "toeplitz")
# sample mean and sd
print_mean_sd(beta_2_toeplitz)
# histogram plot
plot_histogram(beta_2_toeplitz)

# violates homoscedasticity (variance not constant)
beta_2_reducedsd <- random_iterator(x, "reduced_sd")
# b) sample mean and sd
print_mean_sd(beta_2_reducedsd)
# d) histogram plot
plot_histogram(beta_2_reducedsd)








