# 1
# Goal: Analyze variance of bootstrap methods (empirical, residual, and wild) for finding 
# the slope of the linear model of response variable Petal.Width and various covariates
# (Sepal.Length, Sepal.Width, Petal.Length). Using R's built in iris dataset.

# a Linear fits from given data 
fit_iris = lm(Petal.Width ~ Sepal.Length + Sepal.Width + Petal.Length, iris)

# b Bootstrapping linear fits
B = 10000
n = nrow(iris)

# EMPIRICAL
coeff_BT_emp = matrix(NA, nrow=B, ncol = 4)
for (i in 1:B) {
  indices_BT = sample(n, n, replace = T)
  iris_BT = iris[indices_BT, 1:4]
  # Get coefficients from linear model fit
  fit_iris_BT = lm(Petal.Width ~ Sepal.Length + Sepal.Width + Petal.Length, iris_BT)$coefficients
  coeff_BT_emp[i, 1] = fit_iris_BT[1]
  coeff_BT_emp[i, 2] = fit_iris_BT[2]
  coeff_BT_emp[i, 3] = fit_iris_BT[3]
  coeff_BT_emp[i, 4] = fit_iris_BT[4]
}
colnames(coeff_BT_emp) = c("Intercept", "Sepal.Length Slope", 
                       "Sepal.Width Slope", "Petal.Length Slope")
var_coeff_BT_emp = var(coeff_BT_emp)

# RESIDUAL
coeff_BT_res = matrix(NA, nrow=B, ncol = 4)
p0 = predict(fit_iris)
for (i in 1:B){
  # BT indices for samples
  indices_BT = sample(n, n, replace = T)
  # Get residuals 
  res = fit_iris$residuals[indices_BT]
  # BT Y*=predicted Petal.Widt
  y_BT = p0 + res
  # Get coefficients from linear model fit
  fit_BT = lm(y_BT ~ iris$Sepal.Length + iris$Sepal.Width + iris$Petal.Length)$coefficients
  coeff_BT_res[i, 1] = fit_BT[1]
  coeff_BT_res[i, 2] = fit_BT[2]
  coeff_BT_res[i, 3] = fit_BT[3]
  coeff_BT_res[i, 4] = fit_BT[4]
}
colnames(coeff_BT_res) = c("Intercept", "Sepal.Length Slope", 
                           "Sepal.Width Slope", "Petal.Length Slope")
var_coeff_BT_res = var(coeff_BT_res)

# WILD
coeff_BT_wild = matrix(NA, nrow=B, ncol = 4)
for (i in 1:B){
  # BT Y*=predicted Petal.Width 
  y_BT = p0 + fit_iris$residuals*rnorm(n)
  # Get coefficients from linear model fit
  fit_BT = lm(y_BT ~ iris$Sepal.Length + iris$Sepal.Width + iris$Petal.Length)$coefficients
  coeff_BT_wild[i, 1] = fit_BT[1]
  coeff_BT_wild[i, 2] = fit_BT[2]
  coeff_BT_wild[i, 3] = fit_BT[3]
  coeff_BT_wild[i, 4] = fit_BT[4]
}
colnames(coeff_BT_wild) = c("Intercept", "Sepal.Length Slope", 
                           "Sepal.Width Slope", "Petal.Length Slope")
var_coeff_BT_wild = var(coeff_BT_wild)

# Final answer: compare variances of coefficients using the three methods
var_comp = matrix(NA, nrow = 3, ncol = 4)
var_comp[1,] = diag(var_coeff_BT_emp)
var_comp[2,] = diag(var_coeff_BT_res)
var_comp[3,] = diag(var_coeff_BT_wild)
colnames(var_comp) = c("Intercept", "Sepal.Length Slope", 
                       "Sepal.Width Slope", "Petal.Length Slope")
rownames(var_comp) = c("Empirical", "Residual", "Wild")


# c Comparison of variance of bootstrap methods
# Boxplots to compare the intercept of Petal.Width ~ Sepal.Length variance across 3 BT methods
boxplot(coeff_BT_emp[, 1], coeff_BT_res[, 1], coeff_BT_wild[, 1], 
        main = "Variance of LM Intercept Across 3 BT Models",
        names = c("Empirical", "Residual", "Wild"),
        col = c("pink", "red", "orange"))


# -----------------------------------------------------------------------------------------

# 2
# Logistic regression of response variable admission based on covariates gre and gpa

setwd("/Users/nicolesmith/Desktop/School/stat-403-projects")
admissions = read.csv("binary.csv") 
B = 10000
n = nrow(admissions)


# a Parametric bootstrap for 90% CI for slope of GPA
# get estimated model parameter from linear fit
# 1 means admit, 0 means reject (binomial)
lm_GPA = glm(admit ~ gpa, admissions, family="binomial")
slope_sample = lm_GPA$coefficients[2]
p0 = predict(lm_GPA, type = 'response')
slope_GPA_BT = rep(NA, B)
for (i in 1:B){
  # bootstrap sample for admit or not
  Y_BT = rbinom(n, 1, p0)
  # model from bootstrap sample, obtain slope
  lm_BT = glm(Y_BT ~ admissions$gpa, family="binomial")
  slope_GPA_BT[i] = lm_BT$coefficients[2]
}
# histogram of BT sample of model slope
hist(slope_GPA_BT)
# construct 90% CI
quantile(slope_GPA_BT, c(0.05, 0.95))


# b Use both param. and empirical BT to estimate SE of the intercept, slope of admit~GRE, 
# and slope of admit~GPA
# empirical bootstrap
coeff_emp_BT = matrix(NA, B, 3)
for (i in 1:B){
  # get BT samples
  BT_i = sample(n, n, replace = T)
  sample_BT = admissions[BT_i,]
  # fit models from BT samples, obtain intecept, slope of GRE model, slope of GPA model
  lm_GRE_BT = lm(admit ~ gre, sample_BT)
  lm_GPA_BT = lm(admit ~ gpa, sample_BT)
  coeff_emp_BT[i, 1] = lm_GRE_BT$coefficients[1]
  coeff_emp_BT[i, 2] = lm_GRE_BT$coefficients[2]
  coeff_emp_BT[i, 3] = lm_GPA_BT$coefficients[2]
}


# parametric bootstrap
# fit linear models from data 
lm_GRE = glm(admit ~ gre, admissions, family="binomial")
lm_GPA = glm(admit ~ gpa, admissions, family="binomial")
coeff_param_BT = matrix(NA, B, 3)
# estimate probability of admission based on GRE, GPA
p0_GRE = predict(lm_GRE, type = 'response')
p0_GPA = predict(lm_GPA, type = 'response')
for (i in 1:B){
  # bootstrap sample for admit or not
  Y_BT_GRE = rbinom(n, 1, p0_GRE)
  Y_BT_GPA = rbinom(n, 1, p0_GPA)
  # Fit model from BT samples, get coefficients
  lm_GRE_BT = glm(Y_BT_GRE ~ admissions$gre, family="binomial")
  lm_GPA_BT = glm(Y_BT_GPA ~ admissions$gpa, family="binomial")
  coeff_param_BT[i, 1] = lm_GRE_BT$coefficients[1]
  coeff_param_BT[i, 2] = lm_GRE_BT$coefficients[2]
  coeff_param_BT[i, 3] = lm_GPA_BT$coefficients[2]
}

# compare the standard error of the two BT methods for each coefficient estimate
se_emp = c(sd(coeff_emp_BT[,1]), sd(coeff_emp_BT[,2]), sd(coeff_emp_BT[,3]))
se_param = c(sd(coeff_param_BT[,1]), sd(coeff_param_BT[,2]), sd(coeff_param_BT[,3]))
compare_se = matrix(data= c(se_emp, se_param), byrow = T, nrow = 2, ncol = 3)
colnames(compare_se) = c("SE(Intercept)", "SE(GRE Slope)", "SE(GPA Slope)")
rownames(compare_se) = c("Empirical", "Parametric")



# c 90% CI for probability that a student is admitted | gre = 500 and gpa = 3.7
# Empirical Bootstrap for estimates of probabilty p
B = 100000
n = nrow(admissions)
p_BT = rep(NA, B)
for (i in 1:B) {
  BT_i = sample(n, n, replace = T)  
  sample_BT = admissions[BT_i,]
  # fit models from BT samples 
  lm_GRE_BT = lm(admit ~ gre, sample_BT)
  lm_GPA_BT = lm(admit ~ gpa, sample_BT)
  # predict probability of admission for GRE = 500, GPA = 3.7
  p0_GRE = predict(lm_GRE_BT, newdata = 500, type = 'response')
  p0_GPA = predict(lm_GPA_BT, newdata = 3.7, type = 'response')
  
  p_BT[i] = 0
}
# Construct 90% CI using quantile method
CI_lambda = quantile(p_BT, c(0.05, 0.95))


