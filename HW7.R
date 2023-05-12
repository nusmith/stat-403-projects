# 1
# Goal: Analyze variance of bootstrap methods (empirical, residual, and wild) for finding 
# the slope of the linear model of response variable Petal.Width and various covariates
# (Sepal.Length, Sepal.Width, Petal.Length). Using R's built in iris dataset.

# a Linear fits from given data 
fit_slen = lm(Petal.Width ~ Sepal.Length, iris)
fit_swid = lm(Petal.Width ~ Sepal.Width, iris)
fit_plen = lm(Petal.Width ~ Petal.Length, iris)
fit_slen$coefficients
fit_swid$coefficients
fit_plen$coefficients
x = seq(0, 10, 0.1)
plot(iris$Sepal.Length, iris$Petal.Width)
lines(x, fit_slen$coefficients[1] + fit_slen$coefficients[2]*x, col="red")
plot(iris$Sepal.Width, iris$Petal.Width)
lines(x, fit_swid$coefficients[1] + fit_swid$coefficients[2]*x, col="red")
plot(iris$Petal.Length, iris$Petal.Width)
lines(x, fit_plen$coefficients[1] + fit_plen$coefficients[2]*x, col="red")


# b Bootstrapping linear fits
B = 10000
n = nrow(iris)

# EMPIRICAL
coeff_BT_emp = matrix(NA, nrow=B, ncol = 4)
for (i in 1:B) {
  indices_BT = sample(n, n, replace = T)
  iris_BT = iris[indices_BT, 1:4]
  # Get coefficients from linear model fit
  fit_slen_BT = lm(Petal.Width ~ Sepal.Length, iris_BT)$coefficients
  fit_swid_BT = lm(Petal.Width ~ Sepal.Width, iris_BT)$coefficients
  fit_plen_BT = lm(Petal.Width ~ Petal.Length, iris_BT)$coefficients
  # currently using a BT sample for the first intercept (sepal length) only
  coeff_BT_emp[i, 1] = fit_slen_BT[1]
  coeff_BT_emp[i, 2] = fit_slen_BT[2]
  coeff_BT_emp[i, 3] = fit_swid_BT[2]
  coeff_BT_emp[i, 4] = fit_plen_BT[2]
}
colnames(coeff_BT_emp) = c("Sepal.Length Int.", "Sepal.Length Slope", 
                       "Sepal.Width Slope", "Petal.Length Slope")
var_coeff_BT_emp = var(coeff_BT)

# RESIDUAL
coeff_BT_res = matrix(NA, nrow=B, ncol = 4)
slen_pred = predict(fit_slen)
swid_pred = predict(fit_swid)
plen_pred = predict(fit_plen)
for (i in 1:B){
  # BT indices for samples
  indices_BT = sample(n, n, replace = T)
  # Get residuals for each of the three models
  res_slen = fit_slen$residuals[indices_BT]
  res_swid = fit_swid$residuals[indices_BT]
  res_plen = fit_plen$residuals[indices_BT]
  # BT Y*=predicted Petal.Width for each of the three models
  y_slen = slen_pred + res_slen
  y_swid = swid_pred + res_swid
  y_plen = plen_pred + res_plen
  # Get coefficients from linear model fit
  fit_slen_BT = lm(y_slen ~ iris$Sepal.Length)$coefficients
  fit_swid_BT = lm(y_swid ~ iris$Sepal.Width)$coefficients
  fit_plen_BT = lm(y_plen ~ iris$Petal.Length)$coefficients
  # currently using a BT sample for the first intercept (sepal length) only
  coeff_BT_res[i, 1] = fit_slen_BT[1]
  coeff_BT_res[i, 2] = fit_slen_BT[2]
  coeff_BT_res[i, 3] = fit_swid_BT[2]
  coeff_BT_res[i, 4] = fit_plen_BT[2]
}
colnames(coeff_BT_res) = c("Sepal.Length Int.", "Sepal.Length Slope", 
                           "Sepal.Width Slope", "Petal.Length Slope")
var_coeff_BT_res = var(coeff_BT_res)

# WILD
coeff_BT_wild = matrix(NA, nrow=B, ncol = 4)
for (i in 1:B){
  # BT Y*=predicted Petal.Width for each of the three models
  y_slen = slen_pred + fit_slen$residuals*rnorm(n)
  y_swid = swid_pred + fit_swid$residuals*rnorm(n)
  y_plen = plen_pred + fit_plen$residuals*rnorm(n)
  # Get coefficients from linear model fit
  fit_slen_BT = lm(y_slen ~ iris$Sepal.Length)$coefficients
  fit_swid_BT = lm(y_swid ~ iris$Sepal.Width)$coefficients
  fit_plen_BT = lm(y_plen ~ iris$Petal.Length)$coefficients
  # currently using a BT sample for the first intercept (sepal length) only
  coeff_BT_wild[i, 1] = fit_slen_BT[1]
  coeff_BT_wild[i, 2] = fit_slen_BT[2]
  coeff_BT_wild[i, 3] = fit_swid_BT[2]
  coeff_BT_wild[i, 4] = fit_plen_BT[2]
}
colnames(coeff_BT_wild) = c("Sepal.Length Int.", "Sepal.Length Slope", 
                           "Sepal.Width Slope", "Petal.Length Slope")
var_coeff_BT_wild = var(coeff_BT_wild)

# Final answer: compare variances of coefficients using the three methods
var_comp = matrix(NA, nrow = 3, ncol = 4)
var_comp[1,] = diag(var_coeff_BT_emp)
var_comp[2,] = diag(var_coeff_BT_res)
var_comp[3,] = diag(var_coeff_BT_wild)
colnames(var_comp) = c("Sepal.Length Int.", "Sepal.Length Slope", 
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
  lm_BT = glm(Y_BT ~ admissions$gpa, family='binomial')
  slope_GPA_BT[i] = lm_BT$coefficients[2]
}
# histogram of BT sample of model slope
hist(slope_GPA)
# construct 90% CI
quantile(slope_GPA_BT, c(0.05, 0.95))


# b Use both param. and empirical BT to estimate SE of the intercept, slope of admit~GRE, 
# and slope of admit~GPA
# empirical bootstrap
coeff_emp_BT = matrix(NA, B, 3)
for (i in 1:B){
  # get BT samples
  BT_i = sample(n, n, replace = T)
  admit_BT = admissions$admit[BT_i]
  gre_BT = admissions$gre[BT_i]
  gpa_BT = admissions$gpa[BT_i]
  # fit models from BT samples, obtain intecept, slope of GRE model, slope of GPA model
  lm_GRE_BT = lm(admit_BT ~ gre_BT)
  lm_GPA_BT = lm(admit_BT ~ gpa_BT)
  coeff_emp_BT[i, 1] = lm_GRE_BT$coefficients[1]
  coeff_emp_BT[i, 2] = lm_GRE_BT$coefficients[2]
  coeff_emp_BT[i, 3] = lm_GPA_BT$coefficients[2]
}


# parametric bootstrap
# fit linear models from data
lm_GRE = glm(admit ~ gre, admissions, family="binomial")
lm_GPA = glm(admit ~ gpa, admissions, family="binomial")
coeff_param_BT = matrix(NA, B, 3)
# estimate model params
p0_GRE = predict(lm_GRE, type = 'response')
p0_GPA = predict(lm_GPA, type = 'response')
for (i in 1:B){
  # bootstrap sample for admit or not
  Y_BT = rbinom(n, 1, p0)
  # Fit model from BT samples, get coefficients
  lm_GRE_BT = glm(Y_BT ~ admissions$gre, family='binomial')
  lm_GPA_BT = glm(Y_BT ~ admissions$gpa, family='binomial')
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




