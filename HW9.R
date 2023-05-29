library(plyr)
library(ggplot2)

# Problem 1
# rock data: area (response), peri (covariate)
data = data.frame("peri" = rock$peri, "area" = rock$area)
# a-- linear regression 
plot(data$peri, data$area, main="Rock Area vs Peri", xlab="Peri", ylab="Area")
abline(lm(area ~ peri, data),lwd=2, col="blue")
legend("topleft", legend="Linear Regression", col="blue", cex=0.8, lwd=2)
# b-- Gaussian Kernel regression w/ h = 500
h = 500
data_kreg = ksmooth(x=data$peri,y=data$area, kernel = "normal", bandwidth = h)
plot(data$peri, data$area, main="Rock Area vs Peri", xlab="Peri", ylab="Area")
lines(data_kreg,lwd=2, col="magenta")
legend("topleft", legend="Gaussian Kernel Regression", col="magenta", cex=0.8, lwd=2)
# c-- compare Gaussian Kreg across h values
h_seq = c(250, 500, 1000)
colors = c("darkgreen", "turquoise", "dodgerblue")
plot(data$peri, data$area, main="Rock Area vs Peri", xlab="Peri", ylab="Area")
for (i in 1:length(h_seq)){
  data_kreg = ksmooth(x=data$peri,y=data$area, kernel = "normal", bandwidth = h_seq[i])
  lines(data_kreg,lwd=2, col=colors[i])
}
legend("topleft", legend=h_seq, col=colors, cex=0.8, lwd=2, title="h values for Gauss. KReg")
# d-- Cross-validation to find h that minimizes error
n = length(data$area)
#number of iterations
N_cv = 100
# number of folds
k = 3
#determine which of the 3 groups each observation belongs to
cv_lab = sample(n,n,replace=F) %% k
h_seq = seq(250, 1000, 50)
CV_err_h = rep(0,length(h_seq))
for(i_tmp in 1:N_cv){
  CV_err_h_tmp = rep(0, length(h_seq))
  cv_lab = sample(n,n,replace=F) %% k
  for(i in 1:length(h_seq)){
    h0 = h_seq[i]
    CV_err =0 
    for(i_cv in 1:k){
      #take samples in that group
      w_val = which(cv_lab==(i_cv-1))
      #training set
      peri_tr = data$peri[-w_val]
      area_tr = data$area[-w_val]
      #test set
      peri_val = data$peri[w_val]
      area_val = data$area[w_val]
      #fit regression w training set and get fitted values for test data
      kernel_reg = ksmooth(x = peri_tr,y=area_tr,kernel = "normal",bandwidth=h0, x.points=peri_val)
      #get CV_err with test set  (WARNING! The ksmooth() function will order the x.points from smallest to the largest!)
      CV_err = CV_err+mean((area_val[order(peri_val)]-kernel_reg$y)^2,na.rm=T)
    }
    #average CV error for bandwidth i for iteration i_tmp
    CV_err_h_tmp[i] = CV_err/k
  }
  #add to average from other iterations
  CV_err_h = CV_err_h+CV_err_h_tmp
}
#get average CV error for each bandwidth across all iterations
CV_err_h = CV_err_h/N_cv
#look at CV error for different bandwidths
plot(h_seq,CV_err_h, type="b", lwd=4, col="magenta", main= "CV Error vs h", xlab="Smoothing Bandwidth", ylab="CV Error")
h_opt = h_seq[which(CV_err_h==min(CV_err_h))]
# e-- Empirical bootstrap to determine 95% CI for kernel regression (w/ h = 500)
h = 500
B = 10000
# generate x (peri) values and fit to model (predict area)
x_seq = seq(250, 5000, 50)
area_pred = ksmooth(data$peri, data$area, kernel="normal", bandwidth=h, x.points=x_seq)$y
n = length(data$area)
area_kreg_BT = matrix(NA, nrow=B, ncol=length(area_pred))
for (i in 1:B){
  # get sample indices
  w = sample(n, n, replace = T)
  peri_BT = data$peri[w]
  area_BT = data$area[w]
  area_kreg_BT[i,] = ksmooth(peri_BT, area_BT, kernel="normal", bandwidth=h, x.points=x_seq)$y
}
# Get BT SD
area_kreg_SD = sapply(1:length(x_seq), function(x){
  sd(area_kreg_BT[,x], na.rm=T)
})
# Show CI
plot(data$peri,data$area, pch=1, cex=0.5, col="gray", 
  main="Area vs Peri, 95% CI for Kernel Regression",
  xlab = "Peri",
  ylab = "Area")
lines(x_seq, area_pred+qnorm(0.975)*area_kreg_SD, lwd=2, col="blue", lty=2)
lines(x_seq, area_pred-qnorm(0.975)*area_kreg_SD, lwd=2, col="blue", lty=2)
lines(x_seq, area_pred, lwd=2, col="black")
legend("topleft", legend=c("h=500 KReg", "95% CI"), col=c("black", "blue"), cex=0.8, lwd=2)
polygon(x=c(x_seq, rev(x_seq)), y=c(area_pred - qnorm(0.975)*area_kreg_SD, rev(area_pred + qnorm(0.975)*area_kreg_SD)),
        col="lightblue")
lines(x_seq, area_pred, lwd=2, col="black")
legend("topleft", legend=c("h=500 KReg", "95% CI"), col=c("black", "lightblue"), cex=0.8, lwd=2)





# Problem 2
# Still using rock data. Regressogram method
# b -- scatter plot + regressogram curve
plot(data$peri,data$area, pch=1, cex=0.5, col="gray", main="Area vs Peri",
     xlab = "Peri", ylab = "Area")
# m estimates the regression function (regressogram func)
# Using bins B1 = (0, 1000], B2 = (1000, 2000], B3 = (2000, 3000], B4 = (3000, 4000], B5 = (4000, 5000]
m = function(x_arr) {
  mx_arr = rep(NA, length(x_arr))
  for (i in 1:length(x_arr)){
    x = x_arr[i]
    lower = round_any(x, 1000, f=floor)
    upper = round_any(x, 1000, f=ceiling)
    mx = sum(data$area[data$peri < upper & data$peri > lower]) / length(data$area[data$peri < upper & data$peri > lower])
    mx_arr[i] = mx
  }
  mx_arr
}
x = seq(0, 5000, 50)
lines(x, m(x), col="red", lwd=2, type="s")
legend("topleft", legend="Regressogram Linear Regression", col="red", cex=0.8, lwd=2)

