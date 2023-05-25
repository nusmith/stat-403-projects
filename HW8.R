# Problem 2-- Analyzing eruptions from faithful dataset
eruptions = faithful$eruptions
h_values = c(0.1, 0.3, 0.9)

# a
# Using R's density function to form KDEs
kde_1 = density(eruptions, bw = 0.1)
kde_2 = density(eruptions, bw = 0.3)
kde_3 = density(eruptions, bw = 0.9)
# Plot to compare effects of different bandwidths
plot(kde_1, lwd=2, col="red", main = "KDE for varying bandwidths",
     xlab="x", ylab="fh(x)")
lines(kde_2, lwd=2, col="blue")
lines(kde_3, lwd=2, col="green")
legend("top", col=c("red", "blue", "green"), c("h=0.1", "h=0.3", "h=0.9"), lwd=2)

# b
# Compare KDE w/ h=0.3 to histogram
hist(eruptions, main ="Distribution of Eruptions", probability=T, breaks=20)
lines(kde_2, lwd=2, col="blue")
legend("top", col=c("blue"), c("KDE w/ h=0.3"), lwd=2)

# c
# Construct confidence band for KDE w/ h=0.3
h = 0.3
kde = density(eruptions, bw=h, from=0, to=7, n=1000)
n = length(eruptions)
B = 10000
dens_grid_BT = matrix(NA, nrow=B, ncol=length(kde$x))
for (i in 1:B) {
  # BT sample for indices from original sample
  BT_i = sample(n, n, replace = T)
  sample_BT = eruptions[BT_i]
  # Construct KDE from BT sample
  kde_BT = density(sample_BT, bw=h, from=0, to=7, n=1000)
  # Find estimated density at each point in KDE
  dens_grid_BT[i,] = kde_BT$y
}
# Find BT variance and SE for each point in KDE grid
# Note that the diagonal of var(grid) will be the variance of the density at 
# each x point in the grid
dens_var = diag(var(dens_grid_BT))
dens_se = sqrt(dens_var)


plot(kde, lwd=2, col="magenta", main="95% Conf. Band of KDE Using Bootstrap", ylim=c(0,0.6),
     xlab="x", ylab="Density")
polygon(x=c(kde$x, rev(kde$x)), y=c(kde$y - qnorm(0.975)*dens_se, rev(kde$y + qnorm(0.975)*dens_se)),
        col="lightpink", border="lightpink")
# KDE
lines(kde, lwd=3, col="magenta")
# lower line
lines(x = kde$x, y=kde$y - qnorm(0.975)*dens_se, lwd=2, col="pink", lty="dashed")
# upper line
lines(x = kde$x, y=kde$y + qnorm(0.975)*dens_se,lwd=2, col="pink", lty="dashed")
legend("topleft", col=c("magenta", "lightpink"), legend=c("KDE", "95% confidence band"), 
       lwd=2, cex=0.8)



# -----------------------------------------------------------------------------------------------
# Problem 3
# MC simulation to analyze error of KDE
# Target density p: double exponential distribution
p = function(x){
  0.5*exp(-1*abs(x))
}
x = seq(-1,1,0.1)
plot(x, p(x), xlim=c(-1,1), ylim=c(0,1.1))
# CDF of target density P
P = function(x){
  l = length(x)
  y = rep(NA, l)
  for (i in 1:l){
    if (x[i] < 0){
      y[i]=0.5*exp(x[i])
    }
    else {
      y[i]=1-0.5*exp(-1*x[i])
    }
  }
  y
}
lines(x, P(x))
# Inverse of CDF of target density
# takes in a u(0,1) RV and get the +/- from the sign of u-0.5
P_inv = function(x){
  l = length(x)
  y = rep(NA, l)
  for (i in 1:l){
    if (x[i] < 0){
      y[i]= log(2*x[i])
    }
    else {
      y[i]= -log(-2*(x[i]-1))
    }
  }
  y
}
# takes in a u ~ Unif(0,1) sample
P_inv = function(u){
  - sign(u-0.5)*log(1-2*abs(u-0.5))
}
# Unif sample
sample_u = runif(1000, min=0, max=1)
sample_p = P_inv(sample_u)
# Plot the histogram and overlay the true density 
hist(sample_p, breaks=30, probability=T, main = "Histogram of p(x) samples", ylim=c(0,0.5))
x = seq(-10,10,0.1)
lines(x, p(x), col="magenta", lwd=2.5)
# Make KDE from sample
kde_de = density(sample_p, bw=0.2)
plot(kde_de, col="blue", lwd=2, main = "Double Exponential Distribution", ylim=c(0, 0.5))
lines(x, p(x), col="forestgreen", lwd=2)
legend("topleft", col=c("blue", "forestgreen"), legend=c("KDE", "p(x) True"), 
       lwd=2, cex=0.9)
# Find Mean Integrated Square Error (MISE) for h E (0.05, 0.5)
h_seq = seq(0.05, 0.5, 0.05)
# MC Simulation to estimate error
N = 10000
mise_h_seq = rep(NA, length(h_seq))
for (i_h in 1:length(h_seq)) {
  h0 = h_seq[i_h]
  kde_de_MC = matrix(NA, nrow=N, ncol=1000)
  for (i in 1:N) {
    # Get sample from distribution, generate KDE
    sample_u = runif(1000, min=0, max=1)
    sample_p = P_inv(sample_u)
    kde_de = density(sample_p, bw=h0, n=1000, from=-6, to=6)
    kde_de_MC[i,]=kde_de$y
  }
  true_p = p(kde_de$x)
  # Get MSE
  kde_de_mse = colSums((t(t(kde_de_MC)-true_p))^2)/N
  # Save MISE
  mise_h_seq[i_h] = sum(kde_de_mse)
}
# Plot results to see how h-value (bandwidth) affects error in KDE
plot(x=h_seq, y=mise_h_seq, type="l", lwd=3, xlab="h (bandwidth)", ylab ="MISE",
     main= "MISE vs h", col="purple")





