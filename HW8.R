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
hist(eruptions, main ="Distribution of Eruptions", probability=T)
lines(kde_2, lwd=2, col="blue")
legend("topleft", col=c("blue"), c("h=0.3"), lwd=2)

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


