library(archetypal)
library(ggplot2)
# generate random matrix
a <- matrix(runif(100), ncol = 2)
k <- 4
irow <- sample(1:nrow(a), 1)
irows <- FurthestSum(a, irow = irow, k = k)
# create array of length nrow(a) with all 0 except for the k furthest points from indices_furthestPoints which are 1
indices_furthestPoints <- rep(0, nrow(a))
indices_furthestPoints[irows] <- 1


# plot points and signal differently incdices_furthestPoints
ggplot() +
  geom_point(aes(x = a[, 1], y = a[, 2]), color = "blue") +
  geom_point(aes(x = a[irows, 1], y = a[irows, 2]), color = "red") +
  geom_point(aes(x = a[irow, 1], y = a[irow, 2]), color = "green") +
  theme_minimal()


library(archetypes)
# features must be on columns
family <- archetypesFamily(which = "robust", initfn = archetypes:::make.fix.initfn(indizes = irows))
aa <- archetypes(a, k = k, family = family)
