# generate plot with random data, save imagine in output folder
data = data.frame(x = rnorm(100), y = rnorm(100))
png("out/plot.png")
plot(data$x, data$y)
dev.off()
