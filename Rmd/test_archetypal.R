library(archetypal)
data("gallupGPS6")
data("wd25")

df <- gallupGPS6
df <- wd25

maxkappas <- 15
method <- "furthestsum"
nworkers <- 10
df <- gallupGPS6
rseed <- 2024
kappas <- 5:15

opt_kappas <- find_optimal_kappas(df, maxkappas = maxkappas, method = method, nworkers = nworkers)
opt_kappas <- find_optimal_kappas(df, maxkappas = maxkappas, nworkers = nworkers)
kappas <- opt_kappas$optimal_kappas

opt_params <- find_pcha_optimal_parameters(df, kappas = opt_kappas$optimal_kappas, method = method, nworkers = nworkers, rseed = rseed, plot = TRUE)

names(opt_params)
aa <- archetypal(df, opt_kappas$optimal_kappas, method = method, rseed = rseed, save_history = TRUE, nworkers = nworkers)

plot.archetypal(aa)
plot_archs(df, aa)

find_closer_points(df, 3, rseed = rseed)
