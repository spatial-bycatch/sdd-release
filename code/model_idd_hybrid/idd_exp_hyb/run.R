

library(ggplot2)
library(rstan)
library(sdd)

# source run settings
source("../iter.R")

# get model definition
mdl_def <- strsplit(getwd(), "/")[[1]]
mdl_def <- mdl_def[(length(mdl_def) - 2):length(mdl_def)]

mdl_name <- mdl_def[2]
mdl_spp  <- mdl_def[3]

# path to results folders
dat_path <- "../../../../data"
run_path <- "."
stn_path <- "../.."
res_path <- file.path('../../../../results', paste(mdl_def, collapse = "/"))
fig_path <- file.path('../../../../figures', paste(mdl_def, collapse = "/"))

if (!dir.exists(res_path)) dir.create(res_path, recursive = TRUE)
if (!dir.exists(fig_path)) dir.create(fig_path, recursive = TRUE)

save(dat_path, run_path, res_path, fig_path, stn_path, file = "directories.rda")

# load data
if (file.exists(file.path(run_path, 'inputs.rda'))) {
    load(file.path(run_path, 'inputs.rda'))
} else {
    stop("missing 'inputs.rda'")
}    

# compile model
mdl <- stan_model(file = file.path(stn_path, paste0(mdl_name,'.stan')))

# get and save initial values
# from max a posterior estimate
message("Getting initial values")
if (file.exists(file.path(run_path, 'inits.rda'))) {
  load(file.path(run_path, 'inits.rda'))
} else {
  
  mdl.tmp <- optimizing(mdl, mdl.dat, verbose = TRUE, as_vector = FALSE)
  
  mdl.ini <- list()
  mdl.ini$density_log <- mdl.tmp$par$density_log
  mdl.ini$pi_log      <- mdl.tmp$par$pi_log
  mdl.ini$reg_par     <- mdl.tmp$par$reg_par
  mdl.ini$sigma       <- mdl.tmp$par$sigma
  
  save(mdl.ini, file = file.path(run_path, 'inits.rda'))
}       

# MCMC run
pars_exclude <- c("density_hat_log", "omega", "mu_log", "mu_log_nz", "sigma_nz", "density_log_ni", "density_log_is")

mdl.fit <- sampling(mdl, mdl.dat, init = function() mdl.ini, iter = n_iter, cores = n_cores, chains = n_chains, thin = n_thin, pars = pars_exclude, include = FALSE)

mdl.ini <- posterior(mdl.fit, pars = c("density_log", "pi_log", "reg_par", "sigma", "rho", "tau"), fun = "median")
save(mdl.ini, file = file.path(run_path, "inits.rda"))

# number of retained samples
nit <- n_iter_out
        
# save fit
save(mdl, mdl.dat, mdl.ini, mdl.fit, nit, file = file.path(res_path, 'fit.rda'))
    
# trace diganostic plots
gg <- traceplot(mdl.fit, pars = c("density_trace", "catchability_trace", "error_trace"))
gg <- gg + facet_wrap(~parameter, scales = "free_y") + theme_bw(base_size = 12)
gg <- gg + ggtitle("Trace summary statistics")
ggsave(gg, file = file.path(res_path, "trace_summaries.png"), width = 10)


# PLOT POSTERIOR DISTRIBUTIONS
dfr <- posterior(mdl.fit, pars = "pi_log", dim.names = list(list(iter = 1:nit, group = c(X.dimnames$group, "Prior"), label = c("Encounter Rate", "Efficiency"))), melt = TRUE)[[1]]

gg <- ggplot(dfr) + 
    geom_histogram(aes(x = value, y = after_stat(density)), position = "identity", binwidth = 0.5) + 
    labs(x = "", y = "") + theme_bw(base_size = 20) + facet_grid(group~label) + theme(axis.text.y = element_blank())
ggsave(gg, file = file.path(res_path, "posteriors.png"), width = 10)




