
suppressMessages({

library(ggplot2)
library(rstan)

})

# get model definition
mdl_def <- strsplit(getwd(), "/")[[1]]
mdl_def <- mdl_def[(length(mdl_def) - 2):length(mdl_def)]

mdl_name <- mdl_def[2]
mdl_spp  <- mdl_def[3]

# paths to folders
load(file = "directories.rda")

# load data
if (file.exists(file.path(run_path, 'inputs.rda'))) {
    load(file.path(run_path, 'inputs.rda'))
} else {
    stop("missing 'inputs.rda'")
}  

# load initial values
if (file.exists(file.path(run_path, 'inits.rda'))) {
    load(file.path(run_path, 'inits.rda'))
} else {
    stop("missing 'inits.rda'")
} 

# compile model
mdl <- stan_model(file = file.path(stn_path, paste0(mdl_name,'.stan')))

# run model
mdl.tmp <- optimizing(mdl, mdl.dat, init = mdl.ini, verbose = FALSE, as_vector = FALSE)

# check
if (mdl.tmp$return_code > 0) {
	stop("failed ", mdl_spp, ":", mdl_name)
} else {
	message("checked ", mdl_spp, ":", mdl_name)
}