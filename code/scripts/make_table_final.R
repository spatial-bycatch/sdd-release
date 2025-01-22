

library(plyr)
library(dplyr)
library(xtable)
library(reshape2)

spps <- c("SSI", "BBE", "ETB", "SND")

pars <- c("catchability", "density_hat")

wd <- getwd()

aa <- list()
bb <- list()

ff <- function(x) paste0(format(round(quantile(x, 0.5), 2), nsmall = 2), " (", format(round(quantile(x, 0.05), 2), nsmall = 2), " -- ", format(round(quantile(x, 0.95), 2), nsmall = 2), ")")

# SURVEY ONLY MODEL

for (spp in spps) {
  
  def <- "../model_sdd_hybrid/sdd_exp_hyb_reg"
  
  setwd(file.path(def, spp))
  
  load("directories.rda")
  
  load(file.path(run_path, "inputs.rda"))
  load(file.path(res_path, "fit.rda"))
  
  aa[[spp]] <- sdd::posterior(mdl.fit, pars = pars, melt = TRUE)
  
  setwd(wd)
}

aa <- lapply(aa, function(x) { x[["catchability"]] <- data.frame(group = "SURVEY", value = ff(x[["catchability"]][, "value"])); x})
aa <- lapply(aa, function(x) { x[["density_hat"]]  <- ff(x[["density_hat"]][, "value"]); x})

sdd_catchability <- ldply(lapply(aa, function(x) x[["catchability"]]), .id = "Species")
sdd_density      <- ldply(lapply(aa, function(x) x[["density_hat"]]), .id = "Species")
colnames(sdd_density)[2] <- "DENSITY"

sdd_catchability <- data.frame(model = "Survey only", acast(sdd_catchability, Species ~ group))
sdd_catchability$Species <- rownames(sdd_catchability)

sdd_tab <- data.frame(merge(sdd_catchability, sdd_density), "HOKHAKLIN_BT" = NA)


# INTEGRATED MODEL

bb <- list()

for (spp in spps) {
  
  def <- "../model_idd_hybrid/idd_exp_hyb_reg"
  
  setwd(file.path(def, spp))
  
  load("directories.rda")
  
  load(file.path(run_path, "inputs.rda"))
  load(file.path(res_path, "fit.rda"))
  
  bb[[spp]] <- sdd::posterior(mdl.fit, pars = pars, melt = TRUE, dim.names = list(catchability = list(iter = 1:nit, group = X.dimnames$group), density_hat = list(iter = 1:nit, grid = X.dimnames$grid)))
  
  setwd(wd)
}

bb <- lapply(bb, function(x) { x[["catchability"]] <- ddply(x[["catchability"]], "group", summarise, value = ff(value)); x})
bb <- lapply(bb, function(x) { x[["density_hat"]]  <- ff(x[["density_hat"]][, "value"]); x})

idd_catchability <- ldply(lapply(bb, function(x) x[["catchability"]]), .id = "Species")
idd_density      <- ldply(lapply(bb, function(x) x[["density_hat"]]), .id = "Species")
colnames(idd_density)[2] <- "DENSITY"

idd_catchability <- data.frame(model = "Integrated", acast(idd_catchability, Species ~ group))
idd_catchability$Species <- rownames(idd_catchability)

idd_tab <- merge(idd_catchability, idd_density)

catchability <- merge(sdd_catchability, idd_catchability, by = c("Species", "model"), all = TRUE, suffixes = c(".survey", ".integrated"))


# FINAL TABLE

tab <- rbind(sdd_tab, idd_tab)
tab <- tab[order(tab$Species),]

print(xtable(tab[, c(1,2,3,5,4)], caption = "Summary results from final model selections", label = "tab:results"), file = file.path("../../tables", "model_results_draft.tex"), include.rownames = FALSE)


