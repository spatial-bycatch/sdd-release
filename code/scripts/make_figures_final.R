
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(sdd)

spps <- c("SSI", "BBE", "ETB", "SND")

mdl_labels <- data.frame(
    model = c("idd_exp",
              "idd_exp_reg",
              "idd_exp_hyb", 
              "idd_exp_hyb_reg", 
              "idd_wei",
              "idd_wei_reg",
              "idd_wei_hyb",
              "idd_wei_hyb_reg"),
    label = c("E1",
              "E2",
              "E3", 
              "E4", 
              "W1",
              "W2",
              "W3",
              "W4")
)

wd <- getwd()

mdl_use <- "idd_exp_hyb_reg"


#################################
# PLOT REGRESSION RELATIONSHIP  #
#################################

lst1 <- list()
lst2 <- list()

for (spp in spps) {
    
    def <-  mdl_labels$label[match(mdl_use, mdl_labels$model)]
    
    setwd(file.path("..", "model_idd_hybrid", mdl_use, spp))
    
    load("directories.rda")
    
    load(file.path(run_path, "inputs.rda"))
    load(file.path(res_path, "fit.rda"))
    
    message("Plotting regression for ", def, "_", spp)
    
    dfr  <- posterior(mdl.fit, pars = "reg_par")[[1]]
    
    dpts <- log(mdl.dat$dpt)
    dpts <- (dpts - mean(dpts)) / sd(dpts)
    
    if (length(grep("reg", mdl_use)) > 0) { 
        d_log <- sapply(unique(dpts), function(x) apply(dfr, 1, function(y) y[1] + y[2] * x + y[3] * x^2))
    } else {
        d_log <- sapply(unique(dpts), function(x) apply(dfr, 1, function(y) y))
    }
    dimnames(d_log) <- list(iter = 1:nit, dpt = unique(mdl.dat$dpt))
    d_log <- melt(d_log)
    
    dfr <- posterior(mdl.fit, pars = "density_hat", dim.names = list(list(iter = 1:nit, grid = X.dimnames$grid)), fun = "quantiles")[[1]]
    dfr <- adply(dfr, .margins = 2)
    
    dfr2 <- ddply(dat, .(grid), summarise, dpt = mean(depth))
    
    dfr <- merge(dfr, dfr2, by = "grid")
    
    colnames(dfr)[2:4] <- c("med", "low", "upp")
    
    lst1[[spp]] <- d_log
    lst2[[spp]] <- dfr
    
    setwd(wd)
}

dfr1 <- ldply(lst1, .id = "species")
dfr2 <- ldply(lst2, .id = "species")

gg <- ggplot(dfr1) + 
    stat_summary(aes(x = dpt, y = value), geom = "line", col = 4, fun = median, linewidth = 2) +
    stat_summary(aes(x = dpt, y = value), geom = "ribbon", fill = 4, alpha = 0.2, fun = median, fun.max = function(x) quantile(x, 0.95), fun.min = function(x) quantile(x, 0.05)) +
    geom_point(data = dfr2, aes(x = dpt, y = log(med))) +
    geom_errorbar(data = dfr2, aes(x = dpt, ymin = log(low), ymax = log(upp)), alpha = 0.2) +
    theme_bw(base_size = 20) + labs(x = "Depth (metres)", y = "Density (log-scale)") +
    facet_wrap(~species)
fig_path <- "../../figures" 
ggsave(gg, file = file.path(fig_path, "idd_regression.pdf"), width = 10)
#ggsave(gg, file = file.path(fig_path, "idd_regression.png"), width = 10)

##################
# PLOT SPATIALLY #
##################

library(spdep)
source("../../data/grids.R")
library(nzPlot)

dd <- function(x) switch(x, "SND" = 20, "ETB" = 10, "SSI" = 1.5, "BBE" = 2.5)

for (spp in spps) {
    
    def <-  mdl_labels$label[match(mdl_use, mdl_labels$model)]
    
    setwd(file.path("..", "model_idd_hybrid", mdl_use, spp))
    
    load("directories.rda")
    
    load(file.path(run_path, "inputs.rda"))
    load(file.path(res_path, "fit.rda"))
    load(file.path(dat_path, paste0("survey_", spp, "_grid025.rda")))
    
    message("Plotting maps for ", def, "_", spp)
    
    dfr <- posterior(mdl.fit, pars = "density_hat", dim.names = list(list(iter = 1:nit, grid = X.dimnames$grid)), fun = "mean", melt = TRUE)[[1]]
    
    dat <- ddply(subset(dat, group == "SURVEY"), .(grid), summarise, cpua = mean(biomass / sweptArea))
    dat$grid <- as.integer(as.character(dat$grid))
    
    ss <- mean(dfr$value)
    dfr$value <- 1 / (1 + exp(-(dfr$value - ss) / dd(spp)))
    ss <- mean(dat$cpua)
    dat$value <- 1 / (1 + exp(-(dat$cpua - ss) / dd(spp)))
    
    # plot grids
    windows(width = 10, height = 10)
    par(mfrow = c(2, 1), mar = c(0, 0.5, 0, 0.5), oma = c(0, 0, 0, 0))
    nz(xlim = c(min_lon, max_lon), ylim = c(-46, -42), fill.col = "grey", xaxt = "n", yaxt = "n")
    box()
    nz.depth(500, col = "lightblue")
    nz.depth(1000, col = "lightblue")
    nz.depth(1500, col = "lightblue")
    for (i in 1:length(survey_grid_def)) {
        
        j <- names(survey_grid_def)[i]
        
        if (j %in% X.dimnames$grid) {
            
            k <- which(rownames(dfr) == j)
            
            with(survey_grid_def[[i]], nz.polygon(x = lon, y = lat, border = 2, col = grey(1 - dfr$value[k])))  
        }
    }
    mtext("Predicted density", adj = 0, padj = -0.5, cex = 1.5)
    nz(xlim = c(min_lon, max_lon), ylim = c(-46, -42), fill.col = "grey", xaxt = "n", yaxt = "n")
    box()
    nz.depth(500, col = "lightblue")
    nz.depth(1000, col = "lightblue")
    nz.depth(1500, col = "lightblue")
    for (i in 1:length(survey_grid_def)) {
        
        j <- names(survey_grid_def)[i]
        
        if (j %in% X.dimnames$grid) {
            
            k <- which(dat$grid == j)
            
            with(survey_grid_def[[i]], nz.polygon(x = lon, y = lat, border = 2, col = grey(1 - dat$value[k])))  
        }
    }
    mtext("Empirical survey catch rate", adj = 0, padj = -0.5, cex = 1.5)
    fig_path <- "../../../../figures"
    savePlot(filename = file.path(fig_path, paste0("idd_map_density_", tolower(spp),".pdf")), type = "pdf")
    dev.off()

    setwd(wd)
}


#######################
message("End of plots")



