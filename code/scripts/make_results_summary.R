

library(plyr)
library(ggplot2)
library(tidyr)
library(dplyr)

wd <- getwd()

spps <- c("SSI", "BBE", "ETB", "SND")

source("posterior.R")

#############
# SDD MODEL #
#############

mdls <- c("model_sdd/sdd_exp",
          "model_sdd/sdd_exp_reg",
          "model_sdd_hybrid/sdd_exp_hyb", 
          "model_sdd_hybrid/sdd_exp_hyb_reg", 
          "model_sdd/sdd_wei",
          "model_sdd/sdd_wei_reg",
          "model_sdd_hybrid/sdd_wei_hyb",
          "model_sdd_hybrid/sdd_wei_hyb_reg")

mdl_labels <- data.frame(
    model = c("sdd_exp",
              "sdd_exp_reg",
              "sdd_exp_hyb", 
              "sdd_exp_hyb_reg", 
              "sdd_wei",
              "sdd_wei_reg",
              "sdd_wei_hyb",
              "sdd_wei_hyb_reg"),
    label = c("E1",
              "E2",
              "E3", 
              "E4", 
              "W1",
              "W2",
              "W3",
              "W4")
)

aa <- list()

for (i in spps) {
    
    aa[[i]] <- list()
    
    for (j in mdls) {
        
        mdl_label <- mdl_labels$label[match(strsplit(j, "/")[[1]][2], mdl_labels$model)]
        
        message(i, ": SDD ", mdl_label)
        
        aa[[i]][[basename(j)]] <- list()
        
        setwd(file.path("..", j, i))
        
        load("directories.rda")
        
        load(file.path(res_path, "fit.rda"))
        load(file.path(run_path, "inputs.rda"))
        
        x <- posterior(mdl.fit, pars = c("catchability", "cpue_emp", "cpue_hat", "cpue_sim", "pnzero_sim", "pnzero_hat", "pnzero_emp", "cpua_hat", "density_hat", "reg_par"))
        
        nit <- dim(x[["catchability"]])[1]
        
        dat$model  <- "SDD"
        dat$method <- "SURVEY"
        
        mae_1    <- function(x, y) mean(abs(x - y))
        mae_2    <- function(x, y) apply(abs(x - y), 2, mean)
        pvalue_1 <- function(x, y) { p <- y > x; p[y == x] <- 0.5; mean(p)}
        pvalue_2 <- function(x, y) { p <- y > x; p[y == x] <- 0.5; apply(p, 2, mean)}
        
        get_density_mcmc <- function() {
            x <- rstan::extract(mdl.fit, pars = "density_hat", permute = FALSE) 
            x <- apply(x, 1:2, mean)
            x <- adply(x, .margins = 1:2)
            colnames(x)[ncol(x)] <- "value"
            return(x)
        }
        
        get_catchability_mcmc <- function() {
            x <- rstan::extract(mdl.fit, pars = "catchability", permute = FALSE) 
            x <- apply(x, 1:2, mean)
            x <- adply(x, .margins = 1:2)
            colnames(x)[ncol(x)] <- "value"
            return(x)
        }
        
        get_rhat <- function(x) {
            
            pivot_wider(x, names_from = "chains") %>% select(contains("chain")) %>% as.matrix() %>% posterior::rhat_basic()
        }
        
        get_ess <- function(x) {
            
            pivot_wider(x, names_from = "chains") %>% select(contains("chain")) %>% as.matrix() %>% posterior::ess_basic()
        }
        
        dimnames(x[["density_hat"]]) <- list(iter = 1:nit, cell = X.dimnames$grid)
        
        dimnames(x[["cpue_emp"]])   <- list(iter = 1:nit, cell = X.dimnames$grid)
        dimnames(x[["cpue_sim"]])   <- list(iter = 1:nit, cell = X.dimnames$grid)
        dimnames(x[["pnzero_emp"]]) <- list(iter = 1:nit, cell = X.dimnames$grid)
        dimnames(x[["pnzero_sim"]]) <- list(iter = 1:nit, cell = X.dimnames$grid)
        
        # sampling
        aa[[i]][[basename(j)]][["sample_size"]] <- ddply(dat, .(grid, model, method), summarise, effort_n = sum(effort), effort_area = sum(sweptArea))
        
        # results
        aa[[i]][[basename(j)]][["catchability"]]      <- data.frame(iter = 1:nit, model = "SDD", label = mdl_label, method = "SURVEY", value = x[["catchability"]])
        aa[[i]][[basename(j)]][["density_mean"]]      <- data.frame(iter = 1:nit, model = "SDD", label = mdl_label, value = apply(x[["density_hat"]], 1, mean))
        aa[[i]][[basename(j)]][["density_intercept"]] <- data.frame(iter = 1:nit, model = "SDD", label = mdl_label, value = if(length(dim(x[["reg_par"]])) > 1) x[["reg_par"]][,1] else x[["reg_par"]])
        aa[[i]][[basename(j)]][["density_mcmc"]]      <- data.frame(model = "SDD", label = mdl_label, get_density_mcmc())
        aa[[i]][[basename(j)]][["catchability_mcmc"]] <- data.frame(model = "SDD", label = mdl_label, method = "SURVEY", get_catchability_mcmc())
        
        # fits
        aa[[i]][[basename(j)]][["catch_rate"]] <- data.frame(model = "SDD", label = mdl_label, method = "SURVEY", cell = X.dimnames$grid, reshape2::melt(apply(x[["cpue_emp"]],   2, mean), value.name = "value_emp"), reshape2::melt(apply(x[["cpue_hat"]],   2, mean), value.name = "value_sim"), reshape2::melt(apply(x[["cpue_hat"]],   2, median), value.name = "value_hat"))
        aa[[i]][[basename(j)]][["catch_prob"]] <- data.frame(model = "SDD", label = mdl_label, method = "SURVEY", cell = X.dimnames$grid, reshape2::melt(apply(x[["pnzero_emp"]], 2, mean), value.name = "value_emp"), reshape2::melt(apply(x[["pnzero_hat"]], 2, mean), value.name = "value_sim"), reshape2::melt(apply(x[["pnzero_hat"]], 2, median), value.name = "value_hat"))
        aa[[i]][[basename(j)]][["catch_density"]] <- data.frame(cell = X.dimnames$grid, model = "SDD", label = mdl_label, value = apply(x[["density_hat"]], 2, mean))
        
        # distribution
        aa[[i]][[basename(j)]][["distribution"]] <- data.frame(model = "SDD", label = mdl_label, method = "SURVEY", reshape2::melt(x[["density_hat"]], value.name = "density_hat"), reshape2::melt(x[["cpue_emp"]], value.name = "cpue_emp")[,3,drop=FALSE], reshape2::melt(x[["cpue_sim"]], value.name = "cpue_sim")[,3,drop=FALSE])
        
        # diagnostics
        aa[[i]][[basename(j)]][["density_rhat"]]      <- data.frame(model = "SDD", label = mdl_label, value = round(get_rhat(get_density_mcmc()), 3))
        aa[[i]][[basename(j)]][["density_neff"]]      <- data.frame(model = "SDD", label = mdl_label, value = round(get_ess(get_density_mcmc()) / nit, 2))
        aa[[i]][[basename(j)]][["catchability_rhat"]] <- data.frame(model = "SDD", label = mdl_label, method = "SURVEY", value = round(get_rhat(get_catchability_mcmc()), 3))
        aa[[i]][[basename(j)]][["catchability_neff"]] <- data.frame(model = "SDD", label = mdl_label, method = "SURVEY", value = round(get_ess(get_catchability_mcmc()) / nit, 2))
        
        # error values
        aa[[i]][[basename(j)]][["cpue_mae"]]    <- data.frame(model = "SDD", label = mdl_label, method = "SURVEY", value = mae_1(x[["cpue_emp"]], x[["cpue_hat"]]))
        aa[[i]][[basename(j)]][["pnzero_mae"]]  <- data.frame(model = "SDD", label = mdl_label, method = "SURVEY", value = mae_1(x[["pnzero_emp"]], x[["pnzero_hat"]]))
        
        # p-values
        aa[[i]][[basename(j)]][["cpue_pvalue"]]    <- data.frame(model = "SDD", label = mdl_label, method = "SURVEY", value = pvalue_1(x[["cpue_emp"]], x[["cpue_sim"]]))
        aa[[i]][[basename(j)]][["pnzero_pvalue"]]  <- data.frame(model = "SDD", label = mdl_label, method = "SURVEY", value = pvalue_1(x[["pnzero_emp"]], x[["pnzero_sim"]]))
        aa[[i]][[basename(j)]][["density_pvalue"]] <- data.frame(model = "SDD", label = mdl_label, method = "SURVEY", value = pvalue_1(x[["density_hat"]], x[["cpua_hat"]]))
        
        setwd(wd)
    }
}

res_summary_sdd <- aa
save(res_summary_sdd, file = file.path("../../results/res_summary_sdd.rda"))

#############
# IDD MODEL #
#############

mdls <- c("model_idd/idd_exp",
          "model_idd/idd_exp_reg",
          "model_idd_hybrid/idd_exp_hyb", 
          "model_idd_hybrid/idd_exp_hyb_reg", 
          "model_idd/idd_wei",
          "model_idd/idd_wei_reg",
          "model_idd_hybrid/idd_wei_hyb",
          "model_idd_hybrid/idd_wei_hyb_reg")

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

bb <- list()

for (i in spps) {
    
    bb[[i]] <- list()
    
    for (j in mdls) {
        
        mdl_label <- mdl_labels$label[match(strsplit(j, "/")[[1]][2], mdl_labels$model)]
        
        message(i, ": IDD ", mdl_label)
        
        bb[[i]][[basename(j)]] <- list()
        
        setwd(file.path("..", j, i))
        
        load("directories.rda")
        
        load(file.path(res_path, "fit.rda"))
        load(file.path(run_path, "inputs.rda"))
        
        x <- posterior(mdl.fit, pars = c("catchability", "cpue_emp", "cpue_hat", "cpue_sim", "pnzero_sim", "pnzero_hat", "pnzero_emp", "cpua_hat", "density_hat", "reg_par"))
        
        nit       <- dim(x[["catchability"]])[1]
        id_survey <- which(X.dimnames$group == "SURVEY")
        id_comm   <- which(X.dimnames$group != "SURVEY")
        
        dat$model  <- "IDD"
        dat$method <- "SURVEY"
        dat$method[dat$group == id_comm] <- "COMM"
        
        mae_1    <- function(x, y) mean(abs(x - y))
        mae_2    <- function(x, y) apply(abs(x - y), 2, mean)
        pvalue_1 <- function(x, y) { p <- y > x; p[y == x] <- 0.5; mean(p)}
        pvalue_2 <- function(x, y) { p <- y > x; p[y == x] <- 0.5; apply(p, 2, mean)}
        
        get_density_mcmc <- function() {
            x <- rstan::extract(mdl.fit, pars = "density_hat", permute = FALSE) 
            x <- apply(x, 1:2, mean)
            x <- adply(x, .margins = 1:2)
            colnames(x)[ncol(x)] <- "value"
            return(x)
        }
        
        get_catchability_mcmc <- function(id) {
            x <- rstan::extract(mdl.fit, pars = "catchability", permute = FALSE) 
            x <- apply(x[,,id], 1:2, mean)
            x <- adply(x, .margins = 1:2)
            colnames(x)[ncol(x)] <- "value"
            return(x)
        }
        
        get_rhat <- function(x) {
            
            pivot_wider(x, names_from = "chains") %>% select(contains("chain")) %>% as.matrix() %>% posterior::rhat_basic()
        }
        
        get_ess <- function(x) {
            
            pivot_wider(x, names_from = "chains") %>% select(contains("chain")) %>% as.matrix() %>% posterior::ess_basic()
        }
        
        dimnames(x[["density_hat"]])   <- list(iter = 1:nit, cell = X.dimnames$grid)
        
        dimnames(x[["cpue_emp"]])   <- list(iter = 1:nit, group = X.dimnames$group, cell = X.dimnames$grid)
        dimnames(x[["cpue_sim"]])   <- list(iter = 1:nit, group = X.dimnames$group, cell = X.dimnames$grid)
        dimnames(x[["pnzero_emp"]]) <- list(iter = 1:nit, group = X.dimnames$group, cell = X.dimnames$grid)
        dimnames(x[["pnzero_sim"]]) <- list(iter = 1:nit, group = X.dimnames$group, cell = X.dimnames$grid)
        
        # sampling
        bb[[i]][[basename(j)]][["sample_size"]] <- ddply(dat, .(grid, model, method), summarise, effort_n = sum(effort), effort_area = sum(sweptArea))
        
        # results
        bb[[i]][[basename(j)]][["catchability"]]      <- rbind(data.frame(iter = 1:nit, model = "IDD", label = mdl_label, value = x[["catchability"]][,id_survey], method = "SURVEY"), data.frame(iter = 1:nit, model = "IDD", label = mdl_label, value = x[["catchability"]][,id_comm], method = "COMM"))
        bb[[i]][[basename(j)]][["density_mean"]]      <- data.frame(iter = 1:nit, model = "IDD", label = mdl_label, value = apply(x[["density_hat"]], 1, mean))
        bb[[i]][[basename(j)]][["density_intercept"]] <- data.frame(iter = 1:nit, model = "IDD", label = mdl_label, value = if(length(dim(x[["reg_par"]])) > 1) x[["reg_par"]][,1] else x[["reg_par"]])
        bb[[i]][[basename(j)]][["density_mcmc"]]      <- data.frame(model = "IDD", label = mdl_label, get_density_mcmc())
        bb[[i]][[basename(j)]][["catchability_mcmc"]] <- rbind(data.frame(model = "IDD", label = mdl_label, method = "SURVEY", get_catchability_mcmc(id_survey)), data.frame(model = "IDD", label = mdl_label, method = "COMM", get_catchability_mcmc(id_comm)))
        
        # fits
        bb[[i]][[basename(j)]][["catch_rate"]]    <- rbind(data.frame(model = "IDD", label = mdl_label, method = "SURVEY", cell = X.dimnames$grid, reshape2::melt(apply(x[["cpue_emp"]][,id_survey,],   2, mean), value.name = "value_emp"), reshape2::melt(apply(x[["cpue_hat"]][,id_survey,],   2, mean), value.name = "value_sim"), reshape2::melt(apply(x[["cpue_hat"]][,id_survey,],   2, median), value.name = "value_hat")), data.frame(model = "IDD", label = mdl_label, method = "COMM", cell = X.dimnames$grid, reshape2::melt(apply(x[["cpue_emp"]][,id_comm,],   2, mean), value.name = "value_emp"), reshape2::melt(apply(x[["cpue_hat"]][,id_comm,],   2, mean), value.name = "value_sim"), reshape2::melt(apply(x[["cpue_hat"]][,id_comm,],   2, median), value.name = "value_hat")))
        bb[[i]][[basename(j)]][["catch_prob"]]    <- rbind(data.frame(model = "IDD", label = mdl_label, method = "SURVEY", cell = X.dimnames$grid, reshape2::melt(apply(x[["pnzero_emp"]][,id_survey,], 2, mean), value.name = "value_emp"), reshape2::melt(apply(x[["pnzero_hat"]][,id_survey,], 2, mean), value.name = "value_sim"), reshape2::melt(apply(x[["pnzero_hat"]][,id_survey,], 2, median), value.name = "value_hat")), data.frame(model = "IDD", label = mdl_label, method = "COMM", cell = X.dimnames$grid, reshape2::melt(apply(x[["pnzero_emp"]][,id_comm,], 2, mean), value.name = "value_emp"), reshape2::melt(apply(x[["pnzero_hat"]][,id_comm,], 2, mean), value.name = "value_sim"), reshape2::melt(apply(x[["pnzero_hat"]][,id_comm,], 2, median), value.name = "value_hat")))
        bb[[i]][[basename(j)]][["catch_density"]] <- data.frame(cell = X.dimnames$grid, model = "IDD", label = mdl_label, value = apply(x[["density_hat"]], 2, mean))

        # distribution
        bb[[i]][[basename(j)]][["distribution"]] <- rbind(data.frame(model = "IDD", label = mdl_label, method = "SURVEY", reshape2::melt(x[["density_hat"]], value.name = "density_hat"), reshape2::melt(x[["cpue_emp"]][,id_survey,], value.name = "cpue_emp")[,3,drop=FALSE], reshape2::melt(x[["cpue_sim"]][,id_survey,], value.name = "cpue_sim")[,3,drop=FALSE]), data.frame(model = "IDD", label = mdl_label, method = "COMM", reshape2::melt(x[["density_hat"]], value.name = "density_hat"), reshape2::melt(x[["cpue_emp"]][,id_comm,], value.name = "cpue_emp")[,3,drop=FALSE], reshape2::melt(x[["cpue_sim"]][,id_comm,], value.name = "cpue_sim")[,3,drop=FALSE]))
        
        # diagnostics
        bb[[i]][[basename(j)]][["density_rhat"]]      <- data.frame(model = "IDD", label = mdl_label, value = round(get_rhat(get_density_mcmc()), 3))
        bb[[i]][[basename(j)]][["density_neff"]]      <- data.frame(model = "IDD", label = mdl_label, value = round(get_ess(get_density_mcmc()) / nit, 2))
        bb[[i]][[basename(j)]][["catchability_rhat"]] <- rbind(data.frame(model = "IDD", label = mdl_label, method = "SURVEY", value = round(get_rhat(get_catchability_mcmc(id_survey)), 3)), data.frame(model = "IDD", label = mdl_label, method = "COMM", value = round(get_rhat(get_catchability_mcmc(id_comm)), 3)))
        bb[[i]][[basename(j)]][["catchability_neff"]] <- rbind(data.frame(model = "IDD", label = mdl_label, method = "SURVEY", value = round(get_ess(get_catchability_mcmc(id_survey)) / nit, 2)), data.frame(model = "IDD", label = mdl_label, method = "COMM", value = round(get_ess(get_catchability_mcmc(id_comm)) / nit, 2))) 
        
        # error values
        bb[[i]][[basename(j)]][["cpue_mae"]]    <- rbind(data.frame(model = "IDD", label = mdl_label, method = "SURVEY", value = mae_1(x[["cpue_emp"]][,id_survey,], x[["cpue_hat"]][,id_survey,])), data.frame(model = "IDD", label = mdl_label, method = "COMM", value = mae_1(x[["cpue_emp"]][,id_comm,], x[["cpue_hat"]][,id_comm,])))
        bb[[i]][[basename(j)]][["pnzero_mae"]]  <- rbind(data.frame(model = "IDD", label = mdl_label, method = "SURVEY", value = mae_1(x[["pnzero_emp"]][,id_survey,], x[["pnzero_hat"]][,id_survey,])), data.frame(model = "IDD", label = mdl_label, method = "COMM", value = mae_1(x[["pnzero_emp"]][,id_comm,], x[["pnzero_hat"]][,id_comm,])))
        
        # p-values
        bb[[i]][[basename(j)]][["cpue_pvalue"]]    <- rbind(data.frame(model = "IDD", label = mdl_label, method = "SURVEY", value = pvalue_1(x[["cpue_emp"]][,id_survey,], x[["cpue_sim"]][,id_survey,])), data.frame(model = "IDD", label = mdl_label, method = "COMM", value = pvalue_1(x[["cpue_emp"]][,id_comm,], x[["cpue_sim"]][,id_comm,])))
        bb[[i]][[basename(j)]][["pnzero_pvalue"]]  <- rbind(data.frame(model = "IDD", label = mdl_label, method = "SURVEY", value = pvalue_1(x[["pnzero_emp"]][,id_survey,], x[["pnzero_sim"]][,id_survey,])), data.frame(model = "IDD", label = mdl_label, method = "COMM", value = pvalue_1(x[["pnzero_emp"]][,id_comm,], x[["pnzero_sim"]][,id_comm,])))
        bb[[i]][[basename(j)]][["density_pvalue"]] <- rbind(data.frame(model = "IDD", label = mdl_label, method = "SURVEY", value = pvalue_1(x[["density_hat"]], x[["cpua_hat"]][,id_survey,])), data.frame(model = "IDD", label = mdl_label, method = "COMM", value = pvalue_1(x[["density_hat"]], x[["cpua_hat"]][,id_comm,])))
        
        setwd(wd)
    }
}

res_summary_idd <- bb
save(res_summary_idd, file = file.path("../../results/res_summary_idd.rda"))

