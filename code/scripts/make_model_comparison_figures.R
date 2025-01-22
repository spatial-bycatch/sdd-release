

library(ggplot2)
library(plyr)
library(dplyr)
library(sdd)

wd <- getwd()

spps <- c("SSI", "BBE", "ETB", "SND")

#############
# SDD MODEL #
#############

load(file = file.path("../../results/res_summary_sdd.rda"))

#############
# IDD MODEL #
#############

load(file = file.path("../../results/res_summary_idd.rda"))

# ASSIGN

aa <- res_summary_sdd
bb <- res_summary_idd

########
# FITS #
########
dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["catch_rate"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["catch_rate"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd)

dfr$method <- factor(dfr$method, levels = c("SURVEY", "COMM"), labels = c("Survey", "Comm."))
dfr$model  <- factor(dfr$model,  levels = c("SDD", "IDD"), labels = c("Survey only model", "Integrated model"))

for (spp in spps) {
    
    gg1 <- ggplot(filter(dfr, species == spp)) + 
        geom_point(aes(x = value_emp, y = value_hat), alpha = 0.4) +
        facet_grid(method+label~model) + 
        theme_bw(base_size = 12) + scale_y_log10() + scale_x_log10() +
        geom_abline(intercept = 0, slope = 1, col = 'red') +
        labs(x = "Observed", y = "Predicted")
    #ggsave(gg1, file = file.path("../../figures", paste0("fit_catch_rate_per_model_", tolower(spp), ".png")), height = 14)
    ggsave(gg1, file = file.path("../../figures", paste0("fit_catch_rate_per_model_", tolower(spp), ".pdf")), height = 16)
    
}

dfr_final <- dfr %>% filter(label == "E4") %>% filter(species %in% c("SSI", "BBE", "ETB", "SND"))

gg <- ggplot(dfr_final) + 
    geom_point(aes(x = value_emp, y = value_hat), alpha = 0.4) +
    scale_y_log10() + scale_x_log10() +
    facet_grid(species~method+model) + 
    theme_bw(base_size = 12) +
    geom_abline(intercept = 0, slope = 1, col = 'tomato') +
    labs(x = "Observed", y = "Predicted")
#ggsave(gg1, file = file.path("../../figures", paste0("fit_catch_prob_per_model_", tolower(spp), ".png")), height = 14)
ggsave(gg, file = file.path("../../figures", "fit_catch_rate_per_model_final.pdf"), height = 7)


dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["catch_prob"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["catch_prob"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd)

dfr$method <- factor(dfr$method, levels = c("SURVEY", "COMM"), labels = c("Survey", "Comm."))
dfr$model  <- factor(dfr$model,  levels = c("SDD", "IDD"), labels = c("Survey only model", "Integrated model"))

for (spp in spps) {
    
    gg1 <- ggplot(filter(dfr, species == spp)) + 
        geom_point(aes(x = value_emp, y = value_sim), alpha = 0.4) +
        facet_grid(method+label~model) + 
        theme_bw(base_size = 12) +
		scale_x_continuous(breaks = c(0.25, 0.75)) +
		scale_y_continuous(breaks = c(0.25, 0.75)) +
        geom_abline(intercept = 0, slope = 1, col = 'red') +
        labs(x = "Observed", y = "Predicted")
    #ggsave(gg1, file = file.path("../../figures", paste0("fit_catch_prob_per_model_", tolower(spp), ".png")), height = 14)
    ggsave(gg1, file = file.path("../../figures", paste0("fit_catch_prob_per_model_", tolower(spp), ".pdf")), height = 16)
    
}

dfr_final <- dfr %>% filter(label == "E4") %>% filter(species %in% c("SSI", "BBE", "ETB", "SND"))

gg <- ggplot(dfr_final) + 
    geom_point(aes(x = value_emp, y = value_sim), alpha = 0.4) +
    facet_grid(species~method+model) + 
    theme_bw(base_size = 12) +
    scale_x_continuous(breaks = c(0.25, 0.75)) +
    scale_y_continuous(breaks = c(0.25, 0.75)) +
    geom_abline(intercept = 0, slope = 1, col = 'tomato') +
    labs(x = "Observed", y = "Predicted")
#ggsave(gg1, file = file.path("../../figures", paste0("fit_catch_prob_per_model_", tolower(spp), ".png")), height = 14)
ggsave(gg, file = file.path("../../figures", "fit_catch_prob_per_model_final.pdf"), height = 7)


##############
# POSTERIORS #
##############

dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["catchability_mcmc"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["catchability_mcmc"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd) %>% mutate(value_log = log(value))

dfr$method <- factor(dfr$method, levels = c("SURVEY", "COMM"), labels = c("Survey", "Comm."))
dfr$model  <- factor(dfr$model,  levels = c("SDD", "IDD"), labels = c("Survey only model", "Integrated model"))
dfr$iterations <- as.numeric(as.character(dfr$iterations))

dfr_prior <- expand.grid(value_log = log(rbeta(1e5, 1, 1)), method = c("Survey", "Comm."), model = c("Survey only model", "Integrated model"))
dfr_prior$value_log[dfr_prior$method == "Comm." & dfr_prior$model == "Survey only model"] <- NA_real_

for (spp in spps) {
    
    gg1 <- ggplot(filter(dfr, species == spp)) + 
        geom_density(aes(x = value_log, fill = chains), alpha = 0.4) +
        geom_density(data = dfr_prior, aes(x = value_log), col = 'red') +
        facet_grid(method+label~model) + 
        guides(fill = "none") + 
        theme_bw(base_size = 12) + theme(axis.text.y = element_blank()) + 
        scale_x_continuous(limits = c(-5,0)) +
        labs(x = expression(paste("Catchability (", log(pi), ")")), y = "Posterior density")
    #ggsave(gg1, file = file.path("../../figures", paste0("posterior_pi_per_model_", tolower(spp), ".png")), height = 14)
    ggsave(gg1, file = file.path("../../figures", paste0("posterior_pi_per_model_", tolower(spp), ".pdf")), height = 16)
    
    gg2 <- ggplot(filter(dfr, species == spp)) + 
        geom_line(aes(x = iterations, y = value_log, col = chains), alpha = 0.8) + 
        facet_grid(method+label~model) + 
        guides(col = "none") + theme_bw(base_size = 12) +
        labs(y = expression(paste("Catchability (", log(pi), ")")), x = "Posterior sample")
    #ggsave(gg2, file = file.path("../../figures", paste0("trace_pi_per_model_", tolower(spp), ".png")), height = 14)
    ggsave(gg2, file = file.path("../../figures", paste0("trace_pi_per_model_", tolower(spp), ".pdf")), height = 16)
}

dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["density_mcmc"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["density_mcmc"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd) %>% mutate(value_log = log(value))

dfr$model  <- factor(dfr$model,  levels = c("SDD", "IDD"), labels = c("Survey only model", "Integrated model"))
dfr$iterations <- as.numeric(as.character(dfr$iterations))

for (spp in spps) {
    
    gg1 <- ggplot(filter(dfr, species == spp)) + 
        geom_density(aes(x = value_log, fill = chains), alpha = 0.4) +
        facet_grid(label~model) + 
        guides(fill = "none") + 
        theme_bw(base_size = 12) + theme(axis.text.y = element_blank()) + 
        labs(x = expression(paste("Biomass density (", log(E(d)), ")")), y = "Posterior density")
    #ggsave(gg1, file = file.path("../../figures", paste0("posterior_d_per_model_", tolower(spp), ".png")), height = 7)
    ggsave(gg1, file = file.path("../../figures", paste0("posterior_d_per_model_", tolower(spp), ".pdf")), height = 7)
    
    gg2 <- ggplot(filter(dfr, species == spp)) + 
        geom_line(aes(x = iterations, y = value_log, col = chains), alpha = 0.8) + 
        facet_grid(label~model) + 
        guides(col = "none") + theme_bw(base_size = 12) +
        labs(y = expression(paste("Biomass density (", log(E(d)), ")")), x = "Posterior sample")
    #ggsave(gg2, file = file.path("../../figures", paste0("trace_d_per_model_", tolower(spp), ".png")), height = 7)
    ggsave(gg2, file = file.path("../../figures", paste0("trace_d_per_model_", tolower(spp), ".pdf")), height = 7)
}

################
# DISTRIBUTION #
################

ff <- function(x) {
    
    suppressMessages({
        y <- x %>% group_by(model, label, cell) %>% mutate(density_bin = mean(density_hat)) %>% ungroup()
        y <- y %>% 
            group_by(model, label, method, density_bin) %>% 
            summarise(cpue_emp_upp = quantile(cpue_emp, 0.95), cpue_emp_low = quantile(cpue_emp, 0.05), cpue_emp_hat = median(cpue_emp), cpue_emp = mean(cpue_emp), cpue_sim_upp = quantile(cpue_sim, 0.95), cpue_sim_low = quantile(cpue_sim, 0.05), cpue_sim_hat = median(cpue_sim), cpue_sim = mean(cpue_sim))
    })
    return(y)    
}

dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) ff(y[["distribution"]])), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) ff(y[["distribution"]])), .id = "directory")), .id = "species")

dfr <- rbind(dfr_sdd, dfr_idd)

dfr$method <- factor(dfr$method, levels = c("SURVEY", "COMM"), labels = c("Survey", "Comm."))
dfr$model  <- factor(dfr$model,  levels = c("SDD", "IDD"), labels = c("Survey only model", "Integrated model"))

dfr <- dfr %>% group_by(species, model, label, method) %>% mutate(density_order = order(density_bin)) %>% ungroup()
dfr <- dfr %>% group_by(species, model, label, method) %>% mutate(cpue_emp_cs = cumsum(cpue_emp[density_order]) / sum(cpue_emp), cpue_sim_cs = cumsum(cpue_sim_hat[density_order]) / sum(cpue_sim_hat), cpue_sim_cs_upp = cumsum(cpue_sim_upp[density_order]) / sum(cpue_sim_upp), cpue_sim_cs_low = cumsum(cpue_sim_low[density_order]) / sum(cpue_sim_low)) %>% ungroup()

for (spp in spps) {
    
    gg <- ggplot(filter(dfr, species == spp)) + 
        #geom_pointrange(aes(x = density_bin, y = cpue_sim_hat, ymin = cpue_sim_low, ymax = cpue_sim_upp), col = 'tomato', alpha = 0.4) +
        geom_point(aes(x = density_bin, y = cpue_sim_hat), col = 'red', shape = 1) + 
        geom_point(aes(x = density_bin, y = cpue_emp), col = 'darkblue') + 
        scale_y_log10() + scale_x_log10() +
        facet_grid(method+label~model) + theme_bw() +
        labs(x = "Posterior density per cell", y = "Observed and predicted catch rates")
    ggsave(gg, file = file.path("../../figures", paste0("fit_catch_rate_density_per_model_", tolower(spp), ".pdf")), height = 16)
    
    gg <- ggplot(filter(dfr, species == spp)) + 
        geom_line(aes(x = density_order, y = cpue_sim_cs_upp), col = 'tomato') + 
        geom_line(aes(x = density_order, y = cpue_sim_cs_low), col = 'tomato') + 
        geom_line(aes(x = density_order, y = cpue_sim_cs), col = 'red', linewdith = 1.5) + 
        geom_line(aes(x = density_order, y = cpue_emp_cs), col = 'darkblue', linewdith = 1.5) +
        facet_grid(method+label~model) + theme_bw() +
        labs(x = "Ordered grid cell", y = "Cumulative sum of the catch rates")
    ggsave(gg, file = file.path("../../figures", paste0("fit_catch_rate_distribution_per_model_", tolower(spp), ".pdf")), height = 16)
    
}

# RESULTS FIGURE 
dfr_final <- dfr %>% filter(label == "E4") %>% filter(species %in% c("SSI", "BBE", "ETB", "SND"))

gg <- ggplot(filter(dfr_final, cpue_emp > 0 & cpue_sim_hat > 0)) + 
    geom_point(aes(x = density_bin, y = cpue_emp), col = 'darkblue') + 
    geom_point(aes(x = density_bin, y = cpue_sim_hat), col = 'red', shape = 1) + 
    scale_y_log10() + scale_x_log10() +
    facet_grid(species~method+model) + theme_bw() +
    labs(x = "Posterior density per cell", y = "Observed and predicted catch rates")
ggsave(gg, file = file.path("../../figures", "fit_catch_rate_density_per_model_final.pdf"), height = 7)

gg <- ggplot(dfr_final) + 
    geom_line(aes(x = density_order, y = cpue_sim_cs_upp), col = 'tomato') + 
    geom_line(aes(x = density_order, y = cpue_sim_cs_low), col = 'tomato') + 
    geom_line(aes(x = density_order, y = cpue_sim_cs), col = 'red', linewidth = 1) + 
    geom_line(aes(x = density_order, y = cpue_emp_cs), col = 'darkblue', linewidth = 1) +
    scale_y_continuous(breaks = c(0.25, 0.75)) +
    facet_grid(species~method+model) + theme_bw() +
    labs(x = "Ordered grid cell", y = "Cumulative sum of the catch rates")
ggsave(gg, file = file.path("../../figures", "fit_catch_rate_distribution_per_model_final.pdf"), height = 7)

