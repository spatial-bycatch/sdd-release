


library(plyr)
library(dplyr)
library(kableExtra)
library(tidyr)

options(knitr.kable.NA = '--') 

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

###############
# CONVERGENCE #
###############

pal_fnc1 <- colorRamp(c("lightgreen", "lightyellow", "red"))
col_rhat <- function(x) rgb(pal_fnc1(abs(x - 1)), maxColorValue = 255)

pal_fnc2  <- colorRamp(c("lightyellow", "blue"))
col_neff <- function(x) rgb(pal_fnc2(x), maxColorValue = 255)

# get convergence diagnostics
dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["catchability_rhat"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["catchability_rhat"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd[,colnames(dfr_sdd)])

tab1 <- pivot_wider(dfr, values_from = value, names_from = method, names_prefix = "q_rhat_") %>% select(-directory)

dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["catchability_neff"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["catchability_neff"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd[,colnames(dfr_sdd)])

tab2 <- pivot_wider(dfr, values_from = value, names_from = method, names_prefix = "q_neff_") %>% select(-directory)

dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["density_rhat"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["density_rhat"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd[,colnames(dfr_sdd)])

tab3 <- dfr %>% rename(d_rhat = value) %>% select(-directory)

dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["density_neff"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["density_neff"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd[,colnames(dfr_sdd)])

tab4 <- dfr %>% rename(d_neff = value) %>% select(-directory)

tab <- full_join(tab1, tab2) %>% full_join(tab3) %>% full_join(tab4)
tab <- tab %>% select(species, model, label, contains("rhat"), contains("neff"))
tab <- tab %>% mutate(model = factor(model, levels = c("SDD", "IDD"), labels = c("Survey only", "Integrated")))

tab <- tab %>% mutate(across(contains("rhat"), \(x) vapply(x, function(y) {
    z <- ifelse(is.na(y), "--", format(y, nsmall = 3))
    cell_spec(z, background = ifelse(is.na(y), grey(0.99), col_rhat(y)), format = 'latex')
    }, character(1))))

tab <- tab %>% mutate(across(contains("neff"), \(x) vapply(x, function(y) {
    z <- ifelse(is.na(y), "--", format(y, nsmall = 3))
    cell_spec(z, background = ifelse(is.na(y), grey(0.99), col_neff(y)), format = 'latex')
}, character(1))))

tab <- tab %>% split(~species)

lapply(tab, function(x) {
  
      y <- tolower(unique(x$species))
      x <- x %>% select(-species)
      
      x %>% kable(format = 'latex',
                  col.names = c('Model', 'Label', '$\\pi_{\\,\\mathrm{SURVEY}}$', '$\\pi_{\\,\\mathrm{COMM}}$', '$\\E{d_k}$', '$\\pi_{\\,\\mathrm{SURVEY}}$', '$\\pi_{\\,\\mathrm{COMM}}$', '$\\E{d_k}$'),
                  align = c('l', 'l', 'r', 'r', 'r', 'r', 'r', 'r'),
                  label = paste0('rhat_', y),
                  caption = paste0('Convergence diagnostics per model for ', '\\species{', y, '}. Where poor convergence is indicated by $\\hat{R}$, these are shown in red. Where sampling is more efficient, indicated by a higher $N_{eff} / N$, cells are shaded a darker blue.'),
                  escape = FALSE,
                  booktabs = TRUE,
                  linesep = c('','','','','','','','\\addlinespace','','','','','','',''),
                  format.args = list(nsmall = 3)) %>%
          add_header_above(c("", "", "$\\\\hat{R}$" = 3, "$N_{eff} / N$" = 3), escape = FALSE) %>%
          kable_styling(font_size = 7) %>%
          save_kable(file = file.path("../../tables", paste0("tab_rhat_", y, ".tex")))
          
          
})

####################
# PREDICTION ERROR #
####################

# get prediction errors
pal_fnc1   <- colorRamp(c("lightgreen", "red"))
col_pvalue <- function(x) rgb(pal_fnc1(abs(x - 0.5)), maxColorValue = 255)

pal_fnc2  <- colorRamp(c("lightblue", "red"))
col_mae <- function(x) rgb(pal_fnc2(x), maxColorValue = 255)


dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["cpue_mae"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["cpue_mae"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd[,colnames(dfr_sdd)])

tab1 <- pivot_wider(dfr, values_from = value, names_from = method, names_prefix = "cpue_mae_") %>% select(-directory)

dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["pnzero_mae"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["pnzero_mae"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd[,colnames(dfr_sdd)])

tab2 <- pivot_wider(dfr, values_from = value, names_from = method, names_prefix = "pnz_mae_") %>% select(-directory)

dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["cpue_pvalue"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["cpue_pvalue"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd[,colnames(dfr_sdd)])

tab3 <- pivot_wider(dfr, values_from = value, names_from = method, names_prefix = "cpue_pvalue_") %>% select(-directory)

dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["pnzero_pvalue"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["pnzero_pvalue"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd[,colnames(dfr_sdd)])

tab4 <- pivot_wider(dfr, values_from = value, names_from = method, names_prefix = "pnz_pvalue_") %>% select(-directory)

dfr_sdd <- ldply(lapply(aa, function(x) ldply(lapply(x, function(y) y[["density_pvalue"]]), .id = "directory")), .id = "species")
dfr_idd <- ldply(lapply(bb, function(x) ldply(lapply(x, function(y) y[["density_pvalue"]]), .id = "directory")), .id = "species")
dfr     <- rbind(dfr_sdd, dfr_idd[,colnames(dfr_sdd)])

tab5 <- dfr %>% rename(d_pvalue = value) %>% select(-directory)

tab <- full_join(tab1, tab2) %>% full_join(tab3) %>% full_join(tab4)
tab <- tab %>% select(species, model, label, contains("SURVEY"), contains("COMM"))
tab <- tab %>% mutate(model = factor(model, levels = c("SDD", "IDD"), labels = c("Survey only", "Integrated")))

tab <- tab %>% mutate(across(contains("pvalue"), \(x) vapply(x, function(y) {
    z <- ifelse(is.na(y), "--", format(round(y, 2), nsmall = 2))
    cell_spec(z, background = ifelse(is.na(y), grey(0.99), ifelse(y < 0.4 | y > 0.6, "#FF0000", ifelse(y < 0.45 | y > 0.55, "#C77748", "#90EE90"))), format = 'latex')
}, character(1))))

tab <- tab %>% split(~species+model)

for (i in 1:length(tab)) {
    
    tab[[i]] <- tab[[i]] %>% mutate(across(contains("mae"), function(x) { 
        a <- max(x, na.rm = TRUE)
        b <- min(x, na.rm = TRUE)
        vapply(x, function(y) {
            z <- ifelse(is.na(y), "--", format(round(y, 2), nsmall = 2))
            cell_spec(z, background = ifelse(is.na(y), grey(0.99), ifelse(y > (1.5 * b), "#FF0000", ifelse(y > (1.25 * b), "#C597A1", "#ADD8E6"))), format = 'latex')
        }, character(1))
    }))
}

tab <- bind_rows(tab)
tab <- tab %>% split(~species)

lapply(tab, function(x) {
    
    y <- tolower(unique(x$species))
    x <- x %>% select(-species)
    
    x %>% kable(format = 'latex',
                col.names = c('Model', 'Label', rep(c('Catch rate', 'Catch prob.'), times = 4)),
                align = c('l', 'l', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r'),
                label = paste0('mae_', y),
                caption = paste0('Prediction diagnostics per model for ', '\\species{', y, '}. For models with the highest prediction error, or the least satisfactory p-value, the associated cells are shaded brighter red.'),
                escape = FALSE,
                booktabs = TRUE,
                linesep = c('','','','','','','','\\addlinespace','','','','','','',''),
                format.args = list(digits = 2, nsmall = 2)) %>%
        add_header_above(c("", "", "Survey (MAE)" = 2, "Survey (p-value)" = 2, "Comm. (MAE)" = 2, "Comm. (p-value)" = 2)) %>%
        kable_styling(font_size = 7) %>%
        save_kable(file = file.path("../../tables", paste0("tab_mae_", y, ".tex")))
    
    
})




