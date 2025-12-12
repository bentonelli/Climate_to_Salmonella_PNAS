#Script to plot regional masting estimates versus irruptions

library(MCMCvis)
library(ggplot2)
library(dplyr)
library(readr)


mdl_data <- readRDS("data/3_model_in/model_full_data_2024.rds")
mdl_in <- readRDS("data/4_model_out/mdl_fit.rds")

spec_names <- c("Global","Pine siskin","Red-breasted nuthatch","Evening grosbeak",
                "Pine grosbeak","Purple finch","Redpoll","Red crossbill","White-winged crossbill")

par(mfrow=c(2,1))
MCMCplot(mdl_in,params=c("mu_omega1","omega1"),col="#648DE5",
         labels = spec_names,xlim=c(-1.75,1),ci = c(50,89),
         xlab="Effect of Cone Production on Irruption Intensity",ref_ovl = TRUE,sz_ax_txt=1.25)
MCMCplot(mdl_in,params=c("mu_eta1","eta1"),col="#648DE5",
         labels = spec_names,xlim=c(-1.75,1),ci = c(50,89),
         xlab="Lag Effect on Irruption Intensity",ref_ovl = TRUE,sz_ax_txt=1.25)


MCMCplot(mdl_in,params=c("mu_omega2","omega2"),col="#E98A15",
         labels = spec_names,xlim=c(-1.75,1),ci = c(50,89),
         xlab="Effect of Cone Production on Irruption Intensity",ref_ovl = TRUE,sz_ax_txt=1.25)

MCMCplot(mdl_in,params=c("mu_eta2","eta2"),col="#E98A15",
         labels = spec_names,xlim=c(-1.75,1),ci = c(50,89),
         xlab="Lag Effect on Irruption Intensity",ref_ovl = TRUE,sz_ax_txt=1.25)

