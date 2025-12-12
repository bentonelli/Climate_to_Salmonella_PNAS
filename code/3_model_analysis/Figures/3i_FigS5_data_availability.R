# Code to extract the amount of data available

library(dplyr)
library(ggplot2)

# Dummy data
x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)

ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()

data_mdl <- readRDS(paste("data/3_model_in/model_full_data_2024.rds"))

data_mdl$irr_obs_yrs
data_mdl$irr_missing_yrs

dt_dt <- data.frame(xx = data_mdl$mast_years)
dt_dt$yy <- 4
dt_dt$fill[data_mdl$dt_obs_ind] <- 100
dt_dt$fill[data_mdl$dt_miss_ind] <- 0

dt_mast <- data.frame(xx = data_mdl$mast_years)
dt_mast$yy <- 3
dt_mast$fill <- round(data_mdl$N_obs_mast/data_mdl$Ncls,2)*100

dt_irr <- data.frame(xx = data_mdl$mast_years)
dt_irr$yy <- 2
dt_irr$fill <- NA
dt_irr$fill[data_mdl$irr_obs_yrs] <- 100
dt_irr$fill[data_mdl$irr_missing_yrs] <- 0

dt_dis <- data.frame(xx = data_mdl$mast_years)
dt_dis$yy <- 1
dt_dis$fill <- NA
dt_dis$fill[data_mdl$disease_ts_observed] <- 100
dt_dis$fill[data_mdl$disease_ts_unobserved] <- 0

dt_all <- rbind(dt_dt,dt_mast,dt_irr,dt_dis)


dt_all$label <- dt_all$fill
dt_all$label[which(dt_all$fill==100 | dt_all$fill == 0)] <- NA

ggplot(dt_all, aes(y=yy, x=xx, fill=fill)) + 
  geom_tile(color = "black",
            lwd = .5,
            linetype = 1
            ) + theme_minimal() +
  xlab("Year") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  geom_text(aes(label = label), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "orchid4") + coord_fixed()

# Calculate the amount of data by region
mast_region <- data_mdl$yy_obs_mast

w_n_rec <- c()
e_n_rec <- c()
for (nn in 1:nrow(mast_region)){
  w_dt <- mast_region[nn,]
  w_dt <- w_dt[which(data_mdl$er_identity[1,]==1)]
  w_n <- length(which(w_dt != 0.00001))
  w_n_rec <- c(w_n_rec,w_n)
  
  e_dt <- mast_region[nn,]
  e_dt <- e_dt[which(data_mdl$er_identity[2,]==1)]
  e_n <- length(which(e_dt != 0.00001))
  e_n_rec <- c(e_n_rec,e_n)
}

