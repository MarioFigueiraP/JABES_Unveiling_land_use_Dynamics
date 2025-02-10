# =============================================================================
# Article Title: Unveiling land use dynamics: Insights from a hierarchical
#                Bayesian spatio-temporal modelling of Compositional Data
#
# Figure Code: Figure 5. In section: "Downscaling example: LAMASUS data"
# Authors: Mario Figueira Pereira
# 
# Description:
# This script provides the code to obtain the results shown in Figure 5, 
# using the downscaling model for the Compositional Data of 5 categories.
#
# Last Updated: 30/07/2024
# =============================================================================

# Required libraries ----

library(Matrix)
library(INLA)
library(fmesher)
library(compositions)
library(ggplot2)
library(gridExtra)
library(ggtext)
library(dplyr)
library(sf)
library(magick)

# Loading the data ----

NUTS3_sf_complete <- readRDS("./Data/NUTS3_sf_complete.RDS")

# Covariates selected by the stepwise search
selected_cov <- c("Grassland_rent", "Urban_rent", "Forest_rent", "gva", "gva_A_sh", "gva_A", "farm_subsidies_ha", "farm_input_ha", "gdp_pw")

## Generating integration points for dowscaling model ----
Centroids <- NUTS3_sf_complete[NUTS3_sf_complete$year==unique(NUTS3_sf_complete$year)[1],] %>% 
  st_centroid() %>% 
  st_geometry()

NUTS3_sp <- as(NUTS3_sf_complete[NUTS3_sf_complete$year==unique(NUTS3_sf_complete$year)[1],] %>% st_geometry(.), 'Spatial')

spde.puregrid <- as.matrix(expand.grid(x=seq(NUTS3_sp@bbox[1,1], NUTS3_sp@bbox[1,2],length.out=1E2),
                                       y=seq(NUTS3_sp@bbox[2,1], NUTS3_sp@bbox[2,2],length.out=1E2)))
spde.grid <- rbind(spde.puregrid, st_coordinates(Centroids))

Sp_spde.grid <- SpatialPoints(coords=spde.grid)
Sp_spde.grid@proj4string <- NUTS3_sp@proj4string
Spde.grid_over_indx <- over(x=Sp_spde.grid, y=NUTS3_sp)
spde.grid <- spde.grid[!is.na(Spde.grid_over_indx),]
Spde.grid_over_indx <- Spde.grid_over_indx[!is.na(Spde.grid_over_indx)]

## Mesh configuration ----

max_edge <- quantile(as.vector(dist(as.matrix(st_multipoint(st_coordinates(st_convex_hull(st_union(Centroids)))))[,1:2])),
                     p=c(0.03,0.1))
cutoff <- quantile(as.vector(dist(as.matrix(st_multipoint(st_coordinates(st_convex_hull(st_union(Centroids)))))[,1:2])),
                   p=c(0.01))
mesh <- fm_mesh_2d_inla(cutoff=cutoff,
                        boundary=list(fm_nonconvex_hull(x=Centroids, convex=-0.05),
                                      fm_nonconvex_hull(x=Centroids, convex=-0.2)),
                        max.edge=max_edge)

# Build the matrix for one time and then replicate it for the whole time-data inference.
Agrid <- inla.spde.make.A(mesh = mesh, loc = as.matrix(spde.grid))
Ablock_uni <- inla.spde.make.block.A(A = Agrid, block = as.vector(Spde.grid_over_indx), rescale = "count")

Ablock_inf <- Ablock_uni
for(i in 2:length(unique(NUTS3_sf_complete$year))){
  Ablock_inf <- rbind(Ablock_inf, Ablock_uni)
}

spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, prior.range = c(max_edge[1],0.1), prior.sigma=c(1,0.5), constr = TRUE)
spde.index <- inla.spde.make.index(name = "spatial", n.spde = spde$n.spde)
prior.idu <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

# Dowscaling model ----

DF_ALR <- NUTS3_sf_complete[,c("Cropland", "Forest", "Grassland", "Urban")] %>% 
  st_drop_geometry(.) %>% 
  apply(X=., MARGIN=2, FUN=function(x){log(x/NUTS3_sf_complete$OthNatLnd)}) %>%
  st_sf(.,geometry=st_geometry(NUTS3_sf_complete)) %>% 
  mutate(ID=NUTS3_sf_complete$NUTS3, .before="Cropland") %>%
  mutate(year=NUTS3_sf_complete$year, .before="geometry")

prior.idu <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

formula_Downscaling_CoDa = y ~ -1 +
  (beta10 + beta11 + beta12 + beta13 + beta14 + beta15 + beta16 + beta17 + beta18 + beta19 + 
     f(spatial1, model=spde) + 
     f(rw11, model="rw1", constr=TRUE)) + 
  f(u, model="iid", hyper=prior.idu) +
  (beta20 + beta21 + beta22 + beta23 + beta24 + beta25 + beta26 + beta27 + beta28 + beta29 +
     f(spatial2, model=spde) + 
     f(rw12, model="rw1", constr=TRUE)) + f(u, model="iid", hyper=prior.idu) +
  (beta30 + beta31 + beta32 + beta33 + beta34 + beta35 + beta36 + beta37 + beta38 + beta39 +
     f(spatial3, model=spde) + 
     f(rw13, model="rw1", constr=TRUE)) + f(u, model="iid", hyper=prior.idu) +
  (beta40 + beta41 + beta42 + beta43 + beta44 + beta45 + beta46 + beta47 + beta48 + beta49 +
     f(spatial4, model=spde) + 
     f(rw14, model="rw1", constr=TRUE)) +
  f(u, model = "iid", hyper=prior.idu)

inf_spde2_rw1_alry1 <- inla.stack(data = list(y = cbind(DF_ALR$Cropland, NA, NA, NA)),
                                  A = list(Ablock_inf, 1),
                                  effects = list(
                                    list(
                                      spatial1 = spde.index$spatial
                                    ),
                                    list(
                                      beta10 = rep(1, times=nrow(DF_ALR)),
                                      
                                      beta11 = NUTS3_sf_complete$Grassland_rent,
                                      beta12 = NUTS3_sf_complete$Urban_rent,
                                      beta13 = NUTS3_sf_complete$Forest_rent,
                                      beta14 = NUTS3_sf_complete$gva,
                                      beta15 = NUTS3_sf_complete$gva_A_sh,
                                      beta16 = NUTS3_sf_complete$gva_A,
                                      beta17 = NUTS3_sf_complete$farm_subsidies_ha,
                                      beta18 = NUTS3_sf_complete$farm_input_ha,
                                      beta19 = NUTS3_sf_complete$gdp_pw,
                                      
                                      rw11 = rep(1:length(unique(DF_ALR$year)), each=length(unique(DF_ALR$ID))),
                                      u = 1:nrow(DF_ALR)
                                    )
                                  ),
                                  tag = "inf_alry1")

inf_spde2_rw1_alry2 <- inla.stack(data = list(y = cbind(NA, DF_ALR$Forest, NA, NA)),
                                  A = list(Ablock_inf,1),
                                  effects = list(
                                    list(
                                      spatial2 = spde.index$spatial
                                    ),
                                    list(
                                      beta20 = rep(1, times=nrow(DF_ALR)),
                                      
                                      beta21 = NUTS3_sf_complete$Grassland_rent,
                                      beta22 = NUTS3_sf_complete$Urban_rent,
                                      beta23 = NUTS3_sf_complete$Forest_rent,
                                      beta24 = NUTS3_sf_complete$gva,
                                      beta25 = NUTS3_sf_complete$gva_A_sh,
                                      beta26 = NUTS3_sf_complete$gva_A,
                                      beta27 = NUTS3_sf_complete$farm_subsidies_ha,
                                      beta28 = NUTS3_sf_complete$farm_input_ha,
                                      beta29 = NUTS3_sf_complete$gdp_pw,
                                      
                                      rw12 = rep(1:length(unique(DF_ALR$year)), each=length(unique(DF_ALR$ID))),
                                      u = 1:nrow(DF_ALR)
                                    )
                                  ),
                                  tag = "inf_alry2")

inf_spde2_rw1_alry3 <- inla.stack(data = list(y = cbind(NA, NA, DF_ALR$Grassland, NA)),
                                  A = list(Ablock_inf, 1),
                                  effects = list(
                                    list(
                                      spatial3 = spde.index$spatial
                                    ),
                                    list(
                                      beta30 = rep(1, times=nrow(DF_ALR)),
                                      
                                      beta31 = NUTS3_sf_complete$Grassland_rent,
                                      beta32 = NUTS3_sf_complete$Urban_rent,
                                      beta33 = NUTS3_sf_complete$Forest_rent,
                                      beta34 = NUTS3_sf_complete$gva,
                                      beta35 = NUTS3_sf_complete$gva_A_sh,
                                      beta36 = NUTS3_sf_complete$gva_A,
                                      beta37 = NUTS3_sf_complete$farm_subsidies_ha,
                                      beta38 = NUTS3_sf_complete$farm_input_ha,
                                      beta39 = NUTS3_sf_complete$gdp_pw,
                                      
                                      rw13 = rep(1:length(unique(DF_ALR$year)), each=length(unique(DF_ALR$ID))),
                                      u = 1:nrow(DF_ALR)
                                    )
                                  ),
                                  tag="inf_alry3")

inf_spde2_rw1_alry4 <- inla.stack(data = list(y = cbind(NA, NA, NA, DF_ALR$Urban)),
                                  A = list(Ablock_inf, 1),
                                  effects = list(
                                    list(
                                      spatial4 = spde.index$spatial
                                    ),
                                    list(
                                      beta40 = rep(1, times=nrow(DF_ALR)),
                                      
                                      beta41 = NUTS3_sf_complete$Grassland_rent,
                                      beta42 = NUTS3_sf_complete$Urban_rent,
                                      beta43 = NUTS3_sf_complete$Forest_rent,
                                      beta44 = NUTS3_sf_complete$gva,
                                      beta45 = NUTS3_sf_complete$gva_A_sh,
                                      beta46 = NUTS3_sf_complete$gva_A,
                                      beta47 = NUTS3_sf_complete$farm_subsidies_ha,
                                      beta48 = NUTS3_sf_complete$farm_input_ha,
                                      beta49 = NUTS3_sf_complete$gdp_pw,
                                      
                                      rw14 = rep(1:length(unique(DF_ALR$year)), each=length(unique(DF_ALR$ID))),
                                      u = 1:nrow(DF_ALR)
                                    )
                                  ),
                                  tag="inf_alry4")

inf_spde2_rw1 <- inla.stack(inf_spde2_rw1_alry1, inf_spde2_rw1_alry2, inf_spde2_rw1_alry3, inf_spde2_rw1_alry4)

inla_Downscaling_CoDa <- inla(formula = formula_Downscaling_CoDa, data = inla.stack.data(inf_spde2_rw1), family = rep("gaussian", 4),
                              control.predictor = list(A = inla.stack.A(inf_spde2_rw1), compute = TRUE),
                              control.compute = list(config = FALSE, dic = FALSE, waic = FALSE, cpo = FALSE),
                              verbose = FALSE)

Sp_spde.puregrid <- SpatialPoints(coords=spde.puregrid)
Sp_spde.puregrid@proj4string <- NUTS3_sp@proj4string
Spde.puregrid_over_indx <- over(x=Sp_spde.puregrid, y=NUTS3_sp)
spde.puregrid <- spde.puregrid[!is.na(Spde.puregrid_over_indx),]

A_pred <- inla.spde.make.A(mesh = mesh, loc = spde.puregrid)

# "Median", "2.5th Percentile", "97.5th Percentile"

sp1_median <- 
  st_as_sf(data.frame(ID = "Median", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial1$`0.5quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))
sp1_2.5 <- 
  st_as_sf(data.frame(ID = "2.5th Percentile", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial1$`0.025quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))

sp1_97.5 <- 
  st_as_sf(data.frame(ID = "97.5th Percentile", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial1$`0.975quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))

gg_sp1 <- ggplot() + 
  geom_sf(data = sp1_median, mapping = aes(color = z)) + 
  geom_sf(data = sp1_2.5, mapping = aes(color = z)) + 
  geom_sf(data = sp1_97.5, mapping = aes(color = z)) +
  geom_sf(data = NUTS3_sf_complete[NUTS3_sf_complete$NUTS3 %>% unique(.) %>% seq_along(.),], mapping = aes(), fill = NA, color = alpha("grey20",0.2)) +
  scale_color_viridis_c(option = "turbo") + labs(color = "Values", title = "Spatial Effect (ALR 1)", x = "", y = "") + 
  theme_bw() + facet_wrap(~ factor(ID, levels = c("Median", "2.5th Percentile", "97.5th Percentile")), nrow = 1) +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

sp2_median <- 
  st_as_sf(data.frame(ID = "Median", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial2$`0.5quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))
sp2_2.5 <- 
  st_as_sf(data.frame(ID = "2.5th Percentile", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial2$`0.025quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))

sp2_97.5 <- 
  st_as_sf(data.frame(ID = "97.5th Percentile", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial2$`0.975quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))

gg_sp2 <- ggplot() + 
  geom_sf(data = sp2_median, mapping = aes(color = z)) + 
  geom_sf(data = sp2_2.5, mapping = aes(color = z)) + 
  geom_sf(data = sp2_97.5, mapping = aes(color = z)) +
  geom_sf(data = NUTS3_sf_complete[NUTS3_sf_complete$NUTS3 %>% unique(.) %>% seq_along(.),], fill = NA, color = alpha("grey20",0.2)) +
  scale_color_viridis_c(option = "turbo") + labs(color = "Values", title = "Spatial Effect (ALR 2)", x = "", y = "") + 
  theme_bw() + facet_wrap(~ factor(ID, levels = c("Median", "2.5th Percentile", "97.5th Percentile")), nrow = 1) +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

sp3_median <- 
  st_as_sf(data.frame(ID = "Median", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial3$`0.5quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))
sp3_2.5 <- 
  st_as_sf(data.frame(ID = "2.5th Percentile", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial3$`0.025quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))

sp3_97.5 <- 
  st_as_sf(data.frame(ID = "97.5th Percentile", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial3$`0.975quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))

gg_sp3 <- ggplot() + 
  geom_sf(data = sp3_median, mapping = aes(color = z)) + 
  geom_sf(data = sp3_2.5, mapping = aes(color = z)) + 
  geom_sf(data = sp3_97.5, mapping = aes(color = z)) +
  geom_sf(data = NUTS3_sf_complete[NUTS3_sf_complete$NUTS3 %>% unique(.) %>% seq_along(.),], fill = NA, color = alpha("grey20",0.2)) +
  scale_color_viridis_c(option = "turbo") + labs(color = "Values", title = "Spatial Effect (ALR 3)", x = "", y = "") + 
  theme_bw() + facet_wrap(~ factor(ID, levels = c("Median", "2.5th Percentile", "97.5th Percentile")), nrow = 1) +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

sp4_median <- 
  st_as_sf(data.frame(ID = "Median", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial4$`0.5quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))
sp4_2.5 <- 
  st_as_sf(data.frame(ID = "2.5th Percentile", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial4$`0.025quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))

sp4_97.5 <- 
  st_as_sf(data.frame(ID = "97.5th Percentile", spde.puregrid, z = drop(fm_basis(x = mesh, loc = spde.puregrid) %*% inla_Downscaling_CoDa$summary.random$spatial4$`0.975quant`)),
           coords = c("x", "y"), crs = st_crs(NUTS3_sf_complete))

gg_sp4 <- ggplot() + 
  geom_sf(data = sp4_median, mapping = aes(color = z)) + 
  geom_sf(data = sp4_2.5, mapping = aes(color = z)) + 
  geom_sf(data = sp4_97.5, mapping = aes(color = z)) +
  geom_sf(data = NUTS3_sf_complete[NUTS3_sf_complete$NUTS3 %>% unique(.) %>% seq_along(.),], fill = NA, color = alpha("grey20",0.2)) +
  scale_color_viridis_c(option = "turbo") + labs(color = "Values", title = "Spatial Effect (ALR 4)", x = "", y = "") + 
  theme_bw() + facet_wrap(~ factor(ID, levels = c("Median", "2.5th Percentile", "97.5th Percentile")), nrow = 1) + 
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

gg_t1 <- ggplot() +
  geom_ribbon(data = data.frame(inla_Downscaling_CoDa$summary.random$rw11, id = "Temporal Effect") %>%
                rename(., all_of(c(q1 = 'X0.025quant', q3 = 'X0.975quant'))),
              mapping = aes(x = ID + 2006, ymin = q1, ymax = q3), alpha = 0.4) +
  geom_line(data = data.frame(inla_Downscaling_CoDa$summary.random$rw11, id = "Temporal Effect"),
            mapping = aes(x = ID + + 2006, y = mean)) +
  theme_bw() +
  facet_wrap(facets = ~ id, scales = "free") +
  xlab(label = "Year") + ylab(label = "Value") + labs(title = "Temporal Effect (ALR 1)") + 
  scale_x_continuous(breaks = seq(2007, 2018, by = 1)) +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

gg_t2 <- ggplot() +
  geom_ribbon(data = data.frame(inla_Downscaling_CoDa$summary.random$rw12, id = "Temporal Effect") %>%
                rename(., all_of(c(q1 = 'X0.025quant', q3 = 'X0.975quant'))),
              mapping = aes(x = ID + 2006, ymin = q1, ymax = q3), alpha = 0.4) +
  geom_line(data = data.frame(inla_Downscaling_CoDa$summary.random$rw12, id = "Temporal Effect"),
            mapping = aes(x = ID + + 2006, y = mean)) +
  theme_bw() +
  facet_wrap(facets = ~ id, scales = "free") +
  xlab(label = "Year") + ylab(label = "Value") + labs(title = "Temporal Effect (ALR 2)") + 
  scale_x_continuous(breaks = seq(2007, 2018, by = 1)) +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

gg_t3 <- ggplot() +
  geom_ribbon(data = data.frame(inla_Downscaling_CoDa$summary.random$rw13, id = "Temporal Effect") %>%
                rename(., all_of(c(q1 = 'X0.025quant', q3 = 'X0.975quant'))),
              mapping = aes(x = ID + 2006, ymin = q1, ymax = q3), alpha = 0.4) +
  geom_line(data = data.frame(inla_Downscaling_CoDa$summary.random$rw13, id = "Temporal Effect"),
            mapping = aes(x = ID + + 2006, y = mean)) +
  theme_bw() +
  facet_wrap(facets = ~ id, scales = "free") +
  xlab(label = "Year") + ylab(label = "Value") + labs(title = "Temporal Effect (ALR 3)") + 
  scale_x_continuous(breaks = seq(2007, 2018, by = 1)) +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

gg_t4 <- ggplot() +
  geom_ribbon(data = data.frame(inla_Downscaling_CoDa$summary.random$rw14, id = "Temporal Effect") %>%
                rename(., all_of(c(q1 = 'X0.025quant', q3 = 'X0.975quant'))),
              mapping = aes(x = ID + 2006, ymin = q1, ymax = q3), alpha = 0.4) +
  geom_line(data = data.frame(inla_Downscaling_CoDa$summary.random$rw14, id = "Temporal Effect"),
            mapping = aes(x = ID + + 2006, y = mean)) +
  theme_bw() +
  facet_wrap(facets = ~ id, scales = "free") +
  xlab(label = "Year") + ylab(label = "Value") + labs(title = "Temporal Effect (ALR 4)") + 
  scale_x_continuous(breaks = seq(2007, 2018, by = 1)) +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

# Option 1 

gg_sp <- grid.arrange(arrangeGrob(grobs = list(gg_sp1, gg_sp2, gg_sp3, gg_sp4), layout_matrix = matrix(c(1,2,3,4), ncol = 2, byrow = TRUE)))
ggsave(file = "./gg_sp.png", device = "png", plot = gg_sp, width = 48, height = 24, units = "cm", dpi = 300)

gg_temp <- grid.arrange(arrangeGrob(grobs = list(gg_t1, gg_t2, gg_t3, gg_t4), layout_matrix = matrix(c(1,2,3,4), ncol = 2, byrow = TRUE)))
ggsave(file = "./gg_temp.png", device = "png", plot = gg_temp, width = 48, height = 16, units = "cm", dpi = 300)

# Option 2

gg_sp <- grid.arrange(arrangeGrob(grobs = list(gg_sp1, gg_sp2, gg_sp3, gg_sp4), layout_matrix = matrix(c(1,2,3,4), ncol = 1, byrow = TRUE)))
ggsave(file = "./gg_sp.png", device = "png", plot = gg_sp, width = 36, height = 36, units = "cm", dpi = 300)

gg_temp <- grid.arrange(arrangeGrob(grobs = list(gg_t1, gg_t2, gg_t3, gg_t4), layout_matrix = matrix(c(1,2,3,4), ncol = 2, byrow = TRUE)))
ggsave(file = "./gg_temp.png", device = "png", plot = gg_temp, width = 36, height = 12, units = "cm", dpi = 300)

image1 <- image_read("./Figures/gg_sp.png")
image2 <- image_read("./Figures/gg_temp.png")
image_write(image_append(c(image1, image2), stack = TRUE), "./Figures/fig_CoDaExample.png")
