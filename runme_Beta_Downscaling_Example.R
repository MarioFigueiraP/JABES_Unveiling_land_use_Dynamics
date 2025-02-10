# =============================================================================
# Article Title: Unveiling land use dynamics: Insights from a hierarchical
#                Bayesian spatio-temporal modelling of Compositional Data
#
# 
# Description:
# This script provides the code to perform dowscaling models for one category in land use data.
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
  beta10 + beta11 + beta12 + beta13 + beta14 + beta15 + beta16 + beta17 + beta18 + beta19 + 
     f(spatial1, model=spde) + f(rw11, model="rw1", constr=TRUE)

inf_stk <- inla.stack(data = list(y = cbind(NUTS3_sf_complete$Cropland)),
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
                            tag = "inf")

inla_Beta_CoDa <- inla(formula = formula_Downscaling_CoDa, data = inla.stack.data(inf_stk), family = "beta",
                              control.predictor = list(A = inla.stack.A(inf_stk), compute = TRUE),
                              control.compute = list(config = FALSE, dic = FALSE, waic = FALSE, cpo = FALSE),
                              verbose = FALSE)
