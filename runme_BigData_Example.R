# =============================================================================
# Article Title: Unveiling land use dynamics: Insights from a hierarchical
#                Bayesian spatio-temporal modelling of Compositional Data
#
# Figure Code: Figure 6. In section: "Big Data example"
# Authors: Mario Figueira Pereira
# 
# Description:
# Modified script. Properly commented step by step.
# This script provides the code to obtain the results shown in Figure 5,
# using the downscaling model for the Compositional Data of 5 categories.
#
# Last Updated: 27/09/2025
# =============================================================================

remove(list = ls())

# Example Sequential Consensus splitting data.frame ----

library(INLA)
library(inlabru)
library(fmesher)

library(Matrix)
library(compositions)

library(ggplot2)
library(gridExtra)
library(ggnewscale)
library(ggtext)
library(latex2exp)

library(parallel)
library(doParallel) 

library(dplyr)
library(sf)
library(giscoR)

## Custom functions ----

# Function to split a sf_polygon using a Voronoi Diagram

split_poly <- function(seed, sf_poly, n_areas) {
  #' Split polygon
  #' 
  #' @sf_poly Geometry with the polygon to be splitted
  #' @param n_areas The number of resulting areas
  
  if(!missing(seed)){set.seed(seed)}
  
  # Create random points
  points_rnd <- st_sample(sf_poly, size = 1E4)
  # k-means clustering
  points <- do.call(rbind, st_geometry(points_rnd))
  k_means <- kmeans(points, centers = n_areas)
  # Create voronoi polygons
  voronoi_polys <- st_cast(x = st_multipolygon(x = lapply(st_voronoi(x = st_multipoint(x = k_means$centers)), FUN = function(x){x})) %>% st_geometry(obj = .), to = "POLYGON")
  
  # Intersect to set the right boundaries and compute the area for each polygon
  st_crs(voronoi_polys) <- st_crs(sf_poly)
  equal_areas <- st_intersection(voronoi_polys, sf_poly) # sf_poly$geometry
  equal_areas <- st_sf(data.frame(area = st_area(equal_areas)), geometry = equal_areas)
  return(equal_areas)
}

# Function to simulate multivariate gaussian realizations, given the linear predictor 

simulation_MVN <- function(Q, mu){
  # Q: the precision matrix of the Multivariate Gaussian distribution (MVN)
  # mu: matrix of the MVN mean
  R = Matrix::chol(Q)
  
  x <- lapply(X = 1:nrow(mu), FUN = function(i){
    w = rnorm(nrow(Q))
    z = backsolve(R, w)
    z = z + mu[i,]
    return(z)
  }) %>% do.call(what = rbind, .)
  colnames(x) <- paste0("y", 1:nrow(Q))
  return(x)
}

# Function to simulate a matern GP thorugh the SPDE approach

simulation_SPDE <- function(mesh, sigma = 1, range = 1, n_rep = 1, n_cores = 1, seed, big = FALSE, constr = FALSE){
  # mesh: mesh to create the precision matrix of the SPDE
  # sigma: the marginal precision of the CAR prior (the variability between first-order neighbors)
  # range: the spatial range; it means the range at which the correlation falls to 0.1 (in the correl)
  # constr: an argument that, if TRUE, forces the GMRF to have a global mean of zero
  # seed: seed to reproduce the result
  if(!missing(seed)){set.seed(seed)}
  Q = fmesher::fm_matern_precision(x = mesh, alpha = 2, rho = range, sigma = sigma)
  C0 = diag(fmesher::fm_fem(mesh = mesh, Q = Q)$c0)
  if(big){
    # Reordering the Q matrix (or the equivalent Graph) to reduce the fill-in of the Cholesky factor
    idx.irQ <- INLA::inla.qreordering(graph = Q)
    L <- Matrix::chol(Q[idx.irQ$ireordering, idx.irQ$ireordering]) 
    C0 <- C0[idx.irQ$ireordering]
    
    cl <- parallel::makeCluster(spec = n_cores, type='PSOCK')
    clusterEvalQ(cl, { library(Matrix) })
    clusterExport(cl = cl, c("Q", "C0", "L", "idx.irQ", "n_rep", "constr"), envir = environment())
    
    u_sp <- parLapply(cl = cl, X = seq_len(n_rep), fun = function(i){      
      w = rnorm(nrow(Q))
      u_sp = backsolve(L, w)
      if(constr){
        u_sp = u_sp - mean(C0 * u_sp) / C0
      }
      return(u_sp[idx.irQ$reordering])
    }) %>% do.call(what = cbind, .)
    stopCluster(cl)
  } else{
    L = Matrix::chol(Q)
    
    cl <- parallel::makeCluster(spec = n_cores, type='PSOCK')
    clusterEvalQ(cl, { library(Matrix) })
    clusterExport(cl = cl, c("Q", "C0", "L", "n_rep", "constr"), envir = environment())
    
    u_sp <- parLapply(cl = cl, X = seq_len(n_rep), fun = function(i){      
      w = rnorm(nrow(Q))
      u_sp = backsolve(L, w)
      if(constr){
        u_sp = u_sp - mean(C0 * u_sp) / C0
      }
      return(u_sp)
    }) %>% do.call(what = cbind, .)
    stopCluster(cl)
  }
  return(list(u=u_sp, Q=Q))
}

## Simulating data ----

seed <- 1234
set.seed(seed = seed)

### Preparing the Spatial structure for aggregation ----

ES_Polygons <- gisco_get_countries(resolution = "01", country = "ES")$geometry %>% st_cast(x = ., to = "POLYGON")
ES_mainland_boundary <- ES_Polygons[ES_Polygons %>% st_area(.) %>% which.max(.)]

ggplot() + geom_sf(data = ES_mainland_boundary, mapping = aes())

ES_voronoi_1 <- split_poly(sf_poly = ES_mainland_boundary, n_areas = 200)
ES_voronoi_2 <- split_poly(sf_poly = ES_mainland_boundary, n_areas = 250)
ES_voronoi_3 <- split_poly(sf_poly = ES_mainland_boundary, n_areas = 300)
ggplot() + geom_sf(data = ES_voronoi_1, mapping = aes())
ggplot() + geom_sf(data = ES_voronoi_2, mapping = aes())
ggplot() + geom_sf(data = ES_voronoi_3, mapping = aes())

### Simulation of the linear predictor of the 2-MVN ----

coord_ES_mainland_boundary <- (ES_mainland_boundary %>% st_coordinates(x = .))[(ES_mainland_boundary %>% st_coordinates(x = .))[,"L1"]==1,1:2]
simplify_ES_mainland_boundary <- coord_ES_mainland_boundary[inla.simplify.curve(loc = coord_ES_mainland_boundary, idx = 1:nrow(coord_ES_mainland_boundary), eps = 0.01),]

# Internal boundary (mesh)
coord_bnd_int <- st_coordinates(fm_nonconvex_hull(x = simplify_ES_mainland_boundary, convex = -0.03))
boundary_int <- coord_bnd_int[coord_bnd_int[,"L1"]==1,1:2] %>% st_polygon(x = list(.))

# External boundary (mesh)
coord_bnd_ext <- st_coordinates(fm_nonconvex_hull(x = simplify_ES_mainland_boundary, convex = -0.2))
boundary_ext <- coord_bnd_ext[coord_bnd_ext[,"L1"]==1,1:2] %>% st_polygon(x = list(.))

mesh <- fm_mesh_2d_inla(boundary = list(boundary_int, boundary_ext), max.edge = c(0.5,1))
ggplot() + gg(mesh) + geom_sf(data = ES_mainland_boundary, mapping = aes(), alpha = 0.5) + theme_bw()

ggplot() + 
  gg(mesh) + 
  geom_sf(data = ES_voronoi_1, mapping = aes(), color = "red", fill = NA, alpha = 0.5) +
  geom_sf(data = ES_voronoi_2, mapping = aes(), color = "blue", fill = NA, alpha = 0.5) +
  geom_sf(data = ES_voronoi_3, mapping = aes(), color = "forestgreen", fill = NA, alpha = 0.5) +
  theme_bw()

### Simulation first linear predictor ----

nt <- 60 # Change to 600 for the script used for the article results
grid.integration <- st_sample(x = ES_mainland_boundary, size = 1E4, type = "regular")

A_inf_spt <- inla.spde.make.A(mesh = mesh, loc = grid.integration %>% rep(., nt) %>% st_coordinates(.), group = rep(1:nt, each = length(grid.integration))) 
n_sim <- length(grid.integration)*nt

beta01 <- 2
sp_sim <- simulation_SPDE(mesh = mesh, n_rep = nt, sigma = 1, range = 1, constr = TRUE)$u # %>% as.vector(.)

spt_sim <- sp_sim
rho_t <- 0.7 # Autoregressive parameter for the temporal structure
for(i in 2:nt){
  spt_sim[,i] <- rho_t*spt_sim[,i-1] + sqrt(1 - rho_t**2) * sp_sim[,i]
}

# Defined the matrix of constraints. The number of constraints is (mesh$n + nt)
# There is no real need to implent these constraints for this model (no extra spatial or temporal structures)

# vec.i_ind <- c()
# vec.j_ind <- c()
# for(i in 1:(nrow(spt_sim) + ncol(spt_sim))){
#   if(i<=ncol(spt_sim)){
#     vec.i_ind <- c(vec.i_ind, rep(i, time = length((((i-1)*(mesh$n)+1):(i*mesh$n)))))
#     vec.j_ind <- c(vec.j_ind, (((i-1)*(mesh$n)+1):(i*mesh$n)))
#   } else{
#     j <- i - nt
#     vec.i_ind <- c(vec.i_ind, rep(i, time = length(j + ((1:nt)-1)*mesh$n)))
#     vec.j_ind <- c(vec.j_ind, j + ((1:nt)-1)*mesh$n)
#   }
# }
# 
# C <- sparseMatrix(i = vec.i_ind, j = vec.j_ind, x = 1, dims = c(mesh$n+nt, mesh$n*nt))

# Vector of spt effect
spt_sim <- spt_sim %>% as.vector(.)
spt_sim_1 <- spt_sim - mean(spt_sim) # Apply this constraint to vanish the global mean from the spt effect. 

# CCt <- C %*% t(C)
# CCt_inv <- solve(CCt)
# Cconst <- t(C) %*% CCt_inv %*% C
# # Set the constraints
# spt_sim_constr <- spt_sim - Cconst %*% spt_sim
  
# Simulate linear predicor
lin_pred1 <- A_inf_spt %*% spt_sim_1 + beta01

### Simulation second linear predictor ----

beta02 <- -1
sp_sim <- simulation_SPDE(mesh = mesh, n_rep = nt, sigma = 2, range = 4, constr = TRUE)$u # %>% as.vector(.)

spt_sim <- sp_sim
rho_t <- 0.3 # Autoregressive parameter for the temporal structure
for(i in 2:nt){
  spt_sim[,i] <- rho_t*spt_sim[,i-1] + sqrt(1 - rho_t**2) * sp_sim[,i]
}

# Defined the matrix of constraints. The number of constraints is (mesh$n + nt).
# Note: This is not the proper way to introduce the constraints (although it is the standard in MCMC). It will introduce a small bias, but it is faster.
# There is no real need to implent these constraints for this model (no extra spatial or temporal structures)

# vec.i_ind <- c()
# vec.j_ind <- c()
# for(i in 1:(nrow(spt_sim) + ncol(spt_sim))){
#   if(i<=ncol(spt_sim)){
#     vec.i_ind <- c(vec.i_ind, rep(i, time = length((((i-1)*(mesh$n)+1):(i*mesh$n)))))
#     vec.j_ind <- c(vec.j_ind, (((i-1)*(mesh$n)+1):(i*mesh$n)))
#   } else{
#     j <- i - nt
#     vec.i_ind <- c(vec.i_ind, rep(i, time = length(j + ((1:nt)-1)*mesh$n)))
#     vec.j_ind <- c(vec.j_ind, j + ((1:nt)-1)*mesh$n)
#   }
# }
# 
# C <- sparseMatrix(i = vec.i_ind, j = vec.j_ind, x = 1, dims = c(mesh$n+nt, mesh$n*nt))

# Vector of spt effect
spt_sim <- spt_sim %>% as.vector(.)
spt_sim_2 <- spt_sim - mean(spt_sim) # Apply this constraint to vanish the global mean from the spt effect. 

# CCt <- C %*% t(C)
# CCt_inv <- solve(CCt)
# Cconst <- t(C) %*% CCt_inv %*% C
# # Set the constraints
# spt_sim_constr <- spt_sim - Cconst %*% spt_sim

# Simulate linear predicor
lin_pred2 <- A_inf_spt %*% spt_sim_2 + beta02

### Simulation of the 2-MVN data from the linear_predictors ----

lin_pred <- cbind(y1 = lin_pred1, y2 = lin_pred2)

taug1 <- 20; taug2 <- 40; rho12 <- 0.5
Q <- Matrix(data = 0, ncol = 2, nrow = 2)
Q[1,2] <- rho12/sqrt(taug1*taug2)
Q <- t(Q) + Q
diag(Q) <- c(taug1, taug2)

y_MVN <- simulation_MVN(Q = Q, mu = lin_pred)

### Aggregation of the data ----

st_crs(grid.integration) <- st_crs(ES_voronoi_1)
idx.intersects_1 <- st_intersects(x = grid.integration, y = ES_voronoi_1, sparse = TRUE) %>% as.numeric(.)
idx.intersects_2 <- st_intersects(x = grid.integration, y = ES_voronoi_2, sparse = TRUE) %>% as.numeric(.)
idx.intersects_3 <- st_intersects(x = grid.integration, y = ES_voronoi_3, sparse = TRUE) %>% as.numeric(.)

# nt, rep_nt
y1 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(y1) <- c("Group.1", "x")
y2 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(y2) <- c("Group.1", "x")

Support_divgroup_temp <- 1:2 * (nt/3)

for(i in 1:nt){
  if(i<=Support_divgroup_temp[1]){
    agg_y1 <- aggregate(y_MVN[((i-1)*length(grid.integration)+1):(i*length(grid.integration)),1], by = list(idx.intersects_1), FUN = mean) 
    agg_y2 <- aggregate(y_MVN[((i-1)*length(grid.integration)+1):(i*length(grid.integration)),2], by = list(idx.intersects_1), FUN = mean) 
    y1 <- rbind(y1, agg_y1)
    y2 <- rbind(y2, agg_y2)
  } else if(i>Support_divgroup_temp[1] & i<=Support_divgroup_temp[2]){
    agg_y1 <- aggregate(y_MVN[((i-1)*length(grid.integration)+1):(i*length(grid.integration)),1], by = list(idx.intersects_2), FUN = mean) 
    agg_y2 <- aggregate(y_MVN[((i-1)*length(grid.integration)+1):(i*length(grid.integration)),2], by = list(idx.intersects_2), FUN = mean)
    y1 <- rbind(y1, agg_y1)
    y2 <- rbind(y2, agg_y2)
  } else{
    agg_y1 <- aggregate(y_MVN[((i-1)*length(grid.integration)+1):(i*length(grid.integration)),1], by = list(idx.intersects_3), FUN = mean) 
    agg_y2 <- aggregate(y_MVN[((i-1)*length(grid.integration)+1):(i*length(grid.integration)),2], by = list(idx.intersects_3), FUN = mean)
    y1 <- rbind(y1, agg_y1)
    y2 <- rbind(y2, agg_y2)
  }
}

ggplot() + geom_sf(data = ES_voronoi_1, mapping = aes(fill = y1[1:200,2])) + scale_fill_viridis_c(option = "mako")
ggplot() + geom_sf(data = ES_voronoi_1, mapping = aes(fill = y2[1:200,2])) + scale_fill_viridis_c(option = "mako")

### Sequential consensus ----

# The number of temporal related to each different spatial support is nt/3
# To simplify the inferential procedure, simply choose an integer value proportional to (nt/3) for nt_inf
nt_inf <- 5 # The number of temporal nodes for the sequential consensus procedure

if(is.integer((nt/3) / nt_inf)){warning("nt_/3 doesn’t seem to be an integer multiple of nt_inf. Careful consideration of the different spatial support index is required...")}

spde.inf <- inla.spde2.pcmatern(mesh = mesh, prior.range = c(8/5, 0.5), prior.sigma = c(1, 0.5),  constr = TRUE)
spde.idx.inf1 <- inla.spde.make.index(name = "spde.inf1", n.spde = spde.inf$n.spde, n.group = nt_inf)
spde.idx.inf2 <- inla.spde.make.index(name = "spde.inf2", n.spde = spde.inf$n.spde, n.group = nt_inf)

A_inf_spt_inf <- inla.spde.make.A(mesh = mesh, loc = grid.integration %>% rep(., nt_inf) %>% st_coordinates(.), group = rep(1:nt_inf, each = length(grid.integration)))

Sp_supp <- list(ES_voronoi_1, ES_voronoi_2, ES_voronoi_3)
list_idx.intersects <- list(idx.intersects_1, idx.intersects_2, idx.intersects_3)

# Sequential Consensus Alg. 1
# The number of models to evaluate (nt/nt_inf) 
# In case the model is too big, use adaptive memory approach (R does not manage memory efficiently):
# Work around the RAM and the storage memory: writting (tmp files if automated), removing, cleaning and re-reading the require data (storage in tmp files) at each step 

list_sqc_models <- list()
for(i in 1:(nt/nt_inf)){
  cat(paste0('Step ', i, ' of ', (nt/nt_inf),'.\n'))
  
  Spatial_supp_idx <- which(c(
    nt_inf*i<=Support_divgroup_temp[1],
    nt_inf*i>Support_divgroup_temp[1] & nt_inf*i<=Support_divgroup_temp[2],
    nt_inf*i>Support_divgroup_temp[2])
    )
  
  idx <- c()
  for(k in 1:nt_inf){
    idx <- c(idx, list_idx.intersects[[Spatial_supp_idx]]+(k-1)*max(list_idx.intersects[[Spatial_supp_idx]]))
  }
  
  A_inf_spt_block <- inla.spde.make.block.A(A = A_inf_spt_inf, block = idx, rescale = "count")
  
  # Stack for the first ALR comp.
  inf_y1 <- inla.stack(data = list(y = y1[1:(nrow(Sp_supp[[Spatial_supp_idx]])*nt_inf),2]),
                       A = list(A_inf_spt_block, 1),
                       effects = list(
                         spde.idx.inf1,
                         list(
                           beta01 = rep(1, nrow(Sp_supp[[Spatial_supp_idx]])*nt_inf),
                           u = (1):(nrow(Sp_supp[[Spatial_supp_idx]])*nt_inf)   
                         )
                       ),
                       tag = "inf_y1")
  
  # Stack for the second ALR comp.
  inf_y2 <- inla.stack(data = list(y = y2[1:(nrow(Sp_supp[[Spatial_supp_idx]])*nt_inf),2]),
                       A = list(A_inf_spt_block, 1),
                       effects = list(
                         spde.idx.inf2,
                         list(
                           beta02 = rep(1, nrow(Sp_supp[[Spatial_supp_idx]])*nt_inf),
                           u = (nrow(Sp_supp[[Spatial_supp_idx]])*nt_inf+1):(2*nrow(Sp_supp[[Spatial_supp_idx]])*nt_inf)
                         )
                       ),
                       tag = "inf_y2")
  
  inf_yt <- inla.stack(inf_y1, inf_y2)
  
  if(i == 1){
    formula_sqc <- y ~ -1 + 
      beta01 + f(spde.inf1, model = spde.inf, group = spde.inf1.group, control.group = list(model = "ar1")) +
      beta02 + f(spde.inf2, model = spde.inf, group = spde.inf2.group, control.group = list(model = "ar1")) +
      f(u, model = "iid2d", n = 2*(nrow(Sp_supp[[Spatial_supp_idx]])*nt_inf), constr = TRUE)
    
    sqc_model <- inla(formula = formula_sqc, data = inla.stack.data(inf_yt), family = "gaussian",
                      control.predictor = list(A = inla.stack.A(inf_yt)),
                      control.compute = list(config = FALSE), # We don't need to store the conditional posterior as they are different between partitions
                      control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE))),
                      verbose = FALSE)
  } else{
    spde.inf_1 <- inla.spde2.pcmatern(mesh = mesh, prior.range = c(sqc_model$summary.hyperpar$`0.5quant`[1], 0.5), 
                                    prior.sigma = c(sqc_model$summary.hyperpar$`0.5quant`[2], 0.5), constr = TRUE)
    spde.inf_2 <- inla.spde2.pcmatern(mesh = mesh, prior.range = c(sqc_model$summary.hyperpar$`0.5quant`[4], 0.5), 
                                      prior.sigma = c(sqc_model$summary.hyperpar$`0.5quant`[5], 0.5), constr = TRUE)
    
    formula_sqc_mod <- y ~ -1 + 
      beta01 + f(spde.inf1, model = spde.inf_1, group = spde.inf1.group, 
                 control.group = list(model = "ar1",
                                      hyper = list(rho = list(prior = "gaussian", mean = list_sqc_models[[i-1]]$hyper$mean[3], prec = list_sqc_models[[i-1]]$hyper$sd[3]**(-2))))
                 ) +
      beta02 + f(spde.inf2, model = spde.inf_2, group = spde.inf2.group, 
                 control.group = list(model = "ar1",
                                      hyper = list(rho = list(prior = "gaussian", mean = list_sqc_models[[i-1]]$hyper$mean[6], prec = list_sqc_models[[i-1]]$hyper$sd[6]**(-2)))) 
                 ) +
      f(u, model = "iid2d", n = 2*(nrow(Sp_supp[[Spatial_supp_idx]])*nt_inf), constr = TRUE) # Can not be changed, fixed to a Wishart prior
    
    sqc_model <- inla(formula = formula_sqc_mod, data = inla.stack.data(inf_yt), family = "gaussian",
                      control.predictor = list(A = inla.stack.A(inf_yt)),
                      control.mode = list(theta = sqc_model$mode$theta, fixed = FALSE, restart = TRUE), # Init locations, leverage prior results
                      control.fixed = list(
                        mean=list(beta01=list_sqc_models[[i-1]]$fixed.eff$mean[1], beta02=list_sqc_models[[i-1]]$fixed.eff$mean[2], default=0),
                        prec=list(beta01=list_sqc_models[[i-1]]$fixed.eff$sd[1]**(-2), beta02=list_sqc_models[[i-1]]$fixed.eff$sd[2]**(-2), default=0.01)
                      ),
                      control.compute = list(config = FALSE), # We dont need to store the conditional posterior as they are different between partitions
                      control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE))),
                      verbose = FALSE)
  }
  
  list_sqc_models[[i]] <- list(
    fixed.eff = list(mean = sqc_model$summary.fixed$mean, sd = sqc_model$summary.fixed$sd), # Fixed effects
    # In this case we don't need to save the spatio-temporal components
    # spt_1 = list(mean.cor = , 
    #            Qpost = ),
    # spt_2 = list(mean.cor = ,
    #              Qpost = ),
    hyper = list(mean = sqc_model$internal.summary.hyperpar$mean, sd = sqc_model$internal.summary.hyperpar$sd) # Hyperparameters
  )
}

A_pred <- fm_basis(x = mesh, loc = grid.integration)


colsc_fill <- function(...) {
  scale_fill_gradientn(
    colours = viridis::mako(n = 10),
    limits = range(..., na.rm = TRUE)
  )
}

colsc_color <- function(...) {
  scale_color_gradientn(
    colours = viridis::mako(n = 10),
    limits = range(..., na.rm = TRUE)
  )
}

csc_c1 <- colsc_color(c(y1[1:200,2]-beta01, drop(A_pred %*% sqc_model$summary.random$spde.inf1$mean[1:mesh$n])))
csc_f1 <- colsc_fill(c(y1[1:200,2]-beta01, drop(A_pred %*% sqc_model$summary.random$spde.inf1$mean[1:mesh$n])))

gg1.bd <- ggplot() + 
  geom_sf(data = st_sf(grid.integration) %>% mutate(., id = "Spatial effect ALR 1 (downscaled)"), mapping = aes(color = drop(A_pred %*% sqc_model$summary.random$spde.inf1$mean[1:mesh$n]))) +
  geom_sf(data = ES_voronoi_1 %>% mutate(., id = "Spatial effect ALR 1 (aggregated)"), mapping = aes(fill = y1[1:200,2]-beta01)) + 
  labs(color = "Values", title = "Spatial structure for the first temporal node") +
  facet_wrap(facets = ~ id, ncol = 2) + theme_bw() +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

gg1_bd <- gg1.bd + csc_c1 + csc_f1 + guides(fill="none")

csc_c2 <- colsc_color(c(y2[1:200,2]-beta02,drop(A_pred %*% sqc_model$summary.random$spde.inf2$mean[1:mesh$n])))
csc_f2 <- colsc_fill(c(y2[1:200,2]-beta02,drop(A_pred %*% sqc_model$summary.random$spde.inf2$mean[1:mesh$n])))

gg2.bd <- ggplot() + 
  geom_sf(data = st_sf(grid.integration) %>% mutate(., id = "Spatial effect ALR 2 (downscaled)"), mapping = aes(color = drop(A_pred %*% sqc_model$summary.random$spde.inf2$mean[1:mesh$n]))) +
  geom_sf(data = ES_voronoi_1 %>% mutate(., id = "Spatial effect ALR 2 (aggregated)"), mapping = aes(fill = y2[1:200,2]-beta02)) + 
  labs(color = "Values") +
  facet_wrap(facets = ~ id, ncol = 2) + theme_bw() +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

gg2_bd <- gg2.bd + csc_c2 + csc_f2 + guides(fill="none")

gg12_bd <- ggplot() + 
  geom_sf(data = st_sf(grid.integration) %>% mutate(., id = "Spatial effect ALR 1 (downscaled)"), mapping = aes(color = drop(A_pred %*% sqc_model$summary.random$spde.inf1$mean[1:mesh$n]))) +
  geom_sf(data = ES_voronoi_1 %>% mutate(., id = "Spatial effect ALR 1 (aggregated)"), mapping = aes(fill = y1[1:200,2]-beta01)) + 
  labs(color = "Values", title = "Spatial effect (t = 1)") +
  csc_c1 + csc_f1 + guides(fill="none") +
  
  new_scale_color() +
  
  geom_sf(data = st_sf(grid.integration) %>% mutate(., id = "Spatial effect ALR 2 (downscaled)"), mapping = aes(color = drop(A_pred %*% sqc_model$summary.random$spde.inf2$mean[1:mesh$n]))) +
  geom_sf(data = ES_voronoi_1 %>% mutate(., id = "Spatial effect ALR 2 (aggregated)"), mapping = aes(fill = y2[1:200,2]-beta02)) + 
  labs(color = "Values") +
  csc_c2 + csc_f2 + guides(fill="none") +
  
  facet_wrap(facets = ~ id, ncol = 2) + theme_bw() +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

ggplot() + gg(mesh) + geom_sf(data = ES_mainland_boundary, mapping = aes(), alpha = 0.5) + theme_bw()

gg_base <- ggplot() + 
  gg(mesh) +
  geom_sf(data = st_sf(ES_mainland_boundary) %>% mutate(., id = "Mesh"), mapping = aes(), alpha = 0.5) +
  geom_sf(data = ES_voronoi_1 %>% mutate(., id = "Spatial support 1"), mapping = aes()) + 
  geom_sf(data = ES_voronoi_2 %>% mutate(., id = "Spatial support 2"), mapping = aes()) +
  geom_sf(data = ES_voronoi_3 %>% mutate(., id = "Spatial support 3"), mapping = aes()) +
  theme_bw() + labs(title = "Mesh and spatial supports") +
  # labs(title = TeX("\\textbf{Mesh and spatial supports} ($\\log(\\lambda_{\\beta})$)")) +
  facet_wrap(facets = ~ id, ncol = 2) +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

grid.arrange(arrangeGrob(grobs = list(gg_base, gg1_bd, gg2_bd), layout_matrix = matrix(data = c(1,2,1,3), ncol = 2, byrow = TRUE)))
# gg_tot <- grid.arrange(arrangeGrob(grobs = list(gg_base, gg12_bd), layout_matrix = matrix(data = c(1,2), ncol = 2, byrow = TRUE)))
# gg_tot$grobs[[1]]$widths <- c(1,1.29)
# 
# grid::grid.newpage()
# grid::grid.draw(gg_tot)

# ggsave(plot = gg_tot, filename = "./fig_BigData_example_hq.eps", device = "eps", width = 1295, height = 875, units = "px", dpi = 100) # High quality plot

