# =============================================================================
# Article Title: Unveiling land use dynamics: Insights from a hierarchical
#                Bayesian spatio-temporal modelling of Compositional Data
#
# Figure Code: Figure 3. In section: "Dealing with 0's example"
# Authors: Mario Figueira Pereira
# 
# Description:
# This script provides the code to obtain the results shown in Figure 3, 
# from data simulation to data analysis and graphical visualization steps.
#
# Last Updated: 30/07/2024
# =============================================================================

# Dealing with 0's. (Figure 3) 

# Required libraries ----

library(Matrix)
library(INLA)
library(compositions)
library(ggplot2)
library(gridExtra)
library(ggtext)
library(dplyr)
library(sf)
library(giscoR)

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

# Function to simulate a spatial effect with Leroux structure

simulation_Leroux <- function(tau = 1, lambda = 0.65, neigh_structure, constr=TRUE, seed){
  # tau: the marginal precision of the CAR prior (the variability between first-order neighbors)
  # neigh_structure: the adjacency or neighborhood structure (usually expressed as W)
  # lambda: the autoregression coefficient
  # constr: an argument that, if TRUE, forces the GMRF to have a global mean of zero
  # seed: seed to reproduce the result
  if(!missing(seed)){set.seed(seed)}
  W = matrix(0, nrow = length(neigh_structure), ncol = length(neigh_structure))
  for(i in 1:length(neigh_structure)){
    if(length(neigh_structure) > 0){
      W[i, neigh_structure[[i]]] = 1
    }
  }
  D = Diagonal(x = unlist(lapply(neigh_structure, length))) # D is a diagonal matrix defined as: D_ij = n(i~j)
  Q = tau * (lambda * (D - W - Diagonal(x=rep(1, nrow(W)))) + Diagonal(x=rep(1, nrow(W)))); Q = as(object = Q, "sparseMatrix")
  R = Matrix::chol(Q)
  w = rnorm(length(neigh_structure))
  u_sp = backsolve(R, w)
  if(constr){
    u_sp = u_sp - mean(u_sp)
  }
  return(list(u=u_sp, Q=Q))
}

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

# Function to create an adjaccency matrix given an object of class <<nb>>

adj_matrix <- function(nb){
  W <- lapply(X = 1:length(nb), FUN = function(i){
    m0 <- Matrix(data = 0, ncol = length(nb), nrow = 1)
    m0[nb[[i]]] <- 1
    return(m0)
  }) %>% do.call(., what = rbind)
  return(W)
}

# INLA graph from an object of class <<nb>>

inla.graph_nb <- function(nb){
  return(inla.matrix2graph(graph = adj_matrix(nb = nb)))
}

# Function to harmonize different plot scales 

colsc <- function(...) {
  scale_fill_gradientn(
    colours = viridis::mako(n = 11),
    limits = range(..., na.rm = TRUE)
  )
}

# Setting a global seed for replication of results ----

seed <- 123
set.seed(seed = seed)

# Preparing the Spatial structure ----

ES_Polygons <- gisco_get_countries(resolution = "01", country = "ES")$geometry %>% st_cast(x = ., to = "POLYGON")
ES_mainland_boundary <- ES_Polygons[ES_Polygons %>% st_area(.) %>% which.max(.),]

ggplot() + geom_sf(data = ES_mainland_boundary, mapping = aes())

ES_voronoi <- split_poly(sf_poly = ES_mainland_boundary, n_areas = 300)
ggplot() + geom_sf(data = ES_voronoi, mapping = aes())

Es_voronoi_nb <- spdep::poly2nb(pl = ES_voronoi, queen = TRUE) # neighborhood structure
sp_Leroux1 <- simulation_Leroux(tau = 1, lambda = 0.4, neigh_structure = Es_voronoi_nb, constr=TRUE)$u # Leroux spatial effect 1
sp_Leroux2 <- simulation_Leroux(tau = 1.5, lambda = 0.9, neigh_structure = Es_voronoi_nb, constr=TRUE)$u # Leroux spatial effect 2
sp_Leroux3 <- simulation_Leroux(tau = 2, lambda = 0.6, neigh_structure = Es_voronoi_nb, constr=TRUE)$u # Leroux spatial effect 2

DFsim <- st_sf(data.frame(id = 1:length(Es_voronoi_nb), sp1 = sp_Leroux1, sp2 = sp_Leroux2, sp3 = sp_Leroux3), geometry = ES_voronoi$geometry)
ggplot() + geom_sf(data = DFsim, mapping = aes(fill = sp1)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim, mapping = aes(fill = sp2)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim, mapping = aes(fill = sp3)) + scale_fill_viridis_c(option = "mako") + theme_bw()

# CoDa example: presence of 0's in the first CLR (conditional simulation of its values) ----

## Simulation of the linear predictor data: shared between the Bernoulli and MVN related to the first CLR ----

x <- runif(length(Es_voronoi_nb), min = 0, max = 1)

beta01 <- -2; beta11 <- 2
beta02 <- -1; beta12 <- 3
beta03 <- 1; beta13 <- 1
lin_pred1 <- beta01 + beta11*x + sp_Leroux1
lin_pred2 <- beta02 + beta12*x + sp_Leroux2
lin_pred3 <- beta03 + beta13*x + sp_Leroux3

## Simulation of the MVN for the 3 CLR ----

Q <- Matrix(data = 0, nrow = 3, ncol = 3)
tau1 <- 20; tau2 <- 10; tau3 <- 50
rho12 <- 0.2; rho13 <- 0.3; rho23 <- 0.8
Q[1,2] <- rho12/sqrt(tau1*tau2)
Q[1,3] <- rho13/sqrt(tau1*tau3)
Q[2,3] <- rho23/sqrt(tau2*tau3)
Q <- t(Q) + Q
diag(Q) <- c(tau1,tau2,tau3)

y_MVN <- simulation_MVN(Q = Q, mu = cbind(lin_pred1, lin_pred2, lin_pred3))

DFsim1 <- DFsim %>% mutate(.data = ., .before = geometry, lin_pred1, y1 = y_MVN[,1])
DFsim2 <- DFsim %>% mutate(.data = ., .before = geometry, lin_pred2, y2 = y_MVN[,2])
DFsim3 <- DFsim %>% mutate(.data = ., .before = geometry, lin_pred3, y3 = y_MVN[,3])

ggplot() + geom_sf(data = DFsim1, mapping = aes(fill = lin_pred1)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim2, mapping = aes(fill = lin_pred2)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim3, mapping = aes(fill = lin_pred3)) + scale_fill_viridis_c(option = "mako") + theme_bw()

ggplot() + geom_sf(data = DFsim1, mapping = aes(fill = y1)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim2, mapping = aes(fill = y2)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim3, mapping = aes(fill = y3)) + scale_fill_viridis_c(option = "mako") + theme_bw()

## Simulation of 1's and 0's: Bernoulli distribution ----

beta0.ber <- -1
beta0_ber <- beta0.ber - beta01

prob <- exp(beta0_ber + lin_pred1)/(1 + exp(beta0_ber + lin_pred1))
y_ber <- rbinom(n = length(lin_pred1), size = 1, prob = prob)
DFsim.ber <- DFsim %>% mutate(.data = ., .before = geometry, x = x, prob, y.ber = y_ber)

ggplot() + geom_sf(data = DFsim.ber, mapping = aes(fill = prob)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim.ber, mapping = aes(fill = as.factor(y.ber))) + theme_bw()

## Inference without considering the rows with some zero value ----

# spdep::nb2INLA(file = "./Data/ES_voronoi_nb_MVN", nb = Es_voronoi_nb) # writing the file for the adjacency matrix
# g <- inla.read.graph(filename = "./Data/ES_voronoi_nb_MVN") # reading the file related to the adjacency matrix
g <- inla.graph_nb(nb = Es_voronoi_nb)

inf_y1 <- inla.stack(data = list(y = sapply(X = 1:nrow(y_MVN), FUN = function(i){if(DFsim.ber$y.ber[i]==1){y_MVN[i,1]}else{NA}})),
                     A = list(1),
                     effects = list(
                       list(beta01 = rep(1,nrow(y_MVN)),
                            beta11 = x,
                            sp1 = 1:nrow(y_MVN),
                            u = 1:nrow(y_MVN)
                            )
                     ),
                     tag = "inf_y1")

inf_y2 <- inla.stack(data = list(y = sapply(X = 1:nrow(y_MVN), FUN = function(i){if(DFsim.ber$y.ber[i]==1){y_MVN[i,2]}else{NA}})),
                     A = list(1),
                     effects = list(
                       list(beta02 = rep(1,nrow(y_MVN)),
                            beta12 = x,
                            sp2 = 1:nrow(y_MVN),
                            u = (nrow(y_MVN)+1):(2*nrow(y_MVN))
                       )
                     ),
                     tag = "inf_y2")

inf_y3 <- inla.stack(data = list(y = sapply(X = 1:nrow(y_MVN), FUN = function(i){if(DFsim.ber$y.ber[i]==1){y_MVN[i,3]}else{NA}})),
                     A = list(1),
                     effects = list(
                       list(beta03 = rep(1,nrow(y_MVN)),
                            beta13 = x,
                            sp3 = 1:nrow(y_MVN),
                            u = (2*nrow(y_MVN)+1):(3*nrow(y_MVN))
                       )
                     ),
                     tag = "inf_y3")

inf_yMVN <- inla.stack(inf_y1, inf_y2, inf_y3)
formula_yMVN <- y ~ -1 + 
  beta01 + beta11 + f(sp1, model = "besagproper2", graph = g, constr = TRUE) +
  beta02 + beta12 + f(sp2, model = "besagproper2", graph = g, constr = TRUE) +
  beta03 + beta13 + f(sp3, model = "besagproper2", graph = g, constr = TRUE) +
  f(u, model = "iid3d", n = 3*nrow(y_MVN), constr = TRUE)
coda_model <- inla(data = inla.stack.data(inf_yMVN), formula = formula_yMVN, family = "gaussian",
                   control.predictor = list(A = inla.stack.A(inf_yMVN)),
                   control.compute = list(dic = FALSE, waic = FALSE),
                   control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE))),
                   verbose = FALSE)

## Inference considering the rows with zero values (Hurdle like model) ----

inf_y1 <- inla.stack(data = list(y = cbind(sapply(X = 1:nrow(y_MVN), FUN = function(i){if(DFsim.ber$y.ber[i]==1){y_MVN[i,1]}else{NA}}), NA, NA)),
                     A = list(1),
                     effects = list(
                       list(beta01 = rep(1,nrow(y_MVN)),
                            beta11 = x,
                            sp1 = 1:nrow(y_MVN),
                            u = 1:nrow(y_MVN)
                       )
                     ),
                     tag = "inf_y1")

inf_y2 <- inla.stack(data = list(y = cbind(y_MVN[,2], NA, NA)),
                     A = list(1),
                     effects = list(
                       list(beta02 = rep(1,nrow(y_MVN)),
                            beta12 = x,
                            sp2 = 1:nrow(y_MVN),
                            u = (nrow(y_MVN)+1):(2*nrow(y_MVN))
                       )
                     ),
                     tag = "inf_y2")

inf_y3 <- inla.stack(data = list(y = cbind(y_MVN[,3], NA, NA)),
                     A = list(1),
                     effects = list(
                       list(beta03 = rep(1,nrow(y_MVN)),
                            beta13 = x,
                            sp3 = 1:nrow(y_MVN),
                            u = (2*nrow(y_MVN)+1):(3*nrow(y_MVN))
                       )
                     ),
                     tag = "inf_y3")

inf_yber <- inla.stack(data = list(y = cbind(NA, DFsim.ber$y.ber, NA)),
                       A = list(1),
                       effects = list(
                         list(
                           beta0.ber = rep(1, nrow(DFsim.ber)),
                           iid.copied = 1:nrow(DFsim.ber)
                           )
                         ),
                       tag = "inf_yber"
                       )

inf_copy <- inla.stack(data = list(y = cbind(NA, NA, rep(0, nrow(DFsim.ber)))),
                       A = list(1,-1),
                       effects = list(
                         list(beta01 = rep(1,nrow(y_MVN)),
                              beta11 = x,
                              sp1 = 1:nrow(y_MVN)
                              ),
                         list(iid.copy = 1:nrow(DFsim.ber))
                         ),
                       tag = "inf_copy"
                       )

inf_yMVNh <- inla.stack(inf_y1, inf_y2, inf_y3, inf_yber, inf_copy)
formula_yMVNh <- y ~ -1 + 
  beta0.ber + f(iid.copied, copy = "iid.copy", fixed = FALSE) + 
  f(iid.copy, model = "iid", hyper = list(prec = list(initial = -10, fixed = TRUE))) +
  beta01 + beta11 + f(sp1, model = "besagproper2", graph = g, constr = TRUE) +
  beta02 + beta12 + f(sp2, model = "besagproper2", graph = g, constr = TRUE) +
  beta03 + beta13 + f(sp3, model = "besagproper2", graph = g, constr = TRUE) +
  f(u, model = "iid3d", n = 3*nrow(y_MVN), constr = TRUE)
codah_model <- inla(data = inla.stack.data(inf_yMVNh), formula = formula_yMVNh, family = c("gaussian", "binomial", "gaussian"),
                    control.predictor = list(A = inla.stack.A(inf_yMVNh)),
                    control.compute = list(dic = FALSE, waic = FALSE),
                    control.family = list(list(hyper = list(prec = list(initial = 10, fixed = TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed = TRUE)))),
                    verbose = FALSE)

coda_model$summary.fixed
codah_model$summary.fixed

rbind(beta0_ber, beta01, beta11, beta02, beta12, beta03, beta13)

coda_model$summary.hyperpar
codah_model$summary.hyperpar

sqrt(mean((DFsim1$sp1 - coda_model$summary.random$sp1$mean)**2))
sqrt(mean((DFsim1$sp1 - codah_model$summary.random$sp1$mean)**2))

sqrt(mean((DFsim2$sp2 - coda_model$summary.random$sp2$mean)**2))
sqrt(mean((DFsim2$sp2 - codah_model$summary.random$sp2$mean)**2))

sqrt(mean((DFsim3$sp3 - coda_model$summary.random$sp3$mean)**2))
sqrt(mean((DFsim3$sp3 - codah_model$summary.random$sp3$mean)**2))

# saveRDS(file = "./Data/ZerosExample_Models/coda_model.RDS", object = coda_model)
# saveRDS(file = "./Data/ZerosExample_Models/codah_model.RDS", object = codah_model)

# Graphical visualization of the results ----

colsc_sp1 <- colsc(DFsim1$sp1, codah_model$summary.random$sp1$mean)
colsc_sp2 <- colsc(DFsim2$sp2, codah_model$summary.random$sp2$mean)
colsc_sp3 <- colsc(DFsim2$sp3, codah_model$summary.random$sp3$mean)

gg_sp1.sim <- ggplot() + 
  geom_sf(data = DFsim1, mapping = aes(fill = sp1)) + 
  # scale_fill_viridis_c(option = "mako") +
  colsc_sp1 +
  labs(title = "Spatial effect Comp. 1 (simulation)", fill = "Mean") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_sp2.sim <- ggplot() + 
  geom_sf(data = DFsim2, mapping = aes(fill = sp2)) + 
  # scale_fill_viridis_c(option = "mako") + 
  colsc_sp2 +
  labs(title = "Spatial effect Comp. 2 (simulation)", fill = "Mean") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_sp3.sim <- ggplot() + 
  geom_sf(data = DFsim3, mapping = aes(fill = sp3)) + 
  # scale_fill_viridis_c(option = "mako") + 
  colsc_sp3 +
  labs(title = "Spatial effect Comp. 3 (simulation)", fill = "Mean") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_ber.sim <- ggplot() + 
  geom_sf(data = DFsim.ber, mapping = aes(fill = as.factor(y.ber))) + 
  labs(fill = "Bernoulli", title = "Zero and non-zero values") +
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.position = "bottom"
  )

gg_sp1.mean <- ggplot() +
  geom_sf(data = DFsim1 %>% mutate(sp.inf = codah_model$summary.random$sp1$mean, ID = "Hurdle like model"), mapping = aes(fill = sp.inf)) + 
  geom_sf(data = DFsim1 %>% mutate(sp.inf = coda_model$summary.random$sp1$mean, ID = "Standard model"), mapping = aes(fill = sp.inf)) +
  # scale_fill_viridis_c(option = "mako") +
  colsc_sp1 +
  labs(title = "Spatial effect Comp. 1 (inference)", fill = "Mean") +
  facet_wrap(facets = ~ ID, ncol = 2) + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_sp1.sd <- ggplot() +
  geom_sf(data = DFsim1 %>% mutate(sp.inf = codah_model$summary.random$sp1$sd, ID = "Hurdle like model"), mapping = aes(fill = sp.inf)) + 
  geom_sf(data = DFsim1 %>% mutate(sp.inf = coda_model$summary.random$sp1$sd, ID = "Standard model"), mapping = aes(fill = sp.inf)) +
  scale_fill_viridis_c(option = "mako") +
  labs(title = "Spatial effect Comp. 1 (inference)", fill = "Stdev.") +
  facet_wrap(facets = ~ ID, ncol = 2) + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_sp2.mean <- ggplot() +
  geom_sf(data = DFsim1 %>% mutate(sp.inf = codah_model$summary.random$sp2$mean, ID = "Hurdle like model"), mapping = aes(fill = sp.inf)) + 
  geom_sf(data = DFsim1 %>% mutate(sp.inf = coda_model$summary.random$sp2$mean, ID = "Standard model"), mapping = aes(fill = sp.inf)) +
  # scale_fill_viridis_c(option = "mako") +
  colsc_sp2 +
  labs(title = "Spatial effect Comp. 2 (inference)", fill = "Mean") +
  facet_wrap(facets = ~ ID, ncol = 2) + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_sp2.sd <- ggplot() +
  geom_sf(data = DFsim1 %>% mutate(sp.inf = codah_model$summary.random$sp2$sd, ID = "Hurdle like model"), mapping = aes(fill = sp.inf)) + 
  geom_sf(data = DFsim1 %>% mutate(sp.inf = coda_model$summary.random$sp2$sd, ID = "Standard model"), mapping = aes(fill = sp.inf)) +
  scale_fill_viridis_c(option = "mako") +
  labs(title = "Spatial effect Comp. 2 (inference)", fill = "Stdev.") +
  facet_wrap(facets = ~ ID, ncol = 2) + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_sp3.mean <- ggplot() +
  geom_sf(data = DFsim1 %>% mutate(sp.inf = codah_model$summary.random$sp3$mean, ID = "Hurdle like model"), mapping = aes(fill = sp.inf)) + 
  geom_sf(data = DFsim1 %>% mutate(sp.inf = coda_model$summary.random$sp3$mean, ID = "Standard model"), mapping = aes(fill = sp.inf)) +
  # scale_fill_viridis_c(option = "mako") 
  colsc_sp3 +
  labs(title = "Spatial effect Comp. 3 (inference)", fill = "Mean") +
  facet_wrap(facets = ~ ID, ncol = 2) + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_sp3.sd <- ggplot() +
  geom_sf(data = DFsim1 %>% mutate(sp.inf = codah_model$summary.random$sp3$sd, ID = "Hurdle like model"), mapping = aes(fill = sp.inf)) + 
  geom_sf(data = DFsim1 %>% mutate(sp.inf = coda_model$summary.random$sp3$sd, ID = "Standard model"), mapping = aes(fill = sp.inf)) +
  scale_fill_viridis_c(option = "mako") +
  labs(title = "Spatial effect Comp. 3 (inference)", fill = "Stdev.") +
  facet_wrap(facets = ~ ID, ncol = 2) + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_beta0 <- ggplot() + 
  geom_line(data = codah_model$marginals.fixed$beta01, mapping = aes(x = x, y = y, color = "Comp. 1", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.fixed$beta01, mapping = aes(x = x, y = y, color = "Comp. 1", linetype = "Standard")) +
  geom_vline(xintercept = beta01, color = "red") +
  geom_line(data = codah_model$marginals.fixed$beta02, mapping = aes(x = x, y = y, color = "Comp. 2", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.fixed$beta02, mapping = aes(x = x, y = y, color = "Comp. 2", linetype = "Standard")) +
  geom_vline(xintercept = beta02, color = "blue") +
  geom_line(data = codah_model$marginals.fixed$beta03, mapping = aes(x = x, y = y, color = "Comp. 3", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.fixed$beta03, mapping = aes(x = x, y = y, color = "Comp. 3", linetype = "Standard")) +
  geom_vline(xintercept = beta03, color = "forestgreen") +
  labs(title = expression(beta[0])) +
  theme_bw() +
  scale_color_manual(name = 'Composition',
                     breaks = c('Comp. 1', 'Comp. 2', 'Comp. 3'),
                     values = c('Comp. 1' = 'red', 'Comp. 2' = 'blue', 'Comp. 3' = 'forestgreen')) +
  scale_linetype_manual(name = 'Model', breaks = c('Standard', 'Hurdle like'), values = c('Standard' = 'solid', 'Hurdle like' = 'dashed')) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

gg_beta1 <- ggplot() + 
  geom_line(data = codah_model$marginals.fixed$beta11, mapping = aes(x = x, y = y, color = "Comp. 1", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.fixed$beta11, mapping = aes(x = x, y = y, color = "Comp. 1", linetype = "Standard")) +
  geom_vline(xintercept = beta11, color = "red") +
  geom_line(data = codah_model$marginals.fixed$beta12, mapping = aes(x = x, y = y, color = "Comp. 2", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.fixed$beta12, mapping = aes(x = x, y = y, color = "Comp. 2", linetype = "Standard")) +
  geom_vline(xintercept = beta12, color = "blue") +
  geom_line(data = codah_model$marginals.fixed$beta13, mapping = aes(x = x, y = y, color = "Comp. 3", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.fixed$beta13, mapping = aes(x = x, y = y, color = "Comp. 3", linetype = "Standard")) +
  geom_vline(xintercept = beta13, color = "forestgreen") +
  labs(title = expression(beta[1])) +
  theme_bw() +
  scale_color_manual(name = 'Composition',
                     breaks = c('Comp. 1', 'Comp. 2', 'Comp. 3'),
                     values = c('Comp. 1' = 'red', 'Comp. 2' = 'blue', 'Comp. 3' = 'forestgreen')) +
  scale_linetype_manual(name = 'Model', breaks = c('Standard', 'Hurdle like'), values = c('Standard' = 'solid', 'Hurdle like' = 'dashed')) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

gg_tau <- ggplot() + 
  geom_line(data = codah_model$marginals.hyperpar$`Precision for sp1`, mapping = aes(x = x, y = y, color = "Comp. 1", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.hyperpar$`Precision for sp1`, mapping = aes(x = x, y = y, color = "Comp. 1", linetype = "Standard")) +
  geom_vline(xintercept = 1, color = "red") +
  geom_line(data = codah_model$marginals.hyperpar$`Precision for sp2`, mapping = aes(x = x, y = y, color = "Comp. 2", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.hyperpar$`Precision for sp2`, mapping = aes(x = x, y = y, color = "Comp. 2", linetype = "Standard")) +
  geom_vline(xintercept = 1.5, color = "blue") +
  geom_line(data = codah_model$marginals.hyperpar$`Precision for sp3`, mapping = aes(x = x, y = y, color = "Comp. 3", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.hyperpar$`Precision for sp3`, mapping = aes(x = x, y = y, color = "Comp. 3", linetype = "Standard")) +
  geom_vline(xintercept = 2, color = "forestgreen") +
  labs(title = expression(tau)) +
  theme_bw() +
  scale_color_manual(name = 'Composition',
                     breaks = c('Comp. 1', 'Comp. 2', 'Comp. 3'),
                     values = c('Comp. 1' = 'red', 'Comp. 2' = 'blue', 'Comp. 3' = 'forestgreen')) +
  scale_linetype_manual(name = 'Model', breaks = c('Standard', 'Hurdle like'), values = c('Standard' = 'solid', 'Hurdle like' = 'dashed')) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

gg_rho <- ggplot() + 
  geom_line(data = codah_model$marginals.hyperpar$`Lambda for sp1`, mapping = aes(x = x, y = y, color = "Comp. 1", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.hyperpar$`Lambda for sp1`, mapping = aes(x = x, y = y, color = "Comp. 1", linetype = "Standard")) +
  geom_vline(xintercept = 0.4, color = "red") +
  geom_line(data = codah_model$marginals.hyperpar$`Lambda for sp2`, mapping = aes(x = x, y = y, color = "Comp. 2", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.hyperpar$`Lambda for sp2`, mapping = aes(x = x, y = y, color = "Comp. 2", linetype = "Standard")) +
  geom_vline(xintercept = 0.9, color = "blue") +
  geom_line(data = codah_model$marginals.hyperpar$`Lambda for sp3`, mapping = aes(x = x, y = y, color = "Comp. 3", linetype = "Hurdle like")) +
  geom_line(data = coda_model$marginals.hyperpar$`Lambda for sp3`, mapping = aes(x = x, y = y, color = "Comp. 3", linetype = "Standard")) +
  geom_vline(xintercept = 0.6, color = "forestgreen") +
  labs(title = expression(lambda)) +
  theme_bw() +
  scale_color_manual(name = 'Composition',
                     breaks = c('Comp. 1', 'Comp. 2', 'Comp. 3'),
                     values = c('Comp. 1' = 'red', 'Comp. 2' = 'blue', 'Comp. 3' = 'forestgreen')) +
  scale_linetype_manual(name = 'Model', breaks = c('Standard', 'Hurdle like'), values = c('Standard' = 'solid', 'Hurdle like' = 'dashed')) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

gg_tot <- grid.arrange(arrangeGrob(grobs = list(gg_sp1.sim, gg_sp2.sim, gg_sp3.sim, gg_ber.sim, gg_sp1.mean, gg_sp1.sd, gg_sp2.mean, gg_sp2.sd, gg_sp3.mean, gg_sp3.sd, gg_beta0, gg_beta1, gg_tau, gg_rho), 
                                   layout_matrix = matrix(data = c(1,1,5,5,6,6,2,2,7,7,8,8,3,3,9,9,10,10,4,4,11,11,12,12,4,4,13,13,14,14), ncol = 6, byrow = TRUE)))

ggsave(plot = gg_tot, filename = "./Figures/fig_CoDaZeros_example_hq.eps", device = "eps", width = 1597, height = 1195, units = "px", dpi = 100) # High quality plot
ggsave(plot = gg_tot, filename = "./Figures/fig_CoDaZeros_example_lq.png", device = "png", width = 1597, height = 1195, units = "px", dpi = 100) # Low quality plot
