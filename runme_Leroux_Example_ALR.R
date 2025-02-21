# =============================================================================
# Title: ALR example with complex Q_likelihood  
#
# Authors: Mario Figueira Pereira
# 
# Description:
# This script provides the R code to simulate ALR data and fit a multivariate 
# Gaussian model with correlations in the likelihood
# (coding it in the Latent Gaussian Field - a trick for INLA)
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

simulation_MVN <- function(Q, mu, R){
  # Q: the precision matrix of the Multivariate Gaussian distribution (MVN)
  # mu: matrix of the MVN mean
  # R:  upper triangular matrix from Cholesky decomposition. 
  #     It is better to provide it in case we are simulating from the same Q_likelihood.
  if(missing(R)){R = Matrix::chol(Q)}
  
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
diag(Q) <- c(1/tau1,1/tau2,1/tau3)
Q <- solve(Q)

R <- Matrix::chol(Q)
n_repli <- 1
y_MVN <- lapply(X = 1:n_repli, FUN = function(x){simulation_MVN(Q = Q, mu = cbind(lin_pred1, lin_pred2, lin_pred3), R = R)}) %>% do.call(., what = rbind)
y_InvALR_MVN <- compositions::alrInv(z = y_MVN) %>% as.data.frame(.) # Using compositions::alr(x = y_InvALR_MVN) we get back the y_MVN matrix
colnames(y_InvALR_MVN) <- c("inv_y1", "inv_y2", "inv_y3", "inv_y4")
DFsim_invy1 <- st_sf(data.frame(id = 1:length(Es_voronoi_nb), invy1 = y_InvALR_MVN$inv_y1[1:300]), geometry = ES_voronoi$geometry) # 1:300 bcs the number of areas is equal to 300
DFsim_invy2 <- st_sf(data.frame(id = 1:length(Es_voronoi_nb), invy2 = y_InvALR_MVN$inv_y2[1:300]), geometry = ES_voronoi$geometry) # 1:300 bcs the number of areas is equal to 300
DFsim_invy3 <- st_sf(data.frame(id = 1:length(Es_voronoi_nb), invy3 = y_InvALR_MVN$inv_y3[1:300]), geometry = ES_voronoi$geometry) # 1:300 bcs the number of areas is equal to 300
DFsim_invy4 <- st_sf(data.frame(id = 1:length(Es_voronoi_nb), invy4 = y_InvALR_MVN$inv_y4[1:300]), geometry = ES_voronoi$geometry) # 1:300 bcs the number of areas is equal to 300

DFsim1 <- DFsim %>% mutate(.data = ., .before = geometry, lin_pred1, y1 = y_MVN[1:300,1]) # 1:300 bcs the number of areas is equal to 300
DFsim2 <- DFsim %>% mutate(.data = ., .before = geometry, lin_pred2, y2 = y_MVN[1:300,2]) # 1:300 bcs the number of areas is equal to 300
DFsim3 <- DFsim %>% mutate(.data = ., .before = geometry, lin_pred3, y3 = y_MVN[1:300,3]) # 1:300 bcs the number of areas is equal to 300

ggplot() + geom_sf(data = DFsim1, mapping = aes(fill = lin_pred1)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim2, mapping = aes(fill = lin_pred2)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim3, mapping = aes(fill = lin_pred3)) + scale_fill_viridis_c(option = "mako") + theme_bw()

ggplot() + geom_sf(data = DFsim1, mapping = aes(fill = y1)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim2, mapping = aes(fill = y2)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim3, mapping = aes(fill = y3)) + scale_fill_viridis_c(option = "mako") + theme_bw()

ggplot() + geom_sf(data = DFsim_invy1, mapping = aes(fill = invy1)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim_invy2, mapping = aes(fill = invy2)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim_invy3, mapping = aes(fill = invy3)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim_invy4, mapping = aes(fill = invy4)) + scale_fill_viridis_c(option = "mako") + theme_bw()

## Inference without considering the rows with some zero value ----

g <- inla.graph_nb(nb = Es_voronoi_nb)

inf_y1 <- inla.stack(data = list(y = y_MVN[,1]),
                     A = list(1),
                     effects = list(
                       list(beta01 = rep(1,nrow(y_MVN)),
                            beta11 = rep(x, times = n_repli),
                            sp1 = rep(1:length(sp_Leroux1), times = n_repli),
                            u = 1:nrow(y_MVN)
                       )
                     ),
                     tag = "inf_y1")

inf_y2 <- inla.stack(data = list(y = y_MVN[,2]),
                     A = list(1),
                     effects = list(
                       list(beta02 = rep(1,nrow(y_MVN)),
                            beta12 = rep(x, times = n_repli),
                            sp2 = rep(1:length(sp_Leroux2), times = n_repli),
                            u = (nrow(y_MVN)+1):(2*nrow(y_MVN))
                       )
                     ),
                     tag = "inf_y2")

inf_y3 <- inla.stack(data = list(y = y_MVN[,3]),
                     A = list(1),
                     effects = list(
                       list(beta03 = rep(1,nrow(y_MVN)),
                            beta13 = rep(x, times = n_repli),
                            sp3 = rep(1:length(sp_Leroux3), times = n_repli),
                            u = (2*nrow(y_MVN)+1):(3*nrow(y_MVN))
                       )
                     ),
                     tag = "inf_y3")

inf_yMVN <- inla.stack(inf_y1, inf_y2, inf_y3)

formula_yMVN <- y ~ -1 + 
  beta01 + beta11 + f(sp1, model = "besagproper2", graph = g, constr = TRUE) +
  beta02 + beta12 + f(sp2, model = "besagproper2", graph = g, constr = TRUE) +
  beta03 + beta13 + # f(sp3, model = "besagproper2", graph = g, constr = TRUE) +
  f(u, model = "iid3d", n = 3*nrow(y_MVN), constr = TRUE)

coda_model <- inla(data = inla.stack.data(inf_yMVN), formula = formula_yMVN, family = "gaussian",
                   control.predictor = list(A = inla.stack.A(inf_yMVN)),
                   control.compute = list(dic = FALSE, waic = FALSE),
                   control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE))),
                   verbose = FALSE)

