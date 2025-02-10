# =============================================================================
# Article Title: Unveiling land use dynamics: Insights from a hierarchical
#                Bayesian spatio-temporal modelling of Compositional Data
#
# Figure Code: Figure 2. In section: "Dealing with 0's example"
# Authors: Mario Figueira Pereira
# 
# Description:
# This script provides the code to obtain the results shown in Figure 2, 
# from data simulation to data analysis and graphical visualization steps.
#
# Last Updated: 30/07/2024
# =============================================================================

# Dealing with 0's. (Figure 2) 

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

# Function to simulate a spatial effect with the Leroux structure 

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
sp_Leroux <- simulation_Leroux(tau = 1, lambda = 0.7, neigh_structure = Es_voronoi_nb, constr=TRUE)$u # Leroux spatial effect

DFsim <- st_sf(data.frame(id = 1:length(Es_voronoi_nb), sp = sp_Leroux), geometry = ES_voronoi$geometry)
ggplot() + geom_sf(data = DFsim, mapping = aes(fill = sp)) + scale_fill_viridis_c(option = "mako") + theme_bw()

# Beta example: presence of 0's in the data (conditional simulation of Beta values) ----

## Simulation of the linear predictor data: shared between the Bernoulli and Beta distributions ----

x <- runif(length(Es_voronoi_nb), min = 0, max = 1)

beta0 <- 2; beta1 <- -3
lin_pred <- beta0 + beta1*x + sp_Leroux

DFsim <- DFsim %>% mutate(.data = ., .before = geometry, lin_pred)

ggplot() + geom_sf(data = DFsim, mapping = aes(fill = lin_pred)) + scale_fill_viridis_c(option = "mako") + theme_bw()

## Simulation of 1's and 0's: Bernoulli distribution ----

beta0.ber <- 2
beta0_ber <- beta0.ber - beta0

prob <- exp(beta0_ber + lin_pred)/(1 + exp(beta0_ber + lin_pred))
y_ber <- rbinom(n = length(lin_pred), size = 1, prob = prob)
DFsim.ber <- DFsim %>% mutate(.data = ., .before = geometry, x = x, prob, y.ber = y_ber)

ggplot() + geom_sf(data = DFsim.ber, mapping = aes(fill = prob)) + scale_fill_viridis_c(option = "mako") + theme_bw()
ggplot() + geom_sf(data = DFsim.ber, mapping = aes(fill = as.factor(y.ber))) + theme_bw()

## Simulation of proportion values: Beta distribution ----

mean_beta <- exp(lin_pred)/(1 + exp(lin_pred))
var_beta <- 20**(-1) # defining the variance as the inverse of the precision, which is the paratemeter used internalyy by INLA
if(any(var_beta>=min(mean_beta*(1-mean_beta)))){
  var_beta <- min(mean_beta*(1-mean_beta))/1.15
  sprintf("The variance has to be redefined, as it must fulfil that: var_beta < mean_beta*(1-mean_beta), var_beta = %f, prec_beta = %f", var_beta, var_beta**(-1)) %>% cat(.)
}

a_beta <- ((1-mean_beta)/var_beta - 1/mean_beta)*mean_beta**2
b_beta <- a_beta*(1/mean_beta - 1)

cond <- TRUE
while(cond){
  y_beta <- rbeta(n = length(a_beta), shape1 = a_beta, shape2 = b_beta)
  cond <- any(y_beta == 1 | y_beta == 0)
}
DFsim.beta <- DFsim %>% mutate(.data = ., .before = geometry, x = x, mean.beta = mean_beta, y.beta = y_beta)

ggplot() + geom_sf(data = DFsim.beta, mapping = aes(fill = y.beta)) + scale_fill_viridis_c(option = "mako") + theme_bw()

## Inference using porportion only data (single model) ----

# spdep::nb2INLA(file = "./Data/ES_voronoi_nb", nb = Es_voronoi_nb) # writing the file for the adjacency matrix
# g <- inla.read.graph(filename = "./Data/ES_voronoi_nb") # reading the file related to the adjacency matrix
g <- inla.graph_nb(nb = Es_voronoi_nb)

inf_beta <- inla.stack(data = list(y = sapply(X = 1:nrow(DFsim.beta), FUN = function(i){if(DFsim.ber$y.ber[i]==1){DFsim.beta$y.beta[i]}else{NA}})),
                       A = list(1),
                       effects = list(
                         list(beta0 = rep(1, times = nrow(DFsim.beta)), 
                              beta1 = DFsim.beta$x, 
                              sp = 1:nrow(DFsim.beta))
                       ),
                       tag = "inf_beta")

beta_model <- inla(data = inla.stack.data(inf_beta), family = "beta",
                   formula = y ~ -1 + beta0 + beta1 + f(sp, model = "besagproper2", graph = g, constr = TRUE), 
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.fixed = list(prec = list(beta0 = 0.001, beta1 = 0.001)),
                   verbose = FALSE)

## Inference using a Hurdle model like with the Beta and Bernoulli data (joint model) ----

inf_hurdle_beta <- inla.stack(data = list(y = cbind(sapply(X = 1:nrow(DFsim.beta), FUN = function(i){if(DFsim.ber$y.ber[i]==1){DFsim.beta$y.beta[i]}else{NA}}), NA, NA)),
                              A = list(1),
                              effects = list(
                                list(beta0 = rep(1, times = nrow(DFsim.beta)), 
                                     beta1 = DFsim.beta$x, 
                                     sp = 1:nrow(DFsim.beta))
                              ),
                              tag = "inf_beta")

inf_hurdle_ber <- inla.stack(data = list(y = cbind(NA, DFsim.ber$y.ber, NA)),
                             A = list(1),
                             effects = list(
                               list(beta0_ber = rep(1, times = nrow(DFsim.ber)), 
                                    iid.copied = DFsim.ber$id)
                             ),
                             tag = "inf_ber")

inf_hurlde_copy <- inla.stack(data = list(y = cbind(NA, NA, rep(0, times = nrow(DFsim.beta)))),
                              A = list(1, -1),
                              effects = list(
                                list(beta0 = rep(1, times = nrow(DFsim.beta)), 
                                     beta1 = DFsim.beta$x,
                                     sp = 1:nrow(DFsim.ber)),
                                list(iid.copy = 1:nrow(DFsim.beta))
                              ),
                              tag = "inf_copy")

inf_total_stack <- inla.stack(inf_hurdle_beta, inf_hurdle_ber, inf_hurlde_copy)

formula_hurdle <- y ~ -1 + beta0 + beta1 + f(sp, model = "besagproper2", constr = TRUE, graph = g) + beta0_ber + f(iid.copied, copy = "iid.copy", fixed = FALSE) + f(iid.copy, model = "iid", hyper = list(prec = list(initial = -10, fixed = TRUE)))
hurdle_model <- inla(data = inla.stack.data(inf_total_stack), family = c("beta", "binomial", "gaussian"),
                     formula = formula_hurdle,
                     control.predictor = list(A = inla.stack.A(inf_total_stack)),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.fixed = list(prec = list(beta0 = 0.001, beta1 = 0.001, beta0_ber = 0.001)),
                     control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed = TRUE)))),
                     verbose = FALSE)

beta_model$summary.fixed
hurdle_model$summary.fixed

beta_model$summary.hyperpar
hurdle_model$summary.hyperpar

RMSE_beta <- sqrt(mean((DFsim$sp - beta_model$summary.random$sp$mean)**2))
RMSE_hurdle <- sqrt(mean((DFsim$sp - hurdle_model$summary.random$sp$mean)**2))

beta_model$dic$family.dic
hurdle_model$dic$family.dic

# saveRDS(file = "./Data/ZerosExample_Models/beta_model.RDS", object = beta_model)
# saveRDS(file = "./Data/ZerosExample_Models/hurdle_model.RDS", object = hurdle_model)

# Graphical visualization of the results ----

sp.mean_scale <- colsc(DFsim$sp, hurdle_model$summary.random$sp$mean, beta_model$summary.random$sp$mean)

gg_sp.sim <- ggplot() + 
  geom_sf(data = DFsim, mapping = aes(fill = sp)) + 
  # scale_fill_viridis_c(option = "mako") + 
  sp.mean_scale +
  labs(title = "Spatial effect (simulation)", fill = "Mean") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_ber.sim <- ggplot() + 
  geom_sf(data = DFsim.ber, mapping = aes(fill = as.factor(y.ber))) + 
  labs(fill = "Bernoulli", title = "Zero and non-zero values") +
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_sp.mean <- ggplot() +
  geom_sf(data = DFsim.beta %>% mutate(sp.inf = hurdle_model$summary.random$sp$mean, ID = "Hurdle model"), mapping = aes(fill = sp.inf)) + 
  geom_sf(data = DFsim.beta %>% mutate(sp.inf = beta_model$summary.random$sp$mean, ID = "Beta model"), mapping = aes(fill = sp.inf)) +
  # scale_fill_viridis_c(option = "mako") +
  sp.mean_scale +
  labs(title = "Spatial effect (inference)", fill = "Mean") +
  facet_wrap(facets = ~ ID, ncol = 2) + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_sp.mean_diff <- ggplot() + 
  geom_sf(data = DFsim.beta %>% mutate(sp.inf = beta_model$summary.random$sp$mean - hurdle_model$summary.random$sp$mean, ID = "Beta - Hurdle"), mapping = aes(fill = sp.inf)) +
  scale_fill_viridis_c(option = "mako") +
  labs(title = "Spatial difference", fill = "Mean") +
  facet_wrap(facets = ~ ID, ncol = 1) + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )

gg_sp.sd <- ggplot() +
  geom_sf(data = DFsim.beta %>% mutate(sp.inf = hurdle_model$summary.random$sp$sd, ID = "Hurdle model"), mapping = aes(fill = sp.inf)) + 
  geom_sf(data = DFsim.beta %>% mutate(sp.inf = beta_model$summary.random$sp$sd, ID = "Beta model"), mapping = aes(fill = sp.inf)) +
  scale_fill_viridis_c(option = "mako") +
  labs(title = "", fill = "Stdev.") +
  facet_wrap(facets = ~ ID, ncol = 2) + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
    )

gg_sp.sd_diff <- ggplot() + 
  geom_sf(data = DFsim.beta %>% mutate(sp.inf = beta_model$summary.random$sp$sd - hurdle_model$summary.random$sp$sd, ID = "Beta - Hurdle"), mapping = aes(fill = sp.inf)) +
  scale_fill_viridis_c(option = "mako") +
  labs(title = "", fill = "Stdev.") +
  facet_wrap(facets = ~ ID, ncol = 1) + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black")
  )
  

gg_beta0 <- ggplot() + 
  geom_line(data = hurdle_model$marginals.fixed$beta0, mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = beta_model$marginals.fixed$beta0, mapping = aes(x = x, y = y), color = "blue") +
  geom_vline(xintercept = beta0, color = "black") +
  labs(title = expression(beta[0])) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

gg_beta1 <- ggplot() + 
  geom_line(data = hurdle_model$marginals.fixed$beta1, mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = beta_model$marginals.fixed$beta1, mapping = aes(x = x, y = y), color = "blue") +
  geom_vline(xintercept = beta1, color = "black") +
  labs(title = expression(beta[1])) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

gg_tau <- ggplot() + 
  geom_line(data = hurdle_model$marginals.hyperpar$`Precision for sp`, mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = beta_model$marginals.hyperpar$`Precision for sp`, mapping = aes(x = x, y = y), color = "blue") +
  geom_vline(xintercept = 1, color = "black") +
  labs(title = expression(tau)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

gg_rho <- ggplot() + 
  geom_line(data = hurdle_model$marginals.hyperpar$`Lambda for sp`, mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = beta_model$marginals.hyperpar$`Lambda for sp`, mapping = aes(x = x, y = y), color = "blue") +
  geom_vline(xintercept = 0.7, color = "black") +
  labs(title = expression(lambda)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

# gg_tot <- grid.arrange(arrangeGrob(grobs = list(gg_sp.sim, gg_ber.sim, gg_sp.mean, gg_sp.sd, gg_beta0, gg_beta1, gg_tau, gg_rho), layout_matrix = matrix(data = c(1,1,3,3,2,2,4,4,5,6,7,8), nrow = 3, byrow = TRUE))) # Old Figure

gg_tot_space <- arrangeGrob(grobs = list(gg_sp.mean, gg_sp.sd, gg_sp.mean_diff, gg_sp.sd_diff), layout_matrix = matrix(data = c(1,3,2,4), nrow = 2, byrow = TRUE))
gg_tot <- grid.arrange(arrangeGrob(grobs = list(gg_sp.sim, gg_ber.sim, gg_tot_space, gg_beta0, gg_beta1, gg_tau, gg_rho), layout_matrix = matrix(data = c(1,3,3,3,2,3,3,3,4,5,6,7), nrow = 3, byrow = TRUE)))

ggsave(plot = gg_tot, filename = "./Figures/fig_BetaZeros_example_hq.eps", device = "eps", width = 1535, height = 1037, units = "px", dpi = 100) # High quality plot
ggsave(plot = gg_tot, filename = "./Figures/fig_BetaZeros_example_lq.jpg", device = "jpg", width = 1535, height = 1037, units = "px", dpi = 100) # Low quality plot
