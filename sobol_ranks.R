## ----setup, include=FALSE-----------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----load_packages, results="hide", message=FALSE, warning=FALSE--------

# LOAD PACKAGES ---------------------------------------------------------------

# Function to read in all required packages in one go:
loadPackages <- function(x) {
  for(i in x) {
    if(!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}

loadPackages(c("sensobol", "tidyverse", "parallel", "grid", 
               "cowplot", "gridExtra", "sensitivity", "wesanderson", "RColorBrewer", 
               "tikzDevice", "data.table"))

# SET CHECKPOINT --------------------------------------------------------------

dir.create(".checkpoint")

library("checkpoint")

checkpoint("2020-07-02", 
           R.version ="3.6.3", 
           checkpointLocation = getwd())

# CUSTOM FUNCTION TO DEFINE THE PLOT THEMES -----------------------------------

theme_AP <- function() {
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent",
                                         color = NA),
        legend.key = element_rect(fill = "transparent",
                                  color = NA))
}



## ----define_functions, cache=TRUE---------------------------------------

# TEST FUNCTIONS ---------------------------------------------------------------

# Sobol G' function
G_Fun <- function(X) {
  a <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)
  y <- 1
  for (j in 1:8) {
    y <- y * (abs(4 * X[, j] - 2) + a[j])/(1 + a[j])
  }
  return(y)
}

# Special G' function
G_Fun_mod <- function(X) {
  a <- c(0, 0, 0, 0, 0, 0, 0, 0)
  y <- 1
  for (j in 1:8) {
    y <- y * (abs(4 * X[, j] - 2) + a[j])/(1 + a[j])
  }
  return(y)
}

# B1 function (in Kucherenko et al. 2011)
B1 <- function(X) {
  y <- 1
  for(j in 1:ncol(X)) {
    y <- y * (ncol(X) - X[, j]) / (ncol(X) - 0.5)
  }
  return(y)
}

# C1 function (in Kucherenko et al. 2011)
C1 <- function(X) {
  y <- 1
  for (j in 1:ncol(X)) {
    y <- y * (abs(4 * X[, j] - 2))}
  return(y)
}

# The function by Bratley et al. (1992) is defined in the sensobol package
# (Function A1 in Kucherenko et al. 2011)

# The function by Oakley and O'Hagan is defined in the sensobol package

# The function by Morris 1991 is included in the sensitivity package


## ----sample_matrices, cache=TRUE----------------------------------------

# SAMPLE MATRICES ------------------------------------------------------------

CutBySize <- function(m, block.size, nb = ceiling(m / block.size)) {
  int <- m / nb
  upper <- round(1:nb * int)
  lower <- c(1, upper[-nb] + 1)
  size <- c(upper[1], diff(upper))
  cbind(lower, upper)
}

scrambled_sobol <- function(A, B) {
  X <- rbind(A, B)
  for(i in 1:ncol(A)) {
    AB <- A
    AB[, i] <- B[, i]
    X <- rbind(X, AB)
  }
  AB <- X
  return(AB)
}

scrambled_replicas <- function(N, k, version) {
  X <- A <- B <- out <- list()
  df <- randtoolbox::sobol(n = N * version, dim = k * 2)
  indices <- CutBySize(nrow(df), nb = version)
  for(i in 1:nrow(indices)) {
    lower <- indices[i, "lower"]
    upper <- indices[i, "upper"]
    X[[i]] <- df[lower:upper, ]
  }
  for(i in seq_along(X)) {
    A[[i]] <- X[[i]][, 1:k]
    B[[i]] <- X[[i]][, (k + 1) : (k * 2)]
  }
  for(i in seq_along(A)) {
    out[[i]] <- scrambled_sobol(A[[i]], B[[i]])
  }
  return(out)
}


## ----sobol_function, cache=TRUE-----------------------------------------

# SOBOL' INDICES --------------------------------------------------------------

sobol_computeS <- function(Y_A, Y_B, Y_AB) {
  n <- length(Y_A[!is.na(Y_A)])
  f0 <- (1 / (2 * n)) * sum(Y_A + Y_B, na.rm = TRUE)
  VY <- 1 / n * sum((Y_A - f0) ^ 2, na.rm = TRUE)
  Si <- (1 / n) * sum(Y_B * (Y_AB - Y_A), na.rm = TRUE) / VY
  STi <- ((1 / (2 * n)) * sum((Y_A - Y_AB) ^ 2, na.rm = TRUE)) / VY
  return(c(Si, STi))
}

sobol_MapplyS <- function(d) {
  return(mapply(sobol_computeS,
                d[, "Y_A"],
                d[, "Y_B"],
                d[, "Y_AB"]))
}

sobol_compute <- function(Y, params, boot = FALSE, ranks = FALSE) {
  k <- length(params)
  p <- length(1:(length(Y) / (k + 2)))
  Y_A <- Y[1:p]
  Y_B <- Y[(p + 1) : (2 * p)]
  Y_AB <- Y[(2 * p + 1):((length(Y) / (k + 2)) * (k + 2))]
  parameters <- rep(params, each = length(Y_A))
  vec <- cbind(Y_A, Y_B, Y_AB)
  out <- data.table::data.table(vec, parameters)
  if(boot == TRUE) {
    indices <- out[, sample(.I, replace = TRUE), parameters][, V1]
    output <- out[indices, sobol_MapplyS(.SD), by = parameters][
      , sensitivity:= rep(c("Si", "STi"), times = k)]
  }
  if(boot == FALSE) {
    output <- out[, sobol_MapplyS(.SD), by = parameters][
      , sensitivity:= rep(c("Si", "STi"), times = k)]
  }
  if(ranks == FALSE) {
    final <- output[, V1, sensitivity][, V1]
  } 
  if(ranks == TRUE) {
    final <- output[, rank(-V1), sensitivity][, V1]
  }
  return(final)
}


## ----add_diagram, echo=FALSE, fig.align="center", fig.cap="Diagram of the new algorithm", out.width = '100%'----

# INSERT THE ALGORITHM DIAGRAM ------------------------------------------------

knitr::include_graphics("./sobol_ranks_diagram.pdf")


## ----define_algorithm, cache=TRUE, dependson=c("sobol_function", "sample_matrices", "define_functions")----

# DEFINE ALGORITHM -------------------------------------------------------------

# Compute normal
compute_normal <- function(N, params, test_F) {
  A <- sobol_matrices(n = N, k = length(params))
  Y <- test_F(A)
  indices <- sobol_compute(Y, params = params)
  ranks <- sobol_compute(Y, params = params, ranks = TRUE)
  dt <- data.table(cbind(indices, ranks))
  dt[, parameters:= rep(params, times = 2)][
    , sensitivity:= rep(c("Si", "STi"), each = length(params))][
    , method:= "Normal.approach"][
    , N:= N][
    , model.runs:= N * (length(params) + 2)]
  return(dt)
}

# Algorithm to save runs
compute_saving <- function(dt, params, N, test_F, eps = 0.05) {
  dt2 <- copy(dt)
  set <- dt[sensitivity == "STi" & indices <= 0.05][, parameters]
  if(length(set) > 1) {
    Ti.set <- match(set, params)
    Ti.important <- match(setdiff(params, set), params)
    A.important <- sobol_matrices(N, k = length(params), cluster = Ti.important)
    A.set <- sobol_matrices(N, k = length(params), cluster = list(Ti.set))
    A.full <- rbind(A.important, A.set[(2 * N + 1):nrow(A.set), ])
    Y <- test_F(A.full)
    indices <- sobol_compute(Y, params = c(params[Ti.important], "set"))
    ranks <- sobol_compute(Y, params = c(params[Ti.important], "set"), ranks = TRUE)
    dt <- data.table(cbind(indices, ranks))
    dt[, parameters:= rep(c(params[Ti.important], "set"), times = 2)][
      , sensitivity:= rep(c("Si", "STi"), each = length(c(params[Ti.important], "set")))][
      , method:= "New.approach"][
      , N:= N][
      , model.runs:= N * (length(params) - length(set) + 3)]
    final <- rbind(dt2, dt)
  } else {
    final <- dt
  }
  return(final)
}

# Full algorithm
full_algorithm <- function(max.exponent, params, epsilon_1, epsilon_2, test_F) {
  N <- sapply(2:max.exponent, function(x) 2 ^ x)
  dt <- list()
  for(i in seq_along(N)) {
    dt[[i]] <- compute_normal(N[i], params = params, test_F = test_F)
    if(i > 1) {
      a <- dt[[i]][sensitivity == "STi" & method == "Normal.approach"]
      b <- dt[[i-1]][sensitivity == "STi" & method == "Normal.approach"]
      # CHECK CONDITION 1---------------------------------
      condition1 <- identical(a[, ranks], b[, ranks])
      # CHECK CONDITION 2---------------------------------
      ind <- which(!a[, ranks] == b[, ranks])
      condition2 <- condition2 <- all(abs(a[ind][, indices] - b[ind][, indices]) / 
                          (a[ind][, indices] + b[ind][, indices]) < epsilon_1)
      # CHECK CONDITION 3 ---------------------------------
      condition3 <- all(c(a[ind][, indices], b[ind][, indices]) < epsilon_2)
      if(condition1 == TRUE | 
         condition1 == FALSE & condition2 == TRUE |
         condition1 == FALSE & condition2 == FALSE & condition3 == TRUE) {
        dt[[i]] <- compute_saving(dt[[i]], params = params, N = N[i], test_F = test_F)
      } else {
        dt[[i]] <- compute_normal(N[i], params = params, test_F = test_F)
      }
    }
  }
  # Compute total number of model runs for the new algorithm
  final <- rbindlist(dt)
  N_convergence <- final[method == "New.approach" & sensitivity == "STi", min(N)]
  a <- final[N == N_convergence & method == "Normal.approach" 
             & sensitivity == "STi"][, parameters]
  b <- final[N == N_convergence & !parameters == "set" & 
               method == "New.approach" & sensitivity == "STi"][, parameters]
  k_noninfluential <- length(setdiff(a, b))
  final[, model.runs:= ifelse(method %in% "New.approach", 
                              model.runs + N_convergence * k_noninfluential, 
                              model.runs)]
  return(final)
}

# Full algorithm in a function
fun_algorithm <- function(max.exponent, epsilon_1, epsilon_2) {
  test_functions <- c("G_Fun","G_Fun_mod","bratley1992_Fun", 
                      "oakley_Fun", "morris_Fun", "B1", "C1")
  out <- list()
  for(i in test_functions) {
    if(i == "G_Fun") {
      params <- paste("X", 1:8, sep = "")
      test_F <- G_Fun
    } else if(i == "G_Fun_mod") {
      params <- paste("X", 1:8, sep = "")
      test_F <- G_Fun_mod
    } else if(i == "bratley1992_Fun") {
      params <- paste("X", 1:8, sep = "")
      test_F <- bratley1992_Fun
    } else if(i == "oakley_Fun") {
      params <- paste("X", 1:15, sep = "")
      test_F <- oakley_Fun
    } else if(i == "B1") {
      params <- paste("X", 1:8, sep = "")
      test_F <- B1
    } else if(i == "C1") {
      params <- paste("X", 1:8, sep = "")
      test_F <- C1
    } else {
      params <- paste("X", 1:20, sep = "")
      test_F <- sensitivity::morris.fun
    }
    out[[i]] <- full_algorithm(max.exponent = max.exponent, 
                               params = params, 
                               epsilon_1 = epsilon_1,
                               epsilon_2 = epsilon_2,
                               test_F = test_F)
  }  
  final <- rbindlist(out, idcol = "model")
  return(final)
}


## ----run_saving_algorithm, cache=TRUE, dependson="define_algorithm"-----

# RUN THE MODEL --------------------------------------------------------------

# Define epsilons and max exponent
epsilon_1 <- seq(0.05, 0.1, 0.05)
epsilon_2 <- seq(0.01, 0.02, 0.005)
max.exponent <- 16

# Run model
out <- mclapply(epsilon_1, function(x) 
  lapply(epsilon_2, function(y) 
    fun_algorithm(max.exponent = max.exponent, epsilon_1 = x, epsilon_2 = y)), 
  mc.cores = detectCores())


## ----arrange_results_algorithm, cache=TRUE, dependson="run_saving_algorithm"----

# ARRANGE RESULTS -------------------------------------------------------------

test_functions <- c("G_Fun","G_Fun_mod","bratley1992_Fun", 
                    "oakley_Fun", "morris_Fun", "B1", "C1")

# Arrange results
model_names <- c("Sobol' G", "Sobol'G modified", "Bratley et al. 1994", 
                 "Oakley and O'Hagan 2004", "Morris 1991", "B1", "C1")

names(out) <- epsilon_1
for(i in names(out)) {
  names(out[[i]]) <- epsilon_2
}

results <- lapply(out, function(x) lapply(x, function(y) {
  tmp <- y[, parameters:= factor(parameters, 
                            levels = c(paste("X", 1:20, sep = ""), "set"))] %>%
    .[, Type:= ifelse(parameters %in% "set", "Set", "Individual")] %>%
    .[, model:= ifelse(model == "G_Fun", "Sobol G", 
                       ifelse(model == "G_Fun_mod", "Sobol G modified", 
                              ifelse(model == "bratley1992_Fun", "Bratley et al. 1994", 
                                     ifelse(model == "oakley_Fun", "Oakley and O'Hagan 2004", 
                                            ifelse(model == "morris_Fun", "Morris 1991", 
                                                   ifelse(model == "B1", "B1", "C1"))))))]
}))

# Compare number of model runs and savings
tmp <- lapply(out, function(x) lapply(x, function(y) {
  tmp <- y[, .(model.runs = unique(model.runs)), .(model, N, method)] %>%
    dcast(., model + N ~ method, value.var = "model.runs") 
}))

savings.dt <- lapply(tmp, function(x) lapply(x, function(y) {
  tmp <- y[, New.approach:= ifelse(New.approach %in% NA, Normal.approach, New.approach)] %>%
    .[, saving:= 1 - (New.approach / Normal.approach)] %>%
    .[, model:= ifelse(model == "G_Fun", "Sobol G", 
                       ifelse(model == "G_Fun_mod", "Sobol G modified", 
                              ifelse(model == "bratley1992_Fun", "Bratley et al. 1994", 
                                     ifelse(model == "oakley_Fun", "Oakley and O'Hagan 2004", 
                                            ifelse(model == "morris_Fun", "Morris 1991", 
                                                   ifelse(model == "B1", "B1", "C1"))))))]
}))


## ----plot_results_algorithm, cache=TRUE, dependson="arrange_results_algorithm", dev = c("tikz", "jpeg"), dpi=300, fig.height=4, fig.width=5, fig.cap="Evolution of the $T_i$ and the $T_s$ along different sample sizes. The horizontal, red dotted line is at 0.05."----

# PLOT RESULTS ---------------------------------------------------------------

lapply(results, function(x) lapply(x, function(y) 
  y[, Type:= ifelse(parameters %in% "set", "Set", "Individual")] %>%
    .[!sensitivity == "Si" & method == "New.approach"] %>%
    ggplot(., aes(N, indices, 
                  color = Type, 
                  group = parameters)) +
    geom_line() +
    scale_x_continuous(trans="log",
                       breaks = scales::trans_breaks("log2", function(x) 2 ^ x),
                       labels = scales::trans_format("log2", scales::math_format(2^.x))) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("black", "blue"), 
                       labels = c(expression(T[italic(i)]), 
                                  expression(T[italic(s)]))) +
    labs(x = "N", 
         y = "Variance") +
    geom_hline(yintercept = 0.05, 
               lty = 2, 
               color = "red") +
    facet_wrap(~model, 
               ncol = 2, 
               scales = "free") +
    theme_AP() + 
    theme(legend.position = "top"))) 


## ----plot_savings, cache=TRUE, dependson="arrange_results_algorithm", dev = c("tikz", "jpeg"), dpi=300, fig.height=4, fig.width=6.5, fig.cap="Percentage of savings over the total accumulated model runs."----

# PLOT SAVINGS ----------------------------------------------------------------

lapply(savings.dt, function(x) lapply(x, function(y) 
  melt(y, measure.vars = c("New.approach", "Normal.approach")) %>%
    .[, saving:= saving * 100] %>%
    setnames(., "model", "Model")) %>%
    rbindlist(., idcol = "epsilon_2")) %>%
  rbindlist(., idcol = "epsilon_1") %>%
  .[, epsilon_1:= paste("$\\varepsilon_s= ", epsilon_1, "$", sep = "")] %>%
  .[, epsilon_2:= paste("$\\varepsilon_b= ", epsilon_2, "$", sep = "")] %>%
  ggplot(., aes(N, saving, color = Model, shape = Model, group = Model)) +
  geom_point() + 
  geom_line() +
  scale_x_continuous(trans="log",
                     breaks = trans_breaks("log2", function(x) 2 ^ x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(x = "Number of model runs", 
       y = "Saving (\\%)") +
  facet_grid(epsilon_1 ~ epsilon_2) +
  theme_AP() + 
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))


## ----plot_savings_functions2, cache=TRUE, dependson="arrange_results_algorithm", dev = c("tikz", "jpeg"), dpi=300, fig.height=4, fig.width=6.5, fig.cap="Percentage of savings over the total accumulated model runs."----

# PLOT SAVINGS ------------------------------------------------------------------

lapply(savings.dt, function(x) lapply(x, function(y) 
  melt(y, measure.vars = c("New.approach", "Normal.approach")) %>%
    .[, saving:= saving * 100] %>%
    setnames(., "model", "Model")) %>%
    rbindlist(., idcol = "epsilon_2")) %>%
  rbindlist(., idcol = "epsilon_1") %>%
  .[, epsilon_1:= paste("$\\varepsilon_s= ", epsilon_1, "$", sep = "")] %>%
  .[, epsilon_2:= paste("$\\varepsilon_b= ", epsilon_2, "$", sep = "")] %>%
  .[!Model %in% c("C1", "B1", "Sobol G modified")] %>%
  ggplot(., aes(N, saving, color = Model, shape = Model, group = Model)) +
  geom_point() + 
  geom_line() +
  scale_x_continuous(trans="log",
                     breaks = trans_breaks("log2", function(x) 2 ^ x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(x = "Number of model runs", 
       y = "Saving (\\%)") +
  facet_grid(epsilon_1 ~ epsilon_2) +
  theme_AP() + 
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))


## ----plot_ranks, cache=TRUE, dependson="arrange_results_algorithm", dev = c("tikz", "jpeg"), dpi=300, fig.height=8.5, fig.width=4, fig.cap="Evolution of ranks for total-order indices across different base sample sizes"----

# PLOT RANKS ------------------------------------------------------------------

results[[2]][[3]][method == "Normal.approach"] %>%
  .[, sensitivity:= ifelse(sensitivity %in% "Si", "$S_i$", "$T_i$")] %>%
  .[sensitivity == "$T_i$"] %>%
  ggplot(., aes(N, ranks, group = parameters,
                color = parameters)) +
  geom_point(size = 0.5) +
  geom_line() +
  labs(x = "Base sample size",
       y = "Rank") +
  scale_color_discrete(name = "Parameters") +
  scale_x_continuous(trans="log",
                     breaks = trans_breaks("log2", function(x) 2 ^ x),
                     labels = trans_format("log2", math_format(2^.x))) +
  scale_y_reverse() +
  facet_grid(model ~ sensitivity,
             scales = "free_y") +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent",
                                         color = NA),
        legend.key = element_rect(fill = "transparent",
                                  color = NA), 
        strip.text.y = element_text(size = 7))


## ----plot_ranks_functions, cache=TRUE, dependson="arrange_results_algorithm", dev = c("tikz", "jpeg"), dpi=300, fig.height=8, fig.width=4, fig.cap="Evolution of ranks for total-order indices across different base sample sizes."----

# PLOT RANKS ------------------------------------------------------------------

results[[2]][[3]][method == "Normal.approach"] %>%
  .[, sensitivity:= ifelse(sensitivity %in% "Si", "$S_i$", "$T_i$")] %>%
  .[sensitivity == "$T_i$"] %>%
  .[!model %in% c("C1", "B1", "Sobol G modified")] %>%
  ggplot(., aes(N, ranks, group = parameters,
                color = parameters)) +
  geom_point(size = 0.5) +
  geom_line() +
  labs(x = "Base sample size",
       y = "Rank") +
  scale_color_discrete(name = "Parameters") +
  scale_x_continuous(trans="log",
                     breaks = trans_breaks("log2", function(x) 2 ^ x),
                     labels = trans_format("log2", math_format(2^.x))) +
  scale_y_reverse() +
  facet_grid(model ~ sensitivity,
             scales = "free_y") +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent",
                                         color = NA),
        legend.key = element_rect(fill = "transparent",
                                  color = NA))


## ----session_information------------------------------------------------

# SESSION INFORMATION ---------------------------------------------------------

sessionInfo()


