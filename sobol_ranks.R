## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----load_packages, results="hide", message=FALSE, warning=FALSE---------

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

loadPackages(c("data.table", "ggplot2", "sensobol", "scales", "parallel", "grid", 
               "cowplot", "gridExtra", "sensitivity", "wesanderson", "RColorBrewer", 
               "tikzDevice"))

# SET CHECKPOINT --------------------------------------------------------------

#dir.create(".checkpoint")

#library("checkpoint")

#checkpoint("2019-08-28", 
           #R.version ="3.6.1", 
           #checkpointLocation = getwd())

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


## ----define_functions, cache=TRUE----------------------------------------

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

# The function by Bratley et al. (1992) is defined in the sensobol package
# (Function A1 in Kucherenko et al. 2011)

# The function by Oakley and O'Hagan is defined in the sensobol package

# The function by Morris 1991 is included in the sensitivity package


## ----sample_matrices, cache=TRUE-----------------------------------------

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


## ----sobol_function, cache=TRUE------------------------------------------

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
    , method:= "Normal approach"][
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
      , method:= "New approach"][
      , N:= N][
      , model.runs:= N * (length(params) - length(set) + 3)]
    final <- rbind(dt2, dt)
  } else {
    final <- dt
  }
  return(final)
}

# Full algorithm
full_algorithm <- function(max.exponent, params, test_F) {
  N <- sapply(2:max.exponent, function(x) 2 ^ x)
  dt <- list()
  for(i in seq_along(N)) {
    dt[[i]] <- compute_normal(N[i], params = params, test_F = test_F)
    if(i > 1) {
      a <- dt[[i]][sensitivity == "STi" & method == "Normal approach"]
      b <- dt[[i-1]][sensitivity == "STi" & method == "Normal approach"]
      # CHECK CONDITION 1---------------------------------
      condition1 <- identical(a[, ranks], b[, ranks])
      # CHECK CONDITION 2---------------------------------
      ind <- which(!a[, ranks] == b[, ranks])
      condition2 <- condition2 <- all(abs(a[ind][, indices] - b[ind][, indices]) / 
                          (a[ind][, indices] + b[ind][, indices]) < 0.1)
      # CHECK CONDITION 3 ---------------------------------
      condition3 <- all(c(a[ind][, indices], b[ind][, indices]) < 0.02)
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
  N_convergence <- final[method == "New approach" & sensitivity == "STi", min(N)]
  a <- final[N == N_convergence & method == "Normal approach" 
             & sensitivity == "STi"][, parameters]
  b <- final[N == N_convergence & !parameters == "set" & 
               method == "New approach" & sensitivity == "STi"][, parameters]
  k_noninfluential <- length(setdiff(a, b))
  final[, model.runs:= ifelse(method %in% "New approach", 
                              model.runs + N_convergence * k_noninfluential, 
                              model.runs)]
  return(final)
}


## ----run_saving_algorithm, cache=TRUE, dependson="define_algorithm"------

# RUN THE MODEL --------------------------------------------------------------

max.exponent <- 16
test_functions <- c("G_Fun", "bratley1992_Fun", "oakley_Fun", "morris_Fun")
out <- list()

for(i in test_functions) {
  if(i == "G_Fun") {
    params <- paste("X", 1:8, sep = "")
    test_F <- G_Fun
  } else if(i == "bratley1992_Fun") {
    params <- paste("X", 1:8, sep = "")
    test_F <- bratley1992_Fun
  } else if(i == "oakley_Fun") {
    params <- paste("X", 1:15, sep = "")
    test_F <- oakley_Fun
  } else {
    params <- paste("X", 1:20, sep = "")
    test_F <- sensitivity::morris.fun
  }
  out[[i]] <- full_algorithm(max.exponent = max.exponent, 
                             params = params, 
                             test_F = test_F)
}  


## ----arrange_results_algorithm, cache=TRUE, dependson="run_saving_algorithm"----

# ARRANGE RESULTS -------------------------------------------------------------

# Arrange results
model_names <- c("Sobol' G", "Bratley et al. 1994", 
                 "Oakley and O'Hagan 2004", "Morris 1991")
names(out) <- model_names
results <- rbindlist(out, idcol = "model") %>%
  .[, model:= factor(model, levels = model_names)] %>%
  .[, parameters:= factor(parameters, 
                          levels = c(paste("X", 1:20, sep = ""), "set"))] %>%
  .[, method:= ifelse(method %in% "New approach", "New.approach", "Normal.approach")]

# Compare number of model runs and savings
tmp <- results[, .(model.runs = unique(model.runs)), .(model, N, method)]
savings.dt <- dcast(tmp, model + N ~ method, value.var = "model.runs") %>%
  .[, New.approach:= ifelse(New.approach %in% NA, Normal.approach, New.approach)] %>%
  .[, saving:= 1 - (New.approach / Normal.approach)]

# Compute cumulative number of runs
cumulative.runs <- savings.dt[, cumsum(.SD), 
                              .SDcols = c("New.approach", "Normal.approach"), model][
                                , saving:= 1 - (New.approach / Normal.approach)]


## ----export_results, cache=TRUE, dependson="arrange_results_algorithm"----

# EXPORT MODEL RUNS -----------------------------------------------------------

fwrite(results, "results.csv")
fwrite(savings.dt, "savings.dt.csv")


## ----plot_results_algorithm, cache=TRUE, dependson="arrange_results_algorithm", dev = "tikz", fig.height=4, fig.width=5, fig.cap="Evolution of the $T_i$ and the $T_s$ along different sample sizes. The horizontal, red dotted line is at 0.05."----

# PLOT RESULTS ---------------------------------------------------------------

results <- results[, Type:= ifelse(parameters %in% "set", "Set", "Individual")]

formatter <- function(...) function(x) format(round(x, 2), ...)

results[!sensitivity == "Si" & method == "New.approach"] %>%
  ggplot(., aes(N, indices, 
                color = Type, 
                group = parameters)) +
  geom_line() +
  scale_x_continuous(trans="log",
                     breaks = trans_breaks("log2", function(x) 2 ^ x),
                     labels = trans_format("log2", math_format(2^.x))) +
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
  theme(legend.position = "top") 


## ----plot_savings, cache=TRUE, dependson="arrange_results_algorithm", dev = "tikz", fig.height=3, fig.width=6.5, fig.cap="Evolution of ranks for first and total-order indices across different base sample sizes"----

# PLOT SAVINGS ----------------------------------------------------------------

melt(savings.dt, measure.vars = c("New.approach", "Normal.approach")) %>%
  .[, saving:= saving * 100] %>%
  setnames(., "model", "Model") %>%
  ggplot(., aes(N, saving, color = Model, shape = Model, group = Model)) +
  geom_point() + 
  geom_line() +
  scale_x_continuous(trans="log",
                     breaks = trans_breaks("log2", function(x) 2 ^ x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(x = "Number of model runs", 
       y = "Saving (\\%)") +
  theme_AP() + 
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 1,byrow = TRUE))


## ----plot_ranks, cache=TRUE, dependson="arrange_results_algorithm", dev = "tikz", fig.height=8, fig.width=4, fig.cap="Percentage of saving over the total accumulated model runs."----

# PLOT RANKS ------------------------------------------------------------------

results[method == "Normal.approach"] %>%
  .[, sensitivity:= ifelse(sensitivity %in% "Si", "$S_i$", "$T_i$")] %>%
  ggplot(., aes(N, ranks, group = parameters,
                color = parameters)) +
  geom_point(size = 0.5) +
  geom_line() +
  labs(x = "Base sample size",
       y = "Rank") +
  scale_color_discrete(name = "Parameters") +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
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


## ----MAE, cache=TRUE, dependson="arrange_results_algorithm", fig.height = 6.5, fig.width = 3, dev="tikz", fig.cap="Mean Absolute Error for the Oakley and O'Hagan function."----

# COMPUTATION OF MEAN ABSOLUTE ERROR ------------------------------------------

# The analytical values for the Oakley and O'Hagan function
dt.analytics <- data.table(parameters = paste("X", 1:15, sep = ""), 
                           analytical = c(0.059, 0.063, 0.036, 0.055, 0.024, 0.041, 
                                          0.058, 0.082, 0.097, 0.036, 0.151, 0.148, 
                                          0.142, 0.141, 0.155))

# Select rows
AE.test <- results[model == "Oakley and O'Hagan 2004" &
                     sensitivity == "STi" & 
                     N >= 512]

# Compute the Absolute Error (AE) for the non-clustered parameters
# Compute the Absolute Error (AE) for the non-clustered parameters
AE.ns <- merge(AE.test, dt.analytics, by = "parameters", all.x = TRUE) %>%
  .[, AE_ns:= abs(indices - analytical)]

# Extract vector with the clustered parameters
pars <- AE.ns[method =="Normal.approach" &
                indices < 0.05][, unique(parameters)]

# Compute the sum of the analytical values of the clustered parameters
tmp <- AE.ns[parameters %in% pars][, .(tmp = sum(indices)), .(N, model.runs)]

# Compute the AE for the set of clustered parameters (AE.s)
AE.s <- AE.ns[parameters == "set"] %>%
  .[tmp, on = "N"] %>%
  .[, .(AE_s = abs(indices - tmp)), .(N, model.runs)]

# Compute MAE, MAE_s, and MAE average
full.AE <- AE.ns[AE.s, on = "N"] %>%
  .[, .(MAE_S = (sum(AE_ns, na.rm = TRUE) + AE_s) / 15, 
        MAE = mean(abs(indices - analytical), na.rm = TRUE)), .(N, model.runs, method)] %>%
  .[, average.MAE:= (MAE_S + MAE) / 15]

# Plot results
full.AE <- setnames(full.AE, 
                    c("MAE_S", "MAE", "average.MAE"), 
                    c("$T_s$", "$T_i$", "$(T_i + T_s) / k$"))

melt(full.AE, measure.vars = c("$T_s$", "$T_i$", "$(T_i + T_s) / k$")) %>%
  ggplot(., aes(model.runs, value, color = method)) +
  geom_point() +
  geom_line() +
  labs(x = "Number of model runs", 
       y = "MAE") +
  scale_x_continuous(trans="log",
                     breaks = trans_breaks("log2", function(x) 2 ^ x),
                     labels = trans_format("log2", math_format(2^.x))) +
  scale_color_discrete(name = "Method", 
                       labels = c("New approach", "Traditional approach")) +
  facet_wrap(~ variable, 
             ncol = 1,
             scales = "free_y") +
  theme_AP() +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


## ----session_information-------------------------------------------------

# SESSION INFORMATION ---------------------------------------------------------

sessionInfo()


