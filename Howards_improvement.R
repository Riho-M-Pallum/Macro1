###############################################################################
## Howard's Improvement for the Neoclassical Growth Model in R


howard_vfi_track <- function(alpha = 0.36,
                             beta  = 0.95,
                             delta = 0.05,
                             theta = 4,
                             n_k   = 20,
                             tol_value  = 1e-8,
                             tol_policy = 1e-6,
                             max_iter_value  = 10000,
                             max_iter_policy = 1000) {
  u_fun <- function(c) {
    if (theta == 1) return(log(c))
    (c^(1 - theta) - 1) / (1 - theta)
  }
  f_fun <- function(k) k^alpha
  
  # theoretical steady state (from Euler equation)
  k_ss_theory <- (alpha / (1 / beta - 1 + delta))^(1 / (1 - alpha))
  
  # grid and resources
  k_grid <- seq(0.1 * k_ss_theory, 2 * k_ss_theory, length.out = n_k)
  n <- n_k
  fK        <- f_fun(k_grid)
  resources <- fK + (1 - delta) * k_grid
  
  # precompute one-period utility
  util <- matrix(-1e12, nrow = n, ncol = n)
  for (i in 1:n) {
    c_ij <- resources[i] - k_grid
    feasible <- c_ij > 0
    u_ij <- rep(-1e12, n)
    u_ij[feasible] <- u_fun(c_ij[feasible])
    util[i, ] <- u_ij
  }
  
  V <- rep(0, n)
  policy_idx <- rep(1, n)
  
  # to store approximate steady state at each outer iteration
  k_ss_path <- numeric(max_iter_policy)
  it_used <- 0
  
  for (it_pol in 1:max_iter_policy) {
    # inner value iteration with fixed policy
    for (it_val in 1:max_iter_value) {
      V_new <- util[cbind(1:n, policy_idx)] + beta * V[policy_idx]
      diff_V <- max(abs(V_new - V))
      V <- V_new
      if (diff_V < tol_value) break
    }
    
    # policy improvement
    value_candidate <- util + beta * matrix(rep(V, each = n), nrow = n)
    new_policy_idx <- max.col(value_candidate, ties.method = "first")
    
    policy_diff <- max(abs(k_grid[new_policy_idx] - k_grid[policy_idx]))
    policy_idx <- new_policy_idx
    
    # implied policy in levels
    k_policy <- k_grid[policy_idx]
    
    # approximate steady state on the grid: fixed point of g(k) ~ k
    idx_fix <- which.min(abs(k_policy - k_grid))
    k_ss_path[it_pol] <- k_grid[idx_fix]
    
    it_used <- it_pol
    if (policy_diff < tol_policy) break
  }
  
  k_policy <- k_grid[policy_idx]
  c_policy <- resources - k_policy
  
  list(
    k_grid      = k_grid,
    V           = V,
    policy_idx  = policy_idx,
    k_policy    = k_policy,
    c_policy    = c_policy,
    k_ss_theory = k_ss_theory,
    k_ss_path   = k_ss_path[1:it_used]
  )
}


###############################################################################
## Sanity checks with full depreciation and log utility
###############################################################################



# Testing calibration from the problem set:
alpha <- 0.36
beta  <- 0.95
delta <- 1
theta <- 1

res <- howard_vfi_track(alpha = alpha, beta = beta, delta = delta,
                        theta = theta, n_k = 20)

# error per outer iteration
ss_diff <- abs(res$k_ss_path - res$k_ss_theory)
plot(seq_along(ss_diff), ss_diff, type = "b", log = "y",
     xlab = "Howard policy iteration",
     ylab = "|k_ss^approx - k_ss^theory| (log scale)",
     main = "Convergence of Howard algorithm (log error)")


# Bah Ok we converge after one iteration in the testing set this works really well


# Quick sanity checks:
plot(res$k_grid, res$k_policy, type = "l",
     xlab = "k", ylab = "k'(k)", main = "Policy function with Howard's improvement")

###############################################################################
# Now to check the Euler equation on the gridpoints 
################################################################################

# Now to check the Euler equation on the gridpoints 

u_prime <- function(c, theta) {
  if (theta == 1) return(1 / c)
  c^(-theta)
}

# Marginal product of capital
f_prime <- function(k, alpha) {
  alpha * k^(alpha - 1)
}

# Simple linear interpolation of a policy y(x) defined on (x_grid, y_grid)
interp_linear <- function(x, x_grid, y_grid) {
  approx(x = x_grid, y = y_grid, xout = x, rule = 2)$y
}

# Compute Euler equation residuals at midpoints of the capital grid
euler_errors_midpoints <- function(res, alpha, beta, delta, theta) {
  k_grid    <- res$k_grid
  k_policy  <- res$k_policy
  
  n <- length(k_grid)
  
  # Midpoints: Ge = 0.5*(k_i + k_{i+1}), i = 1,...,n-1
  k_mid <- 0.5 * (k_grid[-n] + k_grid[-1])
  
  # Policy at midpoints: k' = g(k_mid)
  k_next_mid <- interp_linear(k_mid, k_grid, k_policy)
  
  # Current consumption at midpoints
  f_k_mid <- k_mid^alpha
  c_mid   <- f_k_mid + (1 - delta) * k_mid - k_next_mid
  
  # Next-period consumption: need k'' = g(k')
  k_next2_mid <- interp_linear(k_next_mid, k_grid, k_policy)
  f_k_next    <- k_next_mid^alpha
  c_next_mid  <- f_k_next + (1 - delta) * k_next_mid - k_next2_mid
  
  # Only keep points with positive c and c_next to avoid NaNs
  good <- (c_mid > 0) & (c_next_mid > 0)
  
  # Euler equation: u'(c) = beta * u'(c') * [f'(k') + 1 - delta]
  lhs <- u_prime(c_mid[good], theta)
  rhs <- beta * u_prime(c_next_mid[good], theta) *
    (f_prime(k_next_mid[good], alpha) + 1 - delta)
  
  # Residual in levels
  ee_resid <- lhs - rhs
  
  # Often people look at log10 of the *relative* error in consumption units:
  # 1 - rhs/lhs  (only where lhs != 0)
  rel_err <- 1 - rhs / lhs
  ee_log10 <- log10(abs(rel_err))
  
  list(
    lhs = lhs,
    rhs = rhs,
    k_mid   = k_mid[good],
    resid   = ee_resid,
    log10err = ee_log10
  )
}

#########################################################
# Non calibration parameters
#########################################################
alpha <- 0.36
beta  <- 0.95
delta <- 0.05
theta <- 4

## n_k = 20
res20 <- howard_vfi_track(
  alpha = alpha,
  beta  = beta,
  delta = delta,
  theta = theta,
  n_k   = 20
)

ee20 <- euler_errors_midpoints(res20, alpha, beta, delta, theta)

cat("n_k = 20\n")
cat("  Max |Euler residual|:", max(abs(ee20$resid)), "\n")
cat("  Max log10 Euler error:", max(ee20$log10err), "\n\n")

## n_k = 30
res30 <- howard_vfi_track(
  alpha = alpha,
  beta  = beta,
  delta = delta,
  theta = theta,
  n_k   = 30
)

plot(ee20$k_mid, ee20$log10err, type = "l",
     xlab = "k (midpoints)", ylab = "log10 Euler error",
     main = "Euler equation errors at midpoints (n_k = 20)")

plot(ee20$k_mid, ee20$lhs, type = "l",
     xlab = "k (midpoints)",
     ylab = "u'(c) and beta u'(c') [f'(k')+1-delta]",
     main = "Euler equation: LHS and RHS over k")
lines(ee20$k_mid, ee20$rhs, lty = 2)

legend("topright", legend = c("LHS = u'(c_t)", "RHS = beta u'(c_{t+1})(...)"),
       lty = c(1, 2), bty = "n")

###########################
# What is up?
###########################
# It generally looks nice, but it has this weird kink on the right for some reason

#########################
# What is up with the euler equation gridpoints
#########################
# So first of all it is not monotonic, I suppose this can suggest that the Euler

ee30 <- euler_errors_midpoints(res30, alpha, beta, delta, theta)

cat("n_k = 30\n")
cat("  Max |Euler residual|:", max(abs(ee30$resid)), "\n")
cat("  Max log10 Euler error:", max(ee30$log10err), "\n")

plot(ee30$k_mid, ee30$log10err, type = "l",
     xlab = "k (midpoints)", ylab = "log10 Euler error",
     main = "Euler equation errors at midpoints (n_k = 30)")
plot(ee30$k_mid, ee30$lhs, type = "l",
     xlab = "k (midpoints)",
     ylab = "u'(c) and beta u'(c') [f'(k')+1-delta]",
     main = "Euler equation: LHS and RHS over k")
lines(ee30$k_mid, ee30$rhs, lty = 2)

legend("topright", legend = c("LHS = u'(c_t)", "RHS = beta u'(c_{t+1})(...)"),
       lty = c(1, 2), bty = "n")

# Same kind persists even with 3 gridpoints, maybe there's numerical instability?

###############################################################################
## Homotopy loops
###############################################################################

# Base parameters (keep alpha, beta as in your PS)
alpha <- 0.36
beta  <- 0.95

# Sequences to explore
theta_vals <- seq(0, 4, by = 1)     # 0,1,2,3,4
delta_vals <- seq(0, 1, by = 0.25)  # 0,0.25,0.5,0.75,1

# Use more grid points so policy functions look smooth-ish
n_k_sweep <- 50

# Set up a grid of panels: rows = theta, cols = delta
par(mfrow = c(length(theta_vals), length(delta_vals)),
    mar = c(3, 3, 2, 1))  # smaller margins

for (th in theta_vals) {
  for (del in delta_vals) {
    
    # Solve with Howard's algorithm for this (theta, delta)
    res_td <- howard_vfi_track(
      alpha = alpha,
      beta  = beta,
      delta = del,
      theta = th,
      n_k   = n_k_sweep
    )
    
    # Plot policy function k'(k)
    plot(res_td$k_grid, res_td$k_policy,
         type = "l",
         xlab = "k",
         ylab = "k'(k)",
         main = bquote(theta == .(th) ~ "," ~ delta == .(del)))
    
    # 45° line for reference
    abline(0, 1, lty = 2)
  }
}


      
      
