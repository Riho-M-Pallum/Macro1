###############################################################################
## Howard's Improvement for the Neoclassical Growth Model in R


howard_vfi <- function(alpha = 0.36,
                       beta  = 0.95,
                       delta = 0.05,
                       theta = 4,
                       n_k   = 20,
                       tol_value  = 1e-8,
                       tol_policy = 1e-6,
                       max_iter_value  = 10000,
                       max_iter_policy = 1000) {
  u_fun <- function(c) {
    # CRRA / iso-elastic utility
    if (theta == 1) return(log(c))
    (c^(1 - theta) - 1) / (1 - theta)
  }
  
  f_fun <- function(k) {
    k^alpha
  }
  
  # Generate an equally spaced capital grid
  # First we calculate the steady state capital stock, under the assumption that f(k) = k^{\alpha}
  k_ss <- (alpha / (1 / beta - 1 + delta))^(1 / (1 - alpha))
  
  # Equally spaced grid for capital
  k_grid <- seq(0.1 * k_ss, 2 * k_ss, length.out = n_k)
  n <- n_k
  
  # Compute output and resources z(k) = f(k) + (1 - δ)k
  fK        <- f_fun(k_grid)
  resources <- fK + (1 - delta) * k_grid
  
  #----------------------------------------------------------
  # 3. Precompute one-period utility u(k, k') over the grid
  #    imposing c >= 0 => k' <= f(k) + (1-δ)k
  #----------------------------------------------------------
  # util[i,j] = u( c(i,j) ) if feasible, else very negative
  util <- matrix(-1e12, nrow = n, ncol = n)  # large negative as -∞ approximation
  
  for (i in 1:n) {
    c_ij <- resources[i] - k_grid                 # c = z_i - k'_j
    feasible <- c_ij > 0                          # c>0; if not, rule out choice
    u_ij <- rep(-1e12, n)
    u_ij[feasible] <- u_fun(c_ij[feasible])
    util[i, ] <- u_ij
  }
  
  #-------------------------------------------------
  # 4. Initialization of value function and policy
  #-------------------------------------------------
  V <- rep(0, n)            # initial guess for V(k)
  policy_idx <- rep(1, n)   # initial policy: choose lowest k' everywhere
  
  #----------------------------------------
  # 5. Howard improvement outer (policy) loop
  #----------------------------------------
  for (it_pol in 1:max_iter_policy) {
    
    #----------------------------------
    # (a) Inner loop: V given policy g
    #     This is Howard's improvement:
    #     assume capital evolution law k' = g(k) fixed
    #----------------------------------
    for (it_val in 1:max_iter_value) {
      # V_new(k_i) = u(k_i, g(k_i)) + β V_old(g(k_i))
      V_new <- util[cbind(1:n, policy_idx)] + beta * V[policy_idx]
      diff_V <- max(abs(V_new - V))
      V <- V_new
      if (diff_V < tol_value) break
    }
    
    #------------------------------------------------
    # (b) Policy improvement step:
    #     Given new V, recompute optimal policy
    #------------------------------------------------
    # For each (i,j): RHS(i,j) = u(i,j) + β V(k'_j)
    # exploing matrix broadcasting
    value_candidate <- util + beta * matrix(rep(V, each = n), nrow = n)
    
    # For each i, find argmax_j value_candidate[i,j]
    new_policy_idx <- max.col(value_candidate, ties.method = "first")
    
    # Check for convergence of policy
    policy_diff <- max(abs(k_grid[new_policy_idx] - k_grid[policy_idx]))
    policy_idx <- new_policy_idx
    
    if (policy_diff < tol_policy) {
      # converged policy
      break
    }
  }
  
  #---------------------------------------
  # 6. Recover policy functions
  #---------------------------------------
  k_policy <- k_grid[policy_idx]        # policy k'(k)
  c_policy <- resources - k_policy      # consumption c(k) = f(k)+(1-δ)k - k'(k)
  
  # Return all objects you might want to compare with your code
  return(list(
    k_grid    = k_grid,
    V         = V,
    policy_idx = policy_idx,
    k_policy  = k_policy,
    c_policy  = c_policy,
    k_ss      = k_ss
  ))
}

###############################################################################
## Example usage (you can tweak parameters to match your PS exactly)
###############################################################################



# Final calibration from the problem set:
res_howard <- howard_vfi(
  alpha = 0.36,
  beta  = 0.95,
  delta = 1,
  theta = 1,
  n_k   = 20,
  max_iter_value = 10
)
res_howard
0.95*0.36*0.243^(0.36)
# Quick sanity checks:
plot(res_howard$k_grid, res_howard$k_policy, type = "l",
      xlab = "k", ylab = "k'(k)", main = "Policy function with Howard's improvement")
      abline(0, 1, lty = 2)  # 45-degree line for reference
