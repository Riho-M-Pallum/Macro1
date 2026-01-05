using LinearAlgebra
using Interpolations
using Roots
using Optim

"""
Time iteration for stochastic growth with detrended resource constraint:
    c = ζ*k^α + (1-δ)*k - G*k'
Euler:
    c^{-θ} = (β/G) * E[ c'^{-θ} * ( ζ' α (k')^{α-1} + (1-δ) ) ]

Inputs
------
kgrid::Vector{Float64}        capital grid (increasing)
svals::Vector{Float64}        shock states (size Ns)
Π::Matrix{Float64}            transition matrix (Ns×Ns), rows sum to 1
β,G,δ,α,θ::Float64            parameters
keyword: maxit, tol, damp

Returns
-------
g::Matrix{Float64}            policy k' on grid (N×Ns)


Calibration
----------
- β = 0.95
- δ = 0.05
- θ = 4
- α = 0.36
- G = 1.04 (not based on previous problem set)
"""
function time_iteration(kgrid, svals, Π; β = 0.95, G = 1.04, δ=0.05, α=0.36, θ = 4,
                        maxit=10000, tol=1e-7, damp=1.0)

    N  = length(kgrid)
    Ns = length(svals)
    @assert size(Π) == (Ns, Ns)

    # Initial policy: keep capital constant
    g_old = repeat(kgrid, 1, Ns)

    # Small number to keep c>0
    ϵ = 1e-12

    for it in 1:maxit
        # Build interpolants for each shock state: g_old(k, s_j)
        itps = Vector{Any}(undef, Ns)
        for j in 1:Ns
            base = interpolate((kgrid,), g_old[:, j], Gridded(Linear())) # Linear interpolation object
            itps[j] = extrapolate(base, Interpolations.Flat())  # For extrapolation outside of grid boundaries, set value to the grid boundary (Flat()).
        end
        # new array of the same size as the old one without copying data.
        g_new = similar(g_old)
        # For each shock
        for j in 1:Ns
            s = svals[j]
            # For each capital grid point
            for i in 1:N
                k = kgrid[i]

                # Feasibility: c = s*k^α + (1-δ)k - G*k'  > 0
                resources = s * k^α + (1 - δ) * k
                kp_max = resources / G - ϵ # Maximum capital next period, ensuring c>0
                kp_min = kgrid[1] # Because the capital grid is an ordered list (see lecture notes)

                # If infeasible at this (k,s), clamp (usually means grid too wide)
                if kp_max <= kp_min + ϵ
                    g_new[i, j] = kp_min
                    continue
                end                # Euler residual as a function of k' (kp), 
                # i.e. LHS - RHS of Euler equation for different values of k' on the grid.
                function resid(kp)
                    c = resources - G * kp
                    c = max(c, ϵ)

                    exp_term = 0.0 #Expected value of capital next period
                    # For each shock value
                    for jp in 1:Ns
                        sp = svals[jp]
                        # k'' from current policy at next state (kp, sp)
                        kpp = itps[jp](kp)
                        # calculate consunption next period by the resource constraint
                        cp = sp * kp^α + (1 - δ) * kp - G * kpp
                        cp = max(cp, ϵ)
                        # Envelope conditoin ∂c/∂k partial derivative of utility w.r.t. k'
                        R = sp * α * kp^(α - 1) + (1 - δ)
                        exp_term += Π[j, jp] * (cp^(-θ)) * R
                    end

                    lhs = c^(-θ)
                    rhs = (β / G) * exp_term
                    return lhs - rhs
                end

                rmin = resid(kp_min)
                rmax = resid(kp_max)

                # We look for the policy that makes the Euler equaiton hold using Brent's method.
                # Note that if there are multiple roots, then we will only find one of them.
                # The structure of the problem is such that there is only one root though.
                # Since the LHS is fixed in kp, while the RHS is monotonically increasing in kp.
                kp_star =
                    if sign(rmin) != sign(rmax)
                        # Bracketed root: robust
                        Roots.find_zero(resid, (kp_min, kp_max), Roots.Brent())
                    else
                        # No sign change: take boundary (handles corners / numerical issues)
                        # Choose the point with smaller absolute residual
                        abs(rmin) < abs(rmax) ? kp_min : kp_max
                    end
                    
                # Optional damping
                g_new[i, j] = damp * kp_star + (1 - damp) * g_old[i, j]
            end
        end

        err = maximum(abs.(g_new .- g_old))
        g_old .= g_new

        if err < tol
            # Converged
            return g_old
        end
    end

    @warn "time_iteration did not converge in maxit iterations"
    return g_old
end
