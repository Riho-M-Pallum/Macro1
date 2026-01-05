using CSV
using DataFrames
using Statistics
using Plots
using StatsBase
using GLM
using QuantEcon  # For Rouwenhorst's approximation
using PrettyTables

include("ar_model.jl")  # Include the AR(1) estimation script
include("vlf.jl")       # Include the time iteration script

function rouwenhorst_to_markov(ρ::Float64, σ::Float64, n::Int=2)
    # Use Rouwenhorst's method to approximate AR(1) as a Markov chain
    mc = QuantEcon.rouwenhorst(n, ρ, σ)
    return mc
end

function simulate_model(kgrid, svals, Π, g, T::Int; β=0.95, δ=0.05, α=0.36, θ=4, G=1.04)
    # Simulate the model for T periods using the policy function g
    N = length(kgrid)
    Ns = length(svals)

    # Initialize variables
    k = kgrid[1]  # Start with the smallest capital
    itps = Vector{Any}(undef, Ns)
    for j in 1:Ns
        base = interpolate((kgrid,), g[:, j], Gridded(Linear())) # Linear interpolation object
        itps[j] = extrapolate(base, Interpolations.Flat())  # For extrapolation outside of grid boundaries, set value to the grid boundary (Flat()).
    end
    
    s_index = 1   # Start with the first shock state
    output = zeros(T)
    consumption = zeros(T)
    investment = zeros(T)

    for t in 1:T
        # Current shock and capital
        s = svals[s_index]

        # Output: y = s * k^α
        output[t] = s * k^α

        # Optimal next-period capital: k' = g(k, s)
        kp= itps[s_index](k)  # Find the index of the current capital
        # Investment: i = G*kp - (1 - δ) * k
        investment[t] = G*kp - (1 - δ) * k

        # Consumption: c = y - i
        consumption[t] = output[t] - investment[t]

        # Update capital for the next period
        k = kp

        # Transition to the next shock state
        s_index = findfirst(rand() .< cumsum(Π[s_index, :]))
    end

    return output, consumption, investment
end

function analyze_results(output, consumption, investment)
    # Calculate variances
    var_output = Statistics.var(output)
    var_consumption = Statistics.var(consumption)
    var_investment = Statistics.var(investment)

    # Calculate covariances
    cor_output_consumption = Statistics.cov(output, consumption)/sqrt(var_output * var_consumption)
    cor_output_investment = Statistics.cov(output, investment)/sqrt(var_output * var_investment)
    cor_consumption_investment = Statistics.cov(consumption, investment)/sqrt(var_consumption * var_investment)

    # Calculate autocorrelations
    ac_output = StatsBase.autocor(output, [1])[1]
    ac_consumption = StatsBase.autocor(consumption, [1])[1]
    ac_investment = StatsBase.autocor(investment, [1])[1]

    # Prepare data for the LaTeX table
    data = (
        Metric = [
            "Variance",
            "Correlation (Output, Consumption)",
            "Correlation (Output, Investment)",
            "Correlation (Consumption, Investment)",
            "Autocorrelation (Output)",
            "Autocorrelation (Consumption)",
            "Autocorrelation (Investment)"
        ],
        Output = [
            var_output,
            cor_output_consumption,
            cor_output_investment,
            "",
            ac_output,
            "",
            ""
        ],
        Consumption = [
            var_consumption,
            "",
            "",
            cor_consumption_investment,
            "",
            ac_consumption,
            ""
        ],
        Investment = [
            var_investment,
            "",
            "",
            "",
            "",
            "",
            ac_investment
        ]
    )

    # Define column headers
    header = ["Metric", "Output", "Consumption", "Investment"]

    # Print the LaTeX table
    println("LaTeX Table:")
    pretty_table(data, backend=:latex, alignment=:lcrr)
end

# Main script
function main()
    print("Estimating and simulating model\n -------------------------------\n")
    # Step 1: Run AR(1) model estimation
    model = process_data("A939RX0Q048SBEA.csv")
    ρ = GLM.coef(model)[2]  # AR(1) coefficient
    σϵ = std(residuals(model); corrected=true)
    σζ = σϵ / sqrt(1 - ρ^2)  # Standard deviation of the shocks
    
    print("Estimated AR(1) parameters: ρ = $ρ, σϵ = $σϵ, σζ = $σζ\n")
    # Step 2: Convert AR(1) to Markov chain using Rouwenhorst's method
    mc = rouwenhorst_to_markov(ρ, σζ, 2)

    svals = exp.(mc.state_values)
    print("The possible shock values: $svals\n")
    Π = mc.p

    # Step 3: Solve the model using time iteration
    kgrid = range(0.1, 10.0, length=100)  # Define a grid for capital
    g = time_iteration(kgrid, svals, Π; β=0.95, G=1.01, δ=0.05, α=0.36, θ=4)
    plot(kgrid, g[:, 1], label="Policy Function (s = $(round(svals[1], digits=3)))", lw=2, xlabel="Capital k", ylabel="Next-period Capital k'", title="Policy Function")
    plot!(kgrid, g[:, 2], label="Policy Function (s = $(round(svals[2], digits=3)))", lw=2)
    savefig("policy_function_plot.png")
    # Step 4: Simulate the model

    burn_in = 1000  # Burn-in periods to discard

    print("Simulating the model for T = 1000\n")
    T = 1_000
    output, consumption, investment = simulate_model(kgrid, svals, Π, g, T+burn_in)
    output = output[burn_in+1:end]
    consumption = consumption[burn_in+1:end]
    investment = investment[burn_in+1:end]
    plot(1:T, output, label="Output", lw=2, xlabel="Time", ylabel="Value", title="Simulation Results (T = 1000)")
    plot!(1:T, consumption, label="Consumption", lw=2)
    plot!(1:T, investment, label="Investment", lw=2)
    savefig("simulation1k_plot.png")
    # Step 5: Analyze results
    analyze_results(output, consumption, investment)

    print("Simulating the model for T = 100000\n")
    T = 100_000
    output, consumption, investment = simulate_model(kgrid, svals, Π, g,T+burn_in)
    output = output[burn_in+1:end]
    consumption = consumption[burn_in+1:end]
    investment = investment[burn_in+1:end]
    analyze_results(output, consumption, investment)
    plot(1:T, output, label="Output", lw=2, xlabel="Time", ylabel="Value", title="Simulation Results (T = 1000)")
    plot!(1:T, consumption, label="Consumption", lw=2)
    plot!(1:T, investment, label="Investment", lw=2)
    savefig("simulation100k_plot.png")

end

main()