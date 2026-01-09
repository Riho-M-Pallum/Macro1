using CSV
using DataFrames
using Statistics
using Plots
using TrendDecomposition
using StatsBase
using GLM

# Load the CSV file
# filepath: .\A939RX0Q048SBEA.csv
function process_data(file_path::String)
    # Step 1: Load the data
    df = CSV.read(file_path, DataFrame)
    output = log.(df.A939RX0Q048SBEA)  # Take log of the output series since otherwise we get stupid large deviations

    # Step 2: Apply the HP filter
    trend = TrendDecomposition.hpFilter(output, 1600)  # λ=1600 is standard for quarterly data
    cycle = output .- trend
    # Step 3: Estimate an AR(1) model
    # Create lagged data
    lagged_cycle = StatsModels.lag(cycle, 1)
    lagged_cycle = lagged_cycle[2:end]  # Remove the first NaN
    cycle = cycle[2:end]  # Align with lagged data

    # Fit AR(1) model using GLM
    data = DataFrame(cycle=cycle, lagged_cycle=lagged_cycle)
    model = lm(@formula(cycle ~ lagged_cycle), data)
    mean_cycle = mean(cycle)
    variance_cylce = var(cycle; corrected=true)
    max_cycle = maximum(cycle)
    min_cycle = minimum(cycle)
    print("Average cycle: $mean_cycle\n")
    print("Variance of cycle: $variance_cylce\n")
    print("maximum of cycle: $max_cycle\n")
    print("minimum of cycle: $min_cycle\n")
    # Print results



    # Plot the original data, trend, and cycle
    plot(1:length(output), output, label="Original Data", lw=2)
    plot!(1:length(trend), trend, label="Trend", lw=2)
    savefig("output_trend.png")
    plot(1:length(cycle), cycle, label="Cycle (HP Filtered)", lw=2)
    savefig("output_cycle.png")
    return model
end

# Call the function with the path to your CSV file
model = process_data("A939RX0Q048SBEA.csv")