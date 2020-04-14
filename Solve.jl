include("Include.jl")
include("flux.jl")
using CSV

#import stoichiometric matrix
S=Matrix(CSV.read("Stoich_Matrix.csv",header=0));

#import bounds
ReactionBounds=Matrix(CSV.read("minmaxBounds.csv",header=0))
SpeciesBounds=Matrix(CSV.read("speciesBounds.csv",header=0))

objective=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag,status_flag = calculate_optimal_flux_distribution(S,ReactionBounds,SpeciesBounds,objective)
print(objective_value)