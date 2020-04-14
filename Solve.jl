#solving Urea Flux with calculated data files
include("Include.jl")
include("flux.jl")
using GLPK

#import stoichiometric matrix
StoicMat=Matrix(CSV.read("Stoich_Matrix.csv",header=0));

#objective: maximize urea production
UMax = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#using Vmax for upper bound
reactionBounds = Matrix(CSV.read("FBounds.csv"))
speciesBounds = Matrix(CSV.read("SPCBounds.csv"))


#Solve problem
objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag,status_flag = calculate_optimal_flux_distribution(StoicMat, reactionBounds, speciesBounds, UMax)
print(objective_value)