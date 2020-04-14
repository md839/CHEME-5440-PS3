#flux.jl solves the Flux balance

using GLPK

function calculate_optimal_flux_distribution(stoichMat::Array{Float64,2},reactionBounds::Array{Float64,2},speciesBounds::Array{Float64,2},UMax::Array{Float64,1})
    
#size of stoichiometric_matrix
    (numSpecies,numFluxes) = size(stoichMat); 
    
#size of reactionBounds
    (numFluxes,numBounds)=size(reactionBounds)
    
#initialize problem and GLPK
    Time_Restart_Lim=60
    lp_problem=GLPK.Prob();GLPK.set_prob_name(lp_problem,"sample")
    GLPK.set_obj_name(lp_problem, "UMax")
    solver_parameters = GLPK.SimplexParam();
    GLPK.init_smcp(solver_parameters);
    solver_parameters.msg_lev = GLPK.MSG_OFF;
    GLPK.add_rows(lp_problem, numSpecies);
    GLPK.add_cols(lp_problem, numFluxes);
    
#always minimize
    GLPK.set_obj_dir(lp_problem, GLPK.MIN);
    
#flux bounds for objective function
    for i=1:numFluxes
        fluxLower=reactionBounds[i,1]
        fluxUpper=reactionBounds[i,2]

        #check if completely constrained
        if (fluxUpper == fluxLower)
    		flux_constraint_type = GLPK.FX
    	else
    		flux_constraint_type = GLPK.DB
    	end

    	# flux symbol? (later use name - for now, fake it)
    	flux_symbol = "R_"*string(i)

    	# Set the bounds in GLPK -
    	GLPK.set_col_name(lp_problem, i, flux_symbol);
    	GLPK.set_col_bnds(lp_problem, i, flux_constraint_type, fluxLower, fluxUpper);
    end

#objective bounds for objective function
    for (j,k) in enumerate(UMax)

    	
# Set the objective function value in GLPK -
    	GLPK.set_obj_coef(lp_problem, j, k);
    end

    
#constraints on metaboliltes
    for l = 1:numSpecies

    	speciesLower = speciesBounds[l,1]
    	speciesUpper= speciesBounds[l,2]
    	species_constraint_type = GLPK.FX

    	if (speciesLower != speciesUpper)
    		species_constraint_type = GLPK.DB
    	end

    	# set the symbol -
    	species_symbol = "x_"*string(l)

    	# Set the species bounds in GLPK -
    	GLPK.set_row_name(lp_problem, l, species_symbol);
    	GLPK.set_row_bnds(lp_problem, l, species_constraint_type, speciesLower, speciesUpper);
    end


 # Setup the stoichiometric array -
 counter = 1;
 row_index_array = zeros(Int,numSpecies*numFluxes);
 col_index_array = zeros(Int,numSpecies*numFluxes);
 species_index_vector = collect(1:numSpecies);
 flux_index_vector = collect(1:numFluxes);
 flat_stoichiometric_array = zeros(Float64,numSpecies*numFluxes);

 for species_index in species_index_vector
     for flux_index in flux_index_vector
         row_index_array[counter] = species_index;
         col_index_array[counter] = flux_index;
         flat_stoichiometric_array[counter] = stoichMat[species_index,flux_index];
         counter+=1;
     end
 end

 GLPK.load_matrix(lp_problem, numSpecies*numFluxes, row_index_array, col_index_array, flat_stoichiometric_array);

 # Call the solver -
 exit_flag = GLPK.simplex(lp_problem, solver_parameters);

 # Get the objective function value -
 objective_value = GLPK.get_obj_val(lp_problem);

 # Get the calculated flux values from GLPK -
 calculated_flux_array = zeros(Float64,numFluxes);
 for flux_index in flux_index_vector
     calculated_flux_array[flux_index] = GLPK.get_col_prim(lp_problem, flux_index);
 end

 # Get the dual values -

 dual_value_array = zeros(Float64,numFluxes);

 for flux_index in flux_index_vector
     dual_value_array[flux_index] = GLPK.get_col_dual(lp_problem, flux_index);
 end

 # is this solution optimal?
 status_flag = GLPK.get_status(lp_problem)

 # Calculate the uptake array -
 uptake_array = stoichMat*calculated_flux_array;

 # Formulate the return tuple -
    return (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag,status_flag);

end





