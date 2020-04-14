# ----------------------------------------------------------------------------------- #
# Copyright (c) 2018 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

function check_unbalanced_boundary_reactions(time_start,time_stop,time_step)

	# load the original data dictionary -
	data_dictionary = maximize_urea_production(time_start,time_stop,time_step)

	# add boundary species to the list of species -
	list_of_metabolite_symbols = data_dictionary["list_of_metabolite_symbols"]
	list_of_metabolite_symbols = [
		list_of_metabolite_symbols ;

		"M_Carbamoyl_phosphate_b"	;	# 1
		"M_L-Aspartate_b"	;			# 2
		"M_Fumarate_b"	;				# 3
		"M_Urea_b"	;					# 4
		"M_ATP_b"	;					# 5
		"M_AMP_b"	;					# 6
		"M_Diphosphate_b"	;			# 7
		"M_Orthophosphate_b"	;		# 8
		"M_Oxygen_b"	;				# 9
		"M_NADPH_b"	;					# 10
		"M_H_b"	;						# 11
		"M_Nitric_oxide_b"	;			# 12
		"M_NADP_b"	;					# 13
		"M_H2O_b"	;					# 14
	]

	# add some boundary species (rows to the stoichiometric array)
	# there will be in the same order as the list_of_metabolite_symbols
	S = data_dictionary["stoichiometric_matrix"]

	# ----------------------------------------------------- #
	# 7,b1,0,[],[],M_Carbamoyl_phosphate_c,0,inf,[]
	# 8,b2,0,[],[],M_L-Aspartate_c,0,inf,[]
	# 9,b3,0,[],M_Fumarate_c,[],0,inf,[]
	# 10,b4,0,[],M_Urea_c,[],0,inf,[]
	# 11,b5,0,[],[],M_ATP_c,0,inf,[]
	# 12,b6,0,[],M_AMP_c,[],0,inf,[]
	# 13,b7,0,[],M_Diphosphate_c,[],0,inf,[]
	# 14,b8,0,[],M_Orthophosphate_c,[],0,inf,[]
	# 15,b9,0,[],[],M_Oxygen_c,0,inf,[]
	# 16,b10,0,[],[],M_NADPH_c,0,inf,[]
	# 17,b11,0,[],[],M_H_c,0,inf,[]
	# 18,b12,0,[],M_Nitric_oxide_c,[],0,inf,[]
	# 19,b13,0,[],M_NADP_c,[],0,inf,[]
	# 20,b14,0,[],M_H2O_c,[],0,inf,[]
	# 21,b14,0,[],[],M_H2O_c,0,inf,[]
	# ------------------------------------------------------ #

	# transfer reaction block -
	transfer_reaction_block = zeros(14,21)
	transfer_reaction_block[1,7] = -1 		# 1 M_Carbamoyl_phosphate_b -> box
	transfer_reaction_block[2,8] = -1 		# 2 M_L-Aspartate_b -> box
	transfer_reaction_block[3,9] =	1		# 3 box -> M_Fumarate_b
	transfer_reaction_block[4,10] = 1		# 4 box -> M_Urea_b
	transfer_reaction_block[5,11] = -1		# 5 M_ATP_b -> box
	transfer_reaction_block[6,12] = 1		# 6 box -> M_AMP_b
	transfer_reaction_block[7,13] = 1		# 7 box -> M_Diphosphate_b
	transfer_reaction_block[8,14] = 1		# 8 box -> M_Orthophosphate_b
	transfer_reaction_block[9,15] = -1		# 9 M_Oxygen_b -> box
	transfer_reaction_block[10,16] = -1		# 10 M_NADPH_b -> box
	transfer_reaction_block[11,17] = -1		# 11 M_H_b -> box
	transfer_reaction_block[12,18] = 1		# 12 box -> M_Nitric_oxide_b
	transfer_reaction_block[13,19] = 1		# 13 box -> M_NADP_b
	transfer_reaction_block[14,20] = 1		# 14 box -> M_H2O_b
	transfer_reaction_block[14,21] = -1		# 15 M_H2O_b -> box
	SAUG = [S ; transfer_reaction_block]

	# repackage -
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["stoichiometric_matrix"] = SAUG

	# return the modified data dictionary -
	return data_dictionary
end

function maximize_urea_production_open(time_start,time_stop,time_step)

	# load the original data dictionary -
	data_dictionary = maximize_urea_production(time_start,time_stop,time_step)

	# lets open up the side products of ec:1.14.13.39 -

	# 2: Update the reaction bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	# Vmax [mmol/gdw-hr] 15	[] --> M_Oxygen_c
	default_flux_bounds_array[15,1] = -10
	default_flux_bounds_array[15,2] = 10

	# Vmax [mmol/gdw-hr] 16	[] --> M_NADPH_c
	default_flux_bounds_array[16,1] = -10
	default_flux_bounds_array[16,2] = 10

	# Vmax [mmol/gdw-hr] 17	[] --> M_H_c
	default_flux_bounds_array[17,1] = -10
	default_flux_bounds_array[17,2] = 10

	# Vmax [mmol/gdw-hr] 18	M_Nitric_oxide_c --> []
	default_flux_bounds_array[18,1] = -10
	default_flux_bounds_array[18,2] = 10

	# Vmax [mmol/gdw-hr] 19	M_NADP_c --> []
	default_flux_bounds_array[19,1] = -10
	default_flux_bounds_array[19,2] = 10

	# Vmax [mmol/gdw-hr] 20	M_H2O_c --> []
	default_flux_bounds_array[20,1] = -10
	default_flux_bounds_array[20,2] = 10

	# Vmax [mmol/gdw-hr] 21	[] --> M_H2O_c
	default_flux_bounds_array[21,1] = -10
	default_flux_bounds_array[21,2] = 10

	# repackage -
	data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array

	# return the updated dictionary -
	return data_dictionary
end

function maximize_urea_production(time_start,time_stop,time_step)

	# load the original data dictionary -
	data_dictionary = DataDictionary(time_start,time_stop,time_step)

	# 1: set the objective function -
	objective_coefficient_array = data_dictionary["objective_coefficient_array"]
	objective_coefficient_array[10] = -1

	# 2: Update the reaction bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	# let all exchanges go from 0,10 mmol/gDW-hr
	range_of_exchange_reactions = collect(7:21)
	for reaction_index in range_of_exchange_reactions
		default_flux_bounds_array[reaction_index,1] = 0.0
		default_flux_bounds_array[reaction_index,2] = 10.0
	end

	# don't allow water exchange -
	default_flux_bounds_array[20,2] = 0
	default_flux_bounds_array[21,2] = 0

	# we have some specific values for v1 -> v5
	E = (0.01)*(1/1000)	# mmol/gDW
	metabolic_vmax_array = [
		203*(3600)*E	;	# v1 ec:6.3.4.5 mmol/gDW-hr
		34.5*(3600)*E	;	# v2 ec:4.3.2.1 mmol/gDW-hr
		249*(3600)*E	;	# v3 ec:3.5.3.1 mmol/gDW-hr
		88.1*(3600)*E	;	# v4 ec:2.1.3.3 mmol/gDW-hr
		13.7*(3600)*E	;	# v5 ec:1.14.13.39 mmol/gDW-hr
		13.7*(3600)*E	;	# v6 ec:1.14.13.39 mmol/gDW-hr
	]
	range_of_cycle_reactions = collect(1:6)
	for reaction_index in range_of_cycle_reactions
		default_flux_bounds_array[reaction_index,1] = 0.0
		default_flux_bounds_array[reaction_index,2] = metabolic_vmax_array[reaction_index]
	end

	# repackage -
	data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array

	# return the updated dictionary -
	return data_dictionary
end


# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2018-03-15T00:00:56.939
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function DataDictionary(time_start,time_stop,time_step)

	# Load the stoichiometric network from disk -
	stoichiometric_matrix = readdlm("Network.dat");

	# What is the system dimension? -
	(number_of_species,number_of_reactions) = size(stoichiometric_matrix)

	# Setup default flux bounds array -
	default_bounds_array = [
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 1	M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 2	M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 3	M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 4	M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 5	2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 6	2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c

		0	2.2148976	;	# Vmax [mmol/gdw-hr] 7	[] --> M_Carbamoyl_phosphate_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 8	[] --> M_L-Aspartate_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 9	M_Fumarate_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 10	M_Urea_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 11	[] --> M_ATP_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 12	M_AMP_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 13	M_Diphosphate_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 14	M_Orthophosphate_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 15	[] --> M_Oxygen_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 16	[] --> M_NADPH_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 17	[] --> M_H_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 18	M_Nitric_oxide_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 19	M_NADP_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 20	M_H2O_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 21	[] --> M_H2O_c
	];


	# Setup default species bounds array -
	species_bounds_array = [
		0.0	0.0	;	# 1 M_AMP_c
		0.0	0.0	;	# 2 M_ATP_c
		0.0	0.0	;	# 3 M_Carbamoyl_phosphate_c
		0.0	0.0	;	# 4 M_Diphosphate_c
		0.0	0.0	;	# 5 M_Fumarate_c
		0.0	0.0	;	# 6 M_H2O_c
		0.0	0.0	;	# 7 M_H_c
		0.0	0.0	;	# 8 M_L-Arginine_c
		0.0	0.0	;	# 9 M_L-Aspartate_c
		0.0	0.0	;	# 10 M_L-Citrulline_c
		0.0	0.0	;	# 11 M_L-Ornithine_c
		0.0	0.0	;	# 12 M_N-(L-Arginino)succinate_c
		0.0	0.0	;	# 13 M_NADPH_c
		0.0	0.0	;	# 14 M_NADP_c
		0.0	0.0	;	# 15 M_Nitric_oxide_c
		0.0	0.0	;	# 16 M_Orthophosphate_c
		0.0	0.0	;	# 17 M_Oxygen_c
		0.0	0.0	;	# 18 M_Urea_c
	];

	# Min/Max flag - default is minimum -
	is_minimum_flag = true

	# Setup the objective coefficient array -
	objective_coefficient_array = [

		0.0	;	# 1 v1::M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
		0.0	;	# 2 v2::M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
		0.0	;	# 3 v3::M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
		0.0	;	# 4 v4::M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
		0.0	;	# 5 v5::2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
		0.0	;	# 6 v5_reverse::2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c
		0.0	;	# 7 b1::[] --> M_Carbamoyl_phosphate_c
		0.0	;	# 8 b2::[] --> M_L-Aspartate_c
		0.0	;	# 9 b3::M_Fumarate_c --> []
		0.0	;	# 10 b4::M_Urea_c --> []
		0.0	;	# 11 b5::[] --> M_ATP_c
		0.0	;	# 12 b6::M_AMP_c --> []
		0.0	;	# 13 b7::M_Diphosphate_c --> []
		0.0	;	# 14 b8::M_Orthophosphate_c --> []
		0.0	;	# 15 b9::[] --> M_Oxygen_c
		0.0	;	# 16 b10::[] --> M_NADPH_c
		0.0	;	# 17 b11::[] --> M_H_c
		0.0	;	# 18 b12::M_Nitric_oxide_c --> []
		0.0	;	# 19 b13::M_NADP_c --> []
		0.0	;	# 20 b14::M_H2O_c --> []
		0.0	;	# 21 b14_reverse::[] --> M_H2O_c

	];

	# List of reation strings - used to write flux report
	list_of_reaction_strings = [
		"v1::M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c"	;	# 1
		"v2::M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c"	;	# 2
		"v3::M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c"	;	# 3
		"v4::M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c"	;	# 4
		"v5::2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c"	;	# 5
		"v5_reverse::2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c"	;	# 6
		"b1::[] --> M_Carbamoyl_phosphate_c"	;	# 7
		"b2::[] --> M_L-Aspartate_c"	;	# 8
		"b3::M_Fumarate_c --> []"	;	# 9
		"b4::M_Urea_c --> []"	;	# 10
		"b5::[] --> M_ATP_c"	;	# 11
		"b6::M_AMP_c --> []"	;	# 12
		"b7::M_Diphosphate_c --> []"	;	# 13
		"b8::M_Orthophosphate_c --> []"	;	# 14
		"b9::[] --> M_Oxygen_c"	;	# 15
		"b10::[] --> M_NADPH_c"	;	# 16
		"b11::[] --> M_H_c"	;	# 17
		"b12::M_Nitric_oxide_c --> []"	;	# 18
		"b13::M_NADP_c --> []"	;	# 19
		"b14::M_H2O_c --> []"	;	# 20
		"b14_reverse::[] --> M_H2O_c"	;	# 21
	];

	# List of metabolite strings - used to write flux report
	list_of_metabolite_symbols = [
		"M_AMP_c"	;	# 1
		"M_ATP_c"	;	# 2
		"M_Carbamoyl_phosphate_c"	;	# 3
		"M_Diphosphate_c"	;	# 4
		"M_Fumarate_c"	;	# 5
		"M_H2O_c"	;	# 6
		"M_H_c"	;	# 7
		"M_L-Arginine_c"	;	# 8
		"M_L-Aspartate_c"	;	# 9
		"M_L-Citrulline_c"	;	# 10
		"M_L-Ornithine_c"	;	# 11
		"M_N-(L-Arginino)succinate_c"	;	# 12
		"M_NADPH_c"	;	# 13
		"M_NADP_c"	;	# 14
		"M_Nitric_oxide_c"	;	# 15
		"M_Orthophosphate_c"	;	# 16
		"M_Oxygen_c"	;	# 17
		"M_Urea_c"	;	# 18
	];


	# Metabolic Vmax array (units: mmol/B-hr) -
	metabolic_vmax_array = [
		2.2148976	;	# Vmax [mmol/gdw-hr] 1	M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 2	M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 3	M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 4	M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 5	2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 6	2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 7	[] --> M_Carbamoyl_phosphate_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 8	[] --> M_L-Aspartate_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 9	M_Fumarate_c --> []
		2.2148976	;	# Vmax [mmol/gdw-hr] 10	M_Urea_c --> []
		2.2148976	;	# Vmax [mmol/gdw-hr] 11	[] --> M_ATP_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 12	M_AMP_c --> []
		2.2148976	;	# Vmax [mmol/gdw-hr] 13	M_Diphosphate_c --> []
		2.2148976	;	# Vmax [mmol/gdw-hr] 14	M_Orthophosphate_c --> []
		2.2148976	;	# Vmax [mmol/gdw-hr] 15	[] --> M_Oxygen_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 16	[] --> M_NADPH_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 17	[] --> M_H_c
		2.2148976	;	# Vmax [mmol/gdw-hr] 18	M_Nitric_oxide_c --> []
		2.2148976	;	# Vmax [mmol/gdw-hr] 19	M_NADP_c --> []
		2.2148976	;	# Vmax [mmol/gdw-hr] 20	M_H2O_c --> []
		2.2148976	;	# Vmax [mmol/gdw-hr] 21	[] --> M_H2O_c
	];

	# Metabolic saturation constant array (units mM) -
	number_of_metabolic_rates = length(metabolic_vmax_array)
	metabolic_saturation_constant_array = 0.130*ones(number_of_metabolic_rates*number_of_species)

	# ------------------------------------------------------------------------------------------#
	# constants (from bionumbers)       units
	# ------------------------------------------------------------------------------------------#
	cell_diameter = 1.1                             # mum
	mass_of_single_cell = 2.8e-13                   # g
	number_of_rnapII = 4600            	            # copies/cells
	number_of_ribosome = 50000         	            # copies/cells
	mRNA_half_life_TF = 0.083                       # hrs
	protein_half_life = 70                          # hrs
	infrastructure_half_life = 300					# hrs
	doubling_time_cell = 0.33                       # hrs
	max_translation_rate = 16.5                     # aa/sec
	max_transcription_rate = 60.0                   # nt/sec
	transcription_initiation_time_contstant = 400  # sec
	average_transcript_length = 1200   	            # nt
	average_protein_length = 400       	            # aa
	fraction_nucleus = 0.0             	            # dimensionless
	av_number = 6.02e23                             # number/mol
	avg_gene_number = 2                             # number of copies of a gene
	polysome_number = 4					            # number of ribsomoses per transcript
	# ------------------------------------------------------------------------------------------#
	#
	# ------------------------------------------------------------------------------------------#
	# Calculate constants using bionumber values
	# ------------------------------------------------------------------------------------------#
	# Calculate the volume (convert to L)
	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)

	# Calculate the rnapII_concentration and ribosome_concentration
	rnapII_concentration = number_of_rnapII*(1/av_number)*(1/mass_of_single_cell)*1e6       # mumol/gdw
	ribosome_concentration = number_of_ribosome*(1/av_number)*(1/mass_of_single_cell)*1e6   # mumol/ddw

	# degrdation rate constants -
	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(0.5)                           # hr^-1
	degradation_constant_protein = -(1/protein_half_life)*log(0.5)                        # hr^-1
	degrdation_constant_infrastructure = -(1/infrastructure_half_life)*log(0.5)			# hr^-1

	# kcats for transcription and translation -
	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)            # hr^-1
	kcat_translation = polysome_number*max_translation_rate*(3600/average_protein_length)   # hr^-1

	# kcat for transcription initiation -
	kcat_transcription_initiation = ((1/3600)*transcription_initiation_time_contstant)^-1   # hr^-1
	kcat_translation_initiation = 10*kcat_transcription_initiation                          # hr^-1

	# Maximum specific growth rate -
	maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                          # hr^-1

	# What is the average gene concentration -
	avg_gene_concentration = avg_gene_number*(1/mass_of_single_cell)*(1/V)*1e9              # nmol/gdw

	# How fast do my cells die?
	death_rate_constant = 0.05*maximum_specific_growth_rate                                 # hr^-1

	# Saturation constants for translation and trascription -
	saturation_transcription = 4600*(1/av_number)*(1/mass_of_single_cell)*1e9               # nmol/gdw
	saturation_translation = 150000*(1/av_number)*(1/mass_of_single_cell)*1e6               # nmol/gdw
	# -------------------------------------------------------------------------------------------#


	# Alias the txtl parameters -
	txtl_parameter_dictionary = Dict{AbstractString,Any}()
	txtl_parameter_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
	txtl_parameter_dictionary["ribosome_concentration"] = ribosome_concentration # muM
	txtl_parameter_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1
	txtl_parameter_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1
	txtl_parameter_dictionary["degrdation_constant_infrastructure"] = degrdation_constant_infrastructure  # hr^-1
	txtl_parameter_dictionary["kcat_transcription"] = kcat_transcription  # hr^-1
	txtl_parameter_dictionary["kcat_translation"] = kcat_translation  # hr^-1
	txtl_parameter_dictionary["maximum_specific_growth_rate"] = maximum_specific_growth_rate  # hr^-1
	txtl_parameter_dictionary["death_rate_constant"] = death_rate_constant
	txtl_parameter_dictionary["avg_gene_concentration"] = avg_gene_concentration
	txtl_parameter_dictionary["saturation_constant_transcription"] = saturation_transcription
	txtl_parameter_dictionary["saturation_constant_translation"] = saturation_translation
	txtl_parameter_dictionary["average_transcript_length"] = average_transcript_length
	txtl_parameter_dictionary["average_protein_length"] = average_protein_length
	txtl_parameter_dictionary["kcat_transcription_initiation"] = kcat_transcription_initiation
	txtl_parameter_dictionary["kcat_translation_initiation"] = kcat_translation_initiation

	# Setup species abundance array -
	species_abundance_array = [
		0.0	;	# 1 M_AMP_c
		0.0	;	# 2 M_ATP_c
		0.0	;	# 3 M_Carbamoyl_phosphate_c
		0.0	;	# 4 M_Diphosphate_c
		0.0	;	# 5 M_Fumarate_c
		0.0	;	# 6 M_H2O_c
		0.0	;	# 7 M_H_c
		0.0	;	# 8 M_L-Arginine_c
		0.0	;	# 9 M_L-Aspartate_c
		0.0	;	# 10 M_L-Citrulline_c
		0.0	;	# 11 M_L-Ornithine_c
		0.0	;	# 12 M_N-(L-Arginino)succinate_c
		0.0	;	# 13 M_NADPH_c
		0.0	;	# 14 M_NADP_c
		0.0	;	# 15 M_Nitric_oxide_c
		0.0	;	# 16 M_Orthophosphate_c
		0.0	;	# 17 M_Oxygen_c
		0.0	;	# 18 M_Urea_c
	];

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["default_flux_bounds_array"] = default_bounds_array;
	data_dictionary["species_abundance_array"] = species_abundance_array
	data_dictionary["species_bounds_array"] = species_bounds_array
	data_dictionary["list_of_reaction_strings"] = list_of_reaction_strings
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["is_minimum_flag"] = is_minimum_flag
	data_dictionary["number_of_species"] = number_of_species
	data_dictionary["number_of_reactions"] = number_of_reactions
	data_dictionary["metabolic_saturation_constant_array"] = metabolic_saturation_constant_array
	data_dictionary["metabolic_vmax_array"] = metabolic_vmax_array
	data_dictionary["characteristic_enzyme_abundance_mM"] = 4.7326871794871794e-5
	data_dictionary["volume_of_cell"] = V
	data_dictionary["mass_of_single_cell"] = mass_of_single_cell
	data_dictionary["txtl_parameter_dictionary"] = txtl_parameter_dictionary
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
