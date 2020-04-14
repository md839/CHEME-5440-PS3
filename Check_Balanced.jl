#If E returns a zero matrix, our stoichiometric matrix is balanced

#ensures that CSV files can be used
using CSV

#checks necessary packages are in place
include("Include.jl")

#import atom matrix
A=Matrix(CSV.read("Atom_Matrix.csv",header=0));

#import stoichiometric matrix
S=Matrix(CSV.read("Stoich_Matrix.csv",header=0));

#First 6 columns should equal 0 if balanced
E=transpose(A)*S  