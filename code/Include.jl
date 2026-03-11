# setup paths -
const _ROOT = @__DIR__
const _PATH_TO_SRC = joinpath(_ROOT, "src");
const _PATH_TO_DATA = joinpath(_ROOT, "data");
const _PATH_TO_FIG = joinpath(_ROOT, "figs");

# check: do we need to download any packages?
using Pkg
if (isfile(joinpath(_ROOT, "Manifest.toml")) == false) # have manifest file, we are good. Otherwise, we need to instantiate the environment
    Pkg.activate("."); Pkg.resolve(); Pkg.instantiate(); Pkg.update();
end

# load the required packages -
using Distributions
using Plots
using StatsPlots
using Colors
using LinearAlgebra
using Statistics
using DataFrames
using PrettyTables
using Random
using CSV
using FileIO
using JLD2
using Dates
using NNlib
using CategoricalArrays
using StatsBase
using MultivariateStats
using Downloads

# set the random seed for reproducibility
Random.seed!(1234); # set the random seed for reproducibility

# include the source code for the project -
include(joinpath(_PATH_TO_SRC, "Data.jl"))
include(joinpath(_PATH_TO_SRC, "Compute.jl"))
include(joinpath(_PATH_TO_SRC, "Utilities.jl"))
include(joinpath(_PATH_TO_SRC, "Protein.jl"))