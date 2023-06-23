""" 
Semiclassical mean value method applied to the kicked rotor system. 

For latest update, see https://github.com/gabrielmlando/MeanVals
"""
module SemiclassicalMeanValsKRS

using LinearAlgebra: norm, I, mul!, det
using UnPack: @unpack
using StaticArrays: SVector, SMatrix

include("config.jl")
include("propagation.jl")
include("get_Lk.jl")
include("sorters.jl")
include("get_means_cla.jl")
include("get_means_osc.jl")
include("get_means.jl")

end
