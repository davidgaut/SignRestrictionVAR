
# using Distributed, JLD

# add path location
# rmprocs(collect(2:150))
# addprocs()
# # rmprocs(collect(2:15))
# println(procs())
# @everywhere 
module SignRestrictionVAR

using Distributed, MAT, SharedArrays, Distributed, JLD, StatsBase
using PyCall, PyPlot, Dates, TimeSeries, ExcelReaders, Parameters
using Statistics, DataFrames, Distributions, LinearAlgebra
using FFTW, MatrixEquations, RecursiveArrayTools, ProgressMeter

default_paths = Dict("save_path" =>  pwd()*"\\","data_path" =>  pwd()*"\\")
default_plot  = Dict("show" => false,"save_figure"=>true,"disp_figure"=>true, "filter" => nothing, "annualize" => false,"hatch_plot" => true, "hatch_map" => ["/", " ",".", " ", "+", " ", "x", " ", "//" , "o" , " "], "precision" => 500, "alphaHPD" => 0.9,"SHK"=>true,"HD"=>true,"IRF"=>true,"VD"=>true,"MM"=>true)
default_var   = Dict("date_sample" => [0000-00-00 9999-00-00],"lags" => 3,"series" => :String, "prior" => "NW", "constant" => true, "IterationsNbr" => 1000, "b_sel" => 1:1000)
default_sr    = Dict("drawcount" => 1e3, "draws" => 300, "DrawLeft" => 1:100, "horizon" => 30 , "alpha" => 1 ,"R" => [], "low_band" => 16, "high_band" => 84, "shocks" => nothing)

@with_kw mutable struct Options
   paths        :: Dict = default_paths
   plot         :: Dict = default_plot
   options_var  :: Dict = default_var
   options_sr   :: Dict = default_sr
end

function unpack(D :: Dict)
    K = collect(keys(D))
    V = collect(values(D))

    for ii in 1:length(K)
        eval(:($(Symbol(K[ii])) = $(V[ii])));
    end
return :($collect(keys(D)))
end

# @everywhere 
export

initialize!,
load_data!,
lag_crit!,
var_base,
bvar,
SR_var,
get_draws!,
load_results,
VAR_VD,
do_figures

include("initialize!.jl")
include("load_data!.jl")
include("lag_crit!.jl")
include("var_base.jl")
include("bvar.jl")
include("SR_var.jl")
include("get_draws!.jl")
include("load_results.jl")
include("set_rgb.jl")
include("VAR_VD.jl")
include("do_figures.jl")

end
